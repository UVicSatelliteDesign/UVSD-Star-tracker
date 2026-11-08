import numpy as np
import cv2
from pathlib import Path
from typing import Callable, List, Tuple
from catalog_parse import parse_Hipparcos, parse_tycho_2
from tumble import camera_rotation_from_ra_dec, create_tumbling_path

class camera:
    def __init__(self, aperture: float, fov: float, aspect_ratio: float, focal_length: float, quantum_efficiency: float, photon_energy: float, angular_spread_parameter: float, transmittance: float, pixel_area: float, dark_current: float, sensor_size: float, width: int, height: int, full_well: int):
        self.aperture = aperture                         # Lens diameter in meters
        self.fov = fov                                   # Vertical field of view in radians
        self.aspect_ratio = aspect_ratio                 # Ratio of width to height
        self.focal_length = focal_length                 # Focal length in meters
        self.quantum_efficiency = quantum_efficiency     # Proportion of photons converted to electrons
        self.photon_energy = photon_energy               # Average energy of a photon
        self.angular_spread_parameter = angular_spread_parameter # 1 std dev of star angular spread
        self.transmittance = transmittance               # Optics transmittance (0 to 1)
        self.pixel_area = pixel_area                     # Physical area of a single pixel (m^2)
        self.dark_current = dark_current                 # Thermal noise in electrons/sec/pixel
        self.sensor_size = sensor_size                   # Height of sensor in meters
        self.width = width                               # Image width in pixels
        self.height = height                             # Image height in pixels
        self.full_well = full_well                       # Electrons required to saturate a pixel

class star_field_star:
    def __init__(self, id: int, direction: np.ndarray, magnitude: float):
        self.id = id
        self.direction = direction                       # 3D unit vector
        self.magnitude = magnitude                       # Visual magnitude

class star_field_centroid:
    def __init__(self, id: int, centroid: np.ndarray,magnitude: float):
        self.id = id                                     # Catalog index / star ID
        self.centroid = centroid                         # 2D pixel coordinates [x, y]
        self.magnitude = magnitude

class star_field_generator:
    def __init__(self, camera: camera, star_count: int, stars: List[star_field_star], path: Callable[[float], np.ndarray]):
        self.camera = camera
        self.star_count = star_count
        self.stars = stars
        self.path = path                                 # Orientation path function

def orientation_to_ra_dec_roll(rot):
    # this treats the camera orientation as a passive body fixed yxz rotation sequence (ra about y, dec about x, roll about z) and extracts the angles
    ra = np.atan2(rot[2, 1], rot[2, 0]) 
    dec = np.pi / 2 - np.arccos(rot[2, 2])
    roll = np.atan2(rot[1, 2], -rot[0, 2]) - np.pi / 2

    return ra, dec, roll

def compute_star_flux(camera: camera, magnitude: float) -> float:
    # flux density of the datum star vega [w/m2]
    f_0 = 3.5e-9
    
    # flux density of this star
    f = f_0 * 10.0 ** (-0.4 * magnitude)

    flux = f * camera.transmittance * 0.25 * np.pi * camera.aperture * camera.aperture
    return flux

def generate_gaussian_image(image_size: int, scale: float, standard_deviation: float, magnitude: float) -> np.ndarray:
    image = np.zeros((image_size, image_size), dtype=np.float32)
    normalization_factor = magnitude / (2.0 * np.pi * standard_deviation * standard_deviation)
    exponent_coefficient = -0.5 / (standard_deviation * standard_deviation)
    
    step = scale / image_size
    x = -scale * 0.5 + 0.5 * step
    
    for i in range(image_size):
        y = -scale * 0.5 + 0.5 * step
        for j in range(image_size):
            squared_dist = x * x + y * y
            intensity = normalization_factor * np.exp(exponent_coefficient * squared_dist)
            image[j, i] = intensity
            y += step
        x += step
        
    return image

def super_impose_image(source: np.ndarray, destination: np.ndarray, offset: np.ndarray, gain: float):
    # get dimensions of the source and destination images
    src_h, src_w = source.shape
    dst_h, dst_w = destination.shape
    
    ox, oy = offset[0], offset[1]
    
    # crop source image to fit in the destination
    x_min = max(0, int(np.floor(ox)))
    y_min = max(0, int(np.floor(oy)))
    x_max = min(dst_w, int(np.ceil(ox + src_w)))
    y_max = min(dst_h, int(np.ceil(oy + src_h)))
    
    for dst_y in range(y_min, y_max):
        for dst_x in range(x_min, x_max):
            src_x_f = dst_x - ox
            src_y_f = dst_y - oy
            
            sx = int(np.floor(src_x_f))
            sy = int(np.floor(src_y_f))
            
            fx = src_x_f - sx
            fy = src_y_f - sy
            
            # interpolation samples the 2x2 grid of source image values nearest to the sub-pixel point of interest (sx, sy). This checks if that 2x2 window extends past the source image
            if sx < 0 or sx + 1 >= src_w or sy < 0 or sy + 1 >= src_h:
                
                # enherit boundary value
                if 0 <= sx < src_w and 0 <= sy < src_h:
                    destination[dst_y, dst_x] += source[sy, sx] * gain
                continue
            
            v00 = source[sy, sx]
            v10 = source[sy, sx + 1]
            v01 = source[sy + 1, sx]
            v11 = source[sy + 1, sx + 1]
            
            interpolated_val = (
                v00 * (1.0 - fx) * (1.0 - fy) +
                v10 * fx * (1.0 - fy) +
                v01 * (1.0 - fx) * fy +
                v11 * fx * fy
            )
            
            destination[dst_y, dst_x] += interpolated_val * gain

def generate_star_field(gen: star_field_generator, start_time: float, end_time: float) -> Tuple[np.ndarray, List[star_field_centroid], np.ndarray, np.ndarray]:
    # !Important!: this is just used to make to make the output images easier to see.
    # set this to 1 when doing actual tests
    debugging_gain = 3


    # maximum distance that gaussians can be placed apart between sub-frames 
    gaussian_spacing = 1.5
    
    start_matrix = gen.path(start_time)
    end_matrix = gen.path(end_time)
    
    # estimate the largest angle a point in the image can sweep out during this frame by finding the maximum angular change of the corners from the start to the end of the frame
    y_offset = np.tan(gen.camera.fov / 2.0)
    x_offset = gen.camera.aspect_ratio * y_offset
    raw_corners = [
        np.array([-x_offset,  y_offset, 1.0]),
        np.array([ x_offset,  y_offset, 1.0]),
        np.array([-x_offset, -y_offset, 1.0]),
        np.array([ x_offset, -y_offset, 1.0])
    ]
    
    # normalize corner vectors:
    corners = [c / np.linalg.norm(c) for c in raw_corners] 

    # find which corner moved the most and use it as a heuristic for the maximmum possible distance a point can move in on frame
    largest_angle = 0.0
    for c in corners:
        rotated_start = start_matrix @ c
        rotated_end = end_matrix @ c
        dot_product = np.dot(rotated_start, rotated_end)
        
        # clamp value to valid range [-1.0, 1.0] to prevent math domain errors in arccos
        clamped_dot = np.clip(dot_product, -1.0, 1.0)
        angle = np.arccos(clamped_dot)
        
        if angle > largest_angle:
            largest_angle = angle
    
    # approximate anglular size of a single pixel
    pixel_angle = gen.camera.fov / gen.camera.height
    
    # 1 standard deviation of the star's 2D gaussian distribution measured in pixels
    sigma_pixels = gen.camera.angular_spread_parameter / pixel_angle

    # 1 standard deviation of the star's 2D gaussian distribution measured in m
    sigma_sensor = gen.camera.sensor_size * gen.camera.angular_spread_parameter / gen.camera.fov

    # minimum number of sub-frames that need to be stitched together this frame to ensure the gaussians are less than one [gaussian_spacing] apart (ensures smooth streaks instead of dotted lines)
    steps = max(1, int(np.ceil(largest_angle / (gaussian_spacing * gen.camera.angular_spread_parameter))))
    time_step = (end_time - start_time) / steps
    
    # this bitmaps is not the actual brightness observed in the image, but rather the average power flux for each pixel 
    flux_bitmap = np.zeros((gen.camera.height, gen.camera.width), dtype=np.float32)
    
    epsilon = np.tan(gen.camera.fov * 0.5)
    pixel_size = gen.camera.sensor_size / gen.camera.height
    
    # make sure the image is large enough to capture nearly the entire normal distribution
    template_size = int(np.ceil(2 * 5 * sigma_pixels))
    if template_size % 2 == 0:
        template_size += 1
    gaussian_template = generate_gaussian_image(template_size, template_size * pixel_size, sigma_sensor, gen.camera.pixel_area)
    
    visible_centroids = []
    start_orientation = gen.path(start_time)
    end_orientation = gen.path(end_time)
    print(f"# of sub-frames needed: {steps}")
    time = start_time 
    for step_idx in range(steps):
        orientation = gen.path(time)

        for j, star in enumerate(gen.stars):
            position = orientation @ star.direction
            
            # check if the star is infront on the camera and in the field of view 
            if position[2] > 0.0 and abs(position[0]) < gen.camera.aspect_ratio * epsilon and abs(position[1]) < epsilon:
                norm_x = -position[0] / position[2]
                norm_y = -position[1] / position[2]
                
                img_x = ((norm_x / (gen.camera.aspect_ratio * epsilon)) + 1.0) * gen.camera.width * 0.5
                img_y = ((norm_y / epsilon) + 1.0) * gen.camera.height * 0.5
                
                sub_frame_energy = time_step * compute_star_flux(gen.camera, star.magnitude)

                # subtract half of the template's size in the x and y directions to get the offset of the top left corner
                offset = np.array([img_x - template_size * 0.5, img_y - template_size * 0.5])
                
                super_impose_image(gaussian_template, flux_bitmap, offset, sub_frame_energy)
                visible_centroids.append(star_field_centroid(id=star.id, centroid=np.array([img_x, img_y]), magnitude=star.magnitude))
                
        time += time_step

        # Calculate completion percentage
        percent = ((step_idx + 1) / steps) * 100
        
        # Print with \r to overwrite the line
        print(f"\rProgress: sub-frame {step_idx + 1}/{steps} [{percent:5.1f}%]", end="", flush=True)

    print("\n")
    # track unique centroid IDs using a set to avoid duplicates    
    unique_centroid_IDs = set()

    # store centroids in a dictionary
    mean_centroids = {}
    instance_counts = {}
    magnitudes = {}
    for vc in visible_centroids:
        if vc.id not in unique_centroid_IDs:
            unique_centroid_IDs.add(vc.id)
            mean_centroids[vc.id] = vc.centroid
            instance_counts[vc.id] = 1
            magnitudes[vc.id] = vc.magnitude
        else:
            mean_centroids[vc.id] = mean_centroids[vc.id] + vc.centroid
            instance_counts[vc.id] = instance_counts[vc.id] + 1
    
    # compute the mean centroid values for each unique star
    compressed_centroids = []
    for ID in unique_centroid_IDs:
        mean_coord = mean_centroids[ID] / instance_counts[ID]
        
        compressed_star = star_field_centroid(id=ID, centroid=mean_coord, magnitude= magnitudes[ID])

        compressed_centroids.append(compressed_star)

    visible_centroids = compressed_centroids

    # this is bitmap contains the actual raw digital value recorded by the camera, accounting for shot noise and thermal fluctuations
    brightness_bitmap = np.zeros((gen.camera.height, gen.camera.width), dtype=np.uint16)
    thermal_noise_mean = (end_time - start_time) * gen.camera.dark_current
    
    rng = np.random.default_rng(0)
    for y in range(gen.camera.height):
        for x in range(gen.camera.width):
            # mean energies
            pixel_energy = flux_bitmap[y, x]
            mean_photons = pixel_energy / gen.camera.photon_energy
            mean_electrons = gen.camera.quantum_efficiency * mean_photons
            
            # since electrons are discrete, arrive randomly and independantly, the number of them reaching the sensor over a given timespan can be found by sampling a poisson distribution 
            thermal_electrons = rng.poisson(thermal_noise_mean)
            signal_electrons = rng.poisson(mean_electrons) if mean_electrons > 0 else 0
            
            total_electrons = thermal_electrons + signal_electrons

            # limit brightness value to saturation threshold 
            digital_val = min(gen.camera.full_well, debugging_gain * total_electrons)
            brightness_bitmap[y, x] = digital_val
            
    return brightness_bitmap, visible_centroids, start_orientation, end_orientation

def save_simulation_output(camera: camera, exposure_time: float, ra: float, dec: float, roll: float, bitmap_img: np.ndarray, centroids: List[star_field_centroid]):
    folder_name = (
        f"exp_{exposure_time}s_"
        f"_RA_{np.degrees(ra):.6f}"
        f"_DEC_{np.degrees(dec):.6f}"
        f"_ROLL_{np.degrees(roll):.6f}"
    )
    
    output_dir = Path("synthetic_output") / folder_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # save PNG Image
    img_path = output_dir / "simulated_star_field.png"
    cv2.imwrite(str(img_path), bitmap_img)
    print(f"Synthetic star field saved to: {img_path}")
    
    # save centroid metadata text file
    txt_path = output_dir / "centroid_metadata.txt"
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write(f"=== Virtual Camera Settings ===\n")
        f.write(f"Fov: {np.degrees(camera.fov):.3f} deg\n")
        f.write(f"Aperture: {camera.aperture}m\n")
        f.write(f"Exposure: {exposure_time}s\n")
        f.write(f"Resolution: {camera.width}x{camera.height}\n")
        
        f.write(f"\n=== Orientation ===\n")
        f.write(f"RA: {np.degrees(ra):.8f} deg\n")
        f.write(f"DEC: {np.degrees(dec):.8f} deg\n")
        f.write(f"ROLL: {np.degrees(roll):.8f} deg\n")

        f.write(f"\n=== Centroids ===\n")
        f.write(f"Number of visible stars: {len(centroids)}\n")
                
        for c in centroids:
            f.write(f"{c.id},{c.centroid[0]:.4f},{c.centroid[1]:.4f}\n")
            
    save_annotated_simulation_output(folder_name, exposure_time=2, bitmap_img=image_bitmap, centroids=centroids, annotation_mag_limit=7.0)

def magnitude_to_heatmap_bgr(magnitude: float, min_mag: float, max_mag: float) -> tuple:
    # normalize magnitude to a 0.0 - 1.0 range, inverted so bright = 1.0
    t = np.clip((max_mag - magnitude) / (max_mag - min_mag), 0.0, 1.0)
    
    # piecewise linear interpolation
    if t < 0.25:
        # black to blue
        u = t * 4.0
        b, g, r = int(127 + 127 * u), 0, 0
    elif t < 0.5:
        # blue to cyan
        u = (t - 0.25) * 4.0
        b, g, r = 255, int(255 * u), 0
    elif t < 0.75:
        # cyan to green/yellow
        u = (t - 0.5) * 4.0
        b, g, r = int(255 * (1.0 - u)), 255, 0
    else:
        # yellow to red
        u = (t - 0.75) * 4.0
        b, g, r = 0, int(255 * (1.0 - u)), 255
        
    return (b, g, r)
def save_annotated_simulation_output(folder, exposure_time: float, bitmap_img: np.ndarray, centroids: List[star_field_centroid], annotation_mag_limit: float = 9.0):
    output_dir = Path("synthetic_output") / folder
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # convert 16-bit or grayscale bitmap to an 8-bit BGR image for OpenCV drawing
    if bitmap_img.dtype == np.uint16:
        img_8bit = cv2.normalize(bitmap_img, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    else:
        img_8bit = bitmap_img.copy()
        
    if len(img_8bit.shape) == 2:
        annotated_img = cv2.cvtColor(img_8bit, cv2.COLOR_GRAY2BGR)
    else:
        annotated_img = img_8bit.copy()
        
    # draw index labels and markers for stars brighter than the magnitude limit
    for c in centroids:
        if c.magnitude <= annotation_mag_limit:
            cx_int = int(round(c.centroid[0]))
            cy_int = int(round(c.centroid[1]))
            
            # draw a red dot at the centroid position
            cv2.circle(annotated_img, (cx_int, cy_int), 1, (0, 0, 255), -1)
            
            # draw the star ID/index text next to the centroid
            # also why not make the label colours fancy
            min_mag = 2
            max_mag = annotation_mag_limit

            label_color = magnitude_to_heatmap_bgr(c.magnitude, min_mag, max_mag)

            label = f"ID:{c.id} (M:{c.magnitude:.1f})"
            cv2.putText(
                annotated_img, 
                label, 
                (cx_int + 5, cy_int + 2), 
                cv2.FONT_HERSHEY_SIMPLEX, 
                0.35, 
                label_color,
                1, 
                cv2.LINE_AA
            )
            
    # save the annotated output image
    annotated_path = output_dir / "simulated_star_field_annotated.png"
    cv2.imwrite(str(annotated_path), annotated_img)
    print(f"Annotated star field image saved to: {annotated_path}")

if __name__ == "__main__":
    magnitude_limit = 10.0

    ##catalog_filename = "hipparcos.txt"
    ##print(f"Loading star catalog from {catalog_filename}...")
    ##ids, mags, vectors = parse_Hipparcos(catalog_filename)
    
    catalog_filename = "tycho_2.dat"
    print(f"Loading star catalog from {catalog_filename}...")
    ids, mags, vectors = parse_tycho_2(catalog_filename, magnitude_limit)

    catalog_stars = []
    for star_id, mag, vec in zip(ids, mags, vectors):
        if mag <= magnitude_limit:
            catalog_stars.append(star_field_star(id = star_id, direction = vec, magnitude = mag))
            
    cam = camera(
        aperture = 0.05,
        fov = np.radians(30.0),
        aspect_ratio = 1.333,
        focal_length = 0.1,
        quantum_efficiency = 0.8,
        photon_energy = 3.5e-19,
        angular_spread_parameter = 0.0002,
        transmittance = 0.9,
        pixel_area = 1.0e-11,
        dark_current = 10.0,
        sensor_size = 0.01,
        width = 1024,
        height = 768,
        full_well = 65535
    )
    
    ## starting orientation:
    initial_rotation = camera_rotation_from_ra_dec(np.radians(56.66), np.radians(24.18), np.radians(15))

    ## angular velocity
    tumbling_omega = [0.0001, 0.0001, 0.0001]  # radians / second
    
    path = create_tumbling_path(omega_vec = tumbling_omega, initial_rotation = initial_rotation)

    generator = star_field_generator(
        camera = cam,
        star_count = len(catalog_stars),
        stars = catalog_stars,
        path = path
    )
    
    print("Generating synthetic star field...")
    image_bitmap, centroids, start_orientation, end_orientation = generate_star_field(generator, start_time = 0.0, end_time = 2.0)
    
    # convert the start and end orientations to RA, DEC, and ROLL angles
    start_ra, start_dec, start_roll = orientation_to_ra_dec_roll(start_orientation)
    end_ra, end_dec, end_roll = orientation_to_ra_dec_roll(end_orientation)

    # compute the average RA, DEC, and ROLL for the exposure
    avg_ra = (start_ra + end_ra) / 2.0
    avg_dec = (start_dec + end_dec) / 2.0
    avg_roll = (start_roll + end_roll) / 2.0

    print(start_orientation)
    print(np.degrees(avg_ra), np.degrees(avg_dec), np.degrees(avg_roll))
    save_simulation_output(cam, exposure_time=2, ra = avg_ra, dec = avg_dec, roll = avg_roll, bitmap_img=image_bitmap, centroids=centroids)