import numpy as np
import cv2
from pathlib import Path
from typing import Callable, List, Tuple
from catalog_parse import parse_Hipparcos, parse_tycho_2
from tumble import camera_rotation_from_ra_dec, create_tumbling_path

class camera:
    def __init__(self, aperture: float, fov: float, aspect_ratio: float, 
                 focal_length: float, quantum_efficiency: float, photon_energy: float, 
                 angular_spread_parameter: float, 
                 transmittance: float, pixel_area: float, dark_current: float, 
                 sensor_size: float, width: int, height: int, full_well: int):
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
    def __init__(self, direction: np.ndarray, magnitude: float):
        self.direction = direction                       # 3D unit vector
        self.magnitude = magnitude                       # Visual magnitude

class star_field_centroid:
    def __init__(self, id: int, centroid: np.ndarray):
        self.id = id                                     # Catalog index / star ID
        self.centroid = centroid                         # 2D pixel coordinates [x, y]

class star_field_generator:
    def __init__(self, camera: camera, star_count: int, stars: List[star_field_star], 
                 exposure_time: int, path: Callable[[float], np.ndarray]):
        self.camera = camera
        self.star_count = star_count
        self.stars = stars
        self.exposure_time = exposure_time               # Exposure time in seconds
        self.path = path                                 # Orientation path function

def compute_star_flux(camera: camera, magnitude: float) -> float:
    # Flux density of the datum star vega [w/m2]
    f_0 = 3.5e-9

    flux = (f_0 * camera.transmittance * 0.25 * np.pi * 
            camera.aperture * camera.aperture * (10.0 ** (-0.4 * magnitude)))
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

def super_impose_image(source: np.ndarray, destination: np.ndarray, offset: np.ndarray, scale: float):
    src_h, src_w = source.shape
    dst_h, dst_w = destination.shape
    
    ox, oy = offset[0], offset[1]
    
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
            
            if sx < 0 or sx + 1 >= src_w or sy < 0 or sy + 1 >= src_h:
                if 0 <= sx < src_w and 0 <= sy < src_h:
                    destination[dst_y, dst_x] += source[sy, sx] * scale
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
            
            destination[dst_y, dst_x] += interpolated_val * scale

def generate_star_field(gen: star_field_generator, start_time: float, end_time: float) -> Tuple[np.ndarray, List[star_field_centroid]]:
    ## Maximum distance that gaussians can be placed apart between sub-frames 
    gaussian_spacing = 1.5
    
    start_matrix = gen.path(start_time)
    end_matrix = gen.path(end_time)
    
    ## Estimate the largest angle a point in the image can sweep out during this frame by finding the maximum angular change of the corners from the start to the end of the frame
    y_offset = np.tan(gen.camera.fov / 2.0)
    x_offset = gen.camera.aspect_ratio * y_offset
    raw_corners = [
        np.array([-x_offset,  y_offset, 1.0]),
        np.array([ x_offset,  y_offset, 1.0]),
        np.array([-x_offset, -y_offset, 1.0]),
        np.array([ x_offset, -y_offset, 1.0])
    ]
    
    ## normalize corner vectors:
    corners = [c / np.linalg.norm(c) for c in raw_corners] 

    ## Find which corner moved the most and use it as a heuristic for the maximmum possible distance a point can move in on frame
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
    
    ## Approximate anglular size out of a single pixel
    pixel_angle = gen.camera.fov / gen.camera.height
    
    ## 1 standard deviation of the star's 2D gaussian distribution measured in pixels
    sigma_pixels = gen.camera.angular_spread_parameter / pixel_angle

    ## 1 standard deviation of the star's 2D gaussian distribution measured in m
    sigma_sensor = gen.camera.sensor_size * gen.camera.angular_spread_parameter / gen.camera.fov

    ## Minimum number of sub-frames that need to stitched together this frame to ensure the gaussians are less than one [gaussian_spacing] apart (ensures smooth streaks instead of dotted lines)
    steps = max(1, int(np.ceil(largest_angle / (gaussian_spacing * gen.camera.angular_spread_parameter))))
    time_step = (end_time - start_time) / steps
    
    ## This bitmaps is not the acutual brightness observed in the image, but rather the average power flux for each pixel 
    flux_bitmap = np.zeros((gen.camera.height, gen.camera.width), dtype=np.float32)
    
    epsilon = np.sin(gen.camera.fov * 0.5)
    pixel_size = gen.camera.sensor_size / gen.camera.height
    template_size = int(np.ceil(2 * 5 * sigma_pixels))
    if template_size % 2 == 0:
        template_size += 1
    gaussian_template = generate_gaussian_image(template_size, template_size * pixel_size, sigma_sensor, gen.camera.pixel_area)
    
    visible_centroids: List[star_field_centroid] = []
    
    print(steps);


    time = start_time
    for step_idx in range(steps):
        orientation = gen.path(time)
        ##print(orientation);
        for j, star in enumerate(gen.stars):
            position = orientation.T @ star.direction
            if position[2] > 0.0 and abs(position[0]) < gen.camera.aspect_ratio * epsilon and abs(position[1]) < epsilon:
                norm_x = position[0] / position[2]
                norm_y = position[1] / position[2]
                img_x = ((norm_x / (gen.camera.aspect_ratio * epsilon)) + 1.0) * gen.camera.width * 0.5
                img_y = ((norm_y / epsilon) + 1.0) * gen.camera.height * 0.5
                
                sub_frame_energy = time_step * compute_star_flux(gen.camera, star.magnitude)
                offset = np.array([img_x - template_size * 0.5, img_y - template_size * 0.5])
                
                super_impose_image(gaussian_template, flux_bitmap, offset, sub_frame_energy)
                visible_centroids.append(star_field_centroid(id=j, centroid=np.array([img_x, img_y])))
                
        time += time_step
        
    ## This is bitmap contains the actual raw digital value recorded by the camera, accounting for shot noise and thermal fluctuations
    brightness_bitmap = np.zeros((gen.camera.height, gen.camera.width), dtype=np.uint16)
    thermal_noise_mean = (end_time - start_time) * gen.camera.dark_current
    
    rng = np.random.default_rng(0)
    for y in range(gen.camera.height):
        for x in range(gen.camera.width):
            pixel_energy = flux_bitmap[y, x]
            mean_photons = pixel_energy / gen.camera.photon_energy
            mean_electrons = gen.camera.quantum_efficiency * mean_photons
            
            thermal_electrons = rng.poisson(thermal_noise_mean)
            signal_electrons = rng.poisson(mean_electrons) if mean_electrons > 0 else 0
            
            total_electrons = thermal_electrons + signal_electrons
            digital_val = min(gen.camera.full_well, 100*total_electrons)
            brightness_bitmap[y, x] = digital_val
            
    return brightness_bitmap, visible_centroids

def save_simulation_output(camera: camera, exposure_time: float, bitmap_img: np.ndarray, centroids: List[star_field_centroid]):
    folder_name = (
        f"fov_{np.degrees(camera.fov):.1f}deg_"
        f"exp_{exposure_time}s_"
        f"ap_{camera.aperture}m_"
        f"foc_{camera.focal_length}m_"
        f"res_{camera.width}x{camera.height}"
    )
    
    output_dir = Path("synthetic_output") / folder_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save PNG Image
    img_path = output_dir / "simulated_star_field.png"
    cv2.imwrite(str(img_path), bitmap_img)
    
    # Save Metadata Text File
    txt_path = output_dir / "star_metadata.txt"
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write(f"=== Star Tracker Synthetic Image Metadata ===\n")
        f.write(f"Total Visible Star Instances: {len(centroids)}\n")
        f.write(f"{'Catalog Index':<15} | {'Centroid X (pixels)':<22} | {'Centroid Y (pixels)'}\n")
        f.write("-" * 65 + "\n")
        
        for c in centroids:
            f.write(f"{c.id:<15} | {c.centroid[0]:<22.4f} | {c.centroid[1]:<22.4f}\n")
            
    print(f"Files saved to: {output_dir}")

if __name__ == "__main__":
    magnitude_limit = 8.0

    ##catalog_filename = "hipparcos.txt"
    ##print(f"Loading star catalog from {catalog_filename}...")
    ##ids, mags, vectors = parse_Hipparcos(catalog_filename)
    
    catalog_filename = "tycho_2.dat"
    print(f"Loading star catalog from {catalog_filename}...")
    ids, mags, vectors = parse_tycho_2(catalog_filename, magnitude_limit)

    catalog_stars = []
    for star_id, mag, vec in zip(ids, mags, vectors):
        if mag <= magnitude_limit:
            catalog_stars.append(star_field_star(direction=vec, magnitude=mag))
            
    cam = camera(
        aperture=0.05,
        fov=np.radians(30.0),
        aspect_ratio=1.333,
        focal_length=0.1,
        quantum_efficiency=0.8,
        photon_energy=3.5e-19,
        angular_spread_parameter=0.0005,
        transmittance=0.9,
        pixel_area=1.0e-11,
        dark_current=10.0,
        sensor_size=0.01,
        width=1024,
        height=768,
        full_well=65535
    )
    
    ## starting orientation:
    initial_rotation = camera_rotation_from_ra_dec(np.radians(266), np.radians(-28.9));

    ## Angular velocity
    tumbling_omega = [0.008, 0.00, 0.002]  # radians / second
    
    path = create_tumbling_path(omega_vec=tumbling_omega, initial_rotation = initial_rotation)
    
    generator = star_field_generator(
        camera=cam,
        star_count=len(catalog_stars),
        stars=catalog_stars,
        exposure_time=2,
        path=path
    )
    
    print("Generating synthetic star field...")
    image_bitmap, centroids = generate_star_field(generator, start_time=0.0, end_time=2.0)
    
    save_simulation_output(cam, exposure_time=2, bitmap_img=image_bitmap, centroids=centroids)