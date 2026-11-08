import cv2
import numpy as np
from pathlib import Path
from centroids import extract_centroids, image_centroid
from feature_extraction import extract_features, save_features_to_csv, load_features_from_csv

class metadata_centroid:
    def __init__(self, id: int, x: float, y: float):
        self.id = id
        self.x = x
        self.y = y

def mark_centroids(folder_path: str, image: np.ndarray, centroids: list[image_centroid]):
    folder = Path(folder_path)

    # prepare image for annotation (convert 16-bit or grayscale to an 8-bit BGR format for OpenCV drawing)
    if image.dtype == np.uint16:
        img_8bit = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    else:
        img_8bit = image.copy()
        
    # check if the image is gray scale (it might no be if we start testing synthetic images with colour)
    # if there are more colour channels, the image will be a 3-D array
    if len(img_8bit.shape) == 2:
        annotated_img = cv2.cvtColor(img_8bit, cv2.COLOR_GRAY2BGR)
    else:
        annotated_img = img_8bit.copy()
        
    # draw bounding boxes and crosshairs for each detected star
    for i, c in enumerate(centroids):
        left, top, width, height = c.bounding_box
        
        # draw green bounding box around the detected star blob
        cv2.rectangle(annotated_img, (left, top), (left + width, top + height), (255, 100, 0), 1)
        
        # draw a red dot at the exact sub-pixel centroid coordinate
        cx = int(round(c.x))
        cy = int(round(c.y))
        cv2.circle(annotated_img, (cx, cy), 1, (0, 0, 255), -1)
        
        # add a lebel for identifying the centroid
        label = f"{i}"
        cv2.putText(
            annotated_img, 
            label, 
            (cx + 5, cy + 2), 
            cv2.FONT_HERSHEY_SIMPLEX, 
            0.35, 
            (0, 255, 0),
            1, 
            cv2.LINE_AA
        )

    # save the annotated output back into the simulation folder
    output_path = folder / "annotated_star_field.png"
    cv2.imwrite(str(output_path), annotated_img)
    print(f"Annotated image successfully saved to: {output_path}")

def load_metadata(metadata_path): 
    print(f"Loading metadata from {metadata_path}...")

    path = Path(metadata_path)
    if not path.exists():
        print("Metadata file not found")
        return

    camera_settings = {}
    orientation = {}
    centroids = []
    
    current_section = None

    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            #skip empty lines
            if not line:
                continue
            
            # determine when one section transitions to another 
            if line.startswith("=== Virtual Camera Settings ==="):
                current_section = "camera"
                continue
            elif line.startswith("=== Orientation ==="):
                current_section = "orientation"
                continue
            elif line.startswith("=== Centroids ==="):
                current_section = "centroids"
                continue
            
            # parse data
            if current_section == "camera":
                if ":" in line:
                    key, val = line.split(":", 1)
                    camera_settings[key.strip()] = val.strip()
                    
            elif current_section == "orientation":
                if ":" in line:
                    key, val = line.split(":", 1)
                    orientation[key.strip()] = val.strip()
                    
            elif current_section == "centroids":
                if line.startswith("Number of visible stars:"):
                    continue
                
                parts = line.split(",")
                if len(parts) == 3:
                    catalog_index = int(parts[0].strip())
                    x_coord = float(parts[1].strip())
                    y_coord = float(parts[2].strip())
                    
                    centroids.append(metadata_centroid(
                        id = catalog_index,
                        x = x_coord,
                        y = y_coord
                    ))

    return camera_settings, orientation, centroids

def match_centroids_to_metadata(true_centroids, extracted_centroids):
    matches = []
    
    # for each extracted centroid, find the closest true centroid based on Euclidean distance
    for c in extracted_centroids:
        closest_dist = np.inf
        closest_candidate = -1
        for true_c in true_centroids:
            dist_sq = (c.x - true_c.x) ** 2 + (c.y - true_c.y) ** 2
            if dist_sq < closest_dist:
                closest_dist = dist_sq
                closest_candidate = true_c.id
        matches.append(closest_candidate)

    return matches

def score_catalog_candidates(observed_angles, catalog_ids, catalog_angles, sigma):
    # angle differences

    # this difference matrix is an MxNmxm matrix storing the differences between each observed angle and each catalog angle for all neighbors

    # note: M is the number of catalog entries, N is the number of observed angles, and m is the number of angles per catalog entry
    diff = catalog_angles[:, None, :] - observed_angles[None, :, None]
    
    # square angle differences and apply exponential
    squared_diff = diff ** 2
    gaussian_weights = np.exp(-squared_diff / (2.0 * (sigma ** 2)))
    
    # sum over catalog neighbours (second axis)
    inner_sum = np.sum(gaussian_weights, axis=2)
    
    # add a tiny value near machine epsilon to prevent log(0) / negative infinity
    inner_sum = np.maximum(inner_sum, 1e-15)
    
    # apply logarithm
    log_terms = np.log(inner_sum)
    
    # sum up contributions from all observed neighbours
    log_probabilities = np.sum(log_terms, axis=1)
    
    # sort indices by highest log-probability in descending order
    sorted_indices = np.argsort(log_probabilities)[::-1]
    
    # convert catalog indices and neighbour distances to ints and float, otherwise when they are printed, they will be formatted np.float64
    sorted_ids = []
    sorted_probabilities = []
    for idx in sorted_indices:
        sorted_ids.append(int(catalog_ids[idx]))
        sorted_probabilities.append(float(log_probabilities[idx]))

    return sorted_ids, sorted_probabilities

def process_simulation_folder(folder_path: str, noise_threshold: float = 50.0):
    folder = Path(folder_path)

    # load the image and metadata
    img_path = folder / "simulated_star_field.png"
    metadata_path = folder / "centroid_metadata.txt"

    # exit prematurely if the target files cannot be located
    if not img_path.exists():
        print(f"Error: Could not find 'simulated_star_field.png' in {folder_path}")
        return

    if not metadata_path.exists():
        print(f"Error: Could not find 'centroid_metadata.txt' in {folder_path}")
        return
        
    print(f"Loading image from {img_path}...")

    # using IMREAD_UNCHANGED here to preserve the full 16-bit dynamic range
    image = cv2.imread(str(img_path), cv2.IMREAD_UNCHANGED)
    
    # exit prematurely if the source image is corrupted
    if image is None:
        print("Error: Failed to load image data.")
        return
    

    camera_settings, orientation, true_centroids = load_metadata(metadata_path)

    print("Extracting centroids...")
    centroids = extract_centroids(image, noise_threshold=noise_threshold)
    print(f"Successfully extracted {len(centroids)} star centroids.")
    
    matches = match_centroids_to_metadata(true_centroids, centroids)

    # separate the vector part of the centroid from the magnitude and ID information
    centroids_vectors = [[c.x, c.y] for c in centroids]

    height, width = image.shape[:2]

    image_angles, image_neighbors = extract_features(
            centroids = centroids_vectors, 
            n_neighbors=10, 
            fov = np.radians(30),
            width = width,
            height = height
        )

    
    mark_centroids(folder_path, image, centroids)
    feature_path = folder / "features.csv"
    # the save_features_to_csv function needs labels for each centroid. Since the extracted centroids are not linked to the catalog, we can just use a range of integers as labels for now. 
    star_ids = list(range(len(centroids)))
    save_features_to_csv(feature_path, star_ids=star_ids, neighbor_indices=image_neighbors, angular_matrix=image_angles)
    
    database_ids, database_neighbors, database_angles = load_features_from_csv("./database/features.csv")

    for i in range(90, 100):
        ids, probabilities = score_catalog_candidates(
            observed_angles=image_angles[i],
            catalog_ids=database_ids,
            catalog_angles=database_angles,
            sigma=0.001
        )
        match = matches[i]
        print(f"{i}, {match}")
        for j in range(10):
            id = ids[j]
            if id == match:
                print(f"  {j}: {ids[j]}, {probabilities[j]:.6f} <-- [MATCH]")
            else:
                print(f"  {j}: {ids[j]}, {probabilities[j]:.6f}")

if __name__ == "__main__":
    path = "./synthetic_output/exp_2s__RA_56.652308_DEC_24.184052_ROLL_14.997422"
    process_simulation_folder(str(path), noise_threshold=14000.0)