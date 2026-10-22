import csv
import numpy as np

def centroids_to_unit_vectors(centroids, fov: float, aspect_ratio: float, width: int, height: int):
    half_sin = np.sin(0.5 * fov)
    half_cos = np.cos(0.5 * fov)
    
    unit_vectors = []
    for center_x_val, center_y_val in centroids:
        v = np.array([
                2 * half_sin * (center_x_val - 0.5 * width) / height, 
                2 * half_sin * (center_y_val - 0.5 * height) / height, 
                -half_cos
            ])
        
        # normalize to get a unit vector
        v /= np.linalg.norm(v)
        unit_vectors.append(v)
        
    return np.array(unit_vectors)

def save_features_to_csv(filename, star_ids, neighbor_indices, angular_matrix):
    with open(filename, mode='w', newline='') as f:
        writer = csv.writer(f)
        
        for i, star_id in enumerate(star_ids):
            row = [star_id]
            
            for j_idx, neighbor_row_idx in enumerate(neighbor_indices[i]):
                neighbor_id = star_ids[neighbor_row_idx]
                distance = angular_matrix[i, j_idx]
                
                row.append(neighbor_id)
                row.append(f"{distance:.6f}")
                
            writer.writerow(row)
            
    print(f"Successfully exported feature descriptors to {filename}")

def extract_features(centroids, n_neighbors, fov: float, aspect_ratio: float, width: int, height: int):
    centroids = np.array(centroids)
    m = len(centroids)
    
    # make sure there are enough star in the image
    if m <= n_neighbors:
            print("Not enough stars")
            return np.array([]), np.array([]) 

    # compute pairwise Euclidean distances in pixel space to find neighbors
    diff = centroids[:, np.newaxis, :] - centroids[np.newaxis, :, :]  # Shape: (m, m, 2)

    pixel_distances_sq = np.sum(diff ** 2, axis=2)  # Shape: (m, m)
    # sort each row by pixel distance and get indices
    sorted_indices = np.argsort(pixel_distances_sq)
    
    # exclude index 0 (the reference star itself) and extract the n closest neighbors
    neighbor_indices = sorted_indices[:, 1 : n_neighbors + 1]  # Shape: (m, n)

    # convert centroids into unit vectors in the camera frame
    unit_vectors = centroids_to_unit_vectors(centroids, fov, aspect_ratio, width, height)

    # compute angular distances for each neighbor entry
    angular_matrix = np.zeros((m, n_neighbors))
    
    for i in range(m):
        v_i = unit_vectors[i]
        for j_idx, neighbor_idx in enumerate(neighbor_indices[i]):
            v_j = unit_vectors[neighbor_idx]
            
            # dot product of unit vectors gives cos(phi)
            cos_phi = np.dot(v_i, v_j)
            
            # calculate angular separation
            phi = np.arccos(cos_phi)
            angular_matrix[i, j_idx] = phi

    return angular_matrix, neighbor_indices