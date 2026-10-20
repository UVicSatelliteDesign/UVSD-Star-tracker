import numpy as np
from typing import List, Tuple

class image_centroid:
    def __init__(self, x: float, y: float, total_intensity: float, pixel_count: int, bounding_box: Tuple[int, int, int, int]):
        self.x = x                               # Image space X coordinate
        self.y = y                               # Image space Y coordinate
        self.total_intensity = total_intensity   # Sum of pixel brightness values in the blob
        self.pixel_count = pixel_count           # Number of pixels belonging to this star blob
        self.bounding_box = bounding_box


def flood_fill(start_y: int, start_x: int, binary_mask: np.ndarray, visited: np.ndarray) -> List[Tuple[int, int]]:
    """
    Performs an iterative stack-based flood fill (8-connectivity) starting from a seed pixel.
    Mirrors the logic of a low-level C implementation using a stack buffer.
    """
    height, width = binary_mask.shape
    component_pixels: List[Tuple[int, int]] = []
    
    # Explicit stack for tracking pixels to process (avoids recursion depth limits in Python)
    stack = [(start_y, start_x)]
    visited[start_y, start_x] = True
    
    # 8-connectivity neighbor offsets (horizontal, vertical, and diagonal)
    neighbors = [
        (-1, -1), (-1, 0), (-1, 1),
        ( 0, -1),          ( 0, 1),
        ( 1, -1), ( 1, 0), ( 1, 1)
    ]
    
    while stack:
        cy, cx = stack.pop()
        component_pixels.append((cy, cx))
        
        for dy, dx in neighbors:
            ny, nx = cy + dy, cx + dx
            
            # Check image boundaries
            if 0 <= ny < height and 0 <= nx < width:
                # If the neighbor is part of the star binary mask and has not been visited yet
                if binary_mask[ny, nx] and not visited[ny, nx]:
                    visited[ny, nx] = True
                    stack.append((ny, nx))
                    
    return component_pixels


def extract_centroids(image_bitmap: np.ndarray, noise_threshold: float) -> List[image_centroid]:
    img_float = image_bitmap.astype(np.float32)
    height, width = img_float.shape
    
    # 1. Threshold the image to create a boolean binary mask of potential star pixels
    binary_mask = (img_float > noise_threshold)
    visited = np.zeros((height, width), dtype=bool)
    
    detected_stars: List[image_centroid] = []
    
    # 2. Raster scan across the entire image grid
    for y in range(height):
        for x in range(width):
            # If we find an unvisited pixel belonging to a star, trigger flood fill
            if binary_mask[y, x] and not visited[y, x]:
                component_pixels = flood_fill(y, x, binary_mask, visited)
                
                pixel_count = len(component_pixels)
                
                # 3. Filter out speckles that are too small
                if pixel_count < 2:
                    continue

                # Compute centroid bounding box
                min_x = min(pixel[1] for pixel in component_pixels)
                max_x = max(pixel[1] for pixel in component_pixels)
                min_y = min(pixel[0] for pixel in component_pixels)
                max_y = max(pixel[0] for pixel in component_pixels)
                
                bbox = (min_x, min_y, (max_x - min_x + 1), (max_y - min_y + 1))

                # perform weighted sum of pixel brightnessses to find centroids
                total_intensity = 0.0
                weighted_x_sum = 0.0
                weighted_y_sum = 0.0
                
                for cy, cx in component_pixels:
                    intensity = img_float[cy, cx]
                    total_intensity += intensity
                    weighted_x_sum += cx * intensity
                    weighted_y_sum += cy * intensity
                
                centroid_x = weighted_x_sum / total_intensity
                centroid_y = weighted_y_sum / total_intensity
                
                detected_stars.append(image_centroid(
                    x=float(centroid_x),
                    y=float(centroid_y),
                    total_intensity=float(total_intensity),
                    pixel_count=pixel_count,
                    bounding_box = bbox
                ))
                
    return detected_stars