import cv2
import numpy as np
from pathlib import Path
from centroids import extract_centroids

def process_simulation_folder(folder_path: str, noise_threshold: float = 50.0):
    folder = Path(folder_path)
    img_path = folder / "simulated_star_field.png"
    
    if not img_path.exists():
        print(f"Error: Could not find 'simulated_star_field.png' in {folder_path}")
        return
        
    print(f"Loading image from {img_path}...")

    # using IMREAD_UNCHANGED here to preserve the full 16-bit dynamic range
    image = cv2.imread(str(img_path), cv2.IMREAD_UNCHANGED)
    
    if image is None:
        print("Error: Failed to load image data.")
        return
        
    print("Extracting centroids...")
    centroids = extract_centroids(image, noise_threshold=noise_threshold)
    print(f"Successfully extracted {len(centroids)} star centroids.")
    
    # prepare image for annotation (convert 16-bit or grayscale to an 8-bit BGR image for OpenCV drawing)
    if image.dtype == np.uint16:
        img_8bit = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    else:
        img_8bit = image.copy()
        
    if len(img_8bit.shape) == 2:
        annotated_img = cv2.cvtColor(img_8bit, cv2.COLOR_GRAY2BGR)
    else:
        annotated_img = img_8bit.copy()
        
    # draw bounding boxes and sub-pixel center crosshairs for each detected star
    for i, c in enumerate(centroids):
        left, top, width, height = c.bounding_box
        
        # draw green bounding box around the detected star blob
        cv2.rectangle(annotated_img, (left, top), (left + width, top + height), (0, 255, 0), 1)
        
        # draw a red dot at the exact sub-pixel centroid coordinate
        cx_int = int(round(c.x))
        cy_int = int(round(c.y))
        cv2.circle(annotated_img, (cx_int, cy_int), 1, (0, 0, 255), -1)
        
    # save the annotated output back into the simulation folder
    output_path = folder / "annotated_star_field.png"
    cv2.imwrite(str(output_path), annotated_img)
    print(f"Annotated image successfully saved to: {output_path}")

if __name__ == "__main__":
    path = "./synthetic_output/fov_30.0deg_exp_2s_ap_0.05m_foc_0.1m_res_1024x768"
    process_simulation_folder(str(path), noise_threshold=6250.0)