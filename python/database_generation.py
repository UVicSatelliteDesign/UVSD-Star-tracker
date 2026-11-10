from catalog_parse import parse_tycho_2
from scipy.spatial import KDTree
import numpy as np
from pathlib import Path
from feature_extraction import save_features_to_csv


def generate_database(folder_path: str, magnitude_limit: float, n_neighbours: int):
    catalog_filename = "tycho_2.dat"
    print(f"Loading star catalog from {catalog_filename}...")
    ids, mags, vectors = parse_tycho_2(catalog_filename, magnitude_limit)

    # this will be way more complicated to implement in c / c++ for the real deal, but for now we can just use this library function:
    tree = KDTree(vectors)

    # query the n nearest neighbours:
    # this returns the distances and indices to the n nearest neighbors. This technically computes neighbour distances as coords through the unit sphere rather than the angles between them, but it still preserves the order so it works. (using n_neighbors + 1 because the first row will all be worthless zeros)
    dists, neighbour_indices = tree.query(vectors, k=n_neighbours + 1)
    
    # crop out the first row:
    dists = dists[:, 1:];
    neighbour_indices = neighbour_indices[:, 1:];

    
    # compute features by converting the coord lengths to the angles between the stars
    angle_cosines = 1.0 - 0.5 * dists ** 2
    angular_matrix = np.arccos(angle_cosines)

    # save to file
    folder = Path(folder_path)
    folder.mkdir(parents=True, exist_ok=True)

    file_path = folder / "features.csv"
    save_features_to_csv(file_path, star_ids=ids, neighbor_indices=neighbour_indices, angular_matrix=angular_matrix)

if __name__ == "__main__":
    generate_database("./database", 6.5, 50);
    