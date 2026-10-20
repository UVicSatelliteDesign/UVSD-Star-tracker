import re
import numpy as np
from typing import Tuple, List

def parse_Hipparcos(file_path):
    # Regular expression to match the catalogue line format:
    pattern = re.compile(
        r"^(\d+):\s*"                               # Star ID
        r"R\.A:\s*(\d+)h\s*(\d+)m\s*(\d+)s\s*"       # RA: hours, minutes, seconds
        r"DEC:\s*([+-]?\d+)°\s*(\d+)m\s*(\d+)s\s*"   # DEC: degrees (with sign), minutes, seconds
        r"brightness:\s*([\d\.]+)$"                  # Brightness / magnitude
    )
    
    ## outputs
    star_ids = []
    magnitudes = []
    unit_vectors = []
    
    with open(file_path, 'r', encoding='cp1252') as f:
        for line in f:
            ## get entry values from the current line
            match = pattern.match(line.strip())
            if not match:
                continue # Skip header lines or malformed rows
                
            star_id, ra_h, ra_m, ra_s, dec_deg, dec_m, dec_s, brightness = match.groups()
            
            # Convert RA from hours/min/sec to radians
            h, m, s = float(ra_h), float(ra_m), float(ra_s)
            alpha = 2 * np.pi * (1.0 / 24.0) * (h + (1.0 / 60.0) * (m + (1.0 / 60.0) * s))
            
            # Convert DEC from deg/min/sec to radians
            d, dm, ds = float(dec_deg), float(dec_m), float(dec_s)
            
            # Ensure minute/second offsets carry the sign of the degree component
            sign = -1.0 if d < 0 else 1.0
            delta = 2 * np.pi * (1.0 / 360.0) * (d + sign * (1.0 / 60.0) * (dm + (1.0 / 60.0) * ds))
            
            # Compute 3D unit vector
            ux = np.cos(alpha) * np.cos(delta)
            uy = np.sin(alpha) * np.cos(delta)
            uz = np.sin(delta)
            
            # Push to the output arrays
            star_ids.append(int(star_id))
            magnitudes.append(float(brightness))
            unit_vectors.append([ux, uy, uz])
            
    return (
        np.array(star_ids, dtype=int),
        np.array(magnitudes, dtype=float),
        np.array(unit_vectors, dtype=float)
    )

def parse_tycho_2(filepath: str, magnitude_limit: float = 11.0) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """
    Parses the Tycho-2 main catalogue (catalog.dat) and returns star identifiers,
    magnitudes, and 3D unit vectors in equatorial coordinates.
    
    Parameters:
    - filepath: Path to the tyc2 catalog.dat file.
    - magnitude_limit: Maximum magnitude threshold to keep file size/density manageable.
    
    Returns:
    - ids: List of string star identifiers (e.g., "TYC 1-13-1").
    - mags: 1D numpy array of visual magnitudes.
    - vectors: Nx3 numpy array of 3D unit vectors pointing to each star.
    """
    ids: List[str] = []
    mags_list: List[float] = []
    vectors_list: List[np.ndarray] = []
    
    print(f"Parsing Tycho-2 catalog from {filepath} (Magnitude limit <= {magnitude_limit})...")
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if len(line) < 130:
                continue
                
            try:
                # Extract Tycho-2 ID components (TYC1, TYC2, TYC3)
                tyc1 = line[0:4].strip()
                tyc2 = line[5:10].strip()
                tyc3 = line[11:12].strip()
                star_id = f"TYC {tyc1}-{tyc2}-{tyc3}"
                
                # Extract RA and Dec in decimal degrees (J2000 epoch)
                ra_str = line[15:27].strip()
                dec_str = line[28:40].strip()
                
                if not ra_str or not dec_str:
                    continue
                    
                ra = float(ra_str)
                dec = float(dec_str)
                
                # Extract VT magnitude (fallback to BT magnitude if VT is missing/blank)
                vt_str = line[123:128].strip()
                bt_str = line[84:90].strip()
                
                if vt_str:
                    mag = float(vt_str)
                elif bt_str:
                    mag = float(bt_str)
                else:
                    continue  # Skip if no magnitude data is present
                    
                # Filter by magnitude limit to optimize simulation performance
                if mag > magnitude_limit:
                    continue
                    
                # Convert RA and Dec to radians
                ra_rad = np.radians(ra)
                dec_rad = np.radians(dec)
                
                # Compute 3D unit vector in equatorial coordinate frame
                x = np.cos(dec_rad) * np.cos(ra_rad)
                y = np.cos(dec_rad) * np.sin(ra_rad)
                z = np.sin(dec_rad)
                vector = np.array([x, y, z], dtype=float)
                
                ids.append(star_id)
                mags_list.append(mag)
                vectors_list.append(vector)
                
            except ValueError:
                # Skip lines with malformed formatting or conversion errors
                continue
                
    mags = np.array(mags_list, dtype=float)
    vectors = np.array(vectors_list, dtype=float)
    
    print(f"Successfully loaded {len(ids)} stars from Tycho-2.")
    return ids, mags, vectors

if __name__ == "__main__":
    # Test parsing with a file named 'catalog.txt' in the same folder
    filename = "catalog.txt"
    print(f"Loading star catalog from {filename}...")
    
    ids, mags, vectors = parse_Hipparcos(filename)
    
    print(f"Successfully loaded {len(ids)} stars.")
    if len(ids) > 0:
        print(f"\nExample Parsed Star (ID: {ids[0]}):")
        print(f"  - Magnitude: {mags[0]}")
        print(f"  - Unit Vector: {vectors[0]}")

        