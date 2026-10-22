import numpy as np
from typing import Callable, Optional, Union

def camera_rotation_from_ra_dec(ra: float, dec: float) -> np.ndarray:
    # convert ra and dec to unit vector components
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)
    target_dir = np.array([x, y, z], dtype=float)
    
    reference_up = np.array([0.0, 0.0, 1.0])
    
    # edge care if pointing directly at a celestial pole
    polar_proximity = np.dot(target_dir, reference_up)
    if np.abs(polar_proximity) == 1.0:
        reference_up = np.array([0.0, np.sign(polar_proximity), 0.0])
        
    # compute camera right (+X) and up (+Y) vectors
    right = np.cross(target_dir, reference_up)
    right /= np.linalg.norm(right)
    
    up = np.cross(right, target_dir)
    up /= np.linalg.norm(up)
    
    # concatenate direction vectors into a rotation matrix
    R = np.column_stack((right, up, target_dir))
    
    return R

def create_tumbling_path(omega_vec:  np.array, initial_rotation:  np.array) -> Callable[[float], np.ndarray]:
    omega_magnitude = np.linalg.norm(omega_vec)
    
    R0 = np.array(initial_rotation, dtype=float)
        
    # if angular velocity is zero, return a static orientation path
    if omega_magnitude == 0.0:
        return lambda t: R0
        
    # unit vector of axis of rotation
    axis = omega_vec / omega_magnitude
    k_x, k_y, k_z = axis
    
    # pre-compute skew-symmetric cross-product matrix K for Rodrigues' formula
    K = np.array([
        [0.0, -k_z, k_y],
        [k_z, 0.0, -k_x],
        [-k_y, k_x, 0.0]
    ])
    
    def path(t: float) -> np.ndarray:
        # total rotation angle at time t
        angle = omega_magnitude * t
        c = np.cos(angle)
        s = np.sin(angle)
        
        # apply rodrigues' formula to incremental the tumble
        R_tumble = np.eye(3) + s * K + (1.0 - c) * np.dot(K, K)
        
        # total orientation is tumble applied relative to the initial orientation
        return R_tumble @ R0
        
    return path