import numpy as np
from typing import Callable, Optional, Union

def camera_rotation_from_ra_dec(ra: float, dec: float, roll: float) -> np.ndarray:
    # performs a zyz rotation sequence to transform vectors from equatorial coordinates to camera coordinates

    # rotation about z-axis by ra
    Rz_ra = np.array([
        [np.cos(ra), np.sin(ra), 0],
        [-np.sin(ra), np.cos(ra), 0],
        [0, 0, 1]
    ])

    # rotation about y-axis by dec
    Ry_dec = np.array([
        [np.cos(np.pi/2 - dec), 0, -np.sin(np.pi/2 - dec)],
        [0, 1, 0],
        [np.sin(np.pi/2 - dec), 0, np.cos(np.pi/2 - dec)]
    ])

    # rotation about z-axis by roll
    Rz_roll = np.array([
        [np.cos(np.pi/2 + roll), np.sin(np.pi/2 + roll), 0],
        [-np.sin(np.pi/2 + roll), np.cos(np.pi/2 + roll), 0],
        [0, 0, 1]
    ])

    # combine all rotations
    R = Rz_roll @ (Ry_dec @ Rz_ra)

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