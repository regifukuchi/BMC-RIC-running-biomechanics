"""Calculates the 3D angles between two cartesian coordinate system using Cardan 
angles sequence."""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"

import numpy as np

def cardanangles(r):
'''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %   This function inputs the ROTATION matrix for a given joint
    # %   at one point in time and uses the cardan angles sequence to
    # %   provide the XYZ joint angles in OUT argument
    # %
    # %   The function uses the following  rotation matrix to calculate angles
    # %      | CzCy-SzSySx  SzCy+CzSySx  -SyCx |
    # %      | -SzCx        CzCx         Sx    |
    # %      | CzSy+SzCySx  SzSy-CzCySx  CyCx  |
    # %  INPUTS
    # %  --------
    # %   R (mat):    A 3x3 rotation matrix for a joint

    # %  OUTPUTS
    # %  -------
    # %  OUT (mat):      The three (1x3) planes of rotation in radians
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%
    # %   angle.L_ankle(i,:) = cardanangles(R.L_ankle(:,:,i));

    # % the use of atan2 increases the stability by avoiding gimble lock in
    # % the physiologicaly posible range of joint angles.
'''
    x = np.arctan2(r[1,2], np.sqrt(r[0,2]**2+r[2,2]**2))
    y = np.arctan2(-r[0,2], r[2,2])
    z = np.arctan2(-r[1,0], r[1,1])
    
    out = np.array([x,y,z])
    
    return out