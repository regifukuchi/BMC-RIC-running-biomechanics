"""Calculates the 3D angles between two cartesian coordinate system using Cardan 
angles sequence."""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"

import numpy as np

def LFTang(bf_L_T):
    '''
    Parameters
    ----------
    bf_L_T : Numpy array
        Versors of the foot segment from anatomical to global.

    Returns
    -------
    None.

    '''
    L_foot_ang = np.empty(bf_L_T.shape[1]) * np.NaN
    # % LEFT FOOT angles
    L_foot_ang[0] = np.arctan(bf_L_T[0,1]/np.sqrt(bf_L_T[1,1]**2 + bf_L_T[2,1]**2))
    # angle of the long axis of the foot about the vertical axis
    L_foot_ang[1] = np.arctan(-bf_L_T[0,0]/np.sqrt(bf_L_T[1,0]**2 + bf_L_T[2,0]**2))
    # and project the long axis into the sagital plane for ID of FOREFOOT
    L_foot_ang[2] = np.arctan2(bf_L_T[1,0], -bf_L_T[2,0])
    
    return L_foot_ang

def RFTang(bf_R_T):
    '''
    Parameters
    ----------
    bf_L_T : Numpy array
        Versors of the foot segment from anatomical to global.

    Returns
    -------
    None.

    '''
    R_foot_ang = np.empty(bf_R_T.shape[1]) * np.NaN
    # % RIGHT FOOT angles
    R_foot_ang[0] = np.arctan(bf_R_T[0,1]/np.sqrt(bf_R_T[1,1]**2 + bf_R_T[2,1]**2))
    # angle of the long axis of the foot about the vertical axis
    R_foot_ang[1] = np.arctan(bf_R_T[0,0]/np.sqrt(bf_R_T[1,0]**2 + bf_R_T[2,0]**2))
    # and project the long axis into the sagital plane for ID of FOREFOOT
    R_foot_ang[2] = np.arctan2(bf_R_T[1,0], -bf_R_T[2,0])
    
    return R_foot_ang
    
def PVang(bpT):
    pelvis_ang = np.empty(bpT.shape[1]) * np.NaN
    # % PELVIS angles
    # % project lateral axis of pelvis
    pelvis_ang[0] = np.arctan2(bpT[0,2],bpT[2,2]) - np.pi/2 #into floor plane
    pelvis_ang[1] = np.arctan(bpT[1,2]/bpT[0,2]) #into frontal plane
    #% and project anterior axis of pelvis
    pelvis_ang[2] = np.arctan2(bpT[1,0],-bpT[2,0])#into sagital plane
    
    return pelvis_ang

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