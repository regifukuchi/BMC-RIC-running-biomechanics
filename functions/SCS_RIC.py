"""Defines the segment coordinate system for the pelvis, thigh, shank and 
# foot segments consistent to the definition adopted by the 3DGaitSystem of the 
# Running Injury Clinic, University of Calgary, Canada.
"""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"

import numpy as np

def pelvisCS(pelvis_x,pelvis_y,jcPelvis):
    '''
    Parameters
    ----------
    pelvis_x : Numpy array 
        Vector representing the A-P axis of the pelvis.
    pelvis_y : Numpy array 
        Vector representing the Vertical axis of the pelvis.
    jcPelvis : Numpy array
        Centroid of the four markers of the pelvis.

    Returns
    -------
    bp :  Numpy array
        Base of the pelvic segment.

    '''
    #%% PELVIS reference system
    # since trochanters are hard to landmark, pelvis will just be orthogonal to the lab
    #anterior axis
    vp_x = (pelvis_x-jcPelvis)/np.linalg.norm(pelvis_x-jcPelvis)
    #long axis pointing up
    vp_y = (pelvis_y-jcPelvis)/np.linalg.norm(pelvis_y-jcPelvis)
    #hinge axis to the subject's right
    vp_z = np.cross(vp_x,vp_y)/np.linalg.norm(np.cross(vp_x,vp_y))

    # combine to create a transformation matrix from anatomical to global
    bp = np.array([vp_x,vp_y,vp_z])
    
    return bp


def thighCS(KNL, KNM, HJC):
    #%% THIGH reference system
    KJC = (KNL+KNM)/2
    # Segment axes
    vt_y=(HJC-KJC)/np.linalg.norm(HJC-KJC); #long axis pointing up
    
    #almost the hinge joing pointing to the right
    if KNL[0] > KNM[0]: # to consider right and left side
        vt_z_temp=(KNL-KNM)/np.linalg.norm(KNL-KNM) 
    else:
        vt_z_temp=(KNM-KNL)/np.linalg.norm(KNM-KNL) 

    #Anterior axis from the cross
    vt_x = np.cross(vt_y,vt_z_temp)/np.linalg.norm(np.cross(vt_y,vt_z_temp))
    #hinge axis, lateral for right
    vt_z = np.cross(vt_x,vt_y)/np.linalg.norm(np.cross(vt_x,vt_y))

    #% foot reference system
    bt = np.array([vt_x,vt_y,vt_z])
    
    return bt

def shankCS(RKNL, RKNM, RMAL, RMAM):
    #%% RIGHT SHANK

    KJC = (RKNL+RKNM)/2 # midpoint of the two knee markers
    AJC = (RMAL+RMAM)/2 # midpoint of the two ankle markers

    # Segment axes
    vs_y=(KJC-AJC)/np.linalg.norm(KJC-AJC); #long axis pointing up
    
    #almost the hinge joing pointing to the right
    if RMAL[0] > RMAM[0]: # to consider right and left side
        vs_z_temp=(RMAL-RMAM)/np.linalg.norm(RMAL-RMAM)
    else:
        vs_z_temp=(RMAM-RMAL)/np.linalg.norm(RMAM-RMAL)
        
    #Anterior axis from the cross
    vs_x = np.cross(vs_y,vs_z_temp)/np.linalg.norm(np.cross(vs_y,vs_z_temp))
    #hinge axis, lateral for right
    vs_z = np.cross(vs_x,vs_y)/np.linalg.norm(np.cross(vs_x,vs_y))

    #% foot reference system
    bs = np.array([vs_x,vs_y,vs_z])
    
    return bs

def footCS(HET, HEB, TOE):
    #%% Foot reference system
    # long axis of the the foot is aligned with the lab AP
    vf_x = (TOE-HEB)/np.linalg.norm(TOE-HEB)
    # SECOND, create a vector from the two left markers (not the lateral one)
    vf_y_temp = (HET-HEB)/np.linalg.norm(HET-HEB)
    # use the temp vertical axis to create the lateral axis
    vf_z = np.cross(vf_x,vf_y_temp)/np.linalg.norm(np.cross(vf_x,vf_y_temp))
    #% and create the 'vertical' axis that provides standing eversion angle
    vf_y = np.cross(vf_z,vf_x)/np.linalg.norm(np.cross(vf_z,vf_x))

    #% foot reference system
    bf = np.array([vf_x,vf_y,vf_z])
    
    return bf