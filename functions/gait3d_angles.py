"""Calculate the 3D joint angles of the hip, knee and ankle using Cardan rotation sequence computation.
It also calculates the segment angle of the pelvis and foot segments. 
"""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"


import numpy as np
from svdt import svdt
import SCS_RIC as scs
import jointangles3d as jang3d

def gait3d_angles(neutral, joints, gait):
    '''
    

    Parameters
    ----------
    neutral : TYPE
        DESCRIPTION.
    joints : TYPE
        DESCRIPTION.
    gait : TYPE
        DESCRIPTION.

    Returns
    -------
    virt_mkrs : TYPE
        DESCRIPTION.
    angles : TYPE
        DESCRIPTION.

    '''
    #%% TRANSFORMATION MATRIX USING SVD ALGORITHM
    # Right Foot segment
    Af_R = neutral.filter(like='R_foot', axis=1).values.mean(axis=0)
    Bf_R = gait.filter(like='R_foot', axis=1).values
    Rf_R, Lf_R, RMSEf_R = svdt(Af_R, Bf_R)
    
    # Left Foot segment
    Af_L = neutral.filter(like='L_foot', axis=1).values.mean(axis=0)
    Bf_L = gait.filter(like='L_foot', axis=1).values
    Rf_L, Lf_L, RMSEf_L = svdt(Af_L, Bf_L)
    
    # Right Shank segment
    As_R = neutral.filter(like='R_shank', axis=1).values.mean(axis=0)
    Bs_R = gait.filter(like='R_shank', axis=1).values
    Rs_R, Ls_R, RMSEs_R = svdt(As_R, Bs_R)
    
    # Left Shank segment
    As_L = neutral.filter(like='L_shank', axis=1).values.mean(axis=0)
    Bs_L = gait.filter(like='L_shank', axis=1).values
    Rs_L, Ls_L, RMSEs_L = svdt(As_L, Bs_L)
    
    # Right Thigh segment
    At_R = neutral.filter(like='R_thigh', axis=1).values.mean(axis=0)
    Bt_R = gait.filter(like='R_thigh', axis=1).values
    Rt_R, Lt_R, RMSEt_R = svdt(At_R, Bt_R)
    
    # Left Thigh segment
    At_L = neutral.filter(like='L_thigh', axis=1).values.mean(axis=0)
    Bt_L = gait.filter(like='L_thigh', axis=1).values
    Rt_L, Lt_L, RMSEt_L = svdt(At_L, Bt_L)
    
    # Pelvic segment
    Ap = neutral.filter(like='pelvis', axis=1).values.mean(axis=0)
    Bp = gait.filter(like='pelvis', axis=1).values
    Rp, Lp, RMSEp = svdt(Ap, Bp)
    
    #%% ANATOMICAL MARKERS OF THE SEGMENTS
    # Pelvis markers
    pelvis_1 = neutral[['pelvis_1_X','pelvis_1_Y','pelvis_1_Z']].values.flatten()
    pelvis_2 = neutral[['pelvis_2_X','pelvis_2_Y','pelvis_2_Z']].values.flatten()
    pelvis_3 = neutral[['pelvis_3_X','pelvis_3_Y','pelvis_3_Z']].values.flatten()
    pelvis_4 = neutral[['pelvis_4_X','pelvis_4_Y','pelvis_4_Z']].values.flatten()
    # the pelvis_jc is simply the average location of the pelvis markers
    jcPelvis = (pelvis_1+pelvis_2+pelvis_3+pelvis_4)/4
    # Pelvic anatomical markers
    pelvis_x_S = jcPelvis + np.array([0,0,-1])
    pelvis_y_S = jcPelvis + np.array([0,1,0])
    # RECONSTRUCT ANATOMICAL MARKERS USING TRANSFORMATION MATRIX
    pelvis_x_R = np.dot(Rp,pelvis_x_S) + Lp
    pelvis_y_R = np.dot(Rp,pelvis_y_S) + Lp
    jcPelvis_R = np.dot(Rp,jcPelvis) + Lp
    
    # Thigh markers
    L_GTR = joints[['L_hip_X','L_hip_Y','L_hip_Z']].values.flatten()
    R_GTR = joints[['R_hip_X','R_hip_Y','R_hip_Z']].values.flatten()
    L_lat_knee = joints[['L_lat_knee_X','L_lat_knee_Y','L_lat_knee_Z']].values.flatten()
    R_lat_knee = joints[['R_lat_knee_X','R_lat_knee_Y','R_lat_knee_Z']].values.flatten()
    L_med_knee = joints[['L_med_knee_X','L_med_knee_Y','L_med_knee_Z']].values.flatten()
    R_med_knee = joints[['R_med_knee_X','R_med_knee_Y','R_med_knee_Z']].values.flatten()
    # HJC
    jcL_hip = L_GTR + (R_GTR-L_GTR)/4
    jcR_hip = R_GTR + (L_GTR-R_GTR)/4
    # RECONSTRUCT ANATOMICAL MARKERS USING TRANSFORMATION MATRIX
    # Right Thigh
    RHJCr = np.dot(Rp,jcR_hip) + Lp
    RKNLr = np.dot(Rt_R,R_lat_knee) + Lt_R
    RKNMr = np.dot(Rt_R,R_med_knee) + Lt_R
    # Left Thigh
    LHJCr = np.dot(Rp,jcL_hip) + Lp
    LKNLr = np.dot(Rt_L,L_lat_knee) + Lt_L
    LKNMr = np.dot(Rt_L,L_med_knee) + Lt_L
    
    # FOOT
    #% FIRST need to identify the heel markers
    #% Combine 3 feet markers into one matrix ... can ignore the created fourth one
    R_foot_S = neutral[['R_foot_1_X','R_foot_1_Y','R_foot_1_Z',
                     'R_foot_2_X','R_foot_2_Y','R_foot_2_Z',
                     'R_foot_3_X','R_foot_3_Y','R_foot_3_Z']].values.reshape((3,3))
    # indices of the sequence of most left to the most right markers
    i_ft_R = list(R_foot_S[:, 0].argsort())
    # Determine what is the location of the 3 foot markers
    if R_foot_S[i_ft_R[0], 1] < R_foot_S[i_ft_R[1], 1]:
        RHEB_lbl = 'R_foot_'+str(i_ft_R[0]+1) # heel bottom
        RHET_lbl = 'R_foot_'+str(i_ft_R[1]+1) # heel top
    else:
        RHEB_lbl = 'R_foot_'+str(i_ft_R[1]+1) # heel bottom
        RHET_lbl = 'R_foot_'+str(i_ft_R[0]+1) # heel top
    #RHEL_lbl = 'R_foot_'+str(i_ft_R[2]+1)
        
    RHEBs = neutral.filter(like=RHEB_lbl, axis=1).values.flatten()
    RHETs = neutral.filter(like=RHET_lbl, axis=1).values.flatten()
    #RHELs = neutral.filter(like=RHEL_lbl, axis=1).values.flatten()
    RTOEs = RHEBs + np.array([0, 0, -1]) # long axis of the the foot is aligned with the lab AP
    
    # Left foot
    #% FIRST need to identify the heel markers
    #% Combine 3 feet markers into one matrix ... can ignore the created fourth one
    L_foot_S = neutral[['L_foot_1_X','L_foot_1_Y','L_foot_1_Z',
                     'L_foot_2_X','L_foot_2_Y','L_foot_2_Z',
                     'L_foot_3_X','L_foot_3_Y','L_foot_3_Z']].values.reshape((3,3))
    # indices of the sequence of most left to the most right markers
    i_ft_L = list(L_foot_S[:, 0].argsort())
    # Determine what is the location of the 3 foot markers
    if L_foot_S[i_ft_L[1], 1] < L_foot_S[i_ft_L[2], 1]:
        LHEB_lbl = 'L_foot_'+str(i_ft_L[1]+1) # heel bottom
        LHET_lbl = 'L_foot_'+str(i_ft_L[2]+1) # heel top
    else:
        LHEB_lbl = 'L_foot_'+str(i_ft_L[2]+1) # heel bottom
        LHET_lbl = 'L_foot_'+str(i_ft_L[1]+1) # heel top
    #LHEL_lbl = 'L_foot_'+str(i_ft_L[0]+1)
        
    LHEBs = neutral.filter(like=LHEB_lbl, axis=1).values.flatten()
    LHETs = neutral.filter(like=LHET_lbl, axis=1).values.flatten()
    #LHELs = neutral.filter(like=LHEL_lbl, axis=1).values.flatten()
    LTOEs = LHEBs + np.array([0, 0, -1]) # long axis of the the foot is aligned with the lab AP
    
    # RECONSTRUCT ANATOMICAL MARKERS USING TRANSFORMATION MATRIX
    # Right Foot
    RHEBr = np.dot(Rf_R,RHEBs)  + Lf_R
    RHETr = np.dot(Rf_R,RHETs)  + Lf_R
    RTOEr = np.dot(Rf_R,RTOEs)  + Lf_R
    # Left Foot
    LHEBr = np.dot(Rf_L,LHEBs)  + Lf_L
    LHETr = np.dot(Rf_L,LHETs)  + Lf_L
    LTOEr = np.dot(Rf_L,LTOEs)  + Lf_L
    
    # SHANK
    # Shank markers
    L_lat_ankle = joints[['L_lat_ankle_X','L_lat_ankle_Y','L_lat_ankle_Z']].values.flatten()
    R_lat_ankle = joints[['R_lat_ankle_X','R_lat_ankle_Y','R_lat_ankle_Z']].values.flatten()
    L_med_ankle = joints[['L_med_ankle_X','L_med_ankle_Y','L_med_ankle_Z']].values.flatten()
    R_med_ankle = joints[['R_med_ankle_X','R_med_ankle_Y','R_med_ankle_Z']].values.flatten()
    
    # RECONSTRUCT ANATOMICAL MARKERS USING TRANSFORMATION MATRIX
    LMALr = np.dot(Rs_L,L_lat_ankle) + Ls_L
    RMALr = np.dot(Rs_R,R_lat_ankle) + Ls_R
    LMAMr = np.dot(Rs_L,L_med_ankle) + Ls_L
    RMAMr = np.dot(Rs_R,R_med_ankle) + Ls_R
    
    # Create dictionary to output results
    virt_mkrs = {'pelvis_x_R':pelvis_x_R, 'pelvis_y_R':pelvis_y_R, 'jcPelvis_R':jcPelvis_R,
                'RHJCr':RHJCr, 'RKNLr':RKNLr, 'RKNMr':RKNMr, 'LHJCr':LHJCr, 'LKNLr':LKNLr,
                'LKNMr':LKNMr,'RHEBr':RHEBr, 'RHETr':RHETr, 'RTOEr':RTOEr, 'LHEBr':LHEBr, 
                 'LHETr':LHETr, 'LTOEr':LTOEr, 'LMALr':LMALr, 'RMALr':RMALr, 'LMAMr':LMAMr, 'RMAMr':RMAMr}
    
    #%% CALCULATE ANGLES
    # preallocate arrays
    L_ankle_ang, R_ankle_ang=np.empty(RHETr.shape) * np.NaN, np.empty(RHETr.shape) * np.NaN
    L_knee_ang, R_knee_ang=np.empty(RHETr.shape) * np.NaN, np.empty(RHETr.shape) * np.NaN
    L_hip_ang, R_hip_ang=np.empty(RHETr.shape) * np.NaN, np.empty(RHETr.shape) * np.NaN
    L_foot_ang, R_foot_ang=np.empty(RHETr.shape) * np.NaN, np.empty(RHETr.shape) * np.NaN
    pelvis_ang  = np.empty(RHETr.shape) * np.NaN
    for i in range(RHETr.shape[0]):
        # Calculate angles
        bf_L = scs.footCS(LHETr[i,:], LHEBr[i,:], LTOEr[i,:])
        bf_R = scs.footCS(RHETr[i,:], RHEBr[i,:], RTOEr[i,:])
        bs_L = scs.shankCS(LKNLr[i,:], LKNMr[i,:], LMALr[i,:], LMAMr[i,:])
        bs_R = scs.shankCS(RKNLr[i,:], RKNMr[i,:], RMALr[i,:], RMAMr[i,:])
        bt_L = scs.thighCS(LKNLr[i,:], LKNMr[i,:], LHJCr[i,:])
        bt_R = scs.thighCS(RKNLr[i,:], RKNMr[i,:], RHJCr[i,:])
        bp = scs.pelvisCS(pelvis_x_R[i,:],pelvis_y_R[i,:],jcPelvis_R[i,:])
        # Global to anatomical
        bf_L_T, bs_L_T, bt_L_T, bpT  = bf_L.T, bs_L.T, bt_L.T, bp.T
        bf_R_T, bs_R_T, bt_R_T  = bf_R.T, bs_R.T, bt_R.T 
    
        # Product between [lab to shank] and [foot to lab]
        L_ankle_ang[i,:] = jang3d.cardanangles(np.dot(bs_L, bf_L_T))
        R_ankle_ang[i,:] = jang3d.cardanangles(np.dot(bs_R, bf_R_T))
        L_knee_ang[i,:]  = jang3d.cardanangles(np.dot(bt_L, bs_L_T))
        R_knee_ang[i,:]  = jang3d.cardanangles(np.dot(bt_R, bs_R_T))
        L_hip_ang[i,:]   = jang3d.cardanangles(np.dot(bp, bt_L_T))
        R_hip_ang[i,:]   = jang3d.cardanangles(np.dot(bp, bt_R_T))
    
        # FOOT angles
        L_foot_ang[i,:] = jang3d.LFTang(bf_L_T)
        R_foot_ang[i,:] = jang3d.RFTang(bf_R_T)    
        # PELVIS angles
        pelvis_ang[i,:] = jang3d.PVang(bpT)
        
    # Dictionary with angles
    angles = {'pelvis_ang':pelvis_ang, 'L_foot_ang':L_foot_ang, 'R_foot_ang':R_foot_ang,
              'L_hip_ang':L_hip_ang, 'R_hip_ang':R_hip_ang, 'L_knee_ang':L_knee_ang,
              'R_knee_ang':R_knee_ang, 'L_ankle_ang':L_ankle_ang, 'R_ankle_ang':R_ankle_ang}
    
    return virt_mkrs, angles