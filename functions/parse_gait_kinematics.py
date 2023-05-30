"""Import files from public data sets, parse them and export as a format accepted by gait_kinematics.py function.

    Parse files from RIC dataset:
      parse_RIC(fname)
    Parse files from RBDS dataset:
      parse_RBDS(fname)

"""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"


import json
from ezc3d import c3d
import numpy as np
import pandas as pd

def parse_RBDS(fn_static, fn_dynamic):
    '''
    Parameters
    ----------
    fn_static : c3d file format
        Containing markers positions during standing calibration trial.
    fn_dynamic : c3d file format
        Containing markers positions during dynamic gait (walking or running) trial.

    Returns
    -------
    None.

    '''
    ##% STATIC TRIAL
    # Import file
    c = c3d(fn_static)
    point_data = c['data']['points']
    mkr_S_labels_RBDS = c['parameters']['POINT']['LABELS']['value']
    # Hardcode marker labels adopted in the RBDS dataset
    mkr_S_labels_RBDS2 = ['L.ASIS','R.ASIS','R.PSIS','L.PSIS','L.Ankle','L.Ankle.Medial',
                          'L.GTR','L.Heel.Bottom','L.Heel.Lateral','L.Heel.Top',
                          'L.Knee', 'L.Knee.Medial','L.Shank.Bottom.Lateral',
                          'L.Shank.Bottom.Medial', 'L.Shank.Top.Lateral',
                          'L.Shank.Top.Medial','L.Thigh.Bottom.Lateral', 
                          'L.Thigh.Bottom.Medial','L.Thigh.Top.Lateral',
                          'L.Thigh.Top.Medial','R.Ankle','R.Ankle.Medial','R.GTR',
                          'R.Heel.Bottom','R.Heel.Lateral','R.Heel.Top','R.Knee',
                          'R.Knee.Medial','R.Shank.Bottom.Lateral','R.Shank.Bottom.Medial',
                          'R.Shank.Top.Lateral','R.Shank.Top.Medial','R.Thigh.Bottom.Lateral',
                          'R.Thigh.Bottom.Medial','R.Thigh.Top.Lateral','R.Thigh.Top.Medial']
    
    # List of corresponding label names of the RIC dataset
    mkr_S_labels_RBDS3 = ['pelvis_1','pelvis_2','pelvis_3','pelvis_4',
                      'L_lat_ankle','L_med_ankle','L_hip','L_foot_2','L_foot_3',
                      'L_foot_1','L_lat_knee','L_med_knee','L_shank_1','L_shank_2',
                      'L_shank_3','L_shank_4','L_thigh_1','L_thigh_2','L_thigh_3',
                      'L_thigh_4','R_lat_ankle','R_med_ankle','R_hip','R_foot_2',
                      'R_foot_3','R_foot_1','R_lat_knee','R_med_knee','R_shank_1',
                      'R_shank_2', 'R_shank_3', 'R_shank_4','R_thigh_1',
                      'R_thigh_2','R_thigh_3','R_thigh_4']
    
    # Indices of mkr labels of the original list
    idx = [mkr_S_labels_RBDS.index(m) for m in mkr_S_labels_RBDS2]
    
    # Array containing markers data
    mkr_S_data_RBDS = np.empty(shape=(point_data.shape[2],3*len(mkr_S_labels_RBDS2)))
    for m, marker in enumerate(mkr_S_labels_RBDS2):
        mkr_S_data_RBDS[:,3*m:3*m+3] = point_data[:3, idx[m], :].T
        
    # List containing columns header of the pandas df
    mkr_S_labels_RBDS4 = list(np.repeat(mkr_S_labels_RBDS3,3))
    mkr_S_labels_RBDS5 = [marker+'_'+xyz for marker, xyz in zip(mkr_S_labels_RBDS4, 
                                                        list('XYZ')*len(mkr_S_labels_RBDS3))]
    
    # Create dataframe with static markers data
    df_S_RBDS = pd.DataFrame(data=mkr_S_data_RBDS, columns=mkr_S_labels_RBDS5)
    df_S_RBDS = df_S_RBDS.mean().to_frame().T
    
    ##% GAIT TRIAL
    c_g = c3d(fn_dynamic)
    point_data_g = c_g['data']['points']
    
    mkr_R_labels_RBDS = c_g['parameters']['POINT']['LABELS']['value']
    
    # Hardcode marker labels adopted in the RBDS dataset
    mkr_R_labels_RBDS2 = ['L.ASIS','R.ASIS','R.PSIS','L.PSIS',
                      'L.Thigh.Bottom.Lateral','L.Thigh.Bottom.Medial','L.Thigh.Top.Lateral','L.Thigh.Top.Medial',
                      'R.Thigh.Bottom.Lateral','R.Thigh.Bottom.Medial','R.Thigh.Top.Lateral','R.Thigh.Top.Medial',
                      'L.Shank.Bottom.Lateral', 'L.Shank.Bottom.Medial', 'L.Shank.Top.Lateral', 'L.Shank.Top.Medial',
                      'R.Shank.Bottom.Lateral','R.Shank.Bottom.Medial','R.Shank.Top.Lateral','R.Shank.Top.Medial',
                      'L.Heel.Bottom','L.Heel.Lateral','L.Heel.Top',
                      'R.Heel.Bottom', 'R.Heel.Lateral', 'R.Heel.Top']
    
    # List of corresponding label names of the RIC dataset
    mkr_R_labels_RBDS3 = ['pelvis_1', 'pelvis_2', 'pelvis_3', 'pelvis_4',
                     'L_thigh_1', 'L_thigh_2', 'L_thigh_3', 'L_thigh_4',
                     'R_thigh_1', 'R_thigh_2', 'R_thigh_3', 'R_thigh_4',
                     'L_shank_1', 'L_shank_2', 'L_shank_3', 'L_shank_4',
                     'R_shank_1', 'R_shank_2', 'R_shank_3', 'R_shank_4',
                     'L_foot_1', 'L_foot_2', 'L_foot_3',
                     'R_foot_1', 'R_foot_2', 'R_foot_3']
    
    # Indices of mkr labels of the original list
    idx = [mkr_R_labels_RBDS.index(m) for m in mkr_R_labels_RBDS2]
    
    # Array containing markers data
    mkr_R_data_RBDS = np.empty(shape=(point_data_g.shape[2],3*len(mkr_R_labels_RBDS2)))
    for m, marker in enumerate(mkr_R_labels_RBDS2):
        mkr_R_data_RBDS[:,3*m:3*m+3] = point_data_g[:3, idx[m], :].T
        
    # List containing columns header of the pandas df
    mkr_R_labels_RBDS4 = list(np.repeat(mkr_R_labels_RBDS3,3))
    mkr_R_labels_RBDS5 = [marker+'_'+xyz for marker, xyz in zip(mkr_R_labels_RBDS4, 
                                                        list('XYZ')*len(mkr_R_labels_RBDS3))]
    
    df_R_RBDS = pd.DataFrame(data=mkr_R_data_RBDS, columns=mkr_R_labels_RBDS5)
    
    ##% APPLY TRANSFORMATION MATRIX TO BE CONSISTENT WITH RIC REFERENCE FRAME CONVENTION
    # Static
    df_Sm = np.empty(df_S_RBDS.shape[1]) * np.NaN
    rot = np.array([[0,0,1],[0,1,0],[-1,0,0]]) # rotate markers 90 deg
    for m in range(len(mkr_S_labels_RBDS3)):
        df_Sm[3*m:3*m+3] = np.dot(rot, df_S_RBDS.values[0,3*m:3*m+3])
    neutral = pd.DataFrame(data=df_Sm[np.newaxis,:], 
                               columns=df_S_RBDS.columns.tolist())
    
    # Add foot 4th marker as a centroid of the other 3 markers
    neutral[['L_foot_4_X','L_foot_4_Y','L_foot_4_Z']] = neutral.filter(like='L_foot', 
                                                                       axis=1).values.reshape((3,
                                                                                               3)).sum(axis=0)/3
    neutral[['R_foot_4_X','R_foot_4_Y','R_foot_4_Z']] = neutral.filter(like='R_foot', 
                                                                       axis=1).values.reshape((3,
                                                                                               3)).sum(axis=0)/3
   # Create joint df
    knee_jc_S  = neutral.filter(like='_knee_', axis=1)
    ankle_jc_S = neutral.filter(like='_ankle_', axis=1)
    hip_jc_S   = neutral.filter(like='_hip_', axis=1)
    
    joints = pd.concat([hip_jc_S, knee_jc_S, ankle_jc_S], axis=1)
    
    # Dynamic trial
    df_Gm = np.empty(shape=df_R_RBDS.shape) * np.NaN
    rot = np.array([[0,0,1],[0,1,0],[-1,0,0]]) # rotate markers 90 deg
    for m in range(len(mkr_R_labels_RBDS3)):
        for i in range(df_R_RBDS.shape[0]):
            df_Gm[i, 3*m:3*m+3] = np.dot(rot, df_R_RBDS.values[i,3*m:3*m+3])
    # Create dataframe with dynamic markers data
    gait = pd.DataFrame(data=df_Gm, 
                               columns=df_R_RBDS.columns.tolist())
    
    # Create foot 4th marker as a centroid of the other 3 markers
    L_foot_4, R_foot_4 = np.empty(shape=(gait.shape[0],3)) * np.NaN, np.empty(shape=(gait.shape[0],3)) * np.NaN
    for i in range(gait.shape[0]):# check better way to do this since it is taking so long
        L_foot_4[i,:] = gait.filter(like='L_foot', axis=1).values[i,:].reshape((3,3)).sum(axis=0)/3
        R_foot_4[i,:] = gait.filter(like='R_foot', axis=1).values[i,:].reshape((3,3)).sum(axis=0)/3
        
    gait[['L_foot_4_X','L_foot_4_Y','L_foot_4_Z']] = L_foot_4
    gait[['R_foot_4_X','R_foot_4_Y','R_foot_4_Z']] = R_foot_4

    # Sampling frequency
    hz = int(c['parameters']['POINT']['RATE']['value'])

    return neutral, joints, gait, hz                      
    
def parse_RIC(filename):
    
    # import JSON file
    with open(filename, 'r') as f:
        data_RIC = json.load(f)
    
    hz = data_RIC['hz_r']
    # Create dataframe column corresponding to the dataset
    neutral_lbls = list(data_RIC['neutral'].keys())
    xyz = list('XYZ')*len(neutral_lbls)
    neutral_lbls = [ele for ele in neutral_lbls for i in range(3)]
    neutral_lbls = [neutral_lbls[i]+'_'+xyz[i] for i in range(len(xyz))]
    
    # Joint marker labels static trial
    joints_lbls = list(data_RIC['joints'].keys())
    xyz = list('XYZ')*len(joints_lbls)
    joints_lbls = [ele for ele in joints_lbls for i in range(3)]
    joints_lbls = [joints_lbls[i]+'_'+xyz[i] for i in range(len(xyz))]
    
    # Marker labels running trial
    gait_lbls = list(data_RIC['running'].keys())
    xyz = list('XYZ')*len(gait_lbls)
    gait_lbls = [ele for ele in gait_lbls for i in range(3)]
    gait_lbls = [gait_lbls[i]+'_'+xyz[i] for i in range(len(xyz))]
    
    # Convert dictionaries into pandas dfs
    neutral = pd.DataFrame.from_dict(data_RIC['neutral']).values.reshape((1,
                                                                          len(neutral_lbls)), 
                                                                         order='F')
    joints  = pd.DataFrame.from_dict(data_RIC['joints']).values.reshape((1,
                                                                         len(joints_lbls)), 
                                                                        order='F')
    
    # Create new pandas dfs with neutral and joint data
    neutral = pd.DataFrame(data=neutral, columns=neutral_lbls)
    joints = pd.DataFrame(data=joints, columns=joints_lbls)
    
    # Array with dynamic markers data
    run_data = np.empty(shape=(5000, len(list(data_RIC['running'].keys()))*3))
    for m, mkr in enumerate(list(data_RIC['running'].keys())):
        run_data[:, 3*m:3*(m+1)] = np.array(data_RIC['running'][mkr])
    
    # Create dataframe with running data
    gait = pd.DataFrame(data = run_data, columns=gait_lbls)
    
    return neutral, joints, gait, hz