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
import numpy as np
import pandas as pd

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