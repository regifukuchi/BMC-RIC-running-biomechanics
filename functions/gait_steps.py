# Prepare Oython environment
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json

from detecta import detect_peaks
#%% Import data parameters
pathname = r'../data'
fn_json=r"C:\Users\Reginaldo\OneDrive - University of Calgary\data\Figshare_SciData\new_unzip\201225\20140515T133244.json"

with open(fn_json, 'r') as f:
    data_RIC = json.load(f)
    
#%% Create dataframe column corresponding to the dataset
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
#%% # LEFT SIDE
# Determine functional measures and gait type (walk vs run) movement speed comes 
# from the A/P position time history of a heel marker so we first need to identify 
# a heel marker.

# % Combine 3 of the foot markers into one matrix (ignore the created fourth)
L_foot = df_S[['L_foot_1_X', 'L_foot_1_Y', 'L_foot_1_Z',
              'L_foot_2_X', 'L_foot_2_Y', 'L_foot_2_Z',
              'L_foot_3_X', 'L_foot_3_Y', 'L_foot_3_Z']].values.reshape((3,3))
# sort the markers from left to right
L_foot = L_foot[L_foot[:, 2].argsort()]

# find the lower of the two medial markers
if L_foot[1,1] < L_foot[2,1]:
    L_heel =  df_R[['L_foot_2_X', 'L_foot_2_Y', 'L_foot_2_Z']].values
else:
    L_heel =  df_R[['L_foot_1_X', 'L_foot_1_Y', 'L_foot_1_Z']].values

# Find peaks location. Signal flipped because of X-axis convention difference.
locs0 = detect_peaks(-np.diff(L_heel[:,0]), mpd=np.round(0.5*data_RIC['hz_r']), 
                    mph=0, show=False)
pks = np.diff(L_heel[:,0])[locs0]

locs = detect_peaks(np.diff(L_heel[:,0]), mpd=np.round(0.5*data_RIC['hz_r']), 
                    mph=0, show=False)

# Gait velocity and cadence
vel    = data_RIC['hz_r']*np.median(-pks)/1000; # gait speed
stRate = 60/(np.median(np.diff(locs))/data_RIC['hz_r']); # cadence
#print('Gait velocity is '+str(vel.round(1))+' m/s')
#print('Stride rate is '+str(stRate.round(1))+' strides/min')

#%% RIGHT SIDE
# % Combine 3 of the foot markers into one matrix (ignore the created fourth)
R_foot = df_S[['R_foot_1_X', 'R_foot_1_Y', 'R_foot_1_Z',
              'R_foot_2_X', 'R_foot_2_Y', 'R_foot_2_Z',
              'R_foot_3_X', 'R_foot_3_Y', 'R_foot_3_Z']].values.reshape((3,3))
# sort the markers from left to right
R_foot = R_foot[R_foot[:, 2].argsort()]

# find the lower of the two medial markers
if R_foot[1,1] < R_foot[2,1]:
    R_heel =  df_R[['R_foot_2_X', 'R_foot_2_Y', 'R_foot_2_Z']].values
else:
    R_heel =  df_R[['R_foot_1_X', 'R_foot_1_Y', 'R_foot_1_Z']].values
#%% Identify gait type using a trained LDA classifier.  This will be more
# robust for shuffle-runners, older adults and speed walkers.  gaitClass
# represents an LDA object which has been trained on 839 test sets of
# walking and running, and validated on ~2000 sets of walking and running.
testSet = np.array([vel, stRate])

import classreg.learning.classif.CompactClassificationDiscriminant
load('gaitClass.mat','gaitClass')


label = predict(gaitClass,testSet);
%label returned as cell
label = label{1};