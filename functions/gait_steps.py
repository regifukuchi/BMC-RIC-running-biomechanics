"""Calculates the gait events based on PCA or FF/FB algorithm.
"""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"

# Prepare environment
import os, sys
import numpy as np
import pandas as pd
#from scipy import signal
import scipy.io as spio
import matlab
import matlab.engine
from scipy.signal import butter, filtfilt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from detecta import detect_peaks
from pca_td import pca_td
from pca_to import pca_to


def gait_steps(neutral, gait, angles, hz):
    '''
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Loads the NEURTAL and DYNAMIC data structures, divides the ANGLES and 
    %   VELOCITIES data structures into the different left and right steps and
    %   provides an index of touchdown and toeoffs (EVENTS/EVENT). The function 
    %   also outputs time normalized angles (NORM_ANG) and velocities 
    %   (NORM_VEL), gait speed (SPEEDOUTPUT), the method of  touchdown/toe off 
    %   detection (EVENTSFLAG), variables of interest (DISCRETE_VARIABLES)
    %   and the determined gait type (LABEL).
    %
    %  INPUTS
    %  --------
    %  NEUTRAL (dict):    Marker shell positions collected as part of 
    %                       the static trial.
    %
    %  GAIT (dict):     Marker shell positions collected as part of
    %                       the dynamic (run/walk) trial.
    %
    %  ANGLES (dict):     Angles (joint angles) structure created as an 
    %                       output from the function: gait_kinematics.
    %
    %  HZ (int):        Data collection sampling frequency.
    %
    %
    %  OUTPUTS
    %  -------
    %
    %  EVENTS (dict):    Matrix of frame numbers for touchdown and toeoffs.
    %
    %  EVENT (dict):     Same matric as EVENTS but also includes midswing in
    %                   order to calculate swing variables
    %
    %  SPEEDOUTPUT (float): Calculated speed of lowest heel marker. 
    %
    %  EVENTSFLAG (mat):Matrix the same size as EVENTS which indicates whether
    %                   PCA event detection was used (1) or if the default FF 
    %                   and FB events were used (0).
    %
    %  LABEL (str):     Returns a string which indicates whether the trial was 
    %                   a 'walk' or 'run' based on the classifier in this 
    %                   function.
    %
    %  LICENSE
    %  -------
    %  See file LICENSE.txt
    %
    % Copyright (C) 2010-2023,  Blayne Hettinga, Sean Osis, Allan Brett and
    %                           The Running Injury Clinic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %for troubleshooting
    
    % neutral = out.neutral;
    % % dynamic = out.walking;
    % dynamic = out.running;
    % angles  = r_angles;
    % velocities = r_velocities;
    % hz = 200;
    % plots = 0;
    '''
   
    #%% Determine functional measures and gait type (walk vs run)
    # % movement speed comes from the A/P position time history of a heel marker
    # % so we first need to identify a heel marker
    # LEFT SIDE
    # Determine functional measures and gait type (walk vs run) movement speed comes 
    # from the A/P position time history of a heel marker so we first need to identify 
    # a heel marker.
    
    # % Combine 3 of the foot markers into one matrix (ignore the created fourth)
    L_foot = neutral[['L_foot_1_X', 'L_foot_1_Y', 'L_foot_1_Z',
                  'L_foot_2_X', 'L_foot_2_Y', 'L_foot_2_Z',
                  'L_foot_3_X', 'L_foot_3_Y', 'L_foot_3_Z']].values.reshape((3,3))
    # sort the markers from left to right
    i_lf   = list(L_foot[:, 0].argsort())
    L_foot = L_foot[L_foot[:, 0].argsort()]
    
    # find the lower of the two medial markers
    if L_foot[1,1] < L_foot[2,1]:
        L_marker = 'L_foot_' + str(i_lf[1]+1)
        L_heel =  gait.filter(like=L_marker).values
    else:
        L_marker = 'L_foot_' + str(i_lf[2]+1)
        L_heel =  gait.filter(like=L_marker).values
        
    # Find peaks location. Signal flipped because of X-axis convention difference.
    locs0 = detect_peaks(np.diff(L_heel[:,2]), mpd=np.round(0.5*hz), 
                        mph=0, show=False)
    pks = np.diff(L_heel[:,2])[locs0]
    
    locs = detect_peaks(-np.diff(L_heel[:,2]), mpd=np.round(0.5*hz), 
                        mph=0, show=False)
    
    # Gait velocity and cadence
    vel    = hz*np.median(pks)/1000; # gait speed
    stRate = 60/(np.median(np.diff(locs))/hz); # cadence
    #print('Gait velocity is '+str(vel)+' m/s')
    #print('Stride rate is '+str(stRate)+' strides/min')
    
    #%% RIGHT SIDE
    # % Combine 3 of the foot markers into one matrix (ignore the created fourth)
    R_foot = neutral[['R_foot_1_X', 'R_foot_1_Y', 'R_foot_1_Z',
                  'R_foot_2_X', 'R_foot_2_Y', 'R_foot_2_Z',
                  'R_foot_3_X', 'R_foot_3_Y', 'R_foot_3_Z']].values.reshape((3,3))
    # sort the markers from left to right
    i_rf   = list(R_foot[:, 0].argsort())
    R_foot = R_foot[R_foot[:, 0].argsort()]
    
    # find the lower of the two medial markers
    if R_foot[0,1] < R_foot[1,1]:
        R_marker = 'R_foot_' + str(i_rf[0]+1)
        R_heel =  gait.filter(like=R_marker).values
    else:
        R_marker = 'R_foot_' + str(i_rf[1]+1)
        R_heel =  gait.filter(like=R_marker).values
        
    # Linear discriminant analysis
    # Import training dataset
    gaitClass = pd.read_csv(os.path.join(r'../data', 'LDA_out.txt'), delimiter='\t', 
                            header=None, names=['Category','Speed','Cadence'], usecols=[0,1,2])
    # Replace numerical by categorical
    gaitClass['Category'] = gaitClass['Category'].replace(1, 'walk')
    gaitClass['Category'] = gaitClass['Category'].replace(2, 'run')
    
    # Input to the model
    X = gaitClass[['Speed','Cadence']].values # training data
    y = gaitClass['Category'].values.tolist() # testing data
    model = LinearDiscriminantAnalysis()# define model
    model.fit(X, y) # Model fit
    
    # make a prediction
    label = model.predict(np.array([vel,stRate]).reshape((1,2)))[0]
    #print('Gait category is '+label)
    
    # Load PCA output from mat file
    # % event_data is a .mat file containing 'coeff' which is the coefficients
    # % from the pre-trained PCA and 'p' which is the list of coefficients of the
    # % linear polynomial relating PCA scores with touchdown timing relative to
    # % the foot acceleration peak.
    event_data_TD = spio.loadmat(os.path.join(r'../data', 'event_data_TD.mat'))
    event_data_TO = spio.loadmat(os.path.join(r'../data', 'event_data_TO.mat'))
    
    #%% Identify Touch Down and Take Off events: Gait Independent
    # % Use PCA touchdown detection based on updated Osis et al. (2014) for
    # % both walking and running.
    # % Use new PCA toeoff detection for both walking and running.
    # % evt variables are NOT rounded
    try:
        evtLtd, evtRtd = pca_td(angles, hz, event_data_TD, label)
        evtLto, evtRto = pca_to(angles, hz, event_data_TO, label)
        
    except Exception as e: 
        #For a small number of people, these functions return errors, or in the
        #case of bad data... default to use FF and FB in these cases
        
        evtRtd = []
        evtRto = []
        
        print('Automated event detection failed, defaulting to foot-forward foot-back')
        print(e)
        
    ## LEFT FOOT EVENTS
    # % when the feet are not tracked very well, discontinuities in the heel
    # % marker can occur causing the findpeaks to pick up additional 'peaks'
    # % for the purposes of simply identifying foot forward and foot back
    # % timing, we can over filter this signal. We do not care about the
    # % magnitude of the signal but only the timing so we can overfit as long
    # % as the filter has a zero phase shift.
    # % Note: signal is now filtered by default.  There is no advantage to not
    # % filtering, as if the signal quality is already good, then the system uses
    # % PCA event detection anyhow, and if the signal is bad, then it has to be
    # % filtered in order to get foot-forward foot-backward events.
    
    # Correct the cutoff frequency for the number of passes in the filter
    b, a = butter(2, 5/(hz/2), btype = 'low')
    # note that Python and Matlab filtfilt behaves slightly different with padding the data
    # see https://mail.python.org/pipermail/scipy-user/2014-April/035646.html
    filtered_L_heel = filtfilt(b, a, L_heel[:,2], padtype='odd')
    
    # Begin by creating a gross estimation of foot forwards and foot backs
    L_FFi = detect_peaks(-filtered_L_heel, mpd=np.round(0.35*hz),
                         show=False)
    
    if label == 'walk':
        # % Use peak foot flexion angle for foot back
        # % To deal with peaks resulting from signal flipping, threshold them
        angSig = angles['L_foot_Z'].values
        angSig[np.abs(angSig) > 90] = np.NaN
        L_FBi = detect_peaks(-angSig, mpd=np.round(0.7*hz),
                         mph=20, show=False)
    else:
        # Use rearmost position of heel marker for foot back
        L_FBi = detect_peaks(filtered_L_heel, mpd=np.round(0.35*hz),
                         show=False)
        
    # %Uncomment block below to enable more aggressive quality control of data
    
    # if (np.nanpercentile(np.abs(angles['foot_Z'].values), 90) > 120) & vel < 4
    #     print('Right ankle values outside of expected ranges, please ensure your shoe markers are properly placed and redo your collection')
    #     sys.exit()
    
    # Remove any leading FB
    L_FBi = L_FBi[L_FBi>L_FFi[0]]
    
    ## find largest chunk of continuous data
    
    #We want to check before and after that there is sufficient data for analysis
    
    if (L_FFi.shape[0] < 2) or (L_FBi.shape[0] < 2):
        print('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
        sys.exit()
        
    # Call Matlab function LARGEST_BLOCK.m from Python
    eng = matlab.engine.start_matlab() # start Matlab engine
    eng.cd(r'../functions', nargout=0) # set path for functions dir
    
    
    L_FFi, L_FBi, L_block_start, block_end = eng.largest_block(matlab.double(list(L_FFi)), 
                                                                 matlab.double(list(L_FBi)), nargout=4)
    
    L_FFi = np.array(L_FFi).flatten().astype(int)
    L_FBi = np.array(L_FBi).flatten().astype(int)
    
    if (L_FFi.shape[0] < 2) or (L_FBi.shape[0] < 2):
        print('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
        sys.exit()
        
    # TOUCHDOWN
    # evtLtd from above
    
    # SELECT SEQUENTIAL STEPS
    # create an ordered set of sequential steps using FFi as guide
    closest = np.abs(np.repeat(L_FFi[:,np.newaxis], evtLtd.shape[0], axis=1)-np.repeat(evtLtd[:,np.newaxis].T, L_FFi.shape[0], axis=0))
    
    mindist, minx = np.nanmin(closest, axis=0), np.nanargmin(closest, axis=0)
    
    for i in np.unique(minx).astype(int):
        if np.sum(np.isin(i,minx)) > 1:
            mindist = mindist[minx!=i]
            evtLtd  = evtLtd[minx!=i]
            minx    = minx[minx!=i]
            
    # Parameter based on the typical frame adjustments observed in 300
    # datasets
    if label=='run':
        maxadj = 0.05*hz
    else:
        maxadj = 0.10*hz
        
    # Preallocate
    L_TD = np.empty(L_FFi.shape[0]) * np.NaN
    evFltd = np.zeros(L_FFi.shape[0])
    
    # Here we replace FF indices with indices from evt where criteria are
    # met... the default is to use FF
    for i in range(L_FFi.shape[0]):
        try:
            if i > np.max(minx):
                break
            elif np.isin(i,minx) and (mindist[minx==i] < maxadj).any():
                #Replace with evtLtd since its more accurate
                L_TD[i] = evtLtd[minx==i][0]
                evFltd[i] = 1
            else:
                #Use FFi since it is more robust
                L_TD[i] = L_FFi[i]
                
        except Exception as e:
            print(e)
            L_TD[i] = L_FFi[i]
            
    # %% TAKEOFF
    
    # % evtLto from above
    
    # % SELECT SEQUENTIAL STEPS
    
    # % Now create an ordered set of sequential steps using FBi as guide
    closest = np.abs(np.repeat(L_FBi[:,np.newaxis], evtLto.shape[0], axis=1)-np.repeat(evtLto[:,np.newaxis].T, L_FBi.shape[0], axis=0))
    
    mindist, minx = np.nanmin(closest, axis=0), np.nanargmin(closest, axis=0)
    for i in np.unique(minx).astype(int):
        if np.sum(np.isin(i,minx)) > 1:
            mindist = mindist[minx!=i]
            evtLtd  = evtLtd[minx!=i]
            minx    = minx[minx!=i]
            
    # Parameter based on the frame adjustment observed from 300 datasets
    maxadj = 0.15*hz
    
    # Preallocate
    L_TO = np.empty(L_FBi.shape[0]) * np.NaN
    evFlto = np.zeros(L_FBi.shape[0])
    
    # Here we replace FB indices with TO from PCA default is to use FB
    for i in range(L_FBi.shape[0]):
        try:
            if i > np.max(minx):
                break
            elif np.isin(i,minx) and (mindist[minx==i] < maxadj).any():
                #Replace with evtLto since its more accurate
                L_TO[i] = evtLto[minx==i][0]
                evFlto[i] = 1
            else:
                #Use FFi since it is more robust
                L_TO[i] = L_FBi[i]
                
        except Exception as e:
            print(e)
            L_TO[i] = L_FBi[i]
    
    # Finally we round to get final indices
    L_TD = L_TD.round()
    L_TO = L_TO.round()
    
    
    # %% RIGHT FOOT EVENTS
    # % the same steps we just took for the left side
    
    # % Begin by creating a gross estimation of foot forwards and foot backs
    # Correct the cutoff frequency for the number of passes in the filter
    b, a = butter(2, 5/(hz/2), btype = 'low')
    # note that Python and Matlab filtfilt behaves slightly different with padding the data
    # see https://mail.python.org/pipermail/scipy-user/2014-April/035646.html
    filtered_R_heel = filtfilt(b, a, R_heel[:,2], padtype='odd')
    
    # Begin by creating a gross estimation of foot forwards and foot backs
    R_FFi = detect_peaks(-filtered_R_heel, mpd=np.round(0.35*hz),
                         show=False)
    
    if label == 'walk':
        # % Use peak foot flexion angle for foot back
        # % To deal with peaks resulting from signal flipping, threshold them
        angSig = angles['R_foot_Z'].values
        angSig[np.abs(angSig) > 90] = np.NaN
        R_FBi = detect_peaks(-angSig, mpd=np.round(0.7*hz),
                         mph=20, show=False)
    else:
        # Use rearmost position of heel marker for foot back
        R_FBi = detect_peaks(filtered_R_heel, mpd=np.round(0.35*hz),
                         show=False)
        
    # %Uncomment block below to enable more aggressive quality control of data
    
    # if (np.nanpercentile(np.abs(angles['R_foot_Z'].values), 90) > 120) & vel < 4
    #     print('Right ankle values outside of expected ranges, please ensure your shoe markers are properly placed and redo your collection')
    #     sys.exit()
    
    # Remove any leading FB
    R_FFi = R_FFi[R_FFi>L_FFi[0]]
    R_FBi = R_FBi[R_FBi>R_FFi[0]]
    
    ## find largest chunk of continuous data
    
    #We want to check before and after that there is sufficient data for analysis
    
    if (R_FFi.shape[0] < 2) or (R_FBi.shape[0] < 2):
        print('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
        sys.exit()
        
        
    # LARGEST_BLOCK
    R_FFi, R_FBi, R_block_start, R_block_end = eng.largest_block(matlab.double(list(R_FFi)), 
                                                                 matlab.double(list(R_FBi)), nargout=4)
    
    R_FFi = np.array(R_FFi).flatten().astype(int)
    R_FBi = np.array(R_FBi).flatten().astype(int)
    
    if (R_FFi.shape[0] < 2) or (R_FBi.shape[0] < 2):
        print('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
        sys.exit()
        
    # %In rare instances a the index will be in incorrect order run below again
    # %in case
    
    # % Remove any leading FF and FB
    R_FFi = R_FFi[R_FFi>L_FFi[0]]
    R_FBi = R_FBi[R_FBi>R_FFi[0]]
    
    # TOUCHDOWN
    # evtRtd from above
    
    # SELECT SEQUENTIAL STEPS
    # create an ordered set of sequential steps using FFi as guide
    closest = np.abs(np.repeat(R_FFi[:,np.newaxis], evtRtd.shape[0], axis=1)-np.repeat(evtRtd[:,np.newaxis].T, R_FFi.shape[0], axis=0))
    
    mindist, minx = np.nanmin(closest, axis=0), np.nanargmin(closest, axis=0)
    
    for i in np.unique(minx).astype(int):
        if np.sum(np.isin(i,minx)) > 1:
            mindist = mindist[minx!=i]
            evtLtd  = evtLtd[minx!=i]
            minx    = minx[minx!=i]
            
    # Parameter based on the typical frame adjustments observed in 300
    # datasets
    if label=='run':
        maxadj = 0.05*hz
    else:
        maxadj = 0.10*hz
        
    # Preallocate
    R_TD = np.empty(R_FFi.shape[0]) * np.NaN
    evFrtd = np.zeros(R_FFi.shape[0])
    
    # Here we replace FF indices with indices from evt where criteria are
    # met... the default is to use FF
    for i in range(R_FFi.shape[0]):
        try:
            if i > np.max(minx):
                break
            elif np.isin(i,minx) and (mindist[minx==i] < maxadj).any():
                #Replace with evtRtd since its more accurate
                R_TD[i] = evtRtd[minx==i][0]
                evFrtd[i] = 1
            else:
                #Use FFi since it is more robust
                R_TD[i] = R_FFi[i]
                
        except Exception as e:
            print(e)
            R_TD[i] = R_FFi[i]
            
            
    # %% TAKEOFF
    
    # % evtRto from above
    
    # % SELECT SEQUENTIAL STEPS
    
    # % Now create an ordered set of sequential steps using FBi as guide
    closest = np.abs(np.repeat(R_FBi[:,np.newaxis], evtRto.shape[0], axis=1)-np.repeat(evtRto[:,np.newaxis].T, R_FBi.shape[0], axis=0))
    
    mindist, minx = np.nanmin(closest, axis=0), np.nanargmin(closest, axis=0)
    for i in np.unique(minx).astype(int):
        if np.sum(np.isin(i,minx)) > 1:
            mindist = mindist[minx!=i]
            evtLtd  = evtLtd[minx!=i]
            minx    = minx[minx!=i]
            
    # Parameter based on the frame adjustment observed from 300 datasets
    maxadj = 0.15*hz
    
    # Preallocate
    R_TO = np.empty(R_FBi.shape[0]) * np.NaN
    evFrto = np.zeros(R_FBi.shape[0])
    
    # Here we replace FB indices with TO from PCA default is to use FB
    for i in range(R_FBi.shape[0]):
        try:
            if i > np.max(minx):
                break
            elif np.isin(i,minx) and (mindist[minx==i] < maxadj).any():
                #Replace with evtLto since its more accurate
                R_TO[i] = evtRto[minx==i][0]
                evFrto[i] = np.float64(1)
            else:
                #Use FFi since it is more robust
                R_TO[i] = R_FBi[i]
                
        except Exception as e:
            print(e)
            R_TO[i] = R_FBi[i]
    
    # Finally we round to get final indices
    R_TD = R_TD.round()
    R_TO = R_TO.round()
    
    # %% if largest chunk of continuous data not at beginning, chop both right and left so they match
    
    
    # %index must begin with left touchdown and end with right toe
    # %off
    
    
    if R_block_start < L_block_start:
        #remove all right indices that occur before the first left touchdown
        
        R_TO = R_TO[(R_TD < L_block_start)!=1]
        R_TD = R_TD[(R_TD < L_block_start)!=1]
        
    flag = 0
    
    if L_block_start < R_block_start:
        #remove left touchdowns more than one touchdown before the first right touchdown
        cut_inds = (L_TD<R_block_start)==1
        cut_inds.astype(int)
        #this loop ensures the first index will be a left touchdown
        for i in range(cut_inds.shape[0]):
            if (cut_inds[i]==1) and (cut_inds[i+1]==0) and flag==0:
                cut_inds[i] = 0
                flag = 1
                
        L_TD = np.delete(L_TD, cut_inds)
        L_TO = np.delete(L_TO, cut_inds)
        
        
    # create an events matrix
    # Remove trailing nans that may have crept in
    evFltd = evFltd[~np.isnan(L_TD)]
    evFlto = evFlto[~np.isnan(L_TO)]
    evFrtd = evFrtd[~np.isnan(R_TD)]
    evFrto = evFrto[~np.isnan(R_TO)]
    
    L_TD = L_TD[~np.isnan(L_TD)]
    L_TO = L_TO[~np.isnan(L_TO)]
    R_TD = R_TD[~np.isnan(R_TD)]
    R_TO = R_TO[~np.isnan(R_TO)]
    
    # Find the closest ordered pairs of L_TO and R_TD to synchronize steps
    closest = np.abs(np.repeat(R_TD[:,np.newaxis], L_TO.shape[0], axis=1)-np.repeat(L_TO[:,np.newaxis].T, R_TD.shape[0], axis=0))
    
    minx = np.nanargmin(closest,axis=0)
    
    #Truncate right stances to match up with left
    evFrtd = evFrtd[np.unique(minx)]
    R_TD = R_TD[np.unique(minx)]
    
    
    testlength = np.min([L_TO.shape[0], R_TD.shape[0]])
    if np.median(L_TO[:testlength]-R_TD[:testlength]) < 0:#Then there is a flight phase
        #Find the closest ordered pairs of R_TD and R_TO to synchronize steps
        closest = np.abs(np.repeat(R_TO[:,np.newaxis], R_TD.shape[0], axis=1)-np.repeat(R_TD[:,np.newaxis].T, R_TO.shape[0], axis=0))
        minx = np.nanargmin(closest,axis=0)
        
    else: # There is no flight phase i.e. grounded running or walking
        #Find the closest ordered pairs of R_TO and L_TD to synchronize steps
        tmp = L_TD[1:]
        closest = np.abs(np.repeat(R_TO[:,np.newaxis], tmp.shape[0], 
                                   axis=1)-np.repeat(tmp[:,np.newaxis].T, R_TO.shape[0], axis=0))
        minx = np.nanargmin(closest,axis=0)
        
    evFrto = evFrto[np.unique(minx)]
    R_TO = R_TO[np.unique(minx)]
    
    events = [L_TD.shape[0], L_TO.shape[0], R_TD.shape[0], R_TO.shape[0]]
    
    # Chop everything to the same length
    L_TD = L_TD[:min(events)]
    L_TO = L_TO[:min(events)]
    R_TD = R_TD[:min(events)]
    R_TO = R_TO[:min(events)]
    
    evFltd = evFltd[:min(events)]
    evFlto = evFlto[:min(events)]
    evFrtd = evFrtd[:min(events)]
    evFrto = evFrto[:min(events)]
    
    # Very rarely, these will wind up empty and assignment doesn't work
    events = np.empty(shape=(L_TD.shape[0],4)) * np.NaN
    events[:,0]=L_TD
    events[:,1]=L_TO
    events[:,2]=R_TD
    events[:,3]=R_TO
    # Very rarely, these will wind up empty and assignment doesn't work
    eventsflag = np.empty(shape=(evFltd.shape[0],4)) * np.NaN
    eventsflag[:,0]=evFltd
    eventsflag[:,1]=evFlto
    eventsflag[:,2]=evFrtd
    eventsflag[:,3]=evFrto
    
    # Remove first row since these will very often be reliant on FF and FB measures
    if events.shape[0] > 1:
        events = np.delete(events,0,axis=0)
        eventsflag = np.delete(eventsflag,0,axis=0)
        
    # %% Occasionally, one stance will drop out, and data becomes
    # % discontinuous...this fix alleviates this by trimming data to largest
    # % continuous block
    try:
        cont = np.array([events[1:,0]>events[:-1,1],events[1:,2]>events[:-1,3]])
        cont = np.hstack((np.zeros((2,1)),cont,np.zeros((2,1)))).T
        F = np.where(np.any(cont==0, axis=1))
        F = np.asarray(F).flatten()
        D = np.diff(F)-2
        M, L = np.max(D), np.argmax(D)
        events = events[F[L]:F[L]+M+1,:]
        eventsflag = eventsflag[F[L]:F[L]+M+1,:]
    except Exception as e:
        print('Could not obtain a continuous block of events')
        events = []
        eventsflag = []
        print(e)
        
    # Worst-case... return to foot forward, foot back detection
    if events.shape[0] < 5:
        print('Automated event detection failed, defaulting to foot-forward foot-back')
        nevents = [L_FFi.shape[0], L_FBi.shape[0], R_FFi.shape[0], R_FBi.shape[0]]
        
        # Chop everything to the same length
        L_FFi = L_FFi[:min(nevents)]
        L_FBi = L_FBi[:min(nevents)]
        R_FFi = R_FFi[:min(nevents)]
        R_FBi = R_FBi[:min(nevents)]
        
        events = np.empty(shape=(min(nevents),4)) * np.NaN
        events[:,0] = L_FFi
        events[:,1] = L_FBi
        events[:,2] = R_FFi
        events[:,3] = R_FBi
        
    # Pull event columns from events so everything is consistent
    L_TD = events[:,0]
    L_TO = events[:,1]
    R_TD = events[:,2]
    R_TO = events[:,3]
    
    return L_TD, L_TO, R_TD, R_TO, eventsflag, label