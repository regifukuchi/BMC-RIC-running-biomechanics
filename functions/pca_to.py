"""Applies a pre-trained algorithm to detect toeoff events"""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"

import sys
import numpy as np
from scipy import signal
from detecta import detect_peaks

def pca_to(angles, hz, event_data, gait_mode):
    """
    %PCA_TO applies a pre-trained algorithm to detect toeoff events.
    %
    %   This function takes in the ANGLES output from a gait_kinematics.m function,
    %   the sampling frequency HZ, and the EVENT_DATA dictionary; and outputs toeoff events as
    %   determined by a trained PCA algorithm.  GAIT is a label for type of
    %   gait (must be either 'walk' or 'run') which determines the model that
    %   will be applied.  EVTL and EVTR are vectors of indices corresponding to
    %   the timing of touchdowns from the original ANGLES signals. NOTE: EVTL
    %   and EVTR are not rounded... rounding is completed in *_steps functions.
    %
    %   Part of this code is based on the programs written by Sean Osis
    """
    
    
    # Input parameters
    # hz = 200
    # gait_mode = yhat
    # angles = angles
    # Load PCA output from mat file
    # % event_data is a .mat file containing 'coeff' which is the coefficients
    # % from the pre-trained PCA and 'p' which is the list of coefficients of the
    # % linear polynomial relating PCA scores with touchdown timing relative to
    # % the foot acceleration peak.
    # event_data = spio.loadmat(os.path.join(pathname, 'event_data_TO.mat'))
    
    # Future use for separate left and right coefficients if needed... right
    # now there appears to be no benefit to having them separate
    coeffL = event_data['coeff']
    coeffR = event_data['coeff']
    
    if gait_mode=='run':
        pL = event_data['pRun'].flatten()
        pR = event_data['pRun'].flatten()
    elif gait_mode=='walk':
        pL = event_data['pWalk'].flatten()
        pR = event_data['pWalk'].flatten()
    
    # Default sampling rate for which model was originally developed
    defaultHz = 200
    
    #Chunk length for PCA... determined to be 35 frames on either side of the
    #foot accel peak which is 35/200 = 0.175
    chnklgth = round(0.175*defaultHz)
    
    # Minimum distance between foot dorsiflexion peaks
    minpkdist = round(0.5*defaultHz)
    
    #Minimum peak height for finding dorsiflexion peaks
    posminpkht = 20
    
    # % Bias added to event timing to compensate for the bias evident from the
    # % methodology used in Osis et al. (2014).  This bias may be due to changes
    # % in timing when differentiating.  Bias currently removed as new model does
    # % not seem to have offset.
    bias = 0*defaultHz
    
    #%% %% Resample signal if needed
    # Original model trained on 200Hz
    
    if hz < 100:
        print("Error! 'Sampling frequency is less than 100 Hz. Event detection may provide inconsistent results at low sampling rates.'!")
        sys.exit()
    
    elif hz!=defaultHz:
        
        # Resample the signals to match 200 Hz for methods below
        temp = angles.values
        temp = signal.resample_poly(temp, defaultHz, hz)
        
        angles.values = temp
        
    # %% Right Side Toeoff Detections
    sig = -angles['foot_Z'].values
    
    # Detect desired positive peaks
    locs = detect_peaks(sig, mpd=minpkdist, 
                        mph=posminpkht, show=False)
    
    #Create acceleration signals for all segments and joints
    foot = np.diff(angles['foot_Z'].values, 2)
    ank  = np.diff(angles['ankle_Z'].values, 2)
    knee = np.diff(angles['knee_Z'].values, 2)
    hip  = np.diff(angles['hip_Z'].values, 2)
    
    signal_ = np.zeros(shape=(locs.shape[0],chnklgth*8+4))
    
    # Fill the signal matrix so that the PCA coefficients can be applied
    for j in list(range(1,locs.shape[0]-1)):#% Skip first and last peaks since we might not have enough data
        # Build signal using - 70 frames of data from each joint
        signal_[j,:] = np.hstack((foot[locs[j]-chnklgth*2:locs[j]+1],
               ank[locs[j]-chnklgth*2:locs[j]+1],
               knee[locs[j]-chnklgth*2:locs[j]+1],
               hip[locs[j]-chnklgth*2:locs[j]+1]))
        
    #Multiply the signal by the 'pre-trained' PCA coefficients
    projected = np.dot(signal_,coeffR)
    
    # Apply the linear polynomial to predict the frame difference
    pred = np.polyval(pR, projected[:,2])
    
    # Apply frame difference to the foot accel peak timing
    evtR = -pred+locs+np.repeat(bias, pred.shape[0])
    
    # Original prediction equation is for 200Hz
    evtR = evtR * (hz/defaultHz)
    
    # Remove first and last events since these are not properly calculated due
    # to lack of data
    evtR = np.delete(evtR, (0,-1))
    
    return evtR