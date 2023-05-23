"""Applies a pre-trained algorithm to detect touchdown events"""

__author__ = "Reginaldo K Fukuchi, https://github.com/regifukuchi"
__version__ = "1.0.1"
__license__ = "MIT"

import sys
import numpy as np
from scipy import signal
from detecta import detect_peaks

def pca_td(angles, hz, event_data, gait_mode):
    """
    %PCA_TD applies a pre-trained algorithm to detect touchdown events.
    %
    %   This function takes in the ANGLES output from a gait_kinematics.py function,
    %   the sampling frequency HZ and EVENT_DATA from PCA output dictionary; and outputs touchdown events as
    %   determined by the Osis et al. (2014) method.  GAIT_MODE is a label for type
    %   of gait (must be either 'walk' or 'run') which determines the model
    %   that will be applied.  EVTL and EVTR are vectors of indices
    %   corresponding to the timing of touchdowns from the original ANGLES
    %   signals. NOTE: EVTL and EVTR are not rounded... rounding is completed
    %   in *_steps functions.
    %
    %   Part of this code is based on the programs written by Sean Osis
    """
    
    
    # Input parameters
    # hz = 200
    # gait_mode = 'run'
    # angles = angles
    # Load PCA output from mat file
    # % event_data is a .mat file containing 'coeff' which is the coefficients
    # % from the pre-trained PCA and 'p' which is the list of coefficients of the
    # % linear polynomial relating PCA scores with touchdown timing relative to
    # % the foot acceleration peak.
    # mat = scipy.io.loadmat(os.path.join(pathname, 'event_data_TD.mat'))
    
    # Future use for separate left and right coefficients if needed... right
    # now there appears to be no benefit to having them separate
    coeffL = event_data['coeff']
    coeffR = event_data['coeff']
    
    pL = event_data['p'].flatten()
    pR = event_data['p'].flatten()
    
    #%% Default sampling rate for which model was originally developed
    defaultHz = 200
    
    # The length of the search window when finding positive peaks in the foot
    # acceleration signal
    srchlgth = round(0.175*defaultHz)
    
    # Chunk length for PCA... determined to be 35 frames on either side of the
    # foot accel peak which is 35/200 = 0.175
    #chnklgth = round(0.175*hz)
    chnklgth = round(0.175*defaultHz)
    
    # Minimum distance between foot acceleration peaks
    if gait_mode=='walk':
        minpkdist = round(0.7*defaultHz)
    else:
        minpkdist = round(0.5*defaultHz)
    
    # Minimum peak height for finding negative acceleration peaks
    negminpkht = 0.1
    
    # Minimum peak height for finding positive acceleration peaks
    posminpkht = 0.01
    
    # Bias added to event timing to compensate for the bias evident from the
    # methodology used in Osis et al. (2014).  This bias may be due to changes
    # in timing when differentiating.  Bias should be added in seconds in order
    # to allow for different frame rates
    bias = 0*defaultHz;
    
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
        
    #%% Right Side Touchdown Detections
    
    #% Flip the foot in sagittal and differentiate
    negsig = -np.diff(angles['foot_Z'].values, 2)
    
    # Find peaks location. Signal flipped because of X-axis convention difference.
    locs = detect_peaks(negsig, mpd=minpkdist, 
                        mph=posminpkht, show=False)
    
    # Create a logical index of search windows to find positive peaks
    loginds = np.zeros(negsig.shape[0])
    
    for i in range(locs.shape[0]):
        loginds[locs[i]:locs[i]+srchlgth] = 1
    
    loginds = loginds[:negsig.shape[0]]
    
    #Floor unwanted signal parts to zero to more easily find positive peaks
    possig = -negsig * loginds
    
    #Detect desired positive peaks
    locs = detect_peaks(possig, mpd=minpkdist, 
                        mph=posminpkht, show=False)
    
    #Create acceleration signals for all segments and joints
    foot = np.diff(angles['foot_Z'].values, 2)
    ank  = np.diff(angles['ankle_Z'].values, 2)
    knee = np.diff(angles['knee_Z'].values, 2)
    hip  = np.diff(angles['hip_Z'].values, 2)
    
    signal_ = np.zeros(shape=(locs.shape[0],chnklgth*8+4))
    
    # Fill the signal matrix so that the PCA coefficients can be applied
    for j in list(range(1,locs.shape[0]-1)):#% Skip first and last peaks since we might not have enough data
        signal_[j,:] = np.hstack((foot[locs[j]-chnklgth:locs[j]+chnklgth+1],
               ank[locs[j]-chnklgth:locs[j]+chnklgth+1],
               knee[locs[j]-chnklgth:locs[j]+chnklgth+1],
               hip[locs[j]-chnklgth:locs[j]+chnklgth+1]))
    
    #Multiply the signal by the 'pre-trained' PCA coefficients
    projected = np.dot(signal_,coeffR)
    
    # Apply the linear polynomial to predict the frame difference
    pred = np.polyval(pR, projected[:,1])
    
    # Apply frame difference to the foot accel peak timing
    evtR = pred+locs+np.repeat(bias, pred.shape[0])
    
    # Original prediction equation is for 200Hz
    evtR = evtR * (hz/defaultHz)
    
    # Remove first and last events since these are not properly calculated due
    # to lack of data
    evtR = np.delete(evtR, (0,-1))
    
    return evtR