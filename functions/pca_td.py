"""
%PCA_TD applies a pre-trained algorithm to detect touchdown events.
%
%   [EVTL,EVTR] = PCA_TD(ANGLES,HZ,GAIT)
%
%   This function takes in the ANGLES output from a *_kinematics.m function
%   and the sampling frequency HZ, and outputs touchdown events as
%   determined by the Osis et al. (2014) method.  GAIT is a label for type
%   of gait (must be either 'walk' or 'run') which determines the model
%   that will be applied.  EVTL and EVTR are vectors of indices
%   corresponding to the timing of touchdowns from the original ANGLES
%   signals. NOTE: EVTL and EVTR are not rounded... rounding is completed
%   in *_steps functions.
%
%   Created By: Sean Osis on March 10, 2014
%
%   Copyright (C) 2014, Sean Osis and The Running Injury Clinic
%
%   Update, July, 2015: Code now implements event detection for both
%   walking and running using model trained for both.  Bias is currently
%   eliminated as it does not drastically improve fit.
%
%   Update, Jan, 2016: Code has been tested for running and walking with
%   full research database.  Roughly 3% of cases demonstrate issues with
%   event detection for running, mostly due to atypical movement patterns
%   and very fast speeds ~4m/s.  Roughly 7% of cases demonstrate issues
%   with detection for walking, mostly due to large variability in walking
%   strides, however, many of these can be screened for contigous blocks or
%   to have bad steps removed.
"""
# Prepare Python environment
import scipy.io
import os

# Input parameters
mat_dir = r'C:\Users\Reginaldo\Documents\Github\RIC_3Dgait\functions'
gait = 'run'
hz = 200
# Load PCA output from mat file
# % event_data is a .mat file containing 'coeff' which is the coefficients
# % from the pre-trained PCA and 'p' which is the list of coefficients of the
# % linear polynomial relating PCA scores with touchdown timing relative to
# % the foot acceleration peak.
mat = scipy.io.loadmat(os.path.join(mat_dir, 'event_data_TD.mat'))

# Future use for separate left and right coefficients if needed... right
# now there appears to be no benefit to having them separate
coeffL = mat['coeff']
coeffR = mat['coeff']

pL = mat['p'].flatten()
pR = mat['p'].flatten()
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
if gait=='walk':
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
       
elif hz != defaultHz:
    
    # Resample the signals to match 200 Hz for methods below
    temp = struct2cell(angles);
    for i = 1:length(temp)
        temp{i} = resample(temp{i},defaultHz,hz);
    end
    
    angles = cell2struct(temp,fieldnames(angles),1);
    
end

## CONTINUE HERE. RESAMPLE SIGNAL IF FREQ IS DIFFERENT