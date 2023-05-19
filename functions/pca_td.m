function [evtL,evtR] = pca_td(angles,hz,gait)

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


%% Load and initialize parameters

% event_data is a .mat file containing 'coeff' which is the coefficients
% from the pre-trained PCA and 'p' which is the list of coefficients of the
% linear polynomial relating PCA scores with touchdown timing relative to
% the foot acceleration peak.
load event_data_TD
% Future use for separate left and right coefficients if needed... right
% now there appears to be no benefit to having them separate
coeffL = coeff;
coeffR = coeff;
pL = p;
pR = p;

% Default sampling rate for which model was originally developed
defaultHz = 200;

% The length of the search window when finding positive peaks in the foot
% acceleration signal
srchlgth = round(0.175*defaultHz);

% Chunk length for PCA... determined to be 35 frames on either side of the
% foot accel peak which is 35/200 = 0.175
%chnklgth = round(0.175*hz);
chnklgth = round(0.175*defaultHz);

% Minimum distance between foot acceleration peaks
if strcmp(gait,'walk')
    minpkdist = round(0.7*defaultHz);
else
    minpkdist = round(0.5*defaultHz);
end

% Minimum peak height for finding negative acceleration peaks
negminpkht = 0.1;

% Minimum peak height for finding positive acceleration peaks
posminpkht = 0.01;

% Bias added to event timing to compensate for the bias evident from the
% methodology used in Osis et al. (2014).  This bias may be due to changes
% in timing when differentiating.  Bias should be added in seconds in order
% to allow for different frame rates
bias = 0*defaultHz;


%% Resample signal if needed
% Original model trained on 200Hz

if hz < 100
    error('Sampling frequency is less than 100 Hz. Event detection may provide inconsistent results at low sampling rates.')
elseif hz ~= defaultHz
    
    % Resample the signals to match 200 Hz for methods below
    temp = struct2cell(angles);
    for i = 1:length(temp)
        temp{i} = resample(temp{i},defaultHz,hz);
    end
    
    angles = cell2struct(temp,fieldnames(angles),1);
    
end


%% Left Side Touchdown Detections

% Flip the foot in sagittal and differentiate
negsig = -diff(angles.L_foot(:,3),2);

% Find negative peaks... this is much more reliable than finding positive
% peaks which is what we need
[~,locs] = findpeaks(negsig,'minpeakdistance',minpkdist,'minpeakheight',negminpkht);

% Create a logical index of search windows to find positive peaks
loginds = zeros(length(negsig),1);
for i = 1:length(locs)
    
    loginds(locs(i):locs(i)+srchlgth) = 1;
    
end
loginds = loginds(1:length(negsig));

% Floor unwanted signal parts to zero to more easily find positive peaks
possig = -negsig.*loginds;

% Detect desired positive peaks
[~,locs] = findpeaks(possig,'minpeakdistance',minpkdist,'minpeakheight',posminpkht);

% Create acceleration signals for all segments and joints
foot = diff(angles.L_foot(:,3),2);
ank = diff(angles.L_ankle(:,3),2);
knee = diff(angles.L_knee(:,3),2);
hip = diff(angles.L_hip(:,3),2);

signal = zeros(length(locs),chnklgth*8+4);

% Fill the signal matrix so that the PCA coefficients can be applied
for j = 2:length(locs)-1  % Skip first and last peaks since we might not have enough data
    % Build signal using +/- 35 frames of data from each joint
    signal(j,:) = [foot(locs(j)-chnklgth:locs(j)+chnklgth)' ank(locs(j)-chnklgth:locs(j)+chnklgth)' knee(locs(j)-chnklgth:locs(j)+chnklgth)' hip(locs(j)-chnklgth:locs(j)+chnklgth)'];
end

% Multiply the signal by the 'pre-trained' PCA coefficients
projected = signal*coeffL;

% Apply the linear polynomial to predict the frame difference
pred = polyval(pL,projected(:,2));

% Apply frame difference to the foot accel peak timing
evtL = pred+locs(:)+repmat(bias,length(pred),1);

% Original prediction equation is for 200Hz
evtL = evtL.*(hz/defaultHz);


%% Right Side Touchdown Detections

% Repeat everything for right side

negsig = -diff(angles.R_foot(:,3),2);

[~,locs] = findpeaks(negsig,'minpeakdistance',minpkdist,'minpeakheight',negminpkht);

loginds = zeros(length(negsig),1);

for i = 1:length(locs)
    loginds(locs(i):locs(i)+srchlgth) = 1;
end

loginds = loginds(1:length(negsig));

possig = -negsig.*loginds;

[~,locs] = findpeaks(possig,'minpeakdistance',minpkdist,'minpeakheight',posminpkht);

foot = diff(angles.R_foot(:,3),2);
ank = diff(angles.R_ankle(:,3),2);
knee = diff(angles.R_knee(:,3),2);
hip = diff(angles.R_hip(:,3),2);

signal = zeros(length(locs),chnklgth*8+4);

for j = 2:length(locs)-1
    
    signal(j,:) = [foot(locs(j)-chnklgth:locs(j)+chnklgth)' ank(locs(j)-chnklgth:locs(j)+chnklgth)' knee(locs(j)-chnklgth:locs(j)+chnklgth)' hip(locs(j)-chnklgth:locs(j)+chnklgth)'];
    
end

projected = signal*coeffR;

pred = polyval(pR,projected(:,2));

evtR = pred+locs(:)+repmat(bias,length(pred),1);

% Original prediction equation is for 200Hz
evtR = evtR.*(hz/defaultHz);

% Remove first and last events since these are not properly calculated due
% to lack of data
evtL([1 end],:) = [];
evtR([1 end],:) = [];