function [evtL,evtR] = pca_to(angles,hz,gait)

%PCA_TO applies a pre-trained algorithm to detect toeoff events.
%
%   [EVTL,EVTR] = PCA_TO(ANGLES,HZ,GAIT)
%
%   This function takes in the ANGLES output from a *_kinematics.m function
%   and the sampling frequency HZ, and outputs touchdown events as
%   determined by a trained PCA algorithm.  GAIT is a label for type of
%   gait (must be either 'walk' or 'run') which determines the model that
%   will be applied.  EVTL and EVTR are vectors of indices corresponding to
%   the timing of touchdowns from the original ANGLES signals. NOTE: EVTL
%   and EVTR are not rounded... rounding is completed in *_steps functions.
%
%   Created By: Sean Osis on July 30, 2015
%
%   Copyright (C) 2014-2015, Sean Osis and The Running Injury Clinic
%
%   Update, Jan, 2016: Code has been tested for running and walking with
%   full research database. Almost zero issues noted for toe-off detection.



%% Load and initialize parameters

% event_data is a .mat file containing 'coeff' which is the coefficients
% from the pre-trained PCA and 'p' which is the list of coefficients of the
% linear polynomial relating PCA scores with touchdown timing relative to
% the foot acceleration peak.
load event_data_TO

% Future use for separate left and right coefficients if needed... right
% now there appears to be no benefit to having them separate
coeffL = coeff;
coeffR = coeff;

switch gait
    case 'run'
        pL = pRun;
        pR = pRun;
    case 'walk'
        pL = pWalk;
        pR = pWalk;
end


% Default sampling rate for which model was originally developed
defaultHz = 200;

% Chunk length for PCA... determined to be 35 frames on either side of the
% foot accel peak which is 35/200 = 0.175
chnklgth = round(0.175*defaultHz);

% Minimum distance between foot dorsiflexion peaks
minpkdist = round(0.5*defaultHz);

% Minimum peak height for finding dorsiflexion peaks
posminpkht = 20;

% Bias added to event timing to compensate for the bias evident from the
% methodology used in Osis et al. (2014).  This bias may be due to changes
% in timing when differentiating.  Bias currently removed as new model does
% not seem to have offset.
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

% Flip the foot in sagittal 
sig = -(angles.L_foot(:,3));

% Detect desired positive peaks
[~,locs] = findpeaks(sig,'minpeakdistance',minpkdist,'minpeakheight',posminpkht);

% Create acceleration signals for all segments and joints
foot = diff(angles.L_foot(:,3),2);
ank = diff(angles.L_ankle(:,3),2);
knee = diff(angles.L_knee(:,3),2);
hip = diff(angles.L_hip(:,3),2);

signal = zeros(length(locs),chnklgth*8+4);

% Fill the signal matrix so that the PCA coefficients can be applied
for j = 2:length(locs)-1  % Skip first and last peaks since we might not have enough data
    % Build signal using - 70 frames of data from each joint
    signal(j,:) = [foot(locs(j)-chnklgth*2:locs(j))' ank(locs(j)-chnklgth*2:locs(j))' knee(locs(j)-chnklgth*2:locs(j))' hip(locs(j)-chnklgth*2:locs(j))'];
end

% Multiply the signal by the 'pre-trained' PCA coefficients
projected = signal*coeffL;

% Apply the linear polynomial to predict the frame difference
pred = polyval(pL,projected(:,3));

% Apply frame difference to the foot dorsiflexion peak timing
evtL = -pred+locs(:)+repmat(bias,length(pred),1);

% Original prediction equation is for 200Hz
evtL = evtL.*(hz/defaultHz);


%% Right Side Touchdown Detections

% Repeat everything for right side

sig = -(angles.R_foot(:,3));

[~,locs] = findpeaks(sig,'minpeakdistance',minpkdist,'minpeakheight',posminpkht);

foot = diff(angles.R_foot(:,3),2);
ank = diff(angles.R_ankle(:,3),2);
knee = diff(angles.R_knee(:,3),2);
hip = diff(angles.R_hip(:,3),2);

signal = zeros(length(locs),chnklgth*8+4);

for j = 2:length(locs)-1
    
    signal(j,:) = [foot(locs(j)-chnklgth*2:locs(j))' ank(locs(j)-chnklgth*2:locs(j))' knee(locs(j)-chnklgth*2:locs(j))' hip(locs(j)-chnklgth*2:locs(j))'];
end 

projected = signal*coeffR;

pred = polyval(pR,projected(:,3));

evtR = -pred+locs(:)+repmat(bias,length(pred),1);

% Original prediction equation is for 200Hz
evtR = evtR.*(hz/defaultHz);


% Remove first and last events since these are not properly calculated due
% to lack of data
evtL([1 end],:) = [];
evtR([1 end],:) = [];


