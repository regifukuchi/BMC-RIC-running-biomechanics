function [norm_ang,norm_vel,events,event,DISCRETE_VARIABLES,speedoutput,eventsflag,label] = gait_steps(neutral,dynamic,angles,velocities,hz,plots)
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
%  NEUTRAL (struct):    Marker shell positions collected as part of 
%                       the static trial.
%
%  DYNAMIC (stuct):     Marker shell positions collected as part of
%                       the dynamic (run/walk) trial.
%
%  ANGLES (struct):     Angles (joint angles) structure created as an 
%                       output from the function: gait_kinematics.
%
%  VELOCITIES (struct): Velocities (joint velocities) structure created as 
%                       an output from the function: gait_kinematics.
%
%  HZ (int):        Data collection sampling frequency.
%
%  PLOTS (bool):    Boolean selected to generate plotted outcomes if 
%                   desired. If no second argument exists, or if PLOTS == 0
%                   , the plotted outputs in this function are suppressed.
%
%
%  OUTPUTS
%  -------
%  NORM_ANG (struct):   Normalized angles from touchown to takeoff across 
%                       all retained steps.
%
%  NORM_VEL (struct):   Normalized velocities from touchown to takeoff 
%                       across  all retained steps.
%
%  EVENTS (mat):    Matrix of frame numbers for touchdown and toeoffs.
%
%  EVENT (mat):     Same matric as EVENTS but also includes midswing in
%                   order to calculate swing variables
%
%  DISCRETE_VARIABLES (mat):    Contains the peaks and values of interest 
%                               for reporting.
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

%% if no indication to plot the data was specified, default to not plot
if (nargin < 7)
    plots = 0;
end

%% Determine functional measures and gait type (walk vs run)
% movement speed comes from the A/P position time history of a heel marker
% so we first need to identify a heel marker

% Combine 3 of the foot markers into one matrix (ignore the created fourth)
L_foot = [neutral.L_foot_1;neutral.L_foot_2;neutral.L_foot_3];
% sort the markers from left to right
[L_foot,i_lf] = sortrows(L_foot,1);
% find the lower of the two medial markers
if L_foot(2,2) < L_foot(3,2)
    L_marker = strcat('L_foot_',num2str(i_lf(2)));
    L_heel = dynamic.(L_marker);
else
    L_marker = strcat('L_foot_',num2str(i_lf(3)));
    L_heel = dynamic.(L_marker);
end

feature accel off
[pks,~] = findpeaks(diff(L_heel(:,3)),'minpeakdistance',round(0.5*hz),'minpeakheight',0);
feature accel on

vel = hz*median(pks)/1000;


feature accel off
[~,locs] = findpeaks(-diff(L_heel(:,3)),'minpeakdistance',round(0.5*hz),'minpeakheight',0);
feature accel on

stRate = 60/(median(diff(locs))/hz);

speedoutput=vel;


% RIGHT SIDE
% Combine 3 of the foot markers into one matrix (ignore the created fourth)
R_foot = [neutral.R_foot_1;neutral.R_foot_2;neutral.R_foot_3];
% sort the markers from left to right
[R_foot,i_rf] = sortrows(R_foot,1);
% find the lower of the two medial markers
if R_foot(1,2) < R_foot(2,2)
    R_marker = strcat('R_foot_',num2str(i_rf(1)));
    R_heel = dynamic.(R_marker);
else
    R_marker = strcat('R_foot_',num2str(i_rf(2)));
    R_heel = dynamic.(R_marker);
end


% Identify gait type using a trained LDA classifier.  This will be more
% robust for shuffle-runners, older adults and speed walkers.  gaitClass
% represents an LDA object which has been trained on 839 test sets of
% walking and running, and validated on ~2000 sets of walking and running.
testSet = [vel stRate];
import classreg.learning.classif.CompactClassificationDiscriminant
load('gaitClass.mat','gaitClass')


label = predict(gaitClass,testSet);
%label returned as cell
label = label{1};


%% Identify Touch Down and Take Off events: Gait Independent

if plots~=0
    if strcmp(label,'walk')
        disp(['The subject was WALKING at ' num2str(vel) 'm/s or ' num2str(vel*3.6/1.6) 'mph' ]); disp(' ');
    else
        disp(['The subject was RUNNING at ' num2str(vel) 'm/s or ' num2str(vel*3.6/1.6) 'mph' ]); disp(' ');
        
    end
end


% Use PCA touchdown detection based on updated Osis et al. (2014) for
% both walking and running.
% Use new PCA toeoff detection for both walking and running.
% evt variables are NOT rounded

feature accel off
try
    [evtLtd,evtRtd] = pca_td(angles,hz,label);
    [evtLto,evtRto] = pca_to(angles,hz,label);
catch ME
    %For a small number of people, these functions return errors, or in the
    %case of bad data... default to use FF and FB in these cases
    
    evtLtd = [];
    evtRtd = [];
    evtLto = [];
    evtRto = [];
    
    disp('Automated event detection failed, defaulting to foot-forward foot-back')
    
    ME.message;
end
feature accel on


%% LEFT FOOT EVENTS

% when the feet are not tracked very well, discontinuities in the heel
% marker can occur causing the findpeaks to pick up additional 'peaks'
% for the purposes of simply identifying foot forward and foot back
% timing, we can over filter this signal. We do not care about the
% magnitude of the signal but only the timing so we can overfit as long
% as the filter has a zero phase shift.
% Note: signal is now filtered by default.  There is no advantage to not
% filtering, as if the signal quality is already good, then the system uses
% PCA event detection anyhow, and if the signal is bad, then it has to be
% filtered in order to get foot-forward foot-backward events.

feature accel off

[B,A] = butter(2, 5/(hz/2), 'low');
filtered_L_heel = filtfilt(B,A,L_heel(:,3));

% Begin by creating a gross estimation of foot forwards and foot backs
[~,L_FFi] = findpeaks(-filtered_L_heel(:,1),'minpeakdistance',round(0.35*hz));

if strcmp(label,'walk')
    % Use peak foot flexion angle for foot back
    % To deal with peaks resulting from signal flipping, threshold them
    angSig = angles.L_foot(:,3);
    angSig(abs(angSig)>90,:) = nan;
    [~,L_FBi] = findpeaks(-angSig,'minpeakdistance',round(0.7*hz),'minpeakheight',20);
else
    % Use rearmost position of heel marker for foot back
    [~,L_FBi] = findpeaks(filtered_L_heel(:,1),'minpeakdistance',round(0.35*hz));
end
feature accel on

%Uncomment block below to enable more aggressive quality control of data

% if prctile(abs(angles.L_foot(:,3)),90) > 120 && vel < 4
%
% error('Left ankle values outside of expected ranges, please ensure your shoe markers are properly placed and redo your collection')
%
% end


% Remove any leading FB
L_FBi(L_FBi<L_FFi(1)) = [];

%% find largest chunk of continuous data

%We want to check before and after that there is sufficient data for
%analysis
if size(L_FFi,1) < 2 || size(L_FBi,1) < 2
    error('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
end


[L_FFi, L_FBi, L_block_start, ~]=largest_block(L_FFi, L_FBi);

if size(L_FFi,1) < 2 || size(L_FBi,1) < 2
    error('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
end

%% TOUCHDOWN
% evtLtd from above


% SELECT SEQUENTIAL STEPS

% create an ordered set of sequential steps using FFi as guide
closest = abs(repmat(L_FFi(:),1,length(evtLtd))-repmat(evtLtd(:)',length(L_FFi),1));
[mindist,minx] = nanmin(closest,[],1);
for i = unique(minx)
    if sum(ismember(minx,i))>1
        mindist(minx==i) = [];
        evtLtd(minx==i) = [];
        minx(minx==i) = [];
    end
end

% Parameter based on the typical frame adjustments observed in 300
% datasets
if strcmp(label,'run')
    maxadj = 0.05*hz;
else
    maxadj = 0.10*hz;
end


% Preallocate
L_TD = nan(length(L_FFi),1);
evFltd = zeros(length(L_FFi),1);

% Here we replace FF indices with indices from evt where criteria are
% met... the default is to use FF
for i = 1:length(L_FFi)
    try
        
        if i > max(minx)
            break
        elseif ismember(i,minx) && mindist(minx==i)<maxadj
            % Replace with evtLtd since its more accurate
            L_TD(i) = evtLtd(minx==i);
            evFltd(i) = 1;
        else
            % Use FFi since it is more robust
            L_TD(i) = L_FFi(i);
        end
        
    catch err
        err.message;
        L_TD(i) = L_FFi(i);
    end
end


% TAKEOFF

% evtLto from above

% SELECT SEQUENTIAL STEPS

% Now create an ordered set of sequential steps using FBi as guide
closest = abs(repmat(L_FBi(:),1,length(evtLto))-repmat(evtLto(:)',length(L_FBi),1));

[mindist,minx] = nanmin(closest,[],1);
for i = unique(minx)
    if sum(ismember(minx,i))>1
        mindist(minx==i) = [];
        evtLto(minx==i) = [];
        minx(minx==i) = [];
    end
end

%Parameter based on the frame adjustment observed from 300 datasets
maxadj = 0.15*hz;

% Preallocate
L_TO = nan(length(L_FBi),1);
evFlto = zeros(length(L_FBi),1);

% Here we replace FB indices with TO from PCA default is to use FB
for i = 1:length(L_FBi)
    try
        
        if i > max(minx)
            break
        elseif ismember(i,minx) && mindist(minx==i)<maxadj
            % Replace with evtLto since its more accurate
            L_TO(i) = evtLto(minx==i);
            evFlto(i) = 1;
        else
            % Use FBi since it is more robust
            L_TO(i) = L_FBi(i);
        end
        
    catch err
        err.message;
        L_TO(i) = L_FBi(i);
    end
end

% Finally we round to get final indices
L_TD = round(L_TD);
L_TO = round(L_TO);



%% RIGHT FOOT EVENTS
% the same steps we just took for the left side

% Begin by creating a gross estimation of foot forwards and foot backs
feature accel off

[B,A] = butter(2, 5/(hz/2), 'low');
filtered_R_heel = filtfilt(B,A,R_heel(:,3));

[~,R_FFi] = findpeaks(-filtered_R_heel(:,1),'minpeakdistance',round(0.35*hz));

if strcmp(label,'walk')
    % To deal with peaks that result from signal flipping, threshold them
    angSig = angles.R_foot(:,3);
    angSig(abs(angSig)>90,:) = nan;
    [~,R_FBi] = findpeaks(-angSig,'minpeakdistance',round(0.7*hz),'minpeakheight',20);
else
    [~,R_FBi] = findpeaks(filtered_R_heel(:,1),'minpeakdistance',round(0.35*hz));
end

feature accel on

%Uncomment block below to enable more aggressive quality control of data

% if prctile(abs(angles.R_foot(:,3)),90) > 120 && vel < 4
%
% error('Right ankle values outside of expected ranges, please ensure your shoe markers are properly placed and redo your collection')
%
% end

% Remove any leading FF and FB
R_FFi(R_FFi<L_FFi(1)) = [];
R_FBi(R_FBi<R_FFi(1)) = [];

%% find largest block of continuous data

%We want to check before and after that there is sufficient data for
%analysis
if size(R_FFi,1) < 2 || size(R_FBi,1) < 2
    error('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
end

[R_FFi, R_FBi, R_block_start, ~]=largest_block(R_FFi, R_FBi);

if size(R_FFi,1) < 2 || size(R_FBi,1) < 2
    error('Automated event detection unable to pull adequate number of strides for analysis. Please redo your data collection.')
end

%In rare instances a the index will be in incorrect order run below again
%in case

% Remove any leading FF and FB
R_FFi(R_FFi<L_FFi(1)) = [];
R_FBi(R_FBi<R_FFi(1)) = [];

%%

% TOUCHDOWN
% evtRtd from above


% SELECT SEQUENTIAL STEPS

% Now create an ordered set of sequential steps using above elements
closest = abs(repmat(R_FFi(:),1,length(evtRtd))-repmat(evtRtd(:)',length(R_FFi),1));
[mindist,minx] = nanmin(closest,[],1);
for i = unique(minx)
    if sum(ismember(minx,i))>1
        mindist(minx==i) = [];
        evtRtd(minx==i) = [];
        minx(minx==i) = [];
    end
end

% Parameter based on the typical frame adjustment observed for 300
% datasets
if strcmp(label,'run')
    maxadj = 0.05*hz;
else
    maxadj = 0.10*hz;
end


R_TD = nan(length(R_FFi),1);
evFrtd = zeros(length(R_FFi),1);

% Here we replace FF indices with evt indices where criteria are met
% default is to use FF
for i = 1:length(R_FFi)
    try
        
        if i > max(minx)
            break
        elseif ismember(i,minx) && mindist(minx==i)<maxadj
            % Replace with evtRtd since its more accurate
            R_TD(i) = evtRtd(minx==i);
            evFrtd(i) = 1;
        else
            % Use FFi since it is more robust
            R_TD(i) = R_FFi(i);
        end
        
    catch err
        err.message;
        R_TD(i) = R_FFi(i);
    end
end


% TAKEOFF

%evtRto from above

% SELECT SEQUENTIAL STEPS

% Now create an ordered set of sequential steps using above elements
closest = abs(repmat(R_FBi(:),1,length(evtRto))-repmat(evtRto(:)',length(R_FBi),1));
[mindist,minx] = nanmin(closest,[],1);
for i = unique(minx)
    if sum(ismember(minx,i))>1
        mindist(minx==i) = [];
        evtRto(minx==i) = [];
        minx(minx==i) = [];
    end
end


%Parameter based on the frame adjustment observed from 300 datasets
maxadj = 0.15*hz;

R_TO = nan(length(R_FBi),1);
evFrto = zeros(length(R_FBi),1);

% Here we replace FB indices with TO from PCA default is to use FB
for i = 1:length(R_FBi)
    try
        
        if i > max(minx)
            break
        elseif ismember(i,minx) && mindist(minx==i)<=maxadj
            % Replace with evtRto since its more accurate
            R_TO(i) = evtRto(minx==i);
            evFrto(i) = 1;
        else
            % Use FBi since it is more robust
            R_TO(i) = R_FBi(i);
        end
        
    catch err
        err.message;
        R_TO(i) = R_FBi(i);
    end
end

% Finally, round to get final indices
R_TD = round(R_TD);
R_TO = round(R_TO);


%% if largest chunk of continuous data not at beginning, chop both right and left so they match


%index must begin with left touchdown and end with right toe
%off

if R_block_start < L_block_start
    %remove all right indices that occur before the first left touchdown
    R_TO((R_TD(:,1)<L_block_start)==1,:)=[];
    R_TD((R_TD(:,1)<L_block_start)==1,:)=[];
end

%end

flag = 0;

if L_block_start < R_block_start
    %remove left touchdowns more than one touchdown before the first right touchdown
    cut_inds = (L_TD(:,1)<R_block_start)==1;
    %this loop ensures the first index will be a left touchdown
    for i = 1:size(cut_inds,1)
        if cut_inds(i) == 1 && cut_inds(i+1) == 0 && flag ==0
            cut_inds(i) = 0;
            flag = 1;
        end
    end
    
    L_TD(cut_inds,:) = [];
    L_TO(cut_inds,:) = [];
    clear cut_inds
end


%% create an events matrix

% Remove trailing nans that may have crept in
evFltd(isnan(L_TD)) = [];
evFlto(isnan(L_TO)) = [];
evFrtd(isnan(R_TD)) = [];
evFrto(isnan(R_TO)) = [];

L_TD(isnan(L_TD)) = [];
L_TO(isnan(L_TO)) = [];
R_TD(isnan(R_TD)) = [];
R_TO(isnan(R_TO)) = [];

% Find the closest ordered pairs of L_TO and R_TD to synchronize steps
closest = abs(repmat(R_TD(:),1,length(L_TO))-repmat(L_TO(:)',length(R_TD),1));
[~,minx] = nanmin(closest,[],1);

% Truncate right stances to match up with left
evFrtd = evFrtd(unique(minx));
R_TD = R_TD(unique(minx));

testlength = min([length(L_TO) length(R_TD)]);
if median(L_TO(1:testlength)-R_TD(1:testlength))<0 % Then there is a flight phase
    
    % Find the closest ordered pairs of R_TD and R_TO to synchronize steps
    closest = abs(repmat(R_TO(:),1,length(R_TD))-repmat(R_TD(:)',length(R_TO),1));
    [~,minx] = nanmin(closest,[],1);
    evFrto = evFrto(unique(minx));
    R_TO = R_TO(unique(minx));
    
else % There is no flight phase i.e. grounded running or walking
    
    % Find the closest ordered pairs of R_TO and L_TD to synchronize steps
    tmp = L_TD(2:end);
    closest = abs(repmat(R_TO(:),1,length(tmp))-repmat(tmp(:)',length(R_TO),1));
    [~,minx] = nanmin(closest,[],1);
    evFrto = evFrto(unique(minx));
    R_TO = R_TO(unique(minx));
    
end


events = [length(L_TD),length(L_TO),length(R_TD),length(R_TO)];

% Chop everything to the same length
L_TD = L_TD(1:min(events));
L_TO = L_TO(1:min(events));
R_TD = R_TD(1:min(events));
R_TO = R_TO(1:min(events));

evFltd = evFltd(1:min(events));
evFlto = evFlto(1:min(events));
evFrtd = evFrtd(1:min(events));
evFrto = evFrto(1:min(events));

% Very rarely, these will wind up empty and assignment doesn't work
if isempty(L_TD) || isempty(L_TO) || isempty(R_TD) || isempty(R_TO)
    % skip
else
    events = nan(min(events),4);
    events(:,1)=L_TD;
    events(:,2)=L_TO;
    events(:,3)=R_TD;
    events(:,4)=R_TO;
end


% Very rarely, these will wind up empty and assignment doesn't work
if isempty(evFltd) || isempty(evFlto) || isempty(evFrtd) || isempty(evFrto)
    % skip
else
    eventsflag(:,1) = evFltd;
    eventsflag(:,2) = evFlto;
    eventsflag(:,3) = evFrtd;
    eventsflag(:,4) = evFrto;
end


% Remove first row since these will very often be reliant on FF and FB
% measures
if size(events,1) > 1
    events(1,:) = [];
    eventsflag(1,:) = [];
end


% Occasionally, one stance will drop out, and data becomes
% discontinuous...this fix alleviates this by trimming data to largest
% continuous block
try
    cont = [events(2:end,1)>events(1:end-1,2) events(2:end,3)>events(1:end-1,4)];  %Touchdowns have to come before Toeoffs for same leg
    F = find(any([0 0;cont;0 0]==0,2));
    D = diff(F)-2;
    [M,L] = max(D);
    events = events(F(L):F(L)+M,:);
    eventsflag = eventsflag(F(L):F(L)+M,:);
catch ME
    disp('Could not obtain a continuous block of events')
    events = [];
    eventsflag = [];
    ME.message;
end


% Worst-case... return to foot forward, foot back detection
if size(events,1) < 5
    disp('Automated event detection failed, defaulting to foot-forward foot-back')
    
    events = [length(L_FFi),length(L_FBi),length(R_FFi),length(R_FBi)];
    
    % Chop everything to the same length
    L_FFi = L_FFi(1:min(events));
    L_FBi = L_FBi(1:min(events));
    R_FFi = R_FFi(1:min(events));
    R_FBi = R_FBi(1:min(events));
    
    events = nan(min(events),4);
    events(:,1)=L_FFi;
    events(:,2)=L_FBi;
    events(:,3)=R_FFi;
    events(:,4)=R_FBi;
    
    eventsflag = zeros(size(events));
    
end


% Pull event columns from events so everything is consistent
L_TD = events(:,1);
L_TO = events(:,2);
R_TD = events(:,3);
R_TO = events(:,4);


%% Normalize the steps 0 to 100 from touchdown to takeoff
% do this by using a cubic spline filling function

% preallocate for speed
norm_ang.L_ankle = zeros(101,length(L_TD),3);
norm_ang.L_knee = norm_ang.L_ankle;
norm_ang.L_hip = norm_ang.L_ankle;
norm_ang.L_foot = norm_ang.L_ankle;
norm_ang.L_pelvis = norm_ang.L_ankle;
norm_vel.L_ankle = norm_ang.L_ankle;
norm_vel.L_knee = norm_ang.L_ankle;
norm_vel.L_hip = norm_ang.L_ankle;
norm_vel.L_pelvis = norm_ang.L_ankle;
norm_pos.L_heel = norm_ang.L_ankle;

norm_ang.R_ankle = zeros(101,length(R_TD),3);
norm_ang.R_knee = norm_ang.R_ankle;
norm_ang.R_hip = norm_ang.R_ankle;
norm_ang.R_foot = norm_ang.R_ankle;
norm_ang.R_pelvis = norm_ang.R_ankle;
norm_vel.R_ankle = norm_ang.R_ankle;
norm_vel.R_knee = norm_ang.R_ankle;
norm_vel.R_hip = norm_ang.R_ankle;
norm_vel.R_pelvis = norm_ang.R_ankle;
norm_pos.R_heel = norm_ang.R_ankle;


for i=1:length(L_TD)
    norm_ang.L_ankle(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),angles.L_ankle(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    norm_ang.L_knee(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),angles.L_knee(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    norm_ang.L_hip(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),angles.L_hip(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    
    norm_ang.L_foot(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),angles.L_foot(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    norm_ang.L_pelvis(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),angles.pelvis(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    
    norm_vel.L_ankle(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),velocities.L_ankle(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    norm_vel.L_knee(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),velocities.L_knee(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    norm_vel.L_hip(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),velocities.L_hip(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    norm_vel.L_pelvis(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),velocities.pelvis(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    
    norm_pos.L_heel(:,i,:) = interp1(0:(L_TO(i)-L_TD(i)),L_heel(L_TD(i):L_TO(i),:),0:(L_TO(i)-L_TD(i))/100:(L_TO(i)-L_TD(i)),'pchip');
    
end
for i=1:length(R_TD)
    norm_ang.R_ankle(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),angles.R_ankle(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    norm_ang.R_knee(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),angles.R_knee(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    norm_ang.R_hip(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),angles.R_hip(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    
    norm_ang.R_foot(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),angles.R_foot(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    norm_ang.R_pelvis(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),angles.pelvis(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    
    norm_vel.R_ankle(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),velocities.R_ankle(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    norm_vel.R_knee(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),velocities.R_knee(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    norm_vel.R_hip(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),velocities.R_hip(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    norm_vel.R_pelvis(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),velocities.pelvis(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
    
    norm_pos.R_heel(:,i,:) = interp1(0:(R_TO(i)-R_TD(i)),R_heel(R_TD(i):R_TO(i),:),0:(R_TO(i)-R_TD(i))/100:(R_TO(i)-R_TD(i)),'pchip');
end

%% in order to identify the heelwhip, we want the foot angle in the global
% projected angle of the long axis of the foot into the floor
% during the swing phase from takeoff to touchdown
L_foot_angle = zeros(101,length(L_TD)-1,3);
for i=1:length(L_TD)-1
    L_foot_angle(:,i,:) = interp1(0:(L_TD(i+1)-L_TO(i)), angles.L_foot(L_TO(i):L_TD(i+1),:), 0:(L_TD(i+1)-L_TO(i))/100:(L_TD(i+1)-L_TO(i)), 'pchip');
end
R_foot_angle = zeros(101,length(R_TD)-1,3);
for i=1:length(R_TD)-1
    R_foot_angle(:,i,:) = interp1(0:(R_TD(i+1)-R_TO(i)), angles.R_foot(R_TO(i):R_TD(i+1),:), 0:(R_TD(i+1)-R_TO(i))/100:(R_TD(i+1)-R_TO(i)), 'pchip');
end

%% use the 'DROP THE BAD' function to seperate the good data

% Specify whether curves are plotted
drop_plot = 1;

[norm_ang.L_ankle,norm_ang.R_ankle]   = drop_the_bad(norm_ang.L_ankle,norm_ang.R_ankle,drop_plot);
[norm_ang.L_knee,norm_ang.R_knee]     = drop_the_bad(norm_ang.L_knee,norm_ang.R_knee,drop_plot);
[norm_ang.L_hip,norm_ang.R_hip]       = drop_the_bad(norm_ang.L_hip,norm_ang.R_hip,drop_plot);
[norm_ang.L_foot,norm_ang.R_foot]     = drop_the_bad(norm_ang.L_foot,norm_ang.R_foot,drop_plot);
[norm_ang.L_pelvis,norm_ang.R_pelvis] = drop_the_bad(norm_ang.L_pelvis,norm_ang.R_pelvis,drop_plot);
[norm_vel.L_ankle,norm_vel.R_ankle]   = drop_the_bad(norm_vel.L_ankle,norm_vel.R_ankle,drop_plot);
[norm_vel.L_knee,norm_vel.R_knee]     = drop_the_bad(norm_vel.L_knee,norm_vel.R_knee,drop_plot);
[norm_vel.L_hip,norm_vel.R_hip]       = drop_the_bad(norm_vel.L_hip,norm_vel.R_hip,drop_plot);
[norm_vel.L_pelvis,norm_vel.R_pelvis] = drop_the_bad(norm_vel.L_pelvis,norm_vel.R_pelvis,drop_plot);
[norm_pos.L_heel,norm_pos.R_heel]     = drop_the_bad(norm_pos.L_heel,norm_pos.R_heel,drop_plot);

close(findobj('tag', 'drop_the_bad_temp_figure'));


%% pick off max and mins for the DISCRETE_VARIABLES automated report file
% MORE ROBUST WITH MEDIANS RATHER THAN MEANS
%
% Note some variables will remain empty. These were used a placeholders 
% incase those variables were determined to be of interest at a later date

DISCRETE_VARIABLES = zeros(77,3);

% 1% 'SIDE'	'Left'	'Right'
% 2% 'STEP WIDTH'
temp_L = sort(norm_pos.L_heel(1,:,1));
temp_R = sort(norm_pos.R_heel(1,:,1));
DISCRETE_VARIABLES(2,2:3) = ( median(temp_R) - median(temp_L) )/1000;

% 3% 'STRIDE RATE'
temp = sort(60*hz./diff(L_TD(:)));
DISCRETE_VARIABLES(3,2)=median(temp);
temp = sort(60*hz./diff(R_TD(:)));
DISCRETE_VARIABLES(3,3)=median(temp);

% 4% 'STRIDE LENGTH'
temp = sort((vel)*diff(L_TD(:))./hz);
DISCRETE_VARIABLES(4,2)=median(temp);
temp = sort((vel)*diff(R_TD(:))./hz);
DISCRETE_VARIABLES(4,3)=median(temp);

% 5% 'SWING TIME'
L_swing_frames = nan(1,length(L_TD)-1);
for i = 1:length(L_TD)-1
    L_swing_frames(i)=L_TD(i+1)-L_TO(i);
end
temp = sort((L_swing_frames)/hz);
DISCRETE_VARIABLES(5,2)=nanmedian(temp);
R_swing_frames = nan(1,length(R_TD)-1);
for i = 1:length(R_TD)-1
    R_swing_frames(i)=R_TD(i+1)-R_TO(i);
end
temp = sort((R_swing_frames)/hz);
DISCRETE_VARIABLES(5,3)=nanmedian(temp);

% 6% 'STANCE TIME'
temp = sort((L_TO-L_TD)/hz);
DISCRETE_VARIABLES(6,2)=nanmedian(temp);
temp = sort((R_TO-R_TD)/hz);
DISCRETE_VARIABLES(6,3)=nanmedian(temp);

% 7% 'PELVIS PEAK DROP ANGLE'
temp = sort(min(norm_ang.L_pelvis(:,:,2)));
DISCRETE_VARIABLES(7,2)= - median(temp);
temp = sort(max(norm_ang.R_pelvis(:,:,2)));
DISCRETE_VARIABLES(7,3)= - median(temp);
% if plots ~=0; figure('tag', 'steps_temp_figure');
%     subplot(211); plot(-norm_ang.L_pelvis(:,:,2),'c');title('LEFT PELVIS PEAK DROP ANGLE')
%     subplot(212); plot(-norm_ang.R_pelvis(:,:,2),'c');title('RIGHT PELVIS PEAK DROP ANGLE')
% for i=1:length(norm_ang.L_pelvis(1,:,2))
%     a = find(norm_ang.L_pelvis(:,i+1,2) == min(norm_ang.L_pelvis(:,i+1,2)));
%     subplot(211); hold on; plot(a,-norm_ang.L_pelvis(a,i+1,2),'b*');
% end
% for i = 1:length(norm_ang.R_pelvis(1,2:end-1,2))
%     a = find(norm_ang.R_pelvis(:,i+1,2) == max(norm_ang.R_pelvis(:,i+1,2)));
%     subplot(212); hold on; plot(a,-norm_ang.R_pelvis(a,i+1,2),'r*');
% end
% end

% 8% 'PELVIS DROP %STANCE'
% 9% 'PELVIS DROP @HS'
% 10% 'PELVIS DROP EXCURSION'
temp = sort((min(norm_ang.L_pelvis(:,:,2)))-(norm_ang.L_pelvis(1,:,2)));
DISCRETE_VARIABLES(10,2)= - median(temp);
temp = sort((max(norm_ang.R_pelvis(:,:,2)))-(norm_ang.R_pelvis(1,:,2)));
DISCRETE_VARIABLES(10,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(-norm_ang.L_pelvis(:,:,2),'c');title('LEFT PELVIS EXCURSION AND MAX DROP ANGLE')
    subplot(212); plot(-norm_ang.R_pelvis(:,:,2),'c');title('RIGHT PELVIS EXCURSION AND MAX DROP ANGLE')
    for i=1:length(norm_ang.L_pelvis(1,:,2))
        a = find(norm_ang.L_pelvis(:,i,2) == min(norm_ang.L_pelvis(:,i,2)));
        subplot(211); hold on; plot(a,-norm_ang.L_pelvis(a,i,2),'b*');
        subplot(211); hold on; plot(1,-norm_ang.L_pelvis(1,i,2),'b*');
    end
    for i = 1:length(norm_ang.R_pelvis(1,2:end-1,2))
        a = find(norm_ang.R_pelvis(:,i,2) == max(norm_ang.R_pelvis(:,i,2)));
        subplot(212); hold on; plot(a,-norm_ang.R_pelvis(a,i,2),'r*');
        subplot(212); hold on; plot(1,-norm_ang.R_pelvis(1,i,2),'r*');
    end
end

% 11% 'ANKLE DF PEAK ANGLE'
temp = sort(min(norm_ang.L_ankle(:,:,3)));
DISCRETE_VARIABLES(11,2)= - median(temp);
temp = sort(min(norm_ang.R_ankle(:,:,3)));
DISCRETE_VARIABLES(11,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(-norm_ang.L_ankle(:,:,3),'c');title('LEFT ANKLE DORSIFLEXION PEAK ANGLE')
    subplot(212); plot(-norm_ang.R_ankle(:,:,3),'c');title('RIGHT ANKLE DORSIFLEXION PEAK ANGLE')
    for i=1:length(norm_ang.L_ankle(1,:,3))
        a = find(norm_ang.L_ankle(:,i,3) == min(norm_ang.L_ankle(:,i,3)));
        subplot(211); hold on; plot(a,-norm_ang.L_ankle(a,i,3),'b*');
    end
    for i=1:length(norm_ang.R_ankle(1,:,3))
        a = find(norm_ang.R_ankle(:,i,3) == min(norm_ang.R_ankle(:,i,3)));
        subplot(212); hold on; plot(a,-norm_ang.R_ankle(a,i,3),'r*');
    end
end

% 12% 'ANKLE DF %STANCE'
% 13% 'ANKLE DF @HS'
% 14% 'ANKLE DF EXCURSION'
% 15% 'ANKLE EVE PEAK ANGLE'
temp = sort(min(norm_ang.L_ankle(:,:,1)));
DISCRETE_VARIABLES(15,2)= median(temp);
temp = sort(max(norm_ang.R_ankle(:,:,1)));
DISCRETE_VARIABLES(15,3)= - median(temp);

% 16% 'ANKLE EVE %STANCE'
for i=1:length(norm_ang.L_ankle(1,:,1))
    a(i) = find(norm_ang.L_ankle(:,i,1) == min(norm_ang.L_ankle(:,i,1)));
end
temp = sort(a)/100;
DISCRETE_VARIABLES(16,2)= median(temp);
for i=1:length(norm_ang.R_ankle(1,:,1))
    a(i) = find(norm_ang.R_ankle(:,i,1) == max(norm_ang.R_ankle(:,i,1)));
end
temp = sort(a)/100;
DISCRETE_VARIABLES(16,3)= median(temp);

% 17% 'ANKLE EVE @HS'
% 18% 'ANKLE EVE EXCURSION'
temp = sort(min(norm_ang.L_ankle(:,:,1))-(norm_ang.L_ankle(1,:,1)));
DISCRETE_VARIABLES(18,2)= median(temp);
temp = sort(max(norm_ang.R_ankle(:,:,1))-(norm_ang.R_ankle(1,:,1)));
DISCRETE_VARIABLES(18,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_ang.L_ankle(:,:,1),'c');title('LEFT ANKLE EVERSION PEAK ANGLE, TIMING, AND EXCURSION')
    subplot(212); plot(-norm_ang.R_ankle(:,:,1),'c');title('RIGHT ANKLE EVERSION PEAK ANGLE, TIMING AND EXCURSION')
    for i=1:length(norm_ang.L_ankle(1,:,1))
        a = find(norm_ang.L_ankle(:,i,1) == min(norm_ang.L_ankle(:,i,1)));
        subplot(211); hold on; plot(a,norm_ang.L_ankle(a,i,1),'b*');
        subplot(211); hold on; plot(1,norm_ang.L_ankle(1,i,1),'b*');
    end
    for i=1:length(norm_ang.R_ankle(1,:,1))
        a = find(norm_ang.R_ankle(:,i,1) == max(norm_ang.R_ankle(:,i,1)));
        subplot(212); hold on; plot(a,-norm_ang.R_ankle(a,i,1),'r*');
        subplot(212); hold on; plot(1,-norm_ang.R_ankle(1,i,1),'r*');
    end
end

% PRONATION ONSET AND OFFSET
% want to find when peak pronation happens and when resupintation happens
L_pros=nan(1,(length(norm_ang.L_ankle(1,:,1))));
L_sups=L_pros;
L_peaks=L_pros;
for i=1:length(norm_ang.L_ankle(1,:,1)) % for each step
    try
        rom = max(norm_ang.L_ankle(:,i,1)) - min(norm_ang.L_ankle(:,i,1));% eversion range of motion
        a=find(norm_ang.L_ankle(:,i,1) < min(norm_ang.L_ankle(:,i,1))+0.2*rom);
        L_pros(i)=a(1);
        L_sups(i)=a(end);
        L_peaks(i)=find(norm_ang.L_ankle(:,i,1) == min(norm_ang.L_ankle(:,i,1)));
    catch ME
        ME.message;
    end
end

R_pros=nan(1,(length(norm_ang.R_ankle(1,:,1))));
R_sups=R_pros;
R_peaks=R_pros;
for i=1:length(norm_ang.R_ankle(1,:,1)) % for each step
    try
        rom = max(norm_ang.R_ankle(:,i,1)) - min(norm_ang.R_ankle(:,i,1));% eversion range of motion
        a=find(-norm_ang.R_ankle(:,i,1) < min(-norm_ang.R_ankle(:,i,1))+0.2*rom);
        R_pros(i)=a(1);
        R_sups(i)=a(end);
        R_peaks(i)=find(norm_ang.R_ankle(:,i,1) == max(norm_ang.R_ankle(:,i,1)));
    catch ME
        ME.message;
    end
end

DISCRETE_VARIABLES(73,2)= round(median(L_pros));
DISCRETE_VARIABLES(73,3)= round(median(R_pros));
DISCRETE_VARIABLES(74,2)= round(median(L_sups));
DISCRETE_VARIABLES(74,3)= round(median(R_sups));

if plots~=0
    L_pro=round(median(L_pros)); L_sup=round(median(L_sups));
    R_pro=round(median(R_pros)); R_sup=round(median(R_sups));
    x=0:100; y=nanmean(norm_ang.L_ankle(:,:,1),2)'; sd=std(norm_ang.L_ankle(:,:,1),1,2)';
    x2=0:100; y2=-nanmean(norm_ang.R_ankle(:,:,1),2)'; sd2=std(norm_ang.R_ankle(:,:,1),1,2)';
    
    figure('tag', 'steps_temp_figure'); hold on;
    
    fill([x,flip(x,2)],[(y+sd),flip((y-sd),2)], ...
        [6 6 7]/8, 'EdgeColor', 'none', 'facealpha', 0.2);
    fill([L_pro-1:L_sup-1,flip(L_pro-1:L_sup-1,2)],...
        [y(L_pro:L_sup)+sd(L_pro:L_sup),flip((y(L_pro:L_sup)-sd(L_pro:L_sup)),2)], ...
        [3 3 7]/8, 'EdgeColor', 'none', 'facealpha', 0.2);
    plot(x,y,'b','LineWidth',1);
    plot(L_pro-1:L_sup-1,nanmean(norm_ang.L_ankle(L_pro:L_sup,:,1),2)','b','LineWidth',2);
    
    fill([x2,flip(x2,2)],[(y2+sd2),flip((y2-sd2),2)], ...
        [7 6 6]/8, 'EdgeColor', 'none', 'facealpha', 0.2);
    fill([R_pro-1:R_sup-1,flip(R_pro-1:R_sup-1,2)],...
        [y2(R_pro:R_sup)+sd2(R_pro:R_sup),flip((y2(R_pro:R_sup)-sd2(R_pro:R_sup)),2)], ...
        [7 3 3]/8, 'EdgeColor', 'none', 'facealpha', 0.2);
    plot(x,y2,'r','LineWidth',1);
    plot(R_pro-1:R_sup-1,nanmean(-norm_ang.R_ankle(R_pro:R_sup,:,1),2),'r','LineWidth',2);
end

% 19% 'ANKLE ROT PEAK ANGLE'
temp = sort(min(norm_ang.L_ankle(:,:,2)));
DISCRETE_VARIABLES(19,2)= median(temp);
temp = sort(max(norm_ang.R_ankle(:,:,2)));
DISCRETE_VARIABLES(19,3)= - median(temp);

% 20% 'ANKLE ROT %STANCE'
% 21% 'ANKLE ROT @HS'
% 22% 'ANKLE ROT EXCURSION'
temp = sort(min(norm_ang.L_ankle(:,:,2))-(norm_ang.L_ankle(1,:,2)));
DISCRETE_VARIABLES(22,2)= median(temp);
temp = sort(max(norm_ang.R_ankle(:,:,2))-(norm_ang.R_ankle(1,:,2)));
DISCRETE_VARIABLES(22,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_ang.L_ankle(:,:,2),'c');title('LEFT ANKLE ROTATION PEAK ANGLE AND EXCURSION')
    subplot(212); plot(-norm_ang.R_ankle(:,:,2),'c');title('RIGHT ANKLE ROTATION PEAK ANGLE AND EXCURSION')
    for i=1:length(norm_ang.L_ankle(1,:,2))
        a = find(norm_ang.L_ankle(:,i,2) == min(norm_ang.L_ankle(:,i,2)));
        subplot(211); hold on; plot(a,norm_ang.L_ankle(a,i,2),'b*');
        subplot(211); hold on; plot(1,norm_ang.L_ankle(1,i,2),'b*');
    end
    for i=1:length(norm_ang.R_ankle(1,:,2))
        a = find(norm_ang.R_ankle(:,i,2) == max(norm_ang.R_ankle(:,i,2)));
        subplot(212); hold on; plot(a,-norm_ang.R_ankle(a,i,2),'r*');
        subplot(212); hold on; plot(1,-norm_ang.R_ankle(1,i,2),'r*');
    end
end

% 23% 'KNEE FLEX PEAK ANGLE'
temp = sort(max(norm_ang.L_knee(1:70,:,3)));
DISCRETE_VARIABLES(23,2)= - median(temp);
temp = sort(max(norm_ang.R_knee(1:70,:,3)));
DISCRETE_VARIABLES(23,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); hold on; plot(-norm_ang.L_knee(:,:,3),'y'); plot(-norm_ang.L_knee(1:70,:,3),'c');title('LEFT KNEE FLEXION PEAK ANGLE')
    subplot(212); hold on; plot(-norm_ang.R_knee(:,:,3),'y'); plot(-norm_ang.R_knee(1:70,:,3),'c');title('RIGHT KNEE FLEXION PEAK ANGLE')
    for i=1:length(norm_ang.L_knee(1,:,3))
        a = find(norm_ang.L_knee(:,i,3) == max(norm_ang.L_knee(1:70,i,3)));
        subplot(211); hold on; plot(a,-norm_ang.L_knee(a,i,3),'b*');
    end
    for i=1:length(norm_ang.R_knee(1,:,3))
        a = find(norm_ang.R_knee(:,i,3) == max(norm_ang.R_knee(1:70,i,3)));
        subplot(212); hold on; plot(a,-norm_ang.R_knee(a,i,3),'r*');
    end
end

% 24% 'KNEE FLEX %STANCE'
% 25% 'KNEE FLEX @HS'
% 26% 'KNEE FLEX EXCURSION'
% 27% 'KNEE ADD PEAK ANGLE'
temp = sort(max(norm_ang.L_knee(:,:,1)));
DISCRETE_VARIABLES(27,2)= median(temp);
temp = sort(min(norm_ang.R_knee(:,:,1)));
DISCRETE_VARIABLES(27,3)= - median(temp);

% 28% 'KNEE ADD %STANCE'
% 29% 'KNEE ADD @HS'
% 30% 'KNEE ADD EXCURSION'
temp = sort(max(norm_ang.L_knee(:,:,1)) - (norm_ang.L_knee(1,:,1)));
DISCRETE_VARIABLES(30,2)= median(temp);
temp = sort(min(norm_ang.R_knee(:,:,1)) - (norm_ang.R_knee(1,:,1)));
DISCRETE_VARIABLES(30,3)= - median(temp);

% 31% 'KNEE ABD PEAK ANGLE'
temp = sort(min(norm_ang.L_knee(:,:,1)));
DISCRETE_VARIABLES(31,2)= median(temp);
temp = sort(max(norm_ang.R_knee(:,:,1)));
DISCRETE_VARIABLES(31,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_ang.L_knee(:,:,1),'c');title('LEFT KNEE ADDUCTION & ABDUCTION PEAK ANGLES AND EXCURSIONS')
    subplot(212); plot(-norm_ang.R_knee(:,:,1),'c');title('RIGHT KNEE ADDUCTION & ABDUCTION PEAK ANGLES AND EXCURSIONS')
    for i=1:length(norm_ang.L_knee(1,:,1))
        a = find(norm_ang.L_knee(:,i,1) == min(norm_ang.L_knee(:,i,1)));
        subplot(211); hold on; plot(a,norm_ang.L_knee(a,i,1),'b*');
        a = find(norm_ang.L_knee(:,i,1) == max(norm_ang.L_knee(:,i,1)));
        subplot(211); hold on; plot(a,norm_ang.L_knee(a,i,1),'b*');
    end
    for i=1:length(norm_ang.R_knee(1,:,1))
        a = find(norm_ang.R_knee(:,i,1) == min(norm_ang.R_knee(:,i,1)));
        subplot(212); hold on; plot(a,-norm_ang.R_knee(a,i,1),'r*');
        a = find(norm_ang.R_knee(:,i,1) == max(norm_ang.R_knee(:,i,1)));
        subplot(212); hold on; plot(a,-norm_ang.R_knee(a,i,1),'r*');
    end
end

% 32% 'KNEE ABD %STANCE'
% 33% 'KNEE ABD @HS'
% 34% 'KNEE ABD EXCURSION'
temp = sort(min(norm_ang.L_knee(:,:,1)) - (norm_ang.L_knee(1,:,1)));
DISCRETE_VARIABLES(34,2)= median(temp);
temp = sort(max(norm_ang.R_knee(:,:,1)) - (norm_ang.R_knee(1,:,1)));
DISCRETE_VARIABLES(34,3)= - median(temp);

% 35% 'KNEE ROT PEAK ANGLE'
temp = sort(max(norm_ang.L_knee(:,:,2)));
DISCRETE_VARIABLES(35,2) = median(temp(1,:));
temp = sort(min(norm_ang.R_knee(:,:,2)));
DISCRETE_VARIABLES(35,3) = - median(temp(1,:));

% 36% 'KNEE ROT %STANCE'
% 37% 'KNEE ROT @HS'
% 38% 'KNEE ROT EXCURSION'
temp = sort(max(norm_ang.L_knee(:,:,2)) - norm_ang.L_knee(1,:,2));
DISCRETE_VARIABLES(38,2) = median(temp(1,:));
temp = sort(min(norm_ang.R_knee(:,:,2)) - norm_ang.R_knee(1,:,2));
DISCRETE_VARIABLES(38,3) = - median(temp(1,:));
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_ang.L_knee(:,:,2),'c');title('LEFT KNEE ROTATION PEAK AND EXCURSION')
    subplot(212); plot(-norm_ang.R_knee(:,:,2),'c');title('RIGHT KNEE ROTATION PEAK AND EXCURSION')
    for i=1:length(norm_ang.L_knee(1,:,2))
        a = find(norm_ang.L_knee(:,i,2) == max(norm_ang.L_knee(:,i,2)));
        subplot(211); hold on; plot(a,norm_ang.L_knee(a,i,2),'b*');
        subplot(211); hold on; plot(norm_ang.L_knee(1,i,2),'b*');
    end
    for i=1:length(norm_ang.R_knee(1,:,2))
        a = find(norm_ang.R_knee(:,i,2) == min(norm_ang.R_knee(:,i,2)));
        subplot(212); hold on; plot(a,-norm_ang.R_knee(a,i,2),'r*');
        subplot(212); hold on; plot(-norm_ang.R_knee(1,i,2),'r*');
    end
end

% 39% 'HIP EXT PEAK ANGLE'
temp = sort(max(norm_ang.L_hip(:,:,3)));
DISCRETE_VARIABLES(39,2)= - median(temp);
temp = sort(max(norm_ang.R_hip(:,:,3)));
DISCRETE_VARIABLES(39,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(-norm_ang.L_hip(:,:,3),'c');title('LEFT HIP EXTENSION PEAK ANGLE')
    subplot(212); plot(-norm_ang.R_hip(:,:,3),'c');title('RIGHT HIP EXTENSION PEAK ANGLE')
    for i=1:length(norm_ang.L_hip(1,:,3))
        a = find(norm_ang.L_hip(:,i,3) == max(norm_ang.L_hip(:,i,3)));
        subplot(211); hold on; plot(a,-norm_ang.L_hip(a,i,3),'b*');
        b=[0,0;100,0];line(b(:,1),b(:,2),'color','k')
    end
    for i=1:length(norm_ang.R_hip(1,:,3))
        a = find(norm_ang.R_hip(:,i,3) == max(norm_ang.R_hip(:,i,3)));
        subplot(212); hold on; plot(a,-norm_ang.R_hip(a,i,3),'r*');
        b=[0,0;100,0];line(b(:,1),b(:,2),'color','k')
    end
end

% 40% 'HIP EXT %STANCE'
% 41% 'HIP EXT @HS'
% 42% 'HIP EXT EXCURSION'
% 43% 'HIP ADD PEAK ANGLE'
temp = sort(max(norm_ang.L_hip(1:60,:,1)));
DISCRETE_VARIABLES(43,2)= median(temp);
temp = sort(min(norm_ang.R_hip(1:60,:,1)));
DISCRETE_VARIABLES(43,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); hold on; plot(norm_ang.L_hip(:,:,1),'y'); plot(norm_ang.L_hip(1:60,:,1),'c');title('LEFT HIP ADDUCTION PEAK ANGLE AND EXCURSION')
    subplot(212); hold on; plot(norm_ang.R_hip(:,:,1),'y'); plot(norm_ang.R_hip(1:60,:,1),'c');title('RIGHT HIP ADDUCTION PEAK ANGLE AND EXCURSION')
    for i=1:length(norm_ang.L_hip(1,:,1))
        a = find(norm_ang.L_hip(:,i,1) == max(norm_ang.L_hip(1:60,i,1)));
        subplot(211); hold on; plot(a,norm_ang.L_hip(a,i,1),'b*');
        subplot(211); hold on; plot(1,norm_ang.L_hip(1,i,1),'b*');
    end
    for i=1:length(norm_ang.R_hip(1,:,1))
        a = find(norm_ang.R_hip(:,i,1) == min(norm_ang.R_hip(1:60,i,1)));
        subplot(212); hold on; plot(a,norm_ang.R_hip(a,i,1),'r*');
        subplot(212); hold on; plot(1,norm_ang.R_hip(1,i,1),'r*');
    end
end

% 44% 'HIP ADD %STANCE'
% 45% 'HIP ADD @HS'
% 46% 'HIP ADD EXCURSION'
temp = sort(max(norm_ang.L_hip(1:60,:,1)) - (norm_ang.L_hip(1,:,1)));
DISCRETE_VARIABLES(46,2)= median(temp);
temp = sort(min(norm_ang.R_hip(1:60,:,1)) - (norm_ang.R_hip(1,:,1)));
DISCRETE_VARIABLES(46,3)= - median(temp);

% 47% 'HIP ROT PEAK ANGLE'
temp = sort(max(norm_ang.L_hip(:,:,2)));
DISCRETE_VARIABLES(47,2)= median(temp);
temp = sort(min(norm_ang.R_hip(:,:,2)));
DISCRETE_VARIABLES(47,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_ang.L_hip(:,:,2),'c');title('LEFT HIP ROTATION PEAK ANGLE AND EXCURSION')
    subplot(212); plot(-norm_ang.R_hip(:,:,2),'c');title('RIGHT HIP ROTATION PEAK ANGLE AND EXCURSION')
    for i=1:length(norm_ang.L_hip(1,:,2))
        a = find(norm_ang.L_hip(:,i,2) == max(norm_ang.L_hip(:,i,2)));
        subplot(211); hold on; plot(a,norm_ang.L_hip(a,i,2),'b*');
    end
    for i=1:length(norm_ang.R_hip(1,:,2))
        a = find(norm_ang.R_hip(:,i,2) == min(norm_ang.R_hip(:,i,2)));
        subplot(212); hold on; plot(a,-norm_ang.R_hip(a,i,2),'r*');
    end
end

% 48% 'HIP ROT %STANCE'
% 49% 'HIP ROT @HS'
% 50% 'HIP ROT EXCURSION'
temp = sort(max(norm_ang.L_hip(:,:,2)) - (norm_ang.L_hip(1,:,2)));
DISCRETE_VARIABLES(50,2)= median(temp);
temp = sort(min(norm_ang.R_hip(:,:,2)) - (norm_ang.R_hip(1,:,2)));
DISCRETE_VARIABLES(50,3)= - median(temp);

% 51% 'FOOT PROG ANGLE'
temp = sort(mean(norm_ang.L_foot(20:50,:,2)));
DISCRETE_VARIABLES(51,2) = - median(temp(1,:));
temp = sort(mean(norm_ang.R_foot(20:50,:,2)));
DISCRETE_VARIABLES(51,3) = - median(temp(1,:));
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); hold on; plot(norm_ang.L_foot(:,:,2),'y'); plot(20:50,norm_ang.L_foot(20:50,:,2),'c');title('LEFT FOOT MEAN PROGRESSION ANGLE')
    subplot(212); hold on; plot(-norm_ang.R_foot(:,:,2),'y'); plot(20:50,-norm_ang.R_foot(20:50,:,2),'c');title('RIGHT FOOT MEAN PROGRESSION ANGLE')
    
    a=[20,50];b=[-DISCRETE_VARIABLES(51,2),-DISCRETE_VARIABLES(51,2)];c=[DISCRETE_VARIABLES(51,3),DISCRETE_VARIABLES(51,3)];
    subplot(211); hold on; line(a,b,'linewidth',3);
    subplot(212); hold on; line(a,c,'linewidth',3,'color','r');
end

% 52% 'FOOT ANG @HS'
DISCRETE_VARIABLES(52,2) = median(norm_ang.L_foot(1,:,3));
DISCRETE_VARIABLES(52,3) = median(norm_ang.R_foot(1,:,3));
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); hold on; plot(norm_ang.L_foot(:,:,3),'c'); title('LEFT FOOT TOUCHDOWN ANGLE')
    subplot(212); hold on; plot(norm_ang.R_foot(:,:,3),'c'); title('RIGHT FOOT TOUCHDOWN ANGLE')
    for i=1:length(norm_ang.L_foot(1,:,1))
        subplot(211);hold on; plot(1,norm_ang.L_foot(1,i,3),'b*');
    end
    for i=1:length(norm_ang.R_foot(1,:,1))
        subplot(212); hold on; plot(1,norm_ang.R_foot(1,i,3),'r*');
    end
end

% 53% 'FOOT ANG @TO'
% 54% 'MED HEEL WHIP PEAK'
% 55% 'MHW %SWING'
% 56% 'MHW EXC FROM TO'
% find the index of when one foot passes the other ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_pass = nan(1,length(L_TD)-1); R_pass = nan(1,length(R_TD)-1); % preallocate
for i=1:length(L_TD)-1
    L_pass(i)= L_TO(i)-1 + find(L_heel(L_TO(i):L_TD(i+1),3) < R_heel(L_TO(i):L_TD(i+1),3),1);
end
for i=1:length(R_TD)-1
    R_pass(i)= R_TO(i)-1 + find(R_heel(R_TO(i):R_TD(i+1),3) < L_heel(R_TO(i):R_TD(i+1),3),1);
end
L_temp = nan(1,length(L_pass)); % preallocate
for i=1:length(L_pass)
    L_temp(i) = sort( max(angles.L_foot(L_TO(i):L_pass(i),2)) - angles.L_foot(L_TO(i),2));
end
DISCRETE_VARIABLES(56,2) = - nanmedian(L_temp);
R_temp = nan(1,length(R_pass)); % preallocate
for i=1:length(R_pass)
    R_temp(i) = sort( max(angles.R_foot(R_TO(i):R_pass(i),2)) - angles.R_foot(R_TO(i),2));
end
DISCRETE_VARIABLES(56,3) = - nanmedian(R_temp);

if plots ~= 0; figure('tag', 'steps_temp_figure')
    subplot(211); hold on; title('LEFT MAX HEEL WHIP EXCERSION FROM TOE OFF')
    for i=1:length(L_pass)
        plot(-angles.L_foot(L_TO(i):L_TD(i+1),2),'y');
        plot(-angles.L_foot(L_TO(i):L_pass(i),2),'c');
        plot(-angles.L_foot(L_TO(i),2),'b*')
        a = find(angles.L_foot(L_TO(i):L_pass(i),2) == max(angles.L_foot(L_TO(i):L_pass(i),2)));
        plot(a,-angles.L_foot(L_TO(i)+a,2),'b*')
    end
    subplot(212); hold on; title('RIGHT MAX HEEL WHIP EXCERSION FROM TOE OFF')
    for i=1:length(R_pass)
        plot(-angles.R_foot(R_TO(i):R_TD(i+1),2),'y');
        plot(-angles.R_foot(R_TO(i):R_pass(i),2),'c');
        plot(-angles.R_foot(R_TO(i),2),'r*')
        a = find(angles.R_foot(R_TO(i):R_pass(i),2) == max(angles.R_foot(R_TO(i):R_pass(i),2)));
        plot(a,-angles.R_foot(R_TO(i)+a,2),'r*')
    end
end


% we want to know the index of the maximal heel whip as well ... create EVENT
% event is a matrix of events used for making the 2d diagrams
% [L_TD, L_Midstance, L_TO, L_heelwhip, R_TD, R_Midstance, R_TO, R_heelwhip]
event=nan(length(events(:,1)),8);
event(:,[1,3,5,7]) = events;
event(:,2)=round((event(:,1)+event(:,3))/2); % midstance = halfway between TD and TO
event(:,6)=round((event(:,5)+event(:,7))/2);% midstance = halfway between TD and TO
% Find index of max heel whip
for i=1:length(L_pass)
    a = find(angles.L_foot(L_TO(i):L_pass(i),2) == max(angles.L_foot(L_TO(i):L_pass(i),2)));
    event(i,4)=L_TO(i)-1+a;
end
for i=1:length(R_pass)
    a = find(angles.R_foot(R_TO(i):R_pass(i),2) == max(angles.R_foot(R_TO(i):R_pass(i),2)));
    event(i,8)=R_TO(i)-1+a;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 57% 'ANKLE DF PEAK VEL'
% 58% 'ANKLE DF VEL %STANCE'
% 59% 'ANKLE EVE PEAK VEL'
temp = sort(min(norm_vel.L_ankle(:,:,1)));
DISCRETE_VARIABLES(59,2)= median(temp);
temp = sort(max(norm_vel.R_ankle(:,:,1)));
DISCRETE_VARIABLES(59,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_vel.L_ankle(:,:,1),'c');title('LEFT ANKLE EVERSION PEAK VELOCITY')
    subplot(212); plot(-norm_vel.R_ankle(:,:,1),'c');title('RIGHT ANKLE EVERSION PEAK VELOCITY')
    for i=1:length(norm_vel.L_ankle(1,:,1))
        a = find(norm_vel.L_ankle(:,i,1) == min(norm_vel.L_ankle(:,i,1)));
        subplot(211); hold on; plot(a,norm_vel.L_ankle(a,i,1),'b*');
    end
    for i=1:length(norm_vel.R_ankle(1,:,1))
        a = find(norm_vel.R_ankle(:,i,1) == max(norm_vel.R_ankle(:,i,1)));
        subplot(212); hold on; plot(a,-norm_vel.R_ankle(a,i,1),'r*');
    end
end

% 60% 'ANKLE EVE VEL %STANCE'
% 61% 'ANKLE ROT PEAK VEL'
temp = sort(min(norm_vel.L_ankle(:,:,2)));
DISCRETE_VARIABLES(61,2)= median(temp);
temp = sort(max(norm_vel.R_ankle(:,:,2)));
DISCRETE_VARIABLES(61,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_vel.L_ankle(:,:,2),'c');title('LEFT ANKLE ROTATION PEAK VELOCITY')
    subplot(212); plot(-norm_vel.R_ankle(:,:,2),'c');title('RIGHT ANKLE ROTATION PEAK VELOCITY')
    for i=1:length(norm_vel.L_ankle(1,:,2))
        a = find(norm_vel.L_ankle(:,i,2) == min(norm_vel.L_ankle(:,i,2)));
        subplot(211); hold on; plot(a,norm_vel.L_ankle(a,i,2),'b*');
    end
    for i=1:length(norm_vel.R_ankle(1,:,2))
        a = find(norm_vel.R_ankle(:,i,2) == max(norm_vel.R_ankle(:,i,2)));
        subplot(212); hold on; plot(a,-norm_vel.R_ankle(a,i,2),'r*');
    end
end

% 62% 'ANKLE ROT VEL %STANCE'
% 63% 'KNEE FLEX PEAK VEL'
% 64% 'KNEE FLEX VEL %STANCE'

% 65% 'KNEE ABD PEAK VEL'

temp = sort(min(norm_vel.L_knee(:,:,1)));
DISCRETE_VARIABLES(65,2)= median(temp);
temp = sort(max(norm_vel.R_knee(:,:,1)));
DISCRETE_VARIABLES(65,3)= - median(temp);

% 66% 'KNEE ABD VEL %STANCE'
% 67% 'KNEE ADD PEAK VEL'
temp = sort(max(norm_vel.L_knee(:,:,1)));
DISCRETE_VARIABLES(67,2)= median(temp);
temp = sort(min(norm_vel.R_knee(:,:,1)));
DISCRETE_VARIABLES(67,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_vel.L_knee(:,:,1),'c');title('LEFT KNEE ADDUCTION & ABDUCTION PEAK VELOCITY')
    subplot(212); plot(-norm_vel.R_knee(:,:,1),'c');title('RIGHT KNEE ADDUCTION & ABDUCTION PEAK VELOCITY')
    for i=1:length(norm_vel.L_knee(1,:,1))
        a = find(norm_vel.L_knee(:,i,1) == min(norm_vel.L_knee(:,i,1)));
        subplot(211); hold on; plot(a,norm_vel.L_knee(a,i,1),'b*');
        a = find(norm_vel.L_knee(:,i,1) == max(norm_vel.L_knee(:,i,1)));
        subplot(211); hold on; plot(a,norm_vel.L_knee(a,i,1),'b*');
    end
    for i=1:length(norm_vel.R_knee(1,:,1))
        a = find(norm_vel.R_knee(:,i,1) == min(norm_vel.R_knee(:,i,1)));
        subplot(212); hold on; plot(a,-norm_vel.R_knee(a,i,1),'r*');
        a = find(norm_vel.R_knee(:,i,1) == max(norm_vel.R_knee(:,i,1)));
        subplot(212); hold on; plot(a,-norm_vel.R_knee(a,i,1),'r*');
    end
end

% 68% 'KNEE ADD VEL %STANCE'
% 69% 'HIP ABD PEAK VEL'
% hip aBduction
temp = sort(min(norm_vel.L_hip(:,:,1)));
DISCRETE_VARIABLES(69,2)= median(temp);
temp = sort(max(norm_vel.R_hip(:,:,1)));
DISCRETE_VARIABLES(69,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_vel.L_hip(:,:,1),'c');title('LEFT HIP ABDUCTION PEAK VELOCITY')
    subplot(212); plot(-norm_vel.R_hip(:,:,1),'c');title('RIGHT HIP ABDUCTION PEAK VELOCITY')
    for i=1:length(norm_vel.L_hip(1,:,1))
        a = find(norm_vel.L_hip(:,i,1) == min(norm_vel.L_hip(:,i,1)));
        subplot(211); hold on; plot(a,norm_vel.L_hip(a,i,1),'b*');
    end
    for i=1:length(norm_vel.R_hip(1,:,1))
        a = find(norm_vel.R_hip(:,i,1) == max(norm_vel.R_hip(:,i,1)));
        subplot(212); hold on; plot(a,-norm_vel.R_hip(a,i,1),'r*');
    end
end


%75 hip aDduction
%L and R should be POSITIVE

temp = sort(max(norm_vel.L_hip(:,:,1)));
DISCRETE_VARIABLES(75,2)= median(temp);
temp = sort(min(norm_vel.R_hip(:,:,1)));
DISCRETE_VARIABLES(75,3)= -median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_vel.L_hip(:,:,1),'c');title('LEFT HIP ADDUCTION PEAK VELOCITY')
    subplot(212); plot(norm_vel.R_hip(:,:,1),'c');title('RIGHT HIP ADDUCTION PEAK VELOCITY')
    for i=1:length(norm_vel.L_hip(1,:,1))
        a = find(norm_vel.L_hip(:,i,1) == max(norm_vel.L_hip(:,i,1)));
        subplot(211); hold on; plot(a,norm_vel.L_hip(a,i,1),'b*');
    end
    for i=1:length(norm_vel.R_hip(1,:,1))
        a = find(norm_vel.R_hip(:,i,1) == min(norm_vel.R_hip(:,i,1)));
        subplot(212); hold on; plot(a,norm_vel.R_hip(a,i,1),'r*');
    end
end

% 70% 'HIP ABD VEL %STANCE'

% 71 - knee rotation velocity
temp = sort(max(norm_vel.L_knee(:,:,2)));
DISCRETE_VARIABLES(71,2)= median(temp);
temp = sort(min(norm_vel.R_knee(:,:,2)));
DISCRETE_VARIABLES(71,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_vel.L_knee(:,:,2),'c');title('LEFT KNEE ROTATION PEAK VELOCITY')
    subplot(212); plot(-norm_vel.R_knee(:,:,2),'c');title('RIGHT KNEE ROTATION PEAK VELOCITY')
    for i=1:length(norm_vel.L_knee(1,:,2))
        a = find(norm_vel.L_knee(:,i,2) == max(norm_vel.L_knee(:,i,2)));
        subplot(211); hold on; plot(a,norm_vel.L_knee(a,i,2),'b*');
    end
    for i=1:length(norm_vel.R_knee(1,:,2))
        a = find(norm_vel.R_knee(:,i,2) == min(norm_vel.R_knee(:,i,2)));
        subplot(212); hold on; plot(a,-norm_vel.R_knee(a,i,2),'r*');
    end
end
% 72 - hip rotation velocity
temp = sort(max(norm_vel.L_hip(:,:,2)));
DISCRETE_VARIABLES(72,2)= median(temp);
temp = sort(min(norm_vel.R_hip(:,:,2)));
DISCRETE_VARIABLES(72,3)= - median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(norm_vel.L_hip(:,:,2),'c');title('LEFT HIP ROTATION PEAK VELOCITY')
    subplot(212); plot(-norm_vel.R_hip(:,:,2),'c');title('RIGHT HIP ROTATION PEAK VELOCITY')
    for i=1:length(norm_vel.L_hip(1,:,2))
        a = find(norm_vel.L_hip(:,i,2) == max(norm_vel.L_hip(:,i,2)));
        subplot(211); hold on; plot(a,norm_vel.L_hip(a,i,2),'b*');
    end
    for i=1:length(norm_vel.R_hip(1,:,2))
        a = find(norm_vel.R_hip(:,i,2) == min(norm_vel.R_hip(:,i,2)));
        subplot(212); hold on; plot(a,-norm_vel.R_hip(a,i,2),'r*');
    end
end

% 76 pelvic drop velocity
% these should be NEGATIVE

temp = sort(min(norm_vel.L_pelvis(:,:,2)));
DISCRETE_VARIABLES(76,2)= median(temp);
temp = sort(max(norm_vel.R_pelvis(:,:,2)));
DISCRETE_VARIABLES(76,3)= -median(temp);
if plots ~=0; figure('tag', 'steps_temp_figure');
    subplot(211); plot(-norm_vel.L_pelvis(:,:,2),'c');title('LEFT PELVIC DROP PEAK VELOCITY')
    subplot(212); plot(-norm_vel.R_pelvis(:,:,2),'c');title('RIGHT PELVIC DROP PEAK VELOCITY')
    for i=1:length(norm_vel.L_pelvis(1,:,2))
        a = find(-norm_vel.L_pelvis(:,i,2) == max(-norm_vel.L_pelvis(:,i,2)));
        subplot(211); hold on; plot(a,-norm_vel.L_pelvis(a,i,2),'b*');
    end
    for i=1:length(norm_vel.R_pelvis(1,:,2))
        a = find(norm_vel.R_pelvis(:,i,2) == max(norm_vel.R_pelvis(:,i,2)));
        subplot(212); hold on; plot(a,-norm_vel.R_pelvis(a,i,2),'r*');
    end
end

%% Vertical Oscillation

filtered_pelvis = filtfilt(B,A,dynamic.pelvis_4);
filtered_L_foot = filtfilt(B,A,dynamic.L_foot_4);
filtered_R_foot = filtfilt(B,A,dynamic.R_foot_4);

try
    
    [vertical_oscillation] = oscillation(filtered_pelvis, L_TD, L_TO, R_TD, R_TO, plots, label);
    
    DISCRETE_VARIABLES(77,2) = median(vertical_oscillation(:,2));
    DISCRETE_VARIABLES(77,3) = median(vertical_oscillation(:,4));
catch ME
    
    error('Error calculating vertical oscillation');
    
end
%% plot the curves that the discrete variables come from


% %% display the outputs
% % if plots ~=0
% disp('VARIABLE                LEFT         RIGHT')
% disp(['STEP WIDTH               ' num2str(DISCRETE_VARIABLES(2,2)) '       ' num2str(DISCRETE_VARIABLES(2,3))])
% disp(['STRIDE RATE              ' num2str(DISCRETE_VARIABLES(3,2)) '        ' num2str(DISCRETE_VARIABLES(3,3))])
% disp(['STRIDE LENGTH            ' num2str(DISCRETE_VARIABLES(4,2)) '         ' num2str(DISCRETE_VARIABLES(4,3))])
% disp(['PELVIS PEAK DROP ANGLE   ' num2str(DISCRETE_VARIABLES(7,2)) '         ' num2str(DISCRETE_VARIABLES(7,3))])
% disp(['PELVIS DROP EXCURSION    ' num2str(DISCRETE_VARIABLES(10,2)) '        ' num2str(DISCRETE_VARIABLES(10,3))])
% disp(['ANKLE DF PEAK ANGLE      ' num2str(DISCRETE_VARIABLES(11,2)) '       ' num2str(DISCRETE_VARIABLES(11,3))])
% disp(['ANKLE EVE PEAK ANGLE     ' num2str(DISCRETE_VARIABLES(15,2)) '       ' num2str(DISCRETE_VARIABLES(15,3))])
% disp(['ANKLE EVE %STANCE        ' num2str(DISCRETE_VARIABLES(16,2)) '          ' num2str(DISCRETE_VARIABLES(16,3))])
% disp(['ANKLE EVE EXCURSION      ' num2str(DISCRETE_VARIABLES(18,2)) '       ' num2str(DISCRETE_VARIABLES(18,3))])
% disp(['ANKLE ROT PEAK ANGLE     ' num2str(DISCRETE_VARIABLES(19,2)) '       ' num2str(DISCRETE_VARIABLES(19,3))])
% disp(['ANKLE ROT EXCURSION      ' num2str(DISCRETE_VARIABLES(22,2)) '       ' num2str(DISCRETE_VARIABLES(22,3))])
% disp(['KNEE FLEX PEAK ANGLE     ' num2str(DISCRETE_VARIABLES(23,2)) '      ' num2str(DISCRETE_VARIABLES(23,3))])
% disp(['KNEE ADD PEAK ANGLE      ' num2str(DISCRETE_VARIABLES(27,2)) '       ' num2str(DISCRETE_VARIABLES(27,3))])
% disp(['KNEE ABD PEAK ANGLE      ' num2str(DISCRETE_VARIABLES(31,2)) '      ' num2str(DISCRETE_VARIABLES(31,3))])
% disp(['KNEE ROT EXCURSION       ' num2str(DISCRETE_VARIABLES(38,2)) '       ' num2str(DISCRETE_VARIABLES(38,3))])
% disp(['HIP EXT PEAK ANGLE       ' num2str(DISCRETE_VARIABLES(39,2)) '        ' num2str(DISCRETE_VARIABLES(39,3))])
% disp(['HIP ADD PEAK ANGLE       ' num2str(DISCRETE_VARIABLES(43,2)) '       ' num2str(DISCRETE_VARIABLES(43,3))])
% disp(['HIP ADD EXCURSION        ' num2str(DISCRETE_VARIABLES(46,2)) '        ' num2str(DISCRETE_VARIABLES(46,3))])
% disp(['HIP ROT PEAK ANGLE       ' num2str(DISCRETE_VARIABLES(47,2)) '       ' num2str(DISCRETE_VARIABLES(47,3))])
% disp(['HIP ROT EXCURSION        ' num2str(DISCRETE_VARIABLES(50,2)) '       ' num2str(DISCRETE_VARIABLES(50,3))])
% disp(['FOOT PROG ANGLE          ' num2str(DISCRETE_VARIABLES(51,2)) '        ' num2str(DISCRETE_VARIABLES(51,3))])
% disp(['MHW EXC FROM TO          ' num2str(DISCRETE_VARIABLES(56,2)) '      ' num2str(DISCRETE_VARIABLES(56,3))])
% disp(['ANKLE EVE PEAK VEL       ' num2str(DISCRETE_VARIABLES(59,2)) '      ' num2str(DISCRETE_VARIABLES(59,3))])
% disp(['ANKLE ROT PEAK VEL       ' num2str(DISCRETE_VARIABLES(61,2)) '     ' num2str(DISCRETE_VARIABLES(61,3))])
% disp(['KNEE ABD PEAK VEL        ' num2str(DISCRETE_VARIABLES(65,2)) '       ' num2str(DISCRETE_VARIABLES(65,3))])
% disp(['KNEE ADD PEAK VEL        ' num2str(DISCRETE_VARIABLES(67,2)) '      ' num2str(DISCRETE_VARIABLES(67,3))])
% disp(['HIP ABD PEAK VEL         ' num2str(DISCRETE_VARIABLES(69,2)) '      ' num2str(DISCRETE_VARIABLES(69,3))])
% disp(['KNEE ROT PEAK VEL        ' num2str(DISCRETE_VARIABLES(71,2)) '       ' num2str(DISCRETE_VARIABLES(71,3))])
% disp(['HIP ROT PEAK VEL         ' num2str(DISCRETE_VARIABLES(72,2)) '       ' num2str(DISCRETE_VARIABLES(72,3))])
% disp(' ');
% % end

close(findobj('tag', 'steps_temp_figure'));
end

function [L_out,R_out] = drop_the_bad(L_data,R_data,plots)

% DROP_THE_BAD   removes the erroneous steps from the LEFT and RIGHT data.
%
%   [ L_OUT , R_OUT ] = DROP_THE_BAD ( L_DATA , R_DATA )
%   loads the LEFT and RIGHT data matricies and determines if
%   any data points are more than THREE standard deviations
%   from the mean value at that normalized time point. The
%   data is simply reassigned using the 'good' steps and
%   ignoring the 'bad' steps.
%
%
% By: Blayne Hettinga
% 10.
% Last modified: 08.17.2010
%   -means and standard deviations are now calculated outside of the loops
%
% for debugging, 'L_data' and 'R_data' need to exist ...
% L_data = norm_ang.L_ankle; R_data = norm_ang.R_ankle;

%% LEFT
% need an initial value that will be removed at the end
L_bad=0; L_good=0;
count = [0,length(L_bad)];
change = 1;
ml = mean(L_data,2);
sdl = std(L_data,1,2);

while change(end) > 0
    
    for j = 1:length(L_data(1,:,1)) % number of steps
        if sum(j == L_bad) == 1 % has the step already been flagged as 'bad'?
            break % if so, skip it
        else % otherwise, carry on
            for k = 1:3 % number of dimensions xyz
                for i = 1:length(L_data(:,1,k)) % normalized time
                    % check if the trial has already been marked as bad and skip it
                    if L_bad(end) == j
                        break
                    else
                        % is the value further than 3SD from the mean (above OR below)
                        if L_data(i,j,k) > ml(i,1,k) + 3*sdl(i,1,k) || L_data(i,j,k) < ml(i,1,k) - 3*sdl(i,1,k)
                            if j > L_bad(end) % if bad, collect the step if not already collected
                                L_bad = [L_bad,j];
                            end
                            break % no need to continue the loop if the step is bad
                        end
                    end
                end
                if k==3
                    if j > L_bad(end) % if the step is not bad, collect the step as good
                        L_good = [L_good,j];
                    end
                end
            end
        end
    end
    
    % remove the leading zero
    L_bad = L_bad(2:end);
    L_good = L_good(2:end);
    
    % an additional catch for the situation where there are two
    % groups of data that are very seperate such that 3*SD doesn't work
    % The data are split into groups using a histogram. If there are only 2
    % groups and they are the first and last, only keep the group with the most
    % data in it and get ride of the smaller group
    for i=1:length(L_data(:,1,1)) % for every normalized time point(1:101)
        for j=1:length(L_data(1,1,:)) % for every dimension xyz
            [n,bin]=hist(L_data(i,L_good,j)); % find the info to make a histogram
            if n(1)>0 && n(10)>0 && sum(n(2:9))< 1 % if there are data in the first, last, but one or less inbetween
                if n(1)>n(10) % and if there are more data in the first bin than the last,
                    [~,b] = find(L_data(i,L_good,j) < bin(2)); % find the trials in the first bin
                    L_good = L_good(b); % and only keep the ones we want
                    [~,d] = find(L_data(i,L_good,j) > bin(2)); % id the trials in the last bin
                    L_bad = [L_bad,d]; % and get rid of them
                else % or if there are more data in the last bin
                    [~,b] = find(L_data(i,L_good,j) > bin(9)); % find the trials in the last bin
                    L_good = L_good(b); % and only keep those ones
                    [~,d] = find(L_data(i,L_good,j) < bin(2)); % id the smaller group
                    L_bad = [L_bad,d]; % and ignore them
                end
                break
            end
        end
    end
    
    % add the leading zeros back in ... had to be a better way to do this ...
    L_bad = [0,L_bad];
    L_good = [0,L_good];
    
    count = [count,length(L_bad)];
    change = diff(count);
    
    ml = mean(L_data(:,L_good(2:end),:),2);
    sdl = std(L_data(:,L_good(2:end),:),1,2);
    
end

% and remove the leading zero
L_bad = L_bad(2:end);
L_good = L_good(2:end);


%%
%% RIGHT
% need an initial value that will be removed at the end
R_bad=0; R_good=0;
count = [0,length(R_bad)];
change = 1;
mr = mean(R_data,2);
sdr = std(R_data,1,2);

while change(end) > 0
    
    for j = 1:length(R_data(1,:,1)) % number of steps
        if sum(j == R_bad) == 1 % has the step already been flagged as 'bad'?
            break % if so, skip it
        else % otherwise, carry on
            for k = 1:3 % number of dimensions xyz
                for i = 1:length(R_data(:,1,k)) % normalized time
                    % check if the trial has already been marked as bad and skip it
                    if R_bad(end) == j
                        break
                    else
                        % is the value further than 3SD from the mean (above OR below)
                        if R_data(i,j,k) > mr(i,1,k) + 3*sdr(i,1,k) || R_data(i,j,k) < mr(i,1,k) - 3*sdr(i,1,k)
                            if j > R_bad(end) % if bad, collect the step if not already collected
                                R_bad = [R_bad,j];
                            end
                            break % no need to continue the loop if the step is bad
                        end
                    end
                end
                if k==3
                    if j > R_bad(end) % if the step is not bad, collect the step as good
                        R_good = [R_good,j];
                    end
                end
            end
        end
    end
    
    % remove the leading zero
    R_bad = R_bad(2:end);
    R_good = R_good(2:end);
    
    % an additional catch for the situation where there are two
    % groups of data that are very seperate such that 3*SD doesn't work
    % The data are split into groups using a histogram. If there are only 2
    % groups and they are the first and last, only keep the group with the most
    % data in it and get ride of the smaller group
    for i=1:length(R_data(:,1,1)) % for every normalized time point(1:101)
        for j=1:length(R_data(1,1,:)) % for every dimension xyz
            [n,bin]=hist(R_data(i,R_good,j)); % find the info to make a histogram
            if n(1)>0 && n(10)>0 && sum(n(2:9))< 1 % if there are data in the first, last, but one or less inbetween
                if n(1)>n(10) % and if there are more data in the first bin than the last,
                    [~,b] = find(R_data(i,R_good,j) < bin(2)); % find the trials in the first bin
                    R_good = R_good(b); % and only keep the ones we want
                    [~,d] = find(R_data(i,R_good,j) > bin(2)); % id the trials in the last bin
                    R_bad = [R_bad,d]; % and get rid of them
                else % or if there are more data in the last bin
                    [~,b] = find(R_data(i,R_good,j) > bin(9)); % find the trials in the last bin
                    R_good = R_good(b); % and only keep those ones
                    [~,d] = find(R_data(i,R_good,j) < bin(2)); % id the smaller group
                    R_bad = [R_bad,d]; % and ignore them
                end
                break
            end
        end
    end
    
    % add the leading zeros back in ... had to be a better way to do this ...
    R_bad = [0,R_bad];
    R_good = [0,R_good];
    
    count = [count,length(R_bad)];
    change = diff(count);
    
    mr = mean(R_data(:,R_good(2:end),:),2);
    sdr = std(R_data(:,R_good(2:end),:),1,2);
    
end

% and remove the leading zero
R_bad = R_bad(2:end);
R_good = R_good(2:end);


%% PLOTS
if isempty(L_good); LG = nan(101,1,3); else LG = L_data(:,L_good,:); end
if isempty(L_bad); LB = nan(101,1,3);  else LB = L_data(:,L_bad,:);  end
if isempty(R_good); RG = nan(101,1,3); else RG = R_data(:,R_good,:); end
if isempty(R_bad); RB = nan(101,1,3);  else RB = R_data(:,R_bad,:);  end

if plots ~= 0
    figure('tag','drop_the_bad_temp_figure'); hold on;
    subplot(321); plot(0:100,LG(:,:,1),'b');hold on;plot(0:100,LB(:,:,1),'r');title('SELECTED LEFT ANGLES - X')
    fill([0:100,flip(0:100,2)],[(ml(:,:,1)+sdl(:,:,1))',(flip((ml(:,:,1)-sdl(:,:,1)),1))'],[4 4 4]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(ml(:,:,1)+2*sdl(:,:,1))',(flip((ml(:,:,1)-2*sdl(:,:,1)),1))'],[5 5 5]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(ml(:,:,1)+3*sdl(:,:,1))',(flip((ml(:,:,1)-3*sdl(:,:,1)),1))'],[6 6 6]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    subplot(323); plot(0:100,LG(:,:,2),'b');hold on;plot(0:100,LB(:,:,2),'r');title('Y')
    fill([0:100,flip(0:100,2)],[(ml(:,:,2)+sdl(:,:,2))',(flip((ml(:,:,2)-sdl(:,:,2)),1))'],[4 4 4]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(ml(:,:,2)+2*sdl(:,:,2))',(flip((ml(:,:,2)-2*sdl(:,:,2)),1))'],[5 5 5]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(ml(:,:,2)+3*sdl(:,:,2))',(flip((ml(:,:,2)-3*sdl(:,:,2)),1))'],[6 6 6]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    subplot(325); plot(0:100,LG(:,:,3),'b');hold on;plot(0:100,LB(:,:,3),'r');title('Z')
    fill([0:100,flip(0:100,2)],[(ml(:,:,3)+sdl(:,:,3))',(flip((ml(:,:,3)-sdl(:,:,3)),1))'],[4 4 4]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(ml(:,:,3)+2*sdl(:,:,3))',(flip((ml(:,:,3)-2*sdl(:,:,3)),1))'],[5 5 5]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(ml(:,:,3)+3*sdl(:,:,3))',(flip((ml(:,:,3)-3*sdl(:,:,3)),1))'],[6 6 6]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    
    subplot(322); plot(0:100,RG(:,:,1),'b');hold on;plot(0:100,RB(:,:,1),'r');title('SELECTED RIGHT ANGLES - X')
    fill([0:100,flip(0:100,2)],[(mr(:,:,1)+sdr(:,:,1))',(flip((mr(:,:,1)-sdr(:,:,1)),1))'],[4 4 4]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(mr(:,:,1)+2*sdr(:,:,1))',(flip((mr(:,:,1)-2*sdr(:,:,1)),1))'],[5 5 5]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(mr(:,:,1)+3*sdr(:,:,1))',(flip((mr(:,:,1)-3*sdr(:,:,1)),1))'],[6 6 6]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    subplot(324); plot(0:100,RG(:,:,2),'b');hold on;plot(0:100,RB(:,:,2),'r');title('Y')
    fill([0:100,flip(0:100,2)],[(mr(:,:,2)+sdr(:,:,2))',(flip((mr(:,:,2)-sdr(:,:,2)),1))'],[4 4 4]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(mr(:,:,2)+2*sdr(:,:,2))',(flip((mr(:,:,2)-2*sdr(:,:,2)),1))'],[5 5 5]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(mr(:,:,2)+3*sdr(:,:,2))',(flip((mr(:,:,2)-3*sdr(:,:,2)),1))'],[6 6 6]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    subplot(326); plot(0:100,RG(:,:,3),'b');hold on;plot(0:100,RB(:,:,3),'r');title('Z')
    fill([0:100,flip(0:100,2)],[(mr(:,:,3)+sdr(:,:,3))',(flip((mr(:,:,3)-sdr(:,:,3)),1))'],[4 4 4]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(mr(:,:,3)+2*sdr(:,:,3))',(flip((mr(:,:,3)-2*sdr(:,:,3)),1))'],[5 5 5]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
    fill([0:100,flip(0:100,2)],[(mr(:,:,3)+3*sdr(:,:,3))',(flip((mr(:,:,3)-3*sdr(:,:,3)),1))'],[6 6 6]/8, 'EdgeColor', 'none', 'facealpha', 0.5);
end

%% MOVE THE CLEANED UP DATA OUT
L_out = L_data(:,L_good,:);
R_out = R_data(:,R_good,:);

end

function [ FFi_mod, FBi_mod, block_start, block_end ] = largest_block( FFi, FBi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function LARGEST_BLOCK loads in  the foot forward/foot back indices 
% (FFi, FBi) and  outputs the longest chunk of continuous events from the
% initial index (FFi_mod, FBi_mod).
% Function also outputs the index of the start/end of the entire continuous
% block (block_start, block_end).
%
% This function is necessary as there are many reason why either a FF or FB
% event could be missed. An imbalanced index (more FFs or more FBs) causes
% many downstream problems with the code. 
%
%
%  INPUTS
%  --------
%  FFI/FBI (mat):   Frame index of foot-forward/foot-back events for
%                   entire length of trial
%
%  OUTPUTS
%  -------
%  FFI_MOD/FBI_MOD (mat):   Abbreviated index of foot-forward/foot-events 
%                           containing only the longest set of continuous 
%                           events.       
%
%  BLOCK_START/BLOCK_END (int): Frame index of when the largest block of 
%                               continuous events begins and ends

% Copyright (C) 2016-2023 Allan Brett and The Running Injury Clinic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Combine and sort the FF and FB
allsort = [FFi(:) zeros(length(FFi),1); FBi(:) ones(length(FBi),1)];
[~,inds] = sort(allsort(:,1));
allsort = allsort(inds,:);

% Remove trailing FF
while allsort(end,2) ~= 1
    allsort(end,:) = [];
end

allsort_bak = allsort;

% Confirm alternating status and remove trailing chunks that aren't
% ordered

i=0;
skip = 0;
longest_length =0;

%series of while loops will search for longest continuous chunk of data
%first while loop adds up the length of a continuous segments adding in the
%indicies that are skipped (because they contain the dicontinuity)
while (sum(longest_length)+skip < size(allsort_bak,1)) && size(allsort,1)>1
    
    idx_skip = 0;
    
    %points, in case of two 1s run below
    k = 1;
    
    %allsort must start with a 0
    while allsort(1,2) == 1
        allsort(1,:) = [];
        %skip is an overall counter for the main while loop in
        %conjuction with longest_length
        skip = skip+1;
        
        %idx_skip keeps track of when values are skipped for the
        %purpose of indexing
        idx_skip = idx_skip+ 1;
        
        %remove discontinuities occuring at the start of allsort
        while allsort(k,2)==allsort(k+1,2) && k+1<size(allsort,1)
            allsort(k:k+1,:) = [];
            skip = skip + 2;
            idx_skip = idx_skip + 2;
        end
    end
    
    k = 2;
    while k <= length(allsort) && mean(allsort(1:2:k,2)) == 0 && mean(allsort(2:2:k,2)) == 1
        k = k + 2;
    end
    
    i=i+1;
    
    %we don't want to use a possibly erroneous point in the data and we
    %must end the sequence on a 1, so when the dicontinuity occurs with two
    %1s in a row, we must roll back by 2
    if size(allsort,1) > k && k > 2
        if allsort(k-1,2)==1 && allsort(k-2,2) == 1
            allsort = allsort(1:k-4,:);
        else
            allsort = allsort(1:k-2,:);
        end
    end
    
    %for the special case where there are two discontinuities of 0s in a row
    if k == 2 && allsort(1,2) == 0
        allsort(1:2,:) = [];
        longest_length(i,1) = 0;
        
        %we want to index one passed the discontinuity
        if i == 1
            %if this occurs for the first index, only includes values
            %skipped
            index(1,1) = idx_skip;
            index(2,1) = idx_skip;
        else
            index(1,i) = index(2,i-1)+idx_skip+1;
            index(2,i) = index(2,i-1)+idx_skip+1;
        end
        skip = skip + 2;
        
    else
        
        %otherwise count as normal
        longest_length(i,1) = size(allsort,1);
        
        
        %create ordered index of where continuous chunks occur
        if i == 1
            index(1,1) = 1 + idx_skip;
            index(2,1) = longest_length(i,1)+idx_skip;
        else
            index(1,i) = index(2,i-1)+3+idx_skip;
            
            %Longest_length can only be 0 when
            %two discontinuities of 1s happen in a row, below accounts that
            %the index end needs to still progress by 1 (but longest_length
            %still needs to be 0 for the main counter)
            if longest_length(i,1) >0
                index(2,i) = index(1,i)+longest_length(i,1)-1;
            else
                index(2,i) = index(1,i)+longest_length(i,1);
            end
        end
        
        
        %reset allsort for next loop iteration to be passed the discontinuity
        allsort = allsort_bak(index(2,i)+3:end,:);
        
        
        %however we want to skip passed the discontinuity to the next footfall.
        %This entails skipping the discontinuity (for example if the
        %discontinuity is two FF, we skip over these two values
        
        skip = skip+ 2;
    end
    
end

%determine which index has the largest continuous block
[~,longest_index]  = max(diff(index));

%reorder allsort to contain only this block
allsort_longest = allsort_bak(index(1,longest_index):index(2,longest_index),:);
allsort = allsort_longest;

% Break back into components
FFi_mod = allsort(allsort(:,2) == 0,1);
FBi_mod = allsort(allsort(:,2) == 1,1);

%want to track frame numbers of when the block on continuous data starts
%and ends
block_start = FFi_mod(1,1);
block_end = FBi_mod(end,1);

end



function [output ] = oscillation(trsegment,L_FFi, L_FBi, R_FFi, R_FBi, plots, label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION oscillation_run reads in a desired cluster and calculates its 
% oscillation in the veritcal plane.
%
%
%  INPUTS
%  --------
%  TRSEGMENT (struct):  Marker segment used to quantify oscillation.
%                       For gait we opt to use the average position of the
%                       3 markers of the PELVIC marker cluster.
%
%  L_FFI, L_FBI (mat):  Touchdown/toe-off  indices of the left feet
%
%  R_FFI, R_FBI (mat):  Touchdown/toe-off  indices of the right feet
%
%  HZ (int):            Sampling frequency.
%
%  PLOTS (bool):        Generate addition plotting for debugging if set
%                       equal to 1.
%
%  LABEL (str):         Gait type (walk/run) which serves to modify the
%                       window where the peak oscillation is calculated.
%
%  OUTPUTS
%  -------
%  OUTPUT (float):      Vertical oscillation (in mm) by stride separated by
%                       left and right stance
%
%
% Copyright (C) 2016-2023 Allan Brett and Running Injury Clinic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
build = 0;

%average out pelvic cluster points in the vertical plane
running_segment = trsegment(:,2);

%backup original curve for debugging
%running_segment_bak = running_segment;

A= [length(L_FFi), length(R_FFi), length(R_FBi), length(L_FBi)];

%% utlize moving window to minimize errors

%ignore first touchdown, too small of an initial window
%causes errors. Similarly, ignore the last touchdown

for k = 1:min(A)-1
    
    %define moving window, window from touchdown to
    %touchdown for each foot, add in 10 frame cushion
    
    if strcmp(label,'run')
        
        l_window = [L_FFi(k);L_FBi(k)+25];
        r_window = [R_FFi(k);R_FBi(k)+25] ;
        
    elseif strcmp(label,'walk')
        
        l_window = [L_FFi(k)-10;R_FFi(k)-5];
        r_window = [R_FFi(k)-10;L_FFi(k+1)-5] ;
        
    end
    
    
    %track peaks/troughs for each approximate stride
    left = running_segment(l_window(1,1):l_window(2,1));
    right = running_segment(r_window(1,1):r_window(2,1));
    
    [left_peak, left_peak_loc ] = findpeaks(left);
    [right_peak, right_peak_loc ] = findpeaks(right);
    [left_trough, left_trough_loc ] = findpeaks(-left);
    [right_trough, right_trough_loc ] = findpeaks(-right);
    
    %Occasionally an error in touchdown or toe-off index can lead to the
    %window being incorrect, this leads to multiple or no peaks. When this
    %occurs, return NaN and move on.
    if size(left_peak,1) ~= 1 || size(right_peak,1) ~= 1 || size(right_trough,1) ~=  1 || size(left_trough,1) ~= 1
        
        build = build + 1;
        peak_locs(build)=NaN;
        trough_locs(build) = NaN;
        peaks(build)=NaN;
        trough(build) = NaN;
        
        build = build+1;
        peak_locs(build) = NaN;
        trough_locs(build) = NaN;
        peaks(build)=NaN;
        trough(build) = NaN;
        
        
        
    else
        
        
        %check for inflection points instead of peaks
        
        for x = 1:length(left_peak)
            check_l_peak(x)  = find(left == left_peak(x));
            if left(check_l_peak(x)) < left(check_l_peak(x) + 1)
                left_peak(x) = 0;
                left_peak_loc(x) = 0;
            end
        end
        
        for x = 1:length(right_peak)
            check_r_peak(x)  = find(right == right_peak(x));
            if right(check_r_peak(x)) < right(check_r_peak(x) + 1)
                right_peak(x) = 0;
                right_peak_loc(x) = 0;
            end
        end
        
        for x = 1:length(left_trough)
            check_l_trough(x)  = find(-left == left_trough(x));
            if left(check_l_trough(x)) > left(check_l_trough(x) + 1)
                left_trough(x) = 0;
                left_trough_loc(x) = 0;
            end
        end
        
        for x = 1:length(right_trough)
            check_r_trough(x)  = find(-right == right_trough(x));
            if right(check_r_trough(x)) > right(check_r_trough(x) + 1)
                right_trough(x) = 0;
                right_trough_loc(x) = 0;
            end
        end
        
        %build peaks and troughs matrices by stacking right and lefts
        %this is necessary for confirming alternating status below
        build = build + 1;
        peak_locs(build)=left_peak_loc + L_FFi(k)-1;
        trough_locs(build) = left_trough_loc + L_FFi(k)-1;
        peaks(build)=left_peak;
        trough(build) = -left_trough;
        
        build = build+1;
        peak_locs(build) = right_peak_loc+R_FFi(k)-1;
        trough_locs(build) = right_trough_loc+R_FFi(k)-1;
        peaks(build)=right_peak;
        trough(build) = -right_trough;
        
    end
end

%set bottom and top edges for plotting
bot = min(running_segment);
top = max(running_segment);

%% calculate oscillation

if peak_locs(1) < trough_locs(1)
    peak_locs(1) = [];
    peaks(1) = [];
end

allsort = [peak_locs(:) zeros(length(peak_locs),1) peaks(:); trough_locs(:) ones(length(trough_locs),1) trough(:)];
[~,inds] = sort(allsort(:,1));
allsort = allsort(inds,:);

%remove trailing trough

while allsort(end,2) ~= 0
    allsort(end,:) = [];
end

k = 2;

%ensure alternating status

while k <= length(allsort) && mean(allsort(1:2:k,2)) == 1 && mean(allsort(2:2:k,2)) == 0
    k = k + 2;
end
allsort = allsort(1:k-2,:);

% Break back into components
peak_locs = allsort(allsort(:,2) == 0,1);
trough_locs = allsort(allsort(:,2) == 1,1);
peaks = allsort(allsort(:,2)==0,3);
trough = allsort(allsort(:,2)==1,3);

%debugging
if plots ~=0
    %plot the curves
    figure
    
    %include offset caused by chopping by first left footfall
    
    for i = 1:length(running_segment)
        if i == 1
            x_axis(i) = L_FFi(1);
        else
            x_axis(i) = x_axis(i-1) + 1;
        end
    end
    
    title('Vertical')
    
    plot (running_segment_bak)
    hold on
    plot(trough_locs,trough, 'r.')
    plot (peak_locs,peaks, 'g.')
    
    %plot touchdown/toeoffs as well
    for p=1:min(A)-1
        fill([L_FFi(p,1),L_FBi(p,1),L_FBi(p,1),L_FFi(p,1)],[bot, bot, top, top], ...
            [3 3 7]/8, 'EdgeColor','none','facealpha',0.2);
        fill([R_FFi(p,1),R_FBi(p,1),R_FBi(p,1),R_FFi(p,1)],[bot, bot, top, top], ...
            [7 3 3]/8, 'EdgeColor','none','facealpha',0.2);
        
    end
    hold off
end


%define oscillation as the distance between peak and trough

vert_oscillation =  peaks-trough;

%vert_oscillation always begin with stance for the left foot
left_stance = vert_oscillation(1:2:end,:);
right_stance =vert_oscillation(2:2:end,:);

%chop to shortest length
min_length = min(length(left_stance), length(right_stance));

left_stance = left_stance(1:min_length,1);
right_stance = right_stance(1:min_length,1);

for i = 1:min_length
    if i == 1
        stride_num(i,1)  = 1;
    else
        stride_num(i,1) = stride_num(i-1,1)+1;
    end
end

stride_num = peak_locs;

output = horzcat(stride_num(1:2:end), left_stance, stride_num(2:2:end), right_stance);

%% quality checks if there are more nans than numbers, something has gone wrong

left_check_nan = sum(isnan(left_stance));
left_check_num = sum(~isnan(left_stance));

right_check_nan = sum(isnan(right_stance));
right_check_num = sum(~isnan(right_stance));

if left_check_nan > left_check_num || right_check_nan > right_check_num
    
    error('Number of rejected strides for vertical oscillation too high')
    
end

clear allsort peak_locs peaks trough trough_locs balance oscillation_up oscillation_down left_ratio right_ratio R_FBi L_FBi L_FFi R_FFi

end

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

end

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


end