function [angles,velocities,jc,R,djc] = gait_kinematics(joints,neutral,dynamic,hz,plots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Loads marker locations for JOINTS centre locations, the segment
%   markers during a standing NEUTRAL trial and a walking or running
%   DYNAMIC trial and uses the data with the joint coordinate
%   system and soderqvist (single value decomposition) to calculate
%   joint angles and joint angular velocities. This data is
%   combined into a structure of joint angles for the hip, knee,
%   and ankle for the x, y, and z axes of rotation. Similarly, the
%   joint angular velocities will be combined into one structure
%   as well. Inetermediairy steps that create the transformation
%   matricies and joint centre locations are also outputted. The
%   fourth input is merely an indication if the user desires for
%   the outputs to be visulized in PLOTS.
%
%
%  INPUTS
%  --------
%  JOINTS (struct): Joint centre locations collected as part of the static
%                   trial
%
%  NEUTRAL (stuct): Marker shell positions collected as part of the static
%                   trial.
%
%  DYNAMIC (struct):    Marker shell positions collected as part of the
%                       dynamic (run/walk) trial.
%
%  HZ (int):    Data collection sampling frequency.
%
%  PLOTS (bool):    Boolean selected to generate plotted outcomes if 
%                   desired. If no second argument exists, or if PLOTS == 0
%                   , the plotted outputs in this function are suppressed.
%
%
%  OUTPUTS
%  -------
%  ANGLES(struct):  Calculated joint angles. 
%
%  VELOCITIES(struct):  Calculated joint velocities.
%
%  JC (struct): Calculated joint centres. 
%
%  R (struct): Rotation matrix used to calculate joint angles.
%
%  DJC (struct): Distance to joint centre from segment centroid.
%
%  LICENSE
%  -------
%  See file LICENSE.txt
%
% Copyright (C) 2010-2023,  Blayne Hettinga, Allan Brett and 
%                           The Running Injury Clinic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE the Bonita Lab Coordinate system
%       X - points to the subject's right
%       Y - points in the direction of walking
%       Z - points vertically upwards
% is switched to the
%% the following Lab Coordinate system
%       X - points to the subject's right
%       Y - points vertically upwards
%       Z - points opposite of the walking direction
%% NOTE the Segment coordinate systems
%       X - Anterior .....................................[AB / AD duction]
%       Y - Vertically upwards ............................[Axial rotation]
%       Z - points to the subject's right side ...[Hinge flexion extension]


%% if no indication to plot the data was specified, the default is to not plot
if (nargin < 5)
    plots = 0;
end


%% calculate joint centre locations from the data in the JOINTS structure

% the pelvis_jc is simply the average location of the pelvis markers
jc.pelvis = (neutral.pelvis_1+neutral.pelvis_2+neutral.pelvis_3+neutral.pelvis_4)/4;
% 25% of the distance between hips
jc.L_hip = joints.L_hip+(joints.R_hip-joints.L_hip)/4;
jc.R_hip = joints.R_hip+(joints.L_hip-joints.R_hip)/4;
% midpoint of the two knee markers
jc.L_knee = (joints.L_lat_knee+joints.L_med_knee)/2;
jc.R_knee = (joints.R_lat_knee+joints.R_med_knee)/2;
% midpoint of the two ankle markers
jc.L_ankle = (joints.L_lat_ankle+joints.L_med_ankle)/2;
jc.R_ankle = (joints.R_med_ankle+joints.R_lat_ankle)/2;

%% plot the markers in xyz to plot the segment coordinate systems

if plots ~= 0
    %%
    figure('tag', 'bonita_kinematics_temp_figure'); hold on; axis equal; %axis([-1 0.2 -0.2 0.6 -0.1 1.2]);
    xlabel('x'); ylabel('y'); zlabel('z')
    rotate3d on
    % plot the neutral marker locations
    nc=struct2cell(neutral);
    for i=1:length(nc)
        plot3(nc{i}(:,1),nc{i}(:,2),nc{i}(:,3),'.b')
    end
    % plot the joint centre locations
    jcc=struct2cell(jc);
    for i=1:length(jcc)
        plot3(jcc{i}(:,1),jcc{i}(:,2),jcc{i}(:,3),'or')
    end
else
end


%% Create the anatomical coordinate system unit vectors for each segment in
% terms of the Global coordinate system

% define how long to make the CS lines if there are being plotted
cs_length = 50;

%% LEFT FOOT
% long axis of the the foot is aligned with the lab
l_foot_x = [0 0 -1];
% vertical axis of the foot is aligned with the two markers on the heel
% first need to identify the two heel markers
% Combine 3 of the foot markers into one matrix (ignore the created fourth)
L_foot = [neutral.L_foot_1;neutral.L_foot_2;neutral.L_foot_3];
% sort the markers from left to right
L_foot = sortrows(L_foot,1);
% second, create a vector from the two right markers (the left is the lateral heel)
l_foot_temp = (L_foot(2,:)-L_foot(3,:))/norm(L_foot(2,:)-L_foot(3,:));
% check that it is pointing up ... and if not, flip it
if l_foot_temp(2) < 0
    l_foot_temp = -l_foot_temp;
end
% use the temp vertical axis to create the lateral axis
l_foot_z = cross(l_foot_x,l_foot_temp)/norm(cross(l_foot_x,l_foot_temp));
% and create the 'vertical' axis that provides standing eversion angle
l_foot_y = cross(l_foot_z,l_foot_x)/norm(cross(l_foot_z,l_foot_x));
% combine to create a transformation matrix from anatomical to global
ag.L_foot=[l_foot_x',l_foot_y',l_foot_z'];
% visualize
if plots ~= 0
    lfx=[jc.L_ankle;(jc.L_ankle+cs_length*l_foot_x)];line(lfx(:,1),lfx(:,2),lfx(:,3),'Color','r')
    lfy=[jc.L_ankle;(jc.L_ankle+cs_length*l_foot_y)];line(lfy(:,1),lfy(:,2),lfy(:,3),'Color','g')
    lfz=[jc.L_ankle;(jc.L_ankle+cs_length*l_foot_z)];line(lfz(:,1),lfz(:,2),lfz(:,3),'Color','b')
else
end
%% RIGHT FOOT
% long axis of the the foot is aligned with the lab
r_foot_x = [0 0 -1];
% vertical axis of the foot is aligned with the two markers on the heel
% first need to identify the heel markers
% Combine 3 feet markers into one matrix ... can ignore the created fourth one
R_foot = [neutral.R_foot_1;neutral.R_foot_2;neutral.R_foot_3];
% sort the points from left to right
R_foot = sortrows(R_foot,1);
% second, create a vector from the two left markers (not the lateral one)
r_foot_temp = (R_foot(1,:)-R_foot(2,:))/norm(R_foot(1,:)-R_foot(2,:));
% check that it is pointing up ... and if not, flip it
if r_foot_temp(2) < 0
    r_foot_temp = -r_foot_temp;
end
% use the temp vertical axis to create the lateral axis
r_foot_z = cross(r_foot_x,r_foot_temp)/norm(cross(r_foot_x,r_foot_temp));
% and create the 'vertical' axis that provides standing eversion angle
r_foot_y = cross(r_foot_z,r_foot_x)/norm(cross(r_foot_z,r_foot_x));
% combine to create a transformation matrix from anatomical to global
ag.R_foot=[r_foot_x',r_foot_y',r_foot_z'];
% visualize
if plots ~= 0
    rfx=[jc.R_ankle;(jc.R_ankle+cs_length*r_foot_x)];line(rfx(:,1),rfx(:,2),rfx(:,3),'Color','r')
    rfy=[jc.R_ankle;(jc.R_ankle+cs_length*r_foot_y)];line(rfy(:,1),rfy(:,2),rfy(:,3),'Color','g')
    rfz=[jc.R_ankle;(jc.R_ankle+cs_length*r_foot_z)];line(rfz(:,1),rfz(:,2),rfz(:,3),'Color','b')
else
end
%% LEFT SHANK
l_shank_y=(jc.L_knee-jc.L_ankle)/norm(jc.L_knee-jc.L_ankle); %long axis pointing up
l_shank_temp=(joints.L_med_ankle-joints.L_lat_ankle)/norm(joints.L_med_ankle-joints.L_lat_ankle); %almost the hinge joing pointing to the right
l_shank_x=(cross(l_shank_y,l_shank_temp))/norm(cross(l_shank_y,l_shank_temp)); %Anterior axis from the cross
l_shank_z=(cross(l_shank_x,l_shank_y))/norm(cross(l_shank_x,l_shank_y)); %hinge axis, points to the right
% combine to create a transformation matrix from anatomical to global
ag.L_shank=[l_shank_x',l_shank_y',l_shank_z'];
% visualize
if plots ~= 0
    rsx=[jc.L_knee;(jc.L_knee+cs_length*l_shank_x)];line(rsx(:,1),rsx(:,2),rsx(:,3),'Color','r')
    rsy=[jc.L_knee;(jc.L_knee+cs_length*l_shank_y)];line(rsy(:,1),rsy(:,2),rsy(:,3),'Color','g')
    rsz=[jc.L_knee;(jc.L_knee+cs_length*l_shank_z)];line(rsz(:,1),rsz(:,2),rsz(:,3),'Color','b')
else
end
%% RIGHT SHANK
r_shank_y=(jc.R_knee-jc.R_ankle)/norm(jc.R_knee-jc.R_ankle); %long axis pointing up
r_shank_temp=(joints.R_lat_ankle-joints.R_med_ankle)/norm(joints.R_lat_ankle-joints.R_med_ankle); %almost the hinge joing pointing to the right
r_shank_x=(cross(r_shank_y,r_shank_temp))/norm(cross(r_shank_y,r_shank_temp)); %Anterior axis from the cross
r_shank_z=(cross(r_shank_x,r_shank_y))/norm(cross(r_shank_x,r_shank_y)); %hinge axis, lateral for right
% combine to create a transformation matrix from anatomical to global
ag.R_shank=[r_shank_x',r_shank_y',r_shank_z'];
% visualize
if plots ~= 0
    rsx=[jc.R_knee;(jc.R_knee+cs_length*r_shank_x)];line(rsx(:,1),rsx(:,2),rsx(:,3),'Color','r')
    rsy=[jc.R_knee;(jc.R_knee+cs_length*r_shank_y)];line(rsy(:,1),rsy(:,2),rsy(:,3),'Color','g')
    rsz=[jc.R_knee;(jc.R_knee+cs_length*r_shank_z)];line(rsz(:,1),rsz(:,2),rsz(:,3),'Color','b')
else
end
%% LEFT THIGH
l_thigh_y=(jc.L_hip-jc.L_knee)/norm(jc.L_hip-jc.L_knee); %long axis pointing up
%almost the hinge joing pointing to the right from the knee joint markers
l_thigh_temp=(joints.L_med_knee-joints.L_lat_knee)/norm(joints.L_med_knee-joints.L_lat_knee);
l_thigh_x=(cross(l_thigh_y,l_thigh_temp))/norm(cross(l_thigh_y,l_thigh_temp)); %anterior axis
l_thigh_z=(cross(l_thigh_x,l_thigh_y))/norm(cross(l_thigh_x,l_thigh_y)); %hinge, to the right
% combine to create a transformation matrix from anatomical to global
ag.L_thigh=[l_thigh_x',l_thigh_y',l_thigh_z'];
% visulize
if plots ~= 0
    rtx=[jc.L_hip;(jc.L_hip+cs_length*l_thigh_x)];line(rtx(:,1),rtx(:,2),rtx(:,3),'Color','r')
    rty=[jc.L_hip;(jc.L_hip+cs_length*l_thigh_y)];line(rty(:,1),rty(:,2),rty(:,3),'Color','g')
    rtz=[jc.L_hip;(jc.L_hip+cs_length*l_thigh_z)];line(rtz(:,1),rtz(:,2),rtz(:,3),'Color','b')
else
end
%% RIGHT THIGH
r_thigh_y=(jc.R_hip-jc.R_knee)/norm(jc.R_hip-jc.R_knee); %long axis pointing up
%almost the hinge joing pointing to the right from the knee joint markers
r_thigh_temp=(joints.R_lat_knee-joints.R_med_knee)/norm(joints.R_lat_knee-joints.R_med_knee);
r_thigh_x=(cross(r_thigh_y,r_thigh_temp))/norm(cross(r_thigh_y,r_thigh_temp)); %anterior axis
r_thigh_z=(cross(r_thigh_x,r_thigh_y))/norm(cross(r_thigh_x,r_thigh_y)); %hinge, to the right
% combine to create a transformation matrix from anatomical to global
ag.R_thigh=[r_thigh_x',r_thigh_y',r_thigh_z'];
% visualize
if plots ~= 0
    rtx=[jc.R_hip;(jc.R_hip+cs_length*r_thigh_x)];line(rtx(:,1),rtx(:,2),rtx(:,3),'Color','r')
    rty=[jc.R_hip;(jc.R_hip+cs_length*r_thigh_y)];line(rty(:,1),rty(:,2),rty(:,3),'Color','g')
    rtz=[jc.R_hip;(jc.R_hip+cs_length*r_thigh_z)];line(rtz(:,1),rtz(:,2),rtz(:,3),'Color','b')
else
end
%% PELVIS
% since trochanters are hard to landmark, pelvis will just be orthogonal to the lab
pelvis_x =[0 0 -1]; %anterior axis
pelvis_y =[0 1 0]; %long axis pointing up
pelvis_z =[1 0 0]; %hinge axis to the subject's right
% combine to create a transformation matrix from anatomical to global
ag.pelvis=[pelvis_x',pelvis_y',pelvis_z'];
% visualize
if plots ~= 0
    px=[jc.pelvis;(jc.pelvis+cs_length*pelvis_x)];line(px(:,1),px(:,2),px(:,3),'Color','r')
    py=[jc.pelvis;(jc.pelvis+cs_length*pelvis_y)];line(py(:,1),py(:,2),py(:,3),'Color','g')
    pz=[jc.pelvis;(jc.pelvis+cs_length*pelvis_z)];line(pz(:,1),pz(:,2),pz(:,3),'Color','b')
else
end

%% use the transformation matrixes we just made to get the segment markers
% in terms of the ANATOMICAL COORDINATE system (ie if the anatomical
% coordinate system was located at the world orgin, where would the markers
% be located). this can be better understood by visualizing the vectors
% defining the marker location.
% the marker locations in this reference frame will be used with the
% soderkvist method as 'before' inputs.

% the foot markers in an ANATOMICAL COORDINATE system ...
ac.L_foot_1=ag.L_foot'*neutral.L_foot_1';
ac.L_foot_2=ag.L_foot'*neutral.L_foot_2';
ac.L_foot_3=ag.L_foot'*neutral.L_foot_3';
ac.L_foot_4=ag.L_foot'*neutral.L_foot_4';
ac.L_foot=[ac.L_foot_1,ac.L_foot_2,ac.L_foot_3,ac.L_foot_4];
% trying to get data to plot the knee joint centre with the dynamic data
% create a vector from a segment marker to the joint centre in the segment
% coordinate system

% distance to joint centre from segment centroid
centroid = mean(ac.L_foot,2);

djc.L_ankle=ag.L_foot'*jc.L_ankle'-centroid;

ac.R_foot_1=ag.R_foot'*neutral.R_foot_1';
ac.R_foot_2=ag.R_foot'*neutral.R_foot_2';
ac.R_foot_3=ag.R_foot'*neutral.R_foot_3';
ac.R_foot_4=ag.R_foot'*neutral.R_foot_4';
ac.R_foot=[ac.R_foot_1,ac.R_foot_2,ac.R_foot_3,ac.R_foot_4];

centroid = mean(ac.R_foot,2);

djc.R_ankle=ag.R_foot'*jc.R_ankle'-centroid;

ac.L_shank_1=ag.L_shank'*neutral.L_shank_1';
ac.L_shank_2=ag.L_shank'*neutral.L_shank_2';
ac.L_shank_3=ag.L_shank'*neutral.L_shank_3';
ac.L_shank_4=ag.L_shank'*neutral.L_shank_4';
ac.L_shank=[ac.L_shank_1,ac.L_shank_2,ac.L_shank_3,ac.L_shank_4];

centroid = mean(ac.L_shank,2);

djc.L_knee=ag.L_shank'*jc.L_knee'-centroid;

ac.R_shank_1=ag.R_shank'*neutral.R_shank_1';
ac.R_shank_2=ag.R_shank'*neutral.R_shank_2';
ac.R_shank_3=ag.R_shank'*neutral.R_shank_3';
ac.R_shank_4=ag.R_shank'*neutral.R_shank_4';
ac.R_shank=[ac.R_shank_1,ac.R_shank_2,ac.R_shank_3,ac.R_shank_4];

centroid = mean(ac.R_shank,2);

djc.R_knee=ag.R_shank'*jc.R_knee'-centroid;
ac.L_thigh_1=ag.L_thigh'*neutral.L_thigh_1';
ac.L_thigh_2=ag.L_thigh'*neutral.L_thigh_2';
ac.L_thigh_3=ag.L_thigh'*neutral.L_thigh_3';
ac.L_thigh_4=ag.L_thigh'*neutral.L_thigh_4';
ac.L_thigh=[ac.L_thigh_1,ac.L_thigh_2,ac.L_thigh_3,ac.L_thigh_4];

centroid = mean(ac.L_thigh,2);

djc.L_hip=ag.L_thigh'*jc.L_hip'-centroid;


ac.R_thigh_1=ag.R_thigh'*neutral.R_thigh_1';
ac.R_thigh_2=ag.R_thigh'*neutral.R_thigh_2';
ac.R_thigh_3=ag.R_thigh'*neutral.R_thigh_3';
ac.R_thigh_4=ag.R_thigh'*neutral.R_thigh_4';
ac.R_thigh=[ac.R_thigh_1,ac.R_thigh_2,ac.R_thigh_3,ac.R_thigh_4];

centroid = mean(ac.R_thigh,2);

djc.R_hip=ag.R_thigh'*jc.R_hip'-centroid;

ac.pelvis_1=ag.pelvis'*neutral.pelvis_1';
ac.pelvis_2=ag.pelvis'*neutral.pelvis_2';
ac.pelvis_3=ag.pelvis'*neutral.pelvis_3';
ac.pelvis_4=ag.pelvis'*neutral.pelvis_4';
ac.pelvis=[ac.pelvis_1,ac.pelvis_2,ac.pelvis_3,ac.pelvis_4];

centroid = mean(ac.pelvis,2);

djc.pelvis=ag.pelvis'*jc.pelvis'-centroid;


%% create the rotation matrix from segment to lab for all the segments for
% every time point. this uses soderkvist. the operations done on the
% neutral trial data only need to be done once and can be done before we
% load the dynamic data and are therefore outside of the loop

% calculate average marker position of the segment markers
avg.ac_L_foot=mean(ac.L_foot,2);
avg.ac_R_foot=mean(ac.R_foot,2);
avg.ac_L_shank=mean(ac.L_shank,2);
avg.ac_R_shank=mean(ac.R_shank,2);
avg.ac_L_thigh=mean(ac.L_thigh,2);
avg.ac_R_thigh=mean(ac.R_thigh,2);
avg.ac_pelvis=mean(ac.pelvis,2);

% calculate the distance from the marker posisiton to the average position
% neutral
dif.ac_L_foot=[ac.L_foot(1,:) - avg.ac_L_foot(1);...
    ac.L_foot(2,:) - avg.ac_L_foot(2);...
    ac.L_foot(3,:) - avg.ac_L_foot(3)];

dif.ac_R_foot=[ac.R_foot(1,:) - avg.ac_R_foot(1);...
    ac.R_foot(2,:) - avg.ac_R_foot(2);...
    ac.R_foot(3,:) - avg.ac_R_foot(3)];

dif.ac_L_shank=[ac.L_shank(1,:) - avg.ac_L_shank(1);...
    ac.L_shank(2,:) - avg.ac_L_shank(2);...
    ac.L_shank(3,:) - avg.ac_L_shank(3)];

dif.ac_R_shank=[ac.R_shank(1,:) - avg.ac_R_shank(1);...
    ac.R_shank(2,:) - avg.ac_R_shank(2);...
    ac.R_shank(3,:) - avg.ac_R_shank(3)];

dif.ac_L_thigh=[ac.L_thigh(1,:) - avg.ac_L_thigh(1);...
    ac.L_thigh(2,:) - avg.ac_L_thigh(2);...
    ac.L_thigh(3,:) - avg.ac_L_thigh(3)];

dif.ac_R_thigh=[ac.R_thigh(1,:) - avg.ac_R_thigh(1);...
    ac.R_thigh(2,:) - avg.ac_R_thigh(2);...
    ac.R_thigh(3,:) - avg.ac_R_thigh(3)];

dif.ac_pelvis=[ac.pelvis(1,:) - avg.ac_pelvis(1);...
    ac.pelvis(2,:) - avg.ac_pelvis(2);...
    ac.pelvis(3,:) - avg.ac_pelvis(3)];

%% now we need the DYNAMIC walking/running data from the input argument

% prealocating the size of these matricies that collect data in the loop
% to increase speed
R.L_foot = zeros(4,4,length(dynamic.L_foot_1(:,1)));
R.R_foot = R.L_foot;
R.L_shank = R.L_foot;
R.R_shank = R.L_foot;
R.L_thigh = R.L_foot;
R.R_thigh = R.L_foot;
R.pelvis = R.L_foot;

R.L_ankle = R.L_foot;
R.R_ankle = R.L_foot;
R.L_knee = R.L_foot;
R.R_knee = R.L_foot;
R.L_hip = R.L_foot;
R.R_hip = R.L_foot;

angle.L_ankle = zeros(size(dynamic.L_foot_1));
angle.R_ankle = angle.L_ankle;
angle.L_knee = angle.L_ankle;
angle.R_knee = angle.L_ankle;
angle.L_hip = angle.L_ankle;
angle.R_hip = angle.L_ankle;
angle.L_foot = angle.L_ankle;
angle.R_foot = angle.L_ankle;
angle.pelvis = angle.L_ankle;

%%
for i=1:length(dynamic.L_foot_1(:,1))
    
    %% need to create matricies from the DYNAMIC data that has the same markers
    % as the " ac.L_foot " in the same order ... matricies for use with soderqvist
    
    d.L_foot=[dynamic.L_foot_1(i,:)',dynamic.L_foot_2(i,:)',dynamic.L_foot_3(i,:)',dynamic.L_foot_4(i,:)'];
    d.R_foot=[dynamic.R_foot_1(i,:)',dynamic.R_foot_2(i,:)',dynamic.R_foot_3(i,:)',dynamic.R_foot_4(i,:)'];
    d.L_shank=[dynamic.L_shank_1(i,:)',dynamic.L_shank_2(i,:)',dynamic.L_shank_3(i,:)',dynamic.L_shank_4(i,:)'];
    d.R_shank=[dynamic.R_shank_1(i,:)',dynamic.R_shank_2(i,:)',dynamic.R_shank_3(i,:)',dynamic.R_shank_4(i,:)'];
    d.L_thigh=[dynamic.L_thigh_1(i,:)',dynamic.L_thigh_2(i,:)',dynamic.L_thigh_3(i,:)',dynamic.L_thigh_4(i,:)'];
    d.R_thigh=[dynamic.R_thigh_1(i,:)',dynamic.R_thigh_2(i,:)',dynamic.R_thigh_3(i,:)',dynamic.R_thigh_4(i,:)'];
    d.pelvis=[dynamic.pelvis_1(i,:)',dynamic.pelvis_2(i,:)',dynamic.pelvis_3(i,:)',dynamic.pelvis_4(i,:)'];
    
    
    %% calculate average position of the segment markers
    
    avg.d_L_foot=mean(d.L_foot,2);
    avg.d_R_foot=mean(d.R_foot,2);
    avg.d_L_shank=mean(d.L_shank,2);
    avg.d_R_shank=mean(d.R_shank,2);
    avg.d_L_thigh=mean(d.L_thigh,2);
    avg.d_R_thigh=mean(d.R_thigh,2);
    avg.d_pelvis=mean(d.pelvis,2);
    
    
    %% caluculate the distance from the marker posisiton to the average position
    
    dif.d_L_foot=[d.L_foot(1,:) - avg.d_L_foot(1);...
        d.L_foot(2,:) - avg.d_L_foot(2);...
        d.L_foot(3,:) - avg.d_L_foot(3)];
    
    dif.d_R_foot=[d.R_foot(1,:) - avg.d_R_foot(1);...
        d.R_foot(2,:) - avg.d_R_foot(2);...
        d.R_foot(3,:) - avg.d_R_foot(3)];
    
    dif.d_L_shank=[d.L_shank(1,:) - avg.d_L_shank(1);...
        d.L_shank(2,:) - avg.d_L_shank(2);...
        d.L_shank(3,:) - avg.d_L_shank(3)];
    
    dif.d_R_shank=[d.R_shank(1,:) - avg.d_R_shank(1);...
        d.R_shank(2,:) - avg.d_R_shank(2);...
        d.R_shank(3,:) - avg.d_R_shank(3)];
    
    dif.d_L_thigh=[d.L_thigh(1,:) - avg.d_L_thigh(1);...
        d.L_thigh(2,:) - avg.d_L_thigh(2);...
        d.L_thigh(3,:) - avg.d_L_thigh(3)];
    
    dif.d_R_thigh=[d.R_thigh(1,:) - avg.d_R_thigh(1);...
        d.R_thigh(2,:) - avg.d_R_thigh(2);...
        d.R_thigh(3,:) - avg.d_R_thigh(3)];
    
    dif.d_pelvis=[d.pelvis(1,:) - avg.d_pelvis(1);...
        d.pelvis(2,:) - avg.d_pelvis(2);...
        d.pelvis(3,:) - avg.d_pelvis(3)];
    
    
    %% step three of soderkvist
    
    C.L_foot = dif.d_L_foot * dif.ac_L_foot';
    C.R_foot = dif.d_R_foot * dif.ac_R_foot';
    C.L_shank = dif.d_L_shank * dif.ac_L_shank';
    C.R_shank = dif.d_R_shank * dif.ac_R_shank';
    C.L_thigh = dif.d_L_thigh * dif.ac_L_thigh';
    C.R_thigh = dif.d_R_thigh * dif.ac_R_thigh';
    C.pelvis = dif.d_pelvis * dif.ac_pelvis';
    
    
    %% step four - do a single value decomp of C
    
    [P.L_foot,T.L_foot,Q.L_foot]=svd(C.L_foot);
    [P.R_foot,T.R_foot,Q.R_foot]=svd(C.R_foot);
    [P.L_shank,T.L_shank,Q.L_shank]=svd(C.L_shank);
    [P.R_shank,T.R_shank,Q.R_shank]=svd(C.R_shank);
    [P.L_thigh,T.L_thigh,Q.L_thigh]=svd(C.L_thigh);
    [P.R_thigh,T.R_thigh,Q.R_thigh]=svd(C.R_thigh);
    [P.pelvis,T.pelvis,Q.pelvis]=svd(C.pelvis);
    
    
    %% step five - calculate a rotation matrix
    
    R_L_foot = P.L_foot * diag([1;1;det(P.L_foot*Q.L_foot')]) * Q.L_foot';
    R_R_foot = P.R_foot * diag([1;1;det(P.R_foot*Q.R_foot')]) * Q.R_foot';
    R_L_shank = P.L_shank * diag([1;1;det(P.L_shank*Q.L_shank')]) * Q.L_shank';
    R_R_shank = P.R_shank * diag([1;1;det(P.R_shank*Q.R_shank')]) * Q.R_shank';
    R_L_thigh = P.L_thigh * diag([1;1;det(P.L_thigh*Q.L_thigh')]) * Q.L_thigh';
    R_R_thigh = P.R_thigh * diag([1;1;det(P.R_thigh*Q.R_thigh')]) * Q.R_thigh';
    R_pelvis = P.pelvis * diag([1;1;det(P.pelvis*Q.pelvis')]) * Q.pelvis';
    
    
    %% step six - calculate the displacement
    
    dis_L_foot = avg.d_L_foot - R_L_foot * avg.ac_L_foot;
    dis_R_foot = avg.d_R_foot - R_R_foot * avg.ac_R_foot;
    dis_L_shank = avg.d_L_shank - R_L_shank * avg.ac_L_shank;
    dis_R_shank = avg.d_R_shank - R_R_shank * avg.ac_R_shank;
    dis_L_thigh = avg.d_L_thigh - R_L_thigh * avg.ac_L_thigh;
    dis_R_thigh = avg.d_R_thigh - R_R_thigh * avg.ac_R_thigh;
    dis_pelvis = avg.d_pelvis - R_pelvis * avg.ac_pelvis;
    
    
    %% combine the rotation matricies in a stack
    
    R.L_foot(:,:,i)=[R_L_foot,dis_L_foot;0,0,0,1];
    R.R_foot(:,:,i)=[R_R_foot,dis_R_foot;0,0,0,1];
    R.L_shank(:,:,i)=[R_L_shank,dis_L_shank;0,0,0,1];
    R.R_shank(:,:,i)=[R_R_shank,dis_R_shank;0,0,0,1];
    R.L_thigh(:,:,i)=[R_L_thigh,dis_L_thigh;0,0,0,1];
    R.R_thigh(:,:,i)=[R_R_thigh,dis_R_thigh;0,0,0,1];
    R.pelvis(:,:,i)=[R_pelvis,dis_pelvis;0,0,0,1];
    
    %% calculate segment angles
    % FOOT
    % the rotation matrix of the foot to lab (R.L_foot) is simply a
    % combination of the XYZ axes of the segment coordinate system unit
    % vectors
    %   footX_labx     footY_labx     footZ_labx
    %   footX_laby     footY_laby     footZ_laby
    %   footX_labz     footY_labz     footZ_labz
    % so in order to calculate foot progression angle / heel whip ...
    % Collect the angle of the long axis of the foot about the
    % vertical axis of the foot.
    % Also by projecting the foot into the sagital plane, we have
    % information to determine heelstrikers or FOREFOOT strikers when
    % identifying events
    
    % angle of the vertical axis of the foot projected into the frontal plane
    %   from a posterior view, 'vertical vector' in the first quadrant is
    %   postitive and second quadrant is negative
    %   . ie. inversion is negative and eversion is positive for the LEFT
    % note that this is projected into the frontal plane so cross talk is
    % present ... dorsi flexion with an abducted foot creates lots of "inversion"
    angle.L_foot(i,1) = atan(R.L_foot(1,2,i)/sqrt(R.L_foot(2,2,i)^2 + R.L_foot(3,2,i)^2));
    % angle of the long axis of the foot about the vertical axis
    angle.L_foot(i,2) = atan(-R.L_foot(1,1,i)/sqrt(R.L_foot(2,1,i)^2 + R.L_foot(3,1,i)^2));
    % and project the long axis into the sagital plane for ID of FOREFOOT
    angle.L_foot(i,3) = atan2(R.L_foot(2,1,i),-R.L_foot(3,1,i));
    
    
    % angle of the vertical axis of the foot projected into the frontal plane
    %   from a posterior view, 'vertical vector' in the first quadrant is
    %   postitive and second quadrant is negative
    %   . ie. inversion is positive and eversion is negative for the RIGHT
    angle.R_foot(i,1) = atan(R.R_foot(1,2,i)/sqrt(R.R_foot(2,2,i)^2 + R.R_foot(3,2,i)^2));
    % angle of the long axis of the foot about the vertical axis
    angle.R_foot(i,2) = atan(R.R_foot(1,1,i)/sqrt(R.R_foot(2,1,i)^2 + R.R_foot(3,1,i)^2));
    % and project the long axis into the sagital plane for ID of FOREFOOT
    angle.R_foot(i,3) = atan2(R.R_foot(2,1,i),-R.R_foot(3,1,i));
    
    
    % PELVIS
    % project lateral axis of pelvis
    angle.pelvis(i,1) = atan2(R.pelvis(1,3,i),R.pelvis(3,3,i)) -pi/2; %into floor plane
    angle.pelvis(i,2) = atan(R.pelvis(2,3,i)/R.pelvis(1,3,i)); %into frontal plane
    % and project anterior axis of pelvis
    angle.pelvis(i,3) = atan2(R.pelvis(2,1,i),-R.pelvis(3,1,i)); %into sagital plane
    
    %% calculate joint angles
    % need rotation matrix from shank to foot ...
    % so multiply [lab to shank] with [foot to lab]
    
    R.L_ankle(:,:,i) = R.L_shank(:,:,i)'*R.L_foot(:,:,i);
    R.R_ankle(:,:,i) = R.R_shank(:,:,i)'*R.R_foot(:,:,i);
    R.L_knee(:,:,i) = R.L_thigh(:,:,i)'*R.L_shank(:,:,i);
    R.R_knee(:,:,i) = R.R_thigh(:,:,i)'*R.R_shank(:,:,i);
    R.L_hip(:,:,i) = R.pelvis(:,:,i)'*R.L_thigh(:,:,i);
    R.R_hip(:,:,i) = R.pelvis(:,:,i)'*R.R_thigh(:,:,i);
    
    % CARDANANGLES uses this rotation matrix to calculate angles
    %   | CzCy-SzSySx  SzCy+CzSySx  -SyCx |
    %   | -SzCx        CzCx         Sx    |
    %   | CzSy+SzCySx  SzSy-CzCySx  CyCx  |
    
    angle.L_ankle(i,:) = cardanangles(R.L_ankle(:,:,i));
    angle.L_knee(i,:)  = cardanangles(R.L_knee(:,:,i));
    angle.L_hip(i,:)   = cardanangles(R.L_hip(:,:,i));
    angle.R_ankle(i,:) = cardanangles(R.R_ankle(:,:,i));
    angle.R_knee(i,:)  = cardanangles(R.R_knee(:,:,i));
    angle.R_hip(i,:)   = cardanangles(R.R_hip(:,:,i));
    
    % start by calculating the x angles which are about the anterior posterior
    % axes. sine of a angle maintains the sign from -90 to +90 which should be
    % well within the physiological range of movement during gait
    
    % Use atan2(a,b) which is stable from -180 to +180. need atan2(sin,cos).
    % Use elements (1,3)/cos(x) and (3,3)/cos(x) for the y axis and
    % (2,1)/cos(x) and (2,2)/cos(x) for the z axis
    
    
end


%% convert angles from radians to degrees
angles.L_ankle = angle.L_ankle * (180/pi);
angles.R_ankle = angle.R_ankle * (180/pi);
angles.L_knee = angle.L_knee * (180/pi);
angles.R_knee = angle.R_knee * (180/pi);
angles.L_hip = angle.L_hip * (180/pi);
angles.R_hip = angle.R_hip * (180/pi);
angles.L_foot = angle.L_foot * (180/pi);
angles.R_foot = angle.R_foot * (180/pi);
angles.pelvis = angle.pelvis * (180/pi);
%% calculate joint angular velocities
velocities.L_ankle=diff(angles.L_ankle)*hz;
velocities.R_ankle=diff(angles.R_ankle)*hz;
velocities.L_knee=diff(angles.L_knee)*hz;
velocities.R_knee=diff(angles.R_knee)*hz;
velocities.L_hip=diff(angles.L_hip)*hz;
velocities.R_hip=diff(angles.R_hip)*hz;
velocities.L_foot=diff(angles.L_foot)*hz;
velocities.R_foot=diff(angles.R_foot)*hz;
velocities.pelvis=diff(angles.pelvis)*hz;


%% plot the angles and velocities

if plots ~= 0
    
    figure('tag', 'bonita_kinematics_temp_figure'); hold on;
    subplot(321);plot(angles.L_ankle(:,1),'r');title('L ankle angle x')
    subplot(323);plot(angles.L_ankle(:,2),'g');title('L ankle angle y')
    subplot(325);plot(angles.L_ankle(:,3),'b');title('L ankle angle z')
    subplot(322);plot(angles.R_ankle(:,1),'r');title('R ankle angle x')
    subplot(324);plot(angles.R_ankle(:,2),'g');title('R ankle angle y')
    subplot(326);plot(angles.R_ankle(:,3),'b');title('R ankle angle z')
    
    figure('tag', 'bonita_kinematics_temp_figure'); hold on;
    subplot(321);plot(angles.L_knee(:,1),'r');title('L knee angle x')
    subplot(323);plot(angles.L_knee(:,2),'g');title('L knee angle y')
    subplot(325);plot(angles.L_knee(:,3),'b');title('L knee angle z')
    subplot(322);plot(angles.R_knee(:,1),'r');title('R knee angle x')
    subplot(324);plot(angles.R_knee(:,2),'g');title('R knee angle y')
    subplot(326);plot(angles.R_knee(:,3),'b');title('R knee angle z')
    
    figure('tag', 'bonita_kinematics_temp_figure'); hold on;
    subplot(321);plot(angles.L_hip(:,1),'r');title('L hip angle x')
    subplot(323);plot(angles.L_hip(:,2),'g');title('L hip angle y')
    subplot(325);plot(angles.L_hip(:,3),'b');title('L hip angle z')
    subplot(322);plot(angles.R_hip(:,1),'r');title('R hip angle x')
    subplot(324);plot(angles.R_hip(:,2),'g');title('R hip angle y')
    subplot(326);plot(angles.R_hip(:,3),'b');title('R hip angle z')
    
    figure('tag', 'bonita_kinematics_temp_figure'); hold on;
    subplot(321);plot(angles.L_foot(:,1),'r');title('L foot angle project long axis onto floor')
    subplot(323);plot(angles.L_foot(:,2),'g');title('L foot angle project long axis onto frontal')
    subplot(325);plot(angles.L_foot(:,3),'b');title('L foot angle project long axis onto sagital')
    subplot(322);plot(angles.R_foot(:,1),'r');title('R foot angle project longs axis onto floor')
    subplot(324);plot(angles.R_foot(:,2),'g');title('R foot angle project long axis onto frontal')
    subplot(326);plot(angles.R_foot(:,3),'b');title('R foot angle project long axis onto sagital')
    
    figure('tag', 'bonita_kinematics_temp_figure'); hold on;
    subplot(311);plot(angles.pelvis(:,1),'r');title('Lateral pevlis axis projected into the floor')
    subplot(312);plot(angles.pelvis(:,2),'g');title('Lateral pevlis axis projected into the frontal plane')
    subplot(313);plot(angles.pelvis(:,3),'b');title('Lateral pevlis axis projected into the sagital plane')
    
    
end


%% plot coordinate systems with markers over time
if plots ~=0
    
    figure('tag', 'bonita_kinematics_temp_figure'); hold on; axis equal;
    
    rotate3d on
    xlabel('z'); ylabel('x'); zlabel('y')
    
    line_length = 100;
    
    for i=1:length(dynamic.L_shank_1(:,1))
        delete(get(gca,'Children'));
        
        view(i/5,20);
        
        plot3(dynamic.pelvis_1(i,3),dynamic.pelvis_1(i,1),dynamic.pelvis_1(i,2),'b.')
        plot3(dynamic.pelvis_2(i,3),dynamic.pelvis_2(i,1),dynamic.pelvis_2(i,2),'r.')
        plot3(dynamic.pelvis_3(i,3),dynamic.pelvis_3(i,1),dynamic.pelvis_3(i,2),'g.')
        plot3(dynamic.pelvis_4(i,3),dynamic.pelvis_4(i,1),dynamic.pelvis_4(i,2),'c.')
        temp=R.pelvis(1:3,1:3,i)*djc.pelvis + dynamic.pelvis_1(i,:)'; %plot3(temp(3),temp(1),temp(2),'*')
        line([temp(3),temp(3)+R.pelvis(3,1,i)*line_length], [temp(1),temp(1)+R.pelvis(1,1,i)*line_length], [temp(2),temp(2)+R.pelvis(2,1,i)*line_length],'color','r')
        line([temp(3),temp(3)+R.pelvis(3,2,i)*line_length], [temp(1),temp(1)+R.pelvis(1,2,i)*line_length], [temp(2),temp(2)+R.pelvis(2,2,i)*line_length],'color','g')
        line([temp(3),temp(3)+R.pelvis(3,3,i)*line_length], [temp(1),temp(1)+R.pelvis(1,3,i)*line_length], [temp(2),temp(2)+R.pelvis(2,3,i)*line_length],'color','b')
        
        plot3(dynamic.L_thigh_1(i,3),dynamic.L_thigh_1(i,1),dynamic.L_thigh_1(i,2),'b.')
        plot3(dynamic.L_thigh_2(i,3),dynamic.L_thigh_2(i,1),dynamic.L_thigh_2(i,2),'r.')
        plot3(dynamic.L_thigh_3(i,3),dynamic.L_thigh_3(i,1),dynamic.L_thigh_3(i,2),'g.')
        plot3(dynamic.L_thigh_4(i,3),dynamic.L_thigh_4(i,1),dynamic.L_thigh_4(i,2),'c.')
        temp=R.L_thigh(1:3,1:3,i)*djc.L_hip + dynamic.L_thigh_1(i,:)'; %plot3(temp(3),temp(1),temp(2),'*')
        line([temp(3),temp(3)+R.L_thigh(3,1,i)*line_length], [temp(1),temp(1)+R.L_thigh(1,1,i)*line_length], [temp(2),temp(2)+R.L_thigh(2,1,i)*line_length],'color','r')
        line([temp(3),temp(3)+R.L_thigh(3,2,i)*line_length], [temp(1),temp(1)+R.L_thigh(1,2,i)*line_length], [temp(2),temp(2)+R.L_thigh(2,2,i)*line_length],'color','g')
        line([temp(3),temp(3)+R.L_thigh(3,3,i)*line_length], [temp(1),temp(1)+R.L_thigh(1,3,i)*line_length], [temp(2),temp(2)+R.L_thigh(2,3,i)*line_length],'color','b')
        
        
        plot3(dynamic.R_thigh_1(i,3),dynamic.R_thigh_1(i,1),dynamic.R_thigh_1(i,2),'b.')
        plot3(dynamic.R_thigh_2(i,3),dynamic.R_thigh_2(i,1),dynamic.R_thigh_2(i,2),'r.')
        plot3(dynamic.R_thigh_3(i,3),dynamic.R_thigh_3(i,1),dynamic.R_thigh_3(i,2),'g.')
        plot3(dynamic.R_thigh_4(i,3),dynamic.R_thigh_4(i,1),dynamic.R_thigh_4(i,2),'c.')
        temp=R.R_thigh(1:3,1:3,i)*djc.R_hip + dynamic.R_thigh_1(i,:)'; %plot3(temp(3),temp(1),temp(2),'*')
        line([temp(3),temp(3)+R.R_thigh(3,1,i)*line_length], [temp(1),temp(1)+R.R_thigh(1,1,i)*line_length], [temp(2),temp(2)+R.R_thigh(2,1,i)*line_length],'color','r')
        line([temp(3),temp(3)+R.R_thigh(3,2,i)*line_length], [temp(1),temp(1)+R.R_thigh(1,2,i)*line_length], [temp(2),temp(2)+R.R_thigh(2,2,i)*line_length],'color','g')
        line([temp(3),temp(3)+R.R_thigh(3,3,i)*line_length], [temp(1),temp(1)+R.R_thigh(1,3,i)*line_length], [temp(2),temp(2)+R.R_thigh(2,3,i)*line_length],'color','b')
        
        
        plot3(dynamic.L_shank_1(i,3),dynamic.L_shank_1(i,1),dynamic.L_shank_1(i,2),'b.')
        plot3(dynamic.L_shank_2(i,3),dynamic.L_shank_2(i,1),dynamic.L_shank_2(i,2),'r.')
        plot3(dynamic.L_shank_3(i,3),dynamic.L_shank_3(i,1),dynamic.L_shank_3(i,2),'g.')
        plot3(dynamic.L_shank_4(i,3),dynamic.L_shank_4(i,1),dynamic.L_shank_4(i,2),'c.')
        temp=R.L_shank(1:3,1:3,i)*djc.L_knee + dynamic.L_shank_1(i,:)'; %plot3(temp(3),temp(1),temp(2),'*')
        line([temp(3),temp(3)+R.L_shank(3,1,i)*line_length], [temp(1),temp(1)+R.L_shank(1,1,i)*line_length], [temp(2),temp(2)+R.L_shank(2,1,i)*line_length],'color','r')
        line([temp(3),temp(3)+R.L_shank(3,2,i)*line_length], [temp(1),temp(1)+R.L_shank(1,2,i)*line_length], [temp(2),temp(2)+R.L_shank(2,2,i)*line_length],'color','g')
        line([temp(3),temp(3)+R.L_shank(3,3,i)*line_length], [temp(1),temp(1)+R.L_shank(1,3,i)*line_length], [temp(2),temp(2)+R.L_shank(2,3,i)*line_length],'color','b')
        
        
        
        
        plot3(dynamic.R_shank_1(i,3),dynamic.R_shank_1(i,1),dynamic.R_shank_1(i,2),'b.')
        plot3(dynamic.R_shank_2(i,3),dynamic.R_shank_2(i,1),dynamic.R_shank_2(i,2),'r.')
        plot3(dynamic.R_shank_3(i,3),dynamic.R_shank_3(i,1),dynamic.R_shank_3(i,2),'g.')
        plot3(dynamic.R_shank_4(i,3),dynamic.R_shank_4(i,1),dynamic.R_shank_4(i,2),'c.')
        temp=R.R_shank(1:3,1:3,i)*djc.R_knee + dynamic.R_shank_1(i,:)'; %plot3(temp(3),temp(1),temp(2),'*')
        line([temp(3),temp(3)+R.R_shank(3,1,i)*line_length], [temp(1),temp(1)+R.R_shank(1,1,i)*line_length], [temp(2),temp(2)+R.R_shank(2,1,i)*line_length],'color','r')
        line([temp(3),temp(3)+R.R_shank(3,2,i)*line_length], [temp(1),temp(1)+R.R_shank(1,2,i)*line_length], [temp(2),temp(2)+R.R_shank(2,2,i)*line_length],'color','g')
        line([temp(3),temp(3)+R.R_shank(3,3,i)*line_length], [temp(1),temp(1)+R.R_shank(1,3,i)*line_length], [temp(2),temp(2)+R.R_shank(2,3,i)*line_length],'color','b')
        
        
        plot3(dynamic.L_foot_1(i,3),dynamic.L_foot_1(i,1),dynamic.L_foot_1(i,2),'b.')
        plot3(dynamic.L_foot_2(i,3),dynamic.L_foot_2(i,1),dynamic.L_foot_2(i,2),'r.')
        plot3(dynamic.L_foot_3(i,3),dynamic.L_foot_3(i,1),dynamic.L_foot_3(i,2),'g.')
        plot3(dynamic.L_foot_4(i,3),dynamic.L_foot_4(i,1),dynamic.L_foot_4(i,2),'c.')
        temp=R.L_foot(1:3,1:3,i)*djc.L_ankle + dynamic.L_foot_1(i,:)'; %plot3(temp(3),temp(1),temp(2),'*')
        line([temp(3),temp(3)+R.L_foot(3,1,i)*line_length], [temp(1),temp(1)+R.L_foot(1,1,i)*line_length], [temp(2),temp(2)+R.L_foot(2,1,i)*line_length],'color','r')
        line([temp(3),temp(3)+R.L_foot(3,2,i)*line_length], [temp(1),temp(1)+R.L_foot(1,2,i)*line_length], [temp(2),temp(2)+R.L_foot(2,2,i)*line_length],'color','g')
        line([temp(3),temp(3)+R.L_foot(3,3,i)*line_length], [temp(1),temp(1)+R.L_foot(1,3,i)*line_length], [temp(2),temp(2)+R.L_foot(2,3,i)*line_length],'color','b')
        
        
        plot3(dynamic.R_foot_1(i,3),dynamic.R_foot_1(i,1),dynamic.R_foot_1(i,2),'b.')
        plot3(dynamic.R_foot_2(i,3),dynamic.R_foot_2(i,1),dynamic.R_foot_2(i,2),'r.')
        plot3(dynamic.R_foot_3(i,3),dynamic.R_foot_3(i,1),dynamic.R_foot_3(i,2),'g.')
        plot3(dynamic.R_foot_4(i,3),dynamic.R_foot_4(i,1),dynamic.R_foot_4(i,2),'c.')
        temp=R.R_foot(1:3,1:3,i)*djc.R_ankle + dynamic.R_foot_1(i,:)'; %plot3(temp(3),temp(1),temp(2),'*')
        line([temp(3),temp(3)+R.R_foot(3,1,i)*line_length], [temp(1),temp(1)+R.R_foot(1,1,i)*line_length], [temp(2),temp(2)+R.R_foot(2,1,i)*line_length],'color','r')
        line([temp(3),temp(3)+R.R_foot(3,2,i)*line_length], [temp(1),temp(1)+R.R_foot(1,2,i)*line_length], [temp(2),temp(2)+R.R_foot(2,2,i)*line_length],'color','g')
        line([temp(3),temp(3)+R.R_foot(3,3,i)*line_length], [temp(1),temp(1)+R.R_foot(1,3,i)*line_length], [temp(2),temp(2)+R.R_foot(2,3,i)*line_length],'color','b')
        
        drawnow
    end
    
end

function out = cardanangles(r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function inputs the ROTATION matrix for a given joint
%   at one point in time and uses the cardan angles sequence to
%   provide the XYZ joint angles in OUT argument
%
%   The function uses the following  rotation matrix to calculate angles
%      | CzCy-SzSySx  SzCy+CzSySx  -SyCx |
%      | -SzCx        CzCx         Sx    |
%      | CzSy+SzCySx  SzSy-CzCySx  CyCx  |
%  INPUTS
%  --------
%   R (mat):    A 3x3 rotation matrix for a joint

%  OUTPUTS
%  -------
%  OUT (mat):      The three (1x3) planes of rotation in radians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%   angle.L_ankle(i,:) = cardanangles(R.L_ankle(:,:,i));

% the use of atan2 increases the stability by avoiding gimble lock in
% the physiologicaly posible range of joint angles.
x = atan2(r(2,3), sqrt(r(1,3)^2+r(3,3)^2));
y = atan2(-r(1,3), r(3,3));
z = atan2(-r(2,1), r(2,2));
out = [x,y,z];
