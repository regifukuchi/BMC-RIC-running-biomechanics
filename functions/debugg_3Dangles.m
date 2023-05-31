%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT processing_code_example.m
% 
%
% Wrapper script for processing data through the Running Injury Clinic
% pipeline
%
% Main functions in wrapper: GAIT_KINEMATICS, GAIT_STEPS. Please open
% these functions for detailed description of INPUTS and OUTPUTS
%
%  LICENSE
%  -------
%  See file LICENSE.txt
%
% Copyright (C) 2022-2023 Allan Brett and the Running Injury Clinic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all
code_folder = 'C:\Users\Reginaldo\Documents\Github\RIC_3Dgait\functions';
%% get user directories

if 0 % RIC dataset
    dataset='RIC';
    %code_folder = uigetdir('C:\','Select folder containing PROCESSING code');
    
    %data_folder = uigetdir('C:\','Select folder containing JSON data');
    data_folder = 'C:\Users\Reginaldo\OneDrive - University of Calgary\data\Figshare_SciData\new_unzip';

    % get fully defined path to json data file
    json_file = fullfile(data_folder, '201225', '20140515T133244.json');

    %load json file
    fid = fopen(json_file);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    out = jsondecode(str);

    % IMPORTANT!%
    %writing to json does not faithfully recreat the structure of
    %out.joints, and out.neutral as stored in original .MAT FILE.
    %Must reformat beforehand for INPUT into bonita_kinematics and
    %bonita_steps

    fields = fieldnames(out.joints);

    for j = 1:size(fields,1)

        out.joints.(fields{j,1}) = transpose(out.joints.(fields{j,1}));

    end

    clear fields
    fields = fieldnames(out.neutral);

    for j = 1:size(fields,1)

        out.neutral.(fields{j,1}) = transpose(out.neutral.(fields{j,1}));

    end
else
    dataset='RBDS';
    data_folder = 'C:\Users\Reginaldo\Documents\data\CNPq\RBDS_v2\Figshare';
    % Neutral
    neutral_lbls = {'pelvis_1', 'pelvis_2', 'pelvis_3', 'pelvis_4', 'L_lat_ankle',...
        'L_med_ankle', 'L_hip', 'L_foot_2', 'L_foot_3', 'L_foot_1', 'L_lat_knee',...
        'L_med_knee', 'L_shank_1', 'L_shank_2', 'L_shank_3', 'L_shank_4',...
        'L_thigh_1', 'L_thigh_2', 'L_thigh_3', 'L_thigh_4', 'R_lat_ankle',...
        'R_med_ankle', 'R_hip', 'R_foot_2', 'R_foot_3', 'R_foot_1',...
        'R_lat_knee', 'R_med_knee', 'R_shank_1', 'R_shank_2', 'R_shank_3', 'R_shank_4',...
        'R_thigh_1', 'R_thigh_2', 'R_thigh_3', 'R_thigh_4', 'L_foot_4', 'R_foot_4'};
    
    x = importdata(fullfile(data_folder,'RBDS001_neutral.csv'));
    x = x.data(2:end);
    
    for i = 1:length(neutral_lbls)
        neutral.(neutral_lbls{i}) = x(3*i-2:3*i);
    end
    % Joints
    joints_lbls = {'L_hip', 'R_hip', 'L_lat_knee', 'L_med_knee', 'R_lat_knee', 'R_med_knee', 'L_lat_ankle', 'L_med_ankle', 'R_lat_ankle', 'R_med_ankle'};
    jt = importdata(fullfile(data_folder,'RBDS001_joints.csv'));
    jt = jt.data(2:end);
    
    for i = 1:length(joints_lbls)
        joints.(joints_lbls{i}) = jt(3*i-2:3*i);
    end
    % Gait
    gait_lbls={'pelvis_1', 'pelvis_2', 'pelvis_3', 'pelvis_4', 'L_thigh_1', 'L_thigh_2', 'L_thigh_3', 'L_thigh_4', 'R_thigh_1', 'R_thigh_2', 'R_thigh_3', 'R_thigh_4', 'L_shank_1', 'L_shank_2', 'L_shank_3', 'L_shank_4', 'R_shank_1', 'R_shank_2', 'R_shank_3', 'R_shank_4', 'L_foot_1', 'L_foot_2', 'L_foot_3', 'R_foot_1', 'R_foot_2', 'R_foot_3', 'L_foot_4', 'R_foot_4'};
    gt = importdata(fullfile(data_folder,'RBDS001_gait.csv'));
    gt = gt.data(:,2:end);
    
    for i = 1:length(gait_lbls)
        gait.(gait_lbls{i}) = gt(:,3*i-2:3*i);
    end
    % Other input params
    out.joints = joints;
    out.neutral = neutral;
    out.running = gait;
    out.hz_r=150;
    plot=0;
end

%% Calculate angles using gait_kinematics.m
[angles,velocities,jc,R,djc] = gait_kinematics(out.joints,out.neutral,out.running,out.hz_r,plot);
%% Calculate gait_steps
[norm_ang,norm_vel,events,event,DISCRETE_VARIABLES,speedoutput,eventsflag,label] = gait_steps(out.neutral,out.running,angles,velocities,out.hz_r,plot);
%% Export data
save(['events_' dataset '.txt'],'events', '-ascii','-tabs')
save(['eventsflag_' dataset '.txt'],'eventsflag', '-ascii','-tabs')
%%
pelvis_ang = angles.pelvis;
save(['pelvis_angle_' dataset '.txt'],'pelvis_ang', '-ascii','-tabs')
foot_ang_R = angles.R_foot;
save(['r_foot_angle_' dataset 'mat.txt'],'foot_ang_R', '-ascii','-tabs')
ankle_ang_R = angles.R_ankle;
save(['r_ankle_angle_' dataset 'mat.txt'],'ankle_ang_R', '-ascii','-tabs')
knee_ang_R = angles.R_knee;
save(['r_knee_angle_' dataset 'mat.txt'],'knee_ang_R', '-ascii','-tabs')
hip_ang_R = angles.R_hip;
save(['r_hip_angle_' dataset 'mat.txt'],'hip_ang_R', '-ascii','-tabs')

foot_ang_L = angles.L_foot;
save(['L_foot_angle_' dataset 'mat.txt'],'foot_ang_L', '-ascii','-tabs')
ankle_ang_L = angles.L_ankle;
save(['L_ankle_angle_' dataset 'mat.txt'],'ankle_ang_L', '-ascii','-tabs')
knee_ang_L = angles.L_knee;
save(['L_knee_angle_' dataset 'mat.txt'],'knee_ang_L', '-ascii','-tabs')
hip_ang_L = angles.L_hip;
save(['L_hip_angle_' dataset 'mat.txt'],'hip_ang_L', '-ascii','-tabs')
%% INPUT PARAMS
% gait_kinematics(out.joints,out.neutral,out.walking,out.hz_w,plots)
joints  = out.joints;
neutral = out.neutral;
dynamic = out.running;
% JOINT CENTRES
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
% CALCULATE FOOT SEGMENT %%
% midpoint of the two ankle markers
jc.L_ankle = (joints.L_lat_ankle+joints.L_med_ankle)/2;
jc.R_ankle = (joints.R_med_ankle+joints.R_lat_ankle)/2;
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
if 0
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
cs_length = 50;
if 0
    % define how long to make the CS lines if there are being plotted
    
    rfx=[jc.R_ankle;(jc.R_ankle+cs_length*r_foot_x)];line(rfx(:,1),rfx(:,2),rfx(:,3),'Color','r')
    rfy=[jc.R_ankle;(jc.R_ankle+cs_length*r_foot_y)];line(rfy(:,1),rfy(:,2),rfy(:,3),'Color','g')
    rfz=[jc.R_ankle;(jc.R_ankle+cs_length*r_foot_z)];line(rfz(:,1),rfz(:,2),rfz(:,3),'Color','b')
end
%% LEFT SHANK
l_shank_y=(jc.L_knee-jc.L_ankle)/norm(jc.L_knee-jc.L_ankle); %long axis pointing up
l_shank_temp=(joints.L_med_ankle-joints.L_lat_ankle)/norm(joints.L_med_ankle-joints.L_lat_ankle); %almost the hinge joing pointing to the right
l_shank_x=(cross(l_shank_y,l_shank_temp))/norm(cross(l_shank_y,l_shank_temp)); %Anterior axis from the cross
l_shank_z=(cross(l_shank_x,l_shank_y))/norm(cross(l_shank_x,l_shank_y)); %hinge axis, points to the right
% combine to create a transformation matrix from anatomical to global
ag.L_shank=[l_shank_x',l_shank_y',l_shank_z'];
% visualize
if 0
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
if 0
    rsx=[jc.R_knee;(jc.R_knee+cs_length*r_shank_x)];line(rsx(:,1),rsx(:,2),rsx(:,3),'Color','r')
    rsy=[jc.R_knee;(jc.R_knee+cs_length*r_shank_y)];line(rsy(:,1),rsy(:,2),rsy(:,3),'Color','g')
    rsz=[jc.R_knee;(jc.R_knee+cs_length*r_shank_z)];line(rsz(:,1),rsz(:,2),rsz(:,3),'Color','b')
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
if 0
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
if 0
    rtx=[jc.R_hip;(jc.R_hip+cs_length*r_thigh_x)];line(rtx(:,1),rtx(:,2),rtx(:,3),'Color','r')
    rty=[jc.R_hip;(jc.R_hip+cs_length*r_thigh_y)];line(rty(:,1),rty(:,2),rty(:,3),'Color','g')
    rtz=[jc.R_hip;(jc.R_hip+cs_length*r_thigh_z)];line(rtz(:,1),rtz(:,2),rtz(:,3),'Color','b')
else
end
% PELVIS
% since trochanters are hard to landmark, pelvis will just be orthogonal to the lab
pelvis_x =[0 0 -1]; %anterior axis
pelvis_y =[0 1 0]; %long axis pointing up
pelvis_z =[1 0 0]; %hinge axis to the subject's right
% combine to create a transformation matrix from anatomical to global
ag.pelvis=[pelvis_x',pelvis_y',pelvis_z'];
% visualize
if 0
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

% the LEFT foot markers in an ANATOMICAL COORDINATE system ...
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

% the RIGHT foot markers in an ANATOMICAL COORDINATE system ...
ac.R_foot_1=ag.R_foot'*neutral.R_foot_1';
ac.R_foot_2=ag.R_foot'*neutral.R_foot_2';
ac.R_foot_3=ag.R_foot'*neutral.R_foot_3';
ac.R_foot_4=ag.R_foot'*neutral.R_foot_4';
ac.R_foot=[ac.R_foot_1,ac.R_foot_2,ac.R_foot_3,ac.R_foot_4];

centroid = mean(ac.R_foot,2);

djc.R_ankle=ag.R_foot'*jc.R_ankle'-centroid;
%% the LEFT shank markers in an ANATOMICAL COORDINATE system ...
ac.L_shank_1=ag.L_shank'*neutral.L_shank_1';
ac.L_shank_2=ag.L_shank'*neutral.L_shank_2';
ac.L_shank_3=ag.L_shank'*neutral.L_shank_3';
ac.L_shank_4=ag.L_shank'*neutral.L_shank_4';
ac.L_shank=[ac.L_shank_1,ac.L_shank_2,ac.L_shank_3,ac.L_shank_4];

centroid = mean(ac.L_shank,2);

djc.L_knee=ag.L_shank'*jc.L_knee'-centroid;

%% the RIGHT shank markers in an ANATOMICAL COORDINATE system ...
ac.R_shank_1=ag.R_shank'*neutral.R_shank_1';
ac.R_shank_2=ag.R_shank'*neutral.R_shank_2';
ac.R_shank_3=ag.R_shank'*neutral.R_shank_3';
ac.R_shank_4=ag.R_shank'*neutral.R_shank_4';
ac.R_shank=[ac.R_shank_1,ac.R_shank_2,ac.R_shank_3,ac.R_shank_4];

centroid = mean(ac.R_shank,2);

djc.R_knee=ag.R_shank'*jc.R_knee'-centroid;
%% the LEFT Thigh markers in an ANATOMICAL COORDINATE system ...
ac.L_thigh_1=ag.L_thigh'*neutral.L_thigh_1';
ac.L_thigh_2=ag.L_thigh'*neutral.L_thigh_2';
ac.L_thigh_3=ag.L_thigh'*neutral.L_thigh_3';
ac.L_thigh_4=ag.L_thigh'*neutral.L_thigh_4';
ac.L_thigh=[ac.L_thigh_1,ac.L_thigh_2,ac.L_thigh_3,ac.L_thigh_4];

centroid = mean(ac.L_thigh,2);

djc.L_hip=ag.L_thigh'*jc.L_hip'-centroid;

%% the RIGHT Thigh markers in an ANATOMICAL COORDINATE system ...
ac.R_thigh_1=ag.R_thigh'*neutral.R_thigh_1';
ac.R_thigh_2=ag.R_thigh'*neutral.R_thigh_2';
ac.R_thigh_3=ag.R_thigh'*neutral.R_thigh_3';
ac.R_thigh_4=ag.R_thigh'*neutral.R_thigh_4';
ac.R_thigh=[ac.R_thigh_1,ac.R_thigh_2,ac.R_thigh_3,ac.R_thigh_4];

centroid = mean(ac.R_thigh,2);

djc.R_hip=ag.R_thigh'*jc.R_hip'-centroid;

% the Pelvis markers in an ANATOMICAL COORDINATE system ...
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
for i=1:length(dynamic.R_foot_1(:,1))
    
    % need to create matricies from the DYNAMIC data that has the same markers
     % as the " ac.L_foot " in the same order ... matricies for use with soderqvist
    
    d.L_foot=[dynamic.L_foot_1(i,:)',dynamic.L_foot_2(i,:)',dynamic.L_foot_3(i,:)',dynamic.L_foot_4(i,:)'];
    d.R_foot=[dynamic.R_foot_1(i,:)',dynamic.R_foot_2(i,:)',dynamic.R_foot_3(i,:)',dynamic.R_foot_4(i,:)'];
    d.L_shank=[dynamic.L_shank_1(i,:)',dynamic.L_shank_2(i,:)',dynamic.L_shank_3(i,:)',dynamic.L_shank_4(i,:)'];
    d.R_shank=[dynamic.R_shank_1(i,:)',dynamic.R_shank_2(i,:)',dynamic.R_shank_3(i,:)',dynamic.R_shank_4(i,:)'];
    d.L_thigh=[dynamic.L_thigh_1(i,:)',dynamic.L_thigh_2(i,:)',dynamic.L_thigh_3(i,:)',dynamic.L_thigh_4(i,:)'];
    d.R_thigh=[dynamic.R_thigh_1(i,:)',dynamic.R_thigh_2(i,:)',dynamic.R_thigh_3(i,:)',dynamic.R_thigh_4(i,:)'];
    d.pelvis=[dynamic.pelvis_1(i,:)',dynamic.pelvis_2(i,:)',dynamic.pelvis_3(i,:)',dynamic.pelvis_4(i,:)'];
    
	% calculate average position of the segment markers
    
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

%% Compare output angles from debugg_3Dangles.m and gait_kinematics.m
% Import files
file_dir = 'C:\Users\Reginaldo\Documents\Github\UFABC_UofC_datasets\data';
ankle_ang = importdata(fullfile(file_dir, 'r_ankle_angle_R.txt'));
knee_ang  = importdata(fullfile(file_dir, 'r_knee_angle_R.txt'));
hip_ang   = importdata(fullfile(file_dir, 'r_hip_angle_R.txt'));
%% Plot comparison
figure % hip angles
subplot(3,1,1)
plot(angles.R_hip(:,1)), hold on
plot(hip_ang(:,1)*(180/pi),'o')
subplot(3,1,2)
plot(angles.R_hip(:,2)), hold on
plot(hip_ang(:,2)*(180/pi),'o')
subplot(3,1,3)
plot(angles.R_hip(:,3)), hold on
plot(hip_ang(:,3)*(180/pi),'o')
%% Knee angles
figure 
subplot(3,1,1)
plot(angles.R_knee(:,1)), hold on
plot(knee_ang(:,1)*(180/pi),'o')
subplot(3,1,2)
plot(angles.R_knee(:,2)), hold on
plot(knee_ang(:,2)*(180/pi),'o')
subplot(3,1,3)
plot(angles.R_knee(:,3)), hold on
plot(knee_ang(:,3)*(180/pi),'o')
%% Ankle angles
figure 
subplot(3,1,1)
plot(angles.R_ankle(:,1)), hold on
plot(ankle_ang(:,1)*(180/pi),'o')
subplot(3,1,2)
plot(angles.R_ankle(:,2)), hold on
plot(ankle_ang(:,2)*(180/pi),'o')
subplot(3,1,3)
plot(angles.R_ankle(:,3)), hold on
plot(ankle_ang(:,3)*(180/pi),'o')