% This script exports processed data (eg joint angles, torques, etc) from
% RBA study as ASCII files
% The ASCII files in Visual 3D were exported separetely for each variables.
% Here the time series (angles, moments, forces and powers) are
% concatednated in a large matrix with 101 rows and 144 columns
% The columns comprises:
% 9 joint angles (3 joints x 3 directions)
% 9 joint moments (3 joints x 3 directions)
% 3 GRFs
% 3 Joint powers
% So there are, 24 variables by 3 gait speeds and 2 legs so 144 curves.


clear all, close all, clc
% Only subjects with good data including 3 gait speeds so far
subjs = [8:10 12:18 20 23:26 27:30 32:40];

txtsave = 0; % Determine whether one wants to save ASCII files. Default should be set to zero

speed = [2.5,3.5,4.5];
side   = {'R','L'};
joints = {'hip','knee','ankle'};
varType = {'Ang','Mom','grf','Pow'};
varListRBA = {'Stride_Length','Strides_Per_Minute','Stride_Width'};

direction = {'X','Y','Z'};

headerProcessRBA = {};
% Loading data
rbaM = load('H:\RBA\Data\rbaM2.mat');

pathDir = 'H:\RBA\Data\RAWall';
%%

for is = 1:length(subjs)
    if is < 10
        subjID = ['RBDS00' num2str(is)];
    else
        subjID = ['RBDS0' num2str(is)];
    end
    
    rbaDataI = [];
    for s = 1:length(side)
        for gs = 1:3
            for ivar = 1:length(varType)
                for j = 1:length(joints)
                    
                    if ivar < 3
                        % Angles and Torques
                        xX = rbaM.([side{s} joints{j} varType{ivar} 'RUN' num2str(speed(gs)*10)])(:,3*is-2:3*is);
                        
                        headerProcessRBA{(3*ivar-3)+j} = [joints{j} varType{ivar}];
                        
                        % Concatenating data
                        rbaDataI = [rbaDataI xX];
                        
                    elseif ivar==3 && j==1
                        % Ground reaction forces
                        xX = rbaM.([side{s} varType{ivar} 'RUN' num2str(speed(gs)*10)])(:,3*is-2:3*is);
                        
                        headerProcessRBA{7} = varType{ivar};
                        
                        % Concatenating data
                        rbaDataI = [rbaDataI xX];
                    elseif ivar==4
                        % Joint Powers
                        xX = rbaM.([side{s} joints{j} varType{ivar} 'RUN' num2str(speed(gs)*10)])(:,3*is-2);
                        
                        headerProcessRBA{7+j} = [joints{j} varType{ivar}];
                        
                        % Concatenating data
                        rbaDataI = [rbaDataI xX];
                    end
                    
                end
            end
            
            for iheader = 1:length(headerProcessRBA)
                if iheader < 8
                    for idir = 1:3
                        % Header for angles, moments and GRF
                        headerProcessRBAside{(72*s-72)+(24*gs-24)+((3*iheader-3)+idir)} = [side{s} headerProcessRBA{iheader} direction{idir} num2str(speed(gs)*10)];
                    end
                else
                    % Header for joint powers
                    headerProcessRBAside{(72*s-72)+(24*gs-24)+14+iheader} = [side{s} headerProcessRBA{iheader} num2str(speed(gs)*10)];
                end
                
            end
        end
    end
    rbaData(:,:,is) = rbaDataI;
    
    
    
    if txtsave
        percGaitCycle = 0:100;
        rbaDataI = [percGaitCycle' rbaDataI];
        fnamP = [pathDir filesep subjID 'processed.txt'];
        % Export processed data
        headerP = [{'PercGcycle'} headerProcessRBAside];
        txtP=sprintf('%s\t',headerP{:});
        txtP(end)='';
        dlmwrite(fnamP,txtP,'');
        dlmwrite(fnamP,rbaDataI,'-append','delimiter','\t','precision',6);
    end
end