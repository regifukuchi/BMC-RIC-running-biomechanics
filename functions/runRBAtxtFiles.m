% This script imports the txt files that were exported from Visual 3D
% and averages the trials for each subject and then saves a struct matlab
% files containing the average pattern of all the analyzed variables for
% each subject in a rba.mat file
clear all, close all, clc
% Only subjects with good data including 3 gait speeds so far
% subjs = [10 12:17 20 23 27];
subjs = [8:10 12:18 20 23:26 27:30 32:40];

pathDir = 'I:\BMClab\RBA\Data';
pathDir = '/Users/regifukuchi/Dropbox/PUBS/Papers/RBA';


speed = [2.5,3.5,4.5];
side   = {'R','L'};
joints = {'pelvis','hip','knee','ankle'};
varListRBA = {'Stride_Length','Strides_Per_Minute','Stride_Width'};
%%
if 0
    
    %
    for sub = 1:length(subjs)
        
        if subjs(sub) < 10
            dirPath = [pathDir '\SUB000' num2str(subjs(sub))];
        else
            dirPath = [pathDir '\SUB00' num2str(subjs(sub))];
        end
        
        
        for gs = 1:3
            % Import Temporal Distance Gait Parameters
            filename = [dirPath filesep 'RUNT' num2str(speed(gs)*10) 'tempSpatGaitParams.txt'];
            xX = importdata(filename);
            t=xX.textdata{2,1}; %assumes your header is in the last row
            h=regexp(t,'\t','split');
            
            % Variable names
            varnames = h(2:end); % 21 variables
            
            % Data
            vdata = xX.data(2,:);
            
            for s = 1:length(side)
                
                % Temporal Distance Gait Parameters
                for ix = 1:length(varListRBA)
                    if ix==3
                        varFullName = [varListRBA{ix} '_Mean'];
                    else
                        if s==1, lado='Right'; else lado = 'Left';end
                        
                        varFullName = [lado '_' varListRBA{ix} '_Mean'];
                    end
                    % Index of the corresponding variables in the big list
                    idx = find(ismember(varnames,varFullName));
                    
                    xData = vdata(idx);
                    % Structure containing the Temp Spatial Variables
                    rba(sub).([side{s} varListRBA{ix} 'RUN' num2str(speed(gs)*10)]) = xData;
                end
                
                
                for j = 1:length(joints)
                    
                    if j==1
                        varType={'Ang'};
                    else
                        varType = {'Ang','Mom','Pow'};
                    end
                    
                    for v = 1:length(varType)
                        
                        
                        % File directory pathname
                        filename = [dirPath filesep 'RUN' num2str(speed(gs)*10) side{s} joints{j} varType{v} '.txt'];
                        
                        X = importdata(filename);
                        
                        rba(sub).([side{s} joints{j} varType{v} 'RUN' num2str(speed(gs)*10)]) = X.data(:,2:end);
                    end
                    
                end
                
                % Including ground reaction forces
                filename = [dirPath filesep 'RUN' num2str(speed(gs)*10) side{s} 'grf.txt'];
                
                X = importdata(filename);
                
                rba(sub).([side{s} 'grf' 'RUN' num2str(speed(gs)*10)]) = X.data(:,2:end);
                
                
            end
        end
        
        % Subjects indices
        rba(sub).idxSubj = ['Subject ' num2str(subjs(sub))];
        
        disp('************************************************************')
        disp(['SUBJECT ' num2str(subjs(sub)) ' DONE!'])
        
    end
    
    % Saving the non-scalar structure variable with all the data (trials and subjects) in a file
    % save([pathDir '\Data\rba2.mat'],'rba') % Including Temporal Spatial Params now
else
    load([pathDir filesep 'rba2.mat'],'rba')
end
%% Calculating average across trials

if 0
    subjs = find(ismember(subjs,[8,9,12,14,17,18,24,28,29,33])); % Non-rearfoot strike subjects
elseif 0
    subjs = find(~ismember(subjs,[8,9,12,14,17,18,24,28,29,33])); % Rearfoot strikers
    
else
    subjs = 1:length(subjs); % everybody
end


direction = {'X','Y','Z'};
joints = {'hip','knee','ankle'};
time = 0:100;
varType = {'Ang','Mom','Pow'};

for sub = 1:length(subjs)
    for gs = 1:3
        for s = 1:length(side)
            
            for v = 1:length(varType)
                for j = 1:length(joints)
                    for di = 1:length(direction)
                        
                        xx = rba(subjs(sub)).([side{s} joints{j} varType{v} 'RUN' num2str(speed(gs)*10)]);
                        rbaM.([side{s} joints{j} varType{v} 'RUN' num2str(speed(gs)*10)])(:,(3*sub-3)+di) = nanmean(xx(:,di:3:end-(3-di)),2);
                        
                        
                        if j==1
                            % Ground Reaction Forces
                            xx1 = rba(subjs(sub)).([side{s} 'grf' 'RUN' num2str(speed(gs)*10)]);
                            rbaM.([side{s} 'grf'  'RUN' num2str(speed(gs)*10)])(:,(3*sub-3)+di) = nanmean(xx1(:,di:3:end-(3-di)),2);
                            
                            % Pelvis Angles
                            xPelvis = rba(subjs(sub)).([side{s} 'pelvisAng' 'RUN' num2str(speed(gs)*10)]);
                            rbaM.([side{s} 'pelvisAng'  'RUN' num2str(speed(gs)*10)])(:,(3*sub-3)+di) = nanmean(xPelvis(:,di:3:end-(3-di)),2);
                        end
                        
                    end
                    % Getting Temp Distance Params for each subject and storing in
                    % the structure
                    xTS = rba(subjs(sub)).([side{s} varListRBA{j} 'RUN' num2str(speed(gs)*10)]);
                    rbaM.([side{s} varListRBA{j} 'RUN' num2str(speed(gs)*10)])(sub) = xTS;
                end
                
            end
        end
    end
end
rbaM.idxSubjALL = subjs;
%% Saving the scalar structure variable with all the data (trials and subjects) in a file
save([pathDir filesep 'rbaM2.mat'],'-struct','rbaM')