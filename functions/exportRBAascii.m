% This script was written to read the ASCII files exported from Cortex
%% Markers data
clc, close all , clear all

pathDir = 'H:\RBA\Data';

subjs = [8:10 12:18 20 23:30 32:40]; % Subjects included in the Data Base
speed = [2.5,3.5,4.5]; % Gait speed
direction = {'x','y','z'};
% Indices of markers that remain in the dynamic trials
iMarkerDyn = [1:28 36 37 46 47];

headerForce = {'Fx','Fy','Fz','COPx','COPy','COPz','Ty'};
%%
for is = 2:length(subjs)
    subjDir = ['SUB00' num2str(subjs(is),'%02i')]; % Subject label
    
    % Detecting filename
    xFile = dir([pathDir filesep subjDir filesep 'RBDS' filesep 'static*.trc']);
    filenameOs = xFile.name;
    
    % Importing static data trc file
    xXkineS = importdata([pathDir filesep subjDir filesep 'RBDS' filesep filenameOs],'\t',5);
    
    % Appending projects and subjects number in the file names
    filenameDs = ['RBDS0' num2str(is,'%02i') filenameOs]; % Subject label
    
    filenameDs(end-4:end)=[]; % Deleting file extension name
    
    % keeping only one second of static trial data
    if size(xXkineS.data,1) > 150
        xXkineS.data(151:end,:) = [];
    end
    
    % Static trial
    markersS = xXkineS.data(:,3:end); % ignore first two columns
    timeS    = xXkineS.data(:,2);
    
    % Marker labels
    xXLabelS = strsplit(xXkineS.textdata{4});
    markerLabels = xXLabelS(3:50);
    
    % Adding marker labels with 3 components x,y,z
    for imarker = 1:length(markerLabels)
        for idir = 1:3
            markerLabelsXYZ{(3*imarker-3)+idir} = [markerLabels{imarker} upper(direction{idir})];
        end
    end
    
    if size(markersS,2)/3 > 48
        % Delete additional markers (eg. virtual markers) exported from Cortex.
        markersS(:,48*3+1:end) = [];
    end
    
    
    markersS = [timeS markersS];
    
    
    fnam = [pathDir filesep subjDir filesep 'RBDS' filesep filenameDs '.txt'];
    
    
    headerS = [{'Time'} markerLabelsXYZ];
    
    txt=sprintf('%s\t',headerS{:});
    txt(end)='';
    dlmwrite(fnam,txt,'');
    dlmwrite(fnam,markersS,'-append','delimiter','\t','precision',6);
    
    
    for gs = 1:3
        % Detecting file name
        xFileD = dir([pathDir filesep subjDir filesep 'RBDS' filesep 'runT' num2str(speed(gs)*10) '*.trc']);
        filenameDd = xFileD.name;
        
        filenameDd(end-3:end)=[]; % Deleting file extension name
        
        xXkineD = importdata([pathDir filesep subjDir filesep 'RBDS' filesep filenameDd '.trc'],'\t',5); % runing trial
        xXforce = importdata([pathDir filesep subjDir filesep 'RBDS' filesep filenameDd '.forces']); % forces
        
        % Appending projects and subjects number in the file names
        filenameD = ['RBDS0' num2str(is,'%02i') filenameDd]; % Subject label
        
        filenameD(end)=[]; % Deleting file extension name
        
        % Deleting data if the duration is more than 30 s
        if size(xXkineD.data,1) > 4500;
            xXkineD.data(4501:end,:) = [];
            xXforce.data(9001:end,:) = [];
        end
        
        markersD= xXkineD.data(:,3:end); % ignore first two columns
        
        % Keeping only markers that will remain in the running trials
        for iMarker = 1:length(iMarkerDyn)
            markersDyn(:,3*iMarker-2:3*iMarker) = markersD(:,3*iMarkerDyn(iMarker)-2:3*iMarkerDyn(iMarker));
            markerLabelsD{iMarker} = markerLabels{iMarkerDyn(iMarker)};
        end
        
        timeM = xXkineD.data(:,2);
        
        markersS = [timeM markersDyn];
        
        forceData = xXforce.data(:,1:8);
        
        fnamM = [pathDir filesep subjDir filesep 'RBDS' filesep filenameD 'markers.txt'];
        fnamF = [pathDir filesep subjDir filesep 'RBDS' filesep filenameD 'forces.txt'];
        
        
        % Adding marker labels with 3 components x,y,z. Have to do for the
        % dynamic because the number of markers are smaller
        for imarker = 1:length(markerLabelsD)
            for idir = 1:3
                markerLabelsXYZd{(3*imarker-3)+idir} = [markerLabelsD{imarker} upper(direction{idir})];
            end
        end
        
        % Markers
        headerM = [{'Time'} markerLabelsXYZd];
        txtM=sprintf('%s\t',headerM{:});
        txtM(end)='';
        dlmwrite(fnamM,txtM,'');
        dlmwrite(fnamM,markersS,'-append','delimiter','\t','precision',6);
        
        % Forces
        headerF = [{'Time'} headerForce];
        txtF=sprintf('%s\t',headerF{:});
        txtF(end)='';
        dlmwrite(fnamF,txtF,'');
        dlmwrite(fnamF,forceData,'-append','delimiter','\t','precision',6);
        
    end
    disp('**********************************************************************')
    disp(['Files subject ' num2str(subjs(is)) ' created!!'])
end
%% Dynamic trials
direction = {'x','y','z'};

speed = [2.5,3.5,4.5];

filePathO = 'H:\RBA\Data\SUB0008\';
filePathD = 'H:\RBA\Data\SUB0008\Temp\';

headerForce = {'Fx','Fy','Fz','Mx','My','Mz','Ty'};

% Marker labels
markerLabels = {'RASIS','LASIS','RPSIS','LPSIS','RILC','LILC','RTTL','RTBL','RTTM','RTBM',...
    'RSTL','RSBL','RSTM','RSBM','RHET','RHEB','RHEL','LTTL','LTBL','LTTM','LTBM',...
    'LSTL','LSBL','LSTM','LSBM','LHET','LHEB','LHEL','RGTR','RKNL','RKNM','RHFB','RTTU','RAKL','RAKM','RMT1','RMT5','RMT2',...
    'LGTR','LKNL','LKNM','LHFB','LTTU','LAKL','LAKM','LMT1','LMT5','LMT2'};



for is = 1:3
    
    
    xXkineD = importdata([filePathO 'runT' num2str(speed(is)*10) '.trc'],'\t',5); % runing trial
    xXforce = importdata([filePathO 'runT' num2str(speed(is)*10) '.forces']); % forces
    
    % Deleting data if the duration is more than 30 s
    if size(xXkineD.data,1)
        xXkineD.data(1:4500,:) = [];
        xXforce.data(1:9000,:) = [];
    end
    
    markersD= xXkineD.data(:,3:end); % ignore first two columns
    
    
    for iMarker = 1:length(iMarkerDyn)
        markersDyn(:,3*iMarker-2:3*iMarker) = markersD(:,3*iMarkerDyn(iMarker)-2:3*iMarkerDyn(iMarker));
        markerLabelsD{iMarker} = markerLabels{iMarkerDyn(iMarker)};
    end
    
    timeM = xXkineD.data(:,2);
    
    markersS = [timeM markersDyn];
    forceData = xXforce.data(:,1:8);
    
    fnamM = [filePathD 'runT' num2str(speed(is)*10) 'markers.txt'];
    fnamF = [filePathD 'runT' num2str(speed(is)*10) 'loads.txt'];
    
    for iS = 1:length(markerLabelsD)
        for iDir = 1:3
            markerLabelsD2{iS} = [markerLabelsD{iS} direction{iDir}];
        end
    end
    
    % Markers
    headerM = [{'Time'} markerLabelsD2];
    txtM=sprintf('%s\t',headerM{:});
    txtM(end)='';
    dlmwrite(fnamM,txtM,'');
    dlmwrite(fnamM,markersS,'-append','delimiter','\t','precision',6);
    
    % Forces
    headerF = [{'Time'} headerForce];
    txtF=sprintf('%s\t',headerF{:});
    txtF(end)='';
    dlmwrite(fnamF,txtF,'');
    dlmwrite(fnamF,forceData,'-append','delimiter','\t','precision',6);
    
end

%% Read markers and forces file
% Same number of frames
fp1     = xXforce.data(:,2:8); % only first force plate

markersD= xXkineD.data(:,3:end); % ignore first two columns

% Indices of markers that remain in the dynamic trials
iMarkerDyn = [1:28 36 37 46 47];

for iMarker = 1:length(iMarkerDyn)
    markersDyn(:,3*iMarker-2:3*iMarker) = markersD(:,3*iMarkerDyn(iMarker)-2:3*iMarkerDyn(iMarker));
    markerLabelsD{iMarker} = markerLabels{iMarkerDyn(iMarker)};
end

% Removing extra markers in case there are more than 48
if size(markersS,2) > length(markerLabels)*3
    markersS(:,length(markerLabels)*3+1:end)=[];
end

framesM = xXkineD.data(:,1);% frames vector markers
timeM   = xXkineD.data(:,2);% time vector markers


framesF = xXforce.data(:,1);% frames vector markers
timeF   = xXforce.data(:,2);% time vector markers
%% Saving and exporting files as ASCII to be uploaded in Figshare
% Force data
headerForce = 'Time Fx Fy Fz Mx My Mz Ty';
headerForce = {'Time','Fx','Fy','Fz','Mx','My','Mz','Ty'};

forceData = xXforce.data(:,1:8);
fnam = 'H:\RBA\Data\SUB0008\Temp\Forces.txt';

txt=sprintf('%s\t',headerForce{:});
txt(end)='';
dlmwrite(fnam,txt,'');
dlmwrite(fnam,forceData,'-append','delimiter','\t','precision',6);
clc, close all , clear all

xX = importdata('H:\RBA\Data\SUB0008\Temp\static.txt');

staticLabel = xX.textdata(2,:); % Check how I got the list markers in a code I did long time ago.

time = xX.data(:,1); % time vector

markersS = xX.data(:,2:end);


markerLabels = {'RASIS','LASIS','RPSIS','LPSIS','RILC','LILC','RTTL','RTBL','RTTM','RTBM',...
    'RSTL','RSBL','RSTM','RSBM','RHET','RHEB','RHEL','LTTL','LTBL','LTTM','LTBM',...
    'LSTL','LSBL','LSTM','LSBM','LHET','LHEB','LHEL','RGTR','RKNL','RKNM','RHFB','RTTU','RAKL','RAKM','RMT1','RMT5','RMT2',...
    'LGTR','LKNL','LKNM','LHFB','LTTU','LAKL','LAKM','LMT1','LMT5','LMT2'};