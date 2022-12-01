clearvars
close all;clear all;clc;

% This code will load Rebecca's ID frome her clustering work, filter out
% the PM click times (ID = 19) and compare it with the click times I have
% for the WAT data
%% User Defined Directories
region = 'WAT'; 
filePrefix = 'WAT_GS';
siteabrev = 'GS';
IDDir = 'G:\GS\RebeccasIDFiles'; %location of Rebecca's ID Files
saveDir = ['I:\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
%% Lod data
%load([saveDir,'\',siteabrev,'_workspaceStep2.mat']); %load my work
%Load Rebecca's ID files
detfn = [filePrefix,'.*','ID1','.mat'];
fileList = cellstr(ls(IDDir)); % Get a list of all the files in the start directory
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0); % Find the file name that matches the filePrefix
concatFiles = fileList(fileMatchIdx); % Concatenate all detections from the same site
%% Concatenate variables
TTall = []; IDall = []; 
for idsk = 1 : length(concatFiles)
    % Load file
    fprintf('Loading %d/%d file %s\n',idsk,length(concatFiles),fullfile(IDDir,concatFiles{idsk}))
    D = load(fullfile(IDDir,concatFiles{idsk}));
   
    MTT = D.zID(:,1); % times
    IDD = D.zID(:,2); % ID labels
    
    % concatenate
    TTall = [TTall; MTT];   % group start times
    IDall = [IDall; IDD]; % group IDs

end
%% After parfor data may not be sorted. Sort all the variables.
[~,sorted] = sort(TTall);
TTall = TTall(sorted);
IDall = IDall(sorted);
%% Index TTall and IDall
idx1 = IDall~=19;
TTall(idx1) = [];
IDall(idx1) = [];
%% Create timetable per click
tbin = datetime(TTall,'ConvertFrom','datenum');
clickData = table(tbin,IDall);
%binidx1 = (clickData.IDall ~= 19);
%clickData(binidx1,:) = [];
clickData = table2timetable(clickData);
clear tbin
%% Convert times to bin vector times
vTT = datevec(TTall);
tbin = datetime([vTT(:,1:4), floor(vTT(:,5)/5)*5, ...
    zeros(length(vTT),1)]);
%% create table and get click counts and max pp per bin
data = timetable(tbin,TTall,IDall);
binData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','IDall');
binData.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binData.Properties.VariableNames{'max_IDall'} = 'ID';
%% Presence absence
binTable = synchronize(binData,binEffort);
binidx1 = (binTable.Count >= 1);
[y,~]=size(binTable);
binTable.PreAbs = zeros(y,1);
binTable.PreAbs(binidx1) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 5min bin
