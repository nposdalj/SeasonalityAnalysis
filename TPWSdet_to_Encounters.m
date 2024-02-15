clearvars
close all
%% Parameters defined by user
tpwsPath = 'L:\Shared drives\MBARC_All\PacFleet_2023_2024\MFAS Per Site\SOCAL_SN\62'; %directory of TPWS files
saveDir = 'L:\Shared drives\MBARC_All\PacFleet_2023_2024\MFAS Per Site\SOCAL_SN\62'; %directory where to save outputs
site = 'SOCAL';
station = 'SN';
deployment = '62';
df = 'df20';
itnum = '2'; %which iteration of TPWS are you looking for
sp = 'MFA'; %species code
encL = 30; %length of time between encounters (min)
%% define subfolder that fit specified iteration
% if itnum > 1
%    for id = 2: str2num(itnum) % iterate id times according to itnum
%        subfolder = ['TPWS',num2str(id)];
%        tpwsPath = (fullfile(tpwsPath,subfolder));
%    end
% end
%% Find all TPWS files that fit your specifications (does not look in subdirectories)
filePrefix = [site,'_',station,'_',deployment]; %concatenate file prefix
% Concatenate parts of file name
if isempty(sp)
    detfn = [filePrefix,'.*','TPWS',itnum,'.mat'];
else
    if isempty(df)
        detfn = [filePrefix,'.*',sp,'.*TPWS',itnum,'.mat'];
    else
        detfn = [filePrefix,'.*',df,'_',sp,'_TPWS',itnum,'.mat'];
    end
end

% Get a list of all the files in the start directory
fileList = cellstr(ls(tpwsPath));

% Find the file name that matches the filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
if isempty(fileMatchIdx)
    % if no matches, throw error
    error('No files matching filePrefix found!')
end
%% Concatenate all detections from the same site
concatFiles = fileList(fileMatchIdx);
%% Concatenate variables
TTall = []; % initialize matrices
LoadingMessage = 'Loading %d/%d file %s\n';
parfor idsk = 1 : length(concatFiles)
    % Load file
    fprintf(LoadingMessage,idsk,length(concatFiles),fullfile(tpwsPath,concatFiles{idsk}))
    D = load(fullfile(tpwsPath,concatFiles{idsk}));
    
    % concatenate
    TTall = [TTall; D.MTT];   % group start times
end
%% After parfor data may not be sorted. Sort all the variables.
[~,sorted] = sort(TTall);
TTall = TTall(sorted);
%% Create timetable per click
tbin = datetime(TTall,'ConvertFrom','datenum');
clickData = timetable(tbin);
clear tbin
%% Difference between click times (the start times eventually become the end 
%times of the encounter and the end times become the start times of the encounters)
clickDataStart = clickData; %seperate table into start times
clickDataEnd = clickData; %seperate table into end times

clickDataStart(end,:) = []; %delete the last click
clickDataEnd(1,:) = []; %delete the first click

timeDiff = clickDataEnd.tbin - clickDataStart.tbin; %take time diff between each click

clickDataStart.Break = timeDiff > duration(00,encL,00); %if time diff is greater than encL, logical = 1
clickDataEnd.Break = timeDiff > duration(00,encL,00); %if time diff is greater than encL, logical = 1
clickDataStart(~clickDataStart.Break,:)=[]; %delete time diffs less than encL
clickDataEnd(~clickDataEnd.Break,:)=[]; %delete time diffs less than encL

EncStart = timetable2table(clickDataEnd); %convert to table
EncStart.Break = []; %delete logical
EncEnd = timetable2table(clickDataStart); %convert to table
EncEnd.Break = []; %delete logical

clickTable = timetable2table(clickData); %conver to table
firstclick = clickTable(1,:); %find time of first click
lastclick = clickTable(end,:); %find time of last click

EncStartTimes = [firstclick;EncStart]; %add first click as start 
EncStartTimes.Properties.VariableNames = {'StartTime'};
EncEndTimes = [EncEnd;lastclick]; %add last click as end 
EncEndTimes.Properties.VariableNames = {'EndTime'};

EncounterTimes = [EncStartTimes EncEndTimes]; %combine start and end times into one table
filenameX = [saveDir,'\',filePrefix,'_',sp,'_EncounterTimes.xlsx'];
filenameM = [saveDir,'\',filePrefix,'_',sp,'_EncounterTimes.mat'];
writetable(EncounterTimes,filenameX)
save(filenameM,'EncounterTimes')