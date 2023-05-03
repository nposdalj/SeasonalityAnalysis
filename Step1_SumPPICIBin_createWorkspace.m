clearvars
close all; clear all;clc;

% Step 1 - concatenates all of the TPWS files and creates a timetable for the click and 5-min bin data
% Saves everything as workspace125 to be used for subsequent steps
% Saves 'dataICIgram' and 'binDataICIgram' as '*_mainicipeak.mat' to be used for ICIgrams

%IMPORTANT OUTPUTS:
% clickData - timetable for click data with peak to peak RL levels, ICI, and Peak Frequency for each click (nothing is excluded)
% binData - timetable for 5-min bin data with max peak to peak RL levls and mean peak frequency (bins with less than...
    % 5 clicks are excluded)
% binEffort - timetale of effort for each 5-min bin (this does not include duty cycles since those are very difficult to account...
    % for on the bin level
% dataICIgram - all clicks that have an ICI wihin the thresholds set --> used for ICIgrams
% binDataICIgram - ICIgram data grouped in 5 min bins
%% Parameters defined by user
%Site names
filePrefix = 'HAT'; %TPWS file names that match
siteabrev = 'HAT'; %abbreviation of site (folder names)
region = 'WAT'; %region
sp = 'Pm'; % your species code
itnum = '2'; % which iteration you are looking for (which TPWS folder)
GDrive = 'G'; %Google Drive

%Other parameters
fft = 512; %length of fft
srate = 200; %sampling rate
pfMIN = 5; %peak frequency min
ClickBinMin = 5; %min number of clicks required in a bin

%Data paths
tpwsPath = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\TPWS_125\',siteabrev]; %directory of TPWS files
effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx'];% specify excel file with effort times
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin'); % get default parameters -- make sure these match for your species
%% define subfolder that fit specified iteration
if itnum > 1
   for id = 2: str2num(itnum) % iterate id times according to itnum
       subfolder = ['TPWS',num2str(id)];
       tpwsPath = (fullfile(tpwsPath,subfolder));
   end
end
%% Find all TPWS files that fit your specifications (does not look in subdirectories)
% Concatenate parts of file name
if isempty(sp)
    detfn = [filePrefix,'.*','TPWS',itnum,'.mat'];
else
    detfn = [filePrefix,'.*',sp,'.*TPWS',itnum,'.mat'];
end

% Get a list of all the files in the start directory
fileList = cellstr(ls(tpwsPath));

% Find the file name that matches the filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
if isempty(fileMatchIdx)
    % if no matches, throw error
    error('No files matching filePrefix found!')
end
%% Get effort times matching prefix file
allEfforts = readtable(effortXls); %read effort table
effTable = allEfforts; %effort is for one site only

% make Variable Names consistent
startVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'Start.*Effort'))>0,1,'first');
endVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'End.*Effort'))>0,1,'first');
effTable.Properties.VariableNames{startVar} = 'Start';
effTable.Properties.VariableNames{endVar} = 'End';

Start = datetime(x2mdate(effTable.Start),'ConvertFrom','datenum');
End = datetime(x2mdate(effTable.End),'ConvertFrom','datenum');

effort = timetable(Start,End);
%% Concatenate all detections from the same site
concatFiles = fileList(fileMatchIdx);
%% Parameters to calculate peak frequency
p.N = fft;
p.frRange = [pfMIN srate/2];
smsp2 = p.N/2; % num fft points
ift = 1:smsp2;
fmsp = ((srate/2)/(smsp2-1))*ift - (srate/2)/(smsp2-1);
fi = find(fmsp > p.frRange(1) & fmsp <= p.frRange(2));
%% Concatenate variables
PPall = []; TTall = []; ICIall = []; PeakFrall = [];% initialize matrices
for idsk = 1 : length(concatFiles)
    % Load file
    fprintf('Loading %d/%d file %s\n',idsk,length(concatFiles),fullfile(tpwsPath,concatFiles{idsk}))
    D = load(fullfile(tpwsPath,concatFiles{idsk}));
    
    % find times outside effort (sometimes there are detections
    % which are from the audio test at the beggining of the wav file)
    within = cell2mat(arrayfun(@(x)sum(isbetween(x,datetime(effort.Start),datetime(effort.End))),datetime(D.MTT,'ConvertFrom','datenum'),'uni',false));
    goodIdx = find(within ~= 0);
    MTT = D.MTT(goodIdx); % only keep the good detections
    MPP = D.MPP(goodIdx);
    MSP = D.MSP(goodIdx,:);
    
    % calculate peak frequency
    [~,im] = max(MSP(:,fi),[],2); % maximum between flow-100kHz       
    peakFr = fmsp(im + fi(1) -1)';
    
    % concatenate all data
    TTall = [TTall; MTT];   % group start times
    PPall = [PPall; MPP];   % group peak-to-peak
    PeakFrall = [PeakFrall; peakFr]; % group peak frequency
    ici = diff(MTT)*24*60*60*1000; % in ms
    ICIall = [ICIall;[ici; nan]];  % group inter-click interval
end
%% If you use parfor data may not be sorted. Sort all the variables.
[~,sorted] = sort(TTall);
TTall = TTall(sorted);
PPall = PPall(sorted);
ICIall = ICIall(sorted);
PeakFrall = PeakFrall(sorted);
%% Create timetable for click data
tbin = datetime(TTall,'ConvertFrom','datenum');
clickData = timetable(tbin,PPall,ICIall,PeakFrall);
clear tbin
%% Convert times to bin vector times in 5 min bins
vTT = datevec(TTall);
tbin = datetime([vTT(:,1:4), floor(vTT(:,5)/p.binDur)*p.binDur, ...
    zeros(length(vTT),1)]);
%% create table and get click counts and max pp per bin
data = timetable(tbin,TTall,PPall,ICIall,PeakFrall);

%Group by max peak to peak recieved level first
binDataMAX = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
binDataMAX.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binDataMAX.Properties.VariableNames{'max_PPall'} = 'maxPP';

%Then group by mean ICI and peak frequency
binDataMEAN = varfun(@mean,data,'GroupingVariable','tbin','InputVariable','PeakFrall');
binDataMEAN.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binDataMEAN.Properties.VariableNames{'mean_PeakFrall'} = 'meanPeakFrall';

binData = synchronize(binDataMAX,binDataMEAN(:,2)); %merge two tables
binData(binData.Count < ClickBinMin,:) = []; %identify any bins with less than 5 clicks and deletes them
%% group effort in bins (this doesn't include any effort lost by duty cycle)
effort.diffSec = seconds(effort.End-effort.Start);
effort.bins = effort.diffSec/(60*p.binDur);
effort.roundbin = round(effort.diffSec/(60*p.binDur));

secMonitEffort = sum(effort.diffSec); %seconds recorded
binMonitEffort = sum(effort.roundbin); %bins recorded

[er,~] = size(effort.Start);

if er > 1
    binEffort = intervalToBinTimetable(effort.Start,effort.End,p); % convert intervals in bins when there is multiple lines of effort
    binEffort.sec = binEffort.bin*(p.binDur*60);
else
    binEffort = intervalToBinTimetable_Only1RowEffort(effort.Start,effort.End,p); % convert intervals in bins when there is only one line of effort
    binEffort.effortSec = binEffort.effortBin*(p.binDur*60);
end
%% get average of detection by effort - excluding bins with less than 5 clicks
positiveCounts = sum(binData.Count); %total # of clicks detected
positiveBins = length(binData.Count); %total # of bins detected 
NktTkt = positiveCounts/secMonitEffort; %detection of clicks according to effort
NktTktbin = positiveBins/binMonitEffort; %detection of bins according to effort
%% Prepare data for ICIgrams
% Exclude ICI above thresholds
idxICI = find(ICIall < p.iciRange(1) | ICIall > p.iciRange(2));
ICISel = ICIall;
ICISel(idxICI) = NaN;

dataICIgram = timetable(tbin,TTall,PPall,ICISel,PeakFrall);
dataICIgram = dataICIgram(~any(ismissing(dataICIgram),2),:); %remove all NaNs from ICIs outside of threshold

dataICIgramPP = varfun(@max,dataICIgram,'GroupingVariable','tbin','InputVariable','PPall');
dataICIgramPP.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
dataICIgramPP.Properties.VariableNames{'max_PPall'} = 'maxPP';

dataICIgramPF = varfun(@mean,dataICIgram,'GroupingVariable','tbin','InputVariable',{'PeakFrall','ICISel'});
dataICIgramPF.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
dataICIgramPF.Properties.VariableNames{'mean_PeakFrall'} = 'avgPeakFr';
dataICIgramPF.Properties.VariableNames{'mean_ICISel'} = 'meanICI';

binDataICIgram = synchronize(dataICIgramPP,dataICIgramPF(:,2:3)); %merge two tables
binDataICIgram(binDataICIgram.Count < ClickBinMin,:) = []; %identify any bins with less than ClickBinMin clicks and delete them

icifn = [siteabrev,'_',p.speName,'_mainicipeak'];
save(fullfile(saveDir,icifn),'dataICIgram','binDataICIgram','data','binEffort');
%% save workspace to avoid running previous parts again
save([saveDir,'\',siteabrev,'_workspace125.mat']);