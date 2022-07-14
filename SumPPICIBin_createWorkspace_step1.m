clearvars
close all

%% Parameters defined by user
filePrefix = 'Palmyra'; % File name to match. 
siteabrev = 'Palmyra'; %abbreviation of site.
region = 'CentralPac'; %region
sp = 'Pm'; % your species code
itnum = '2'; % which iteration you are looking for (which TPWS folder)
GDrive = 'I'; %Google Drive

tpwsPath = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\TPWS_125\',siteabrev]; %directory of TPWS files
effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx'];% specify excel file with effort times
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
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
%% get default parameters
p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin');
%% Parameters to calculate peak frequency
% fix parameters that should be given from sp_setting_default (now it is
% just hard coded)
p.N = 512;
srate = 200;
p.frRange = [5 srate/2];

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
    
    % concatenate
    TTall = [TTall; MTT];   % group start times
    PPall = [PPall; MPP];   % group peak-to-peak
    PeakFrall = [PeakFrall; peakFr]; % group peak frequency

    
    % Inter-Click Interval
    ici = diff(MTT)*24*60*60*1000; % in ms
    ICIall = [ICIall;[ici; nan]];  % group inter-click interval
end
%% After parfor data may not be sorted. Sort all the variables.
[~,sorted] = sort(TTall);
TTall = TTall(sorted);
PPall = PPall(sorted);
ICIall = ICIall(sorted);
%% Create timetable per click
tbin = datetime(TTall,'ConvertFrom','datenum');
clickData = timetable(tbin,PPall,ICIall);
clear tbin
%% Convert times to bin vector times in 5 min bins
vTT = datevec(TTall);
tbin = datetime([vTT(:,1:4), floor(vTT(:,5)/p.binDur)*p.binDur, ...
    zeros(length(vTT),1)]);
%% create table and get click counts and max pp per bin
data = timetable(tbin,TTall,PPall);
binData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
binData.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binData.Properties.VariableNames{'max_PPall'} = 'maxPP';

positiveCounts = sum(binData.Count);
positiveBins = length(binData.Count);
%% group effort in bins
effort.diffSec = seconds(effort.End-effort.Start);
effort.bins = effort.diffSec/(60*p.binDur);
effort.roundbin = round(effort.diffSec/(60*p.binDur));

secMonitEffort = sum(effort.diffSec);
binMonitEffort = sum(effort.roundbin);

[er,~] = size(effort.Start);

if er > 1
    binEffort = intervalToBinTimetable(effort.Start,effort.End,p); % convert intervals in bins when there is multiple lines of effort
    binEffort.sec = binEffort.bin*(p.binDur*60);
else
    binEffort = intervalToBinTimetable_Only1RowEffort(effort.Start,effort.End,p); % convert intervals in bins when there is only one line of effort
    binEffort.sec = binEffort.bin*(p.binDur*60);
end
%% get average of detection by effort
NktTkt = positiveCounts/secMonitEffort;
NktTktbin = positiveBins/binMonitEffort;
%% ICIgram information
%  Exclude ICI above thresholds
idxICI = find(ICIall < p.iciRange(1) | ICIall > p.iciRange(2));
ICISel = ICIall;
ICISel(idxICI) = NaN;

dataICIgram = timetable(tbin,TTall,PPall,ICISel,PeakFrall);
binDataICIgram = varfun(@max,dataICIgram,'GroupingVariable','tbin','InputVariable','PPall');
binDataICIgram.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binDataICIgram.Properties.VariableNames{'max_PPall'} = 'maxPP';

binData2ICIgram = varfun(@mean,dataICIgram,'GroupingVariable','tbin','InputVariable','PeakFrall');
binData2ICIgram.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binData2ICIgram.Properties.VariableNames{'mean_PeakFrall'} = 'avgPeakFr';

binDataICIgram = synchronize(binDataICIgram,binData2ICIgram(:,2));

icifn = [filePrefix,'_',p.speName,'_mainicipeak'];
save(fullfile(saveDir,icifn),'dataICIgram','binDataICIgram');
%% save workspace to avoid running previous parts again
save([saveDir,'\',siteabrev,'_workspace125.mat']);