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
SITE = {'BD','KS'};
region = 'GofAK'; %region
sp = 'Pm'; % your species code
GDrive = 'I'; %Google Drive

for i = 1:length(SITE)
siteabrev = SITE{i};
filePrevix = SITE{i};
saveDirr = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
%% load workspace
GDrive_correct = GDrive; % Preserve correct GDrive as it was entered above
%% Sex Data
Genderfn = [saveDirr,'\',region,'_',siteabrev,'_Pm_gender.mat'];
if isfile(Genderfn)
    load(Genderfn);
else
    load([saveDirr,'\',siteabrev,'_Pm_gender.mat']);
end
gender_binData = binData;
load([saveDirr,'\',siteabrev,'_workspace125.mat']);

p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin'); % get default parameters -- make sure these match for your species
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
binDataICIgram(binDataICIgram.Count < 5,:) = []; %identify any bins with less than ClickBinMin clicks and delete them

%merge original bin data and cleaned up dataICIgram table
binData = synchronize(binDataICIgram(:,1:4),gender_binData(:,4:8));
binData(isnan(binData.Count),:) = []; %remove NaNs
%% Change column names
binData.OtherB = binData.Other;
binData.Properties.VariableNames("Female") = "SocialGroup";
binData.Properties.VariableNames("Juvenile") = "MidSize";
binData.Properties.VariableNames("Other") = "OtherA";
binData.BigMale = [];

icifn = [siteabrev,'_',p.speName,'_gender'];
save(fullfile(saveDir,icifn),'binData');
end