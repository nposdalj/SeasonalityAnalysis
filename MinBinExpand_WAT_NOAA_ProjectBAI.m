close all;clear all;clc;
%% Load data
siteabrev = 'OC';
saveDir = ['E:\',siteabrev,'\SeasonalityAnalysis\']; %specify directory to save files
%saveDir = 'I:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis\GS';
load([saveDir,'\',siteabrev,'_workspace125.mat']);
%% Find start and end times
startTime = binTable.tbin(1);
endTime = binTable.tbin(end);
missing = (startTime:minutes(1):endTime)';
%Make a new table
[p,q]=size(missing);
x = zeros(p, 1);
allBins = timetable(missing,x);
binSynch = synchronize(binTable,allBins);
binSynch.Count(isnan(binSynch.Count))=0;
binSynch.Effort_Bin(isnan(binSynch.Effort_Bin))=0;
binSynch.Effort_Sec(isnan(binSynch.Effort_Sec))=0;
binSynch(:,'x') = [];

%% Save table
writetable(timetable2table(binSynch), [saveDir,'\',siteabrev,'_1minBinnedData.csv']); %table with the mean for each day of the year



