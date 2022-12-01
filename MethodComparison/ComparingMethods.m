clearvars
close all

%% Parameters defined by user
filePrefix = 'GS'; % File name to match. 
siteabrev = 'GS'; %abbreviation of site.
region = 'WAT'; %region
sp = 'Pm'; % your species code
itnum = '3'; % which iteration you are looking for
srate = 200; % sample rate
GDrive = 'I'; %Google Drive

effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx']; % specify excel file with effort times

saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
%% load workspace
load([saveDir,'\',siteabrev,'_workspace125.mat']);
%% group data by 5min bins, hourly, days, weeks, and seasons 
%group data by 5 minute bins
binTable = synchronize(binData,binEffort);
binTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
binTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
binTable.maxPP = [];
binidx1 = (binTable.Count >= 1);
[y,~]=size(binTable);
binTable.PreAbs = zeros(y,1);
binTable.PreAbs(binidx1) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 5min bin
%no effort bins are excluded 

%hourly
Click = retime(binData(:,1),'hourly','sum'); % #click per day
Bin = retime(binData(:,1),'hourly','count'); % #bin per day

hourData = synchronize(Click,Bin);
hourlyEffort = retime(binEffort,'hourly','sum');
hourlyTab = synchronize(hourData,hourlyEffort);
hourlyTab.Properties.VariableNames{'bin'} = 'Effort_Bin';
hourlyTab.Properties.VariableNames{'sec'} = 'Effort_Sec';
hourlyTab(~hourlyTab.Effort_Bin,:)=[]; %removes days with no effort, NOT days with no presence
binidx_hourly = (hourlyTab.Count_Bin >= 1);
[y,~]=size(hourlyTab);
hourlyTab.PreAbs = zeros(y,1);
hourlyTab.PreAbs(binidx_hourly) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

%daily
Click = retime(binData(:,1),'daily','sum'); % #click per day
Bin = retime(binData(:,1),'daily','count'); % #bin per day

dayData = synchronize(Click,Bin);
dayEffort = retime(binEffort,'daily','sum');
dayTab = synchronize(dayData,dayEffort);
dayTable = synchronize(dayData,dayEffort);
dayTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
dayTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
dayTableZeros = dayTable;
dayTable(~dayTable.Effort_Bin,:)=[]; %removes days with no effort, NOT days with no presence
%% day table with days grouped together (summed and averaged) ** USE THIS **
[MD,~] = findgroups(dayTable.day);

if length(MD) < 365
    meantab365 = table(dayTable.day(:), dayTable.HoursProp(:));
    meantab365.Properties.VariableNames = {'Day' 'HoursProp'};
    sumtab365 = meantab365;
    meantab365.HoursProp(isnan(meantab365.HoursProp)) = 0;
    dayTable.HoursProp(isnan(dayTable.HoursProp)) = 0; 
else
dayTable.day = categorical(dayTable.day);
%mean
[mean, sem, std, var, range] = grpstats(dayTable.HoursProp, dayTable.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
meantable = array2table(mean);
semtable = array2table(sem);
stdtable = array2table(std);
vartable = array2table(var);
rangetable = array2table(range);
newcol_mean = (1:length(mean))';
meanarray365 = [newcol_mean mean sem std var range];
meantab365 = array2table(meanarray365);
meantab365.Properties.VariableNames = {'Day' 'HoursProp' 'SEM' 'Std' 'Var' 'Range'};
end

[pp,~]=size(meantab365);
meantab365.Season = zeros(pp,1);
meantab365.month = month(meantab365.Day);

%Winter starts on January 1st
summeridxD = (meantab365.month == 7  | meantab365.month == 8 | meantab365.month == 9);
fallidxD = (meantab365.month == 10  | meantab365.month == 11 | meantab365.month == 12);
winteridxD = (meantab365.month == 1  | meantab365.month == 2 | meantab365.month == 3);
springidxD = (meantab365.month == 4  | meantab365.month == 5 | meantab365.month == 6);

%adds the season according to the month the data was collected
meantab365.Season(summeridxD) = 1;
meantab365.Season(fallidxD) = 2;
meantab365.Season(winteridxD) = 3;
meantab365.Season(springidxD) = 4;