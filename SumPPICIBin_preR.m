clearvars
close all

%% Parameters defined by user
filePrefix = 'GofAK_CB'; % File name to match. 
siteabrev = 'CB'; %abbreviation of site.
sp = 'Pm'; % your species code
itnum = '2'; % which iteration you are looking for
srate = 200; % sample rate
tpwsPath = 'E:\Project_Sites\CB\TPWS_125'; %directory of TPWS files
effortXls = 'E:\Project_Sites\CB\Pm_Effort_CB.xlsx'; % specify excel file with effort times
saveDir = 'E:\Project_Sites\CB\Seasonality'; %specify directory to save files
%% load workspace
load([saveDir,'\',siteabrev,'_workspace125.mat']);
%% group data by 5min bins, days, weeks, and seasons
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
dayTable.Minutes = dayTable.Count_Bin * 5; %convert bins to minutes
dayTable.Hours = (dayTable.Count_Bin * 5) ./ 60; %convert the number of bins sperm whales were detected in to hours per day
dayTable.HoursProp = dayTable.Hours./(dayTable.Effort_Sec ./ (60 * 60)); %proportion of hours per day w/clicks
%% accounting for the duty cycle and effort
[p,~]=size(dayTable);
dayTable.MaxEffort_Bin = ones(p,1)*(288); %total number of bins possible in one day
dayTable.MaxEffort_Sec = ones(p,1)*(86400); %seconds in one day

%dealing with duty cycled data
if strcmp(siteabrev,'CB');
dayTable.Effort_Bin(222:507) = 127;%for CB02 ONLY - only .44 of each hour is recorded...
%so effort of 5 min bins for each day is 127 bins
    else
if strcmp(siteabrev,'BD');
dayTable.Effort_Bin(274:end) = 96; %for ALEUT03BD ONLY - only 0.33 of each hour is recorded...
%so effort of 5 min bins for each day is 96
    else
dayTable.MaxEffort_Bin = ones(p,1)*(288);
end
end

dayTable.NormEffort_Bin = dayTable.Effort_Bin./dayTable.MaxEffort_Bin; %what proportion of the day was there effort
dayTable.NormEffort_Sec = dayTable.Effort_Sec./dayTable.MaxEffort_Sec; %what proportion of the day was there effort
dayTable.NormBin = dayTable.Count_Bin ./ dayTable.NormEffort_Bin; %what would the normalized bin count be given the amount of effort
dayTable.NormClick = dayTable.Count_Click ./ dayTable.NormEffort_Sec; %what would be the normalized click count given the amount of effort
dayTable.HoursNorm = (dayTable.NormBin * 5) ./ 60; %convert the number of 5-min bins per day to hours

dayTable.Season = zeros(p,1);
dayTable.month = month(dayTable.tbin);
summeridxD = (dayTable.month == 6  | dayTable.month == 7 | dayTable.month == 8);
fallidxD = (dayTable.month == 9  | dayTable.month == 10 | dayTable.month == 11);
winteridxD = (dayTable.month == 12  | dayTable.month == 1 | dayTable.month == 2);
springidxD = (dayTable.month == 3  | dayTable.month == 4 | dayTable.month == 5);

%adds the season according to the month the data was collected
dayTable.Season(summeridxD) = 1;
dayTable.Season(fallidxD) = 2;
dayTable.Season(winteridxD) = 3;
dayTable.Season(springidxD) = 4;

%add year to data
dayTable.Year = year(dayTable.tbin); 

NANidx = ismissing(dayTable(:,{'NormBin'}));
dayTable{:,{'NormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero

%save day table to csv for R
dayBinTAB = timetable2table(dayTable);
writetable(dayBinTAB,[saveDir,'\',siteabrev,'_dayData_forGLMR125.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
%% day table with days grouped together (summed and averaged) ** USE THIS **
dayTable.day = day(dayTable.tbin,'dayofyear');
[MD,~] = findgroups(dayTable.day);
dayTable.day = categorical(dayTable.day);

if length(MD) < 365
    meantab365 = table(dayTable.day(:), dayTable.HoursProp(:));
    meantab365.Properties.VariableNames = {'Day' 'HoursProp'};
    sumtab365 = meantab365;
else
%mean
meanarray = grpstats(dayTable.HoursProp, dayTable.day, @mean); %takes the mean of each day of the year
meantable = array2table(meanarray);
newcol_mean = (1:length(meanarray))';
meanarray365 = [newcol_mean meanarray];
meantab365 = array2table(meanarray365);
meantab365.Properties.VariableNames = {'Day' 'HoursProp'};
end

writetable(meantab365, [saveDir,'\',siteabrev,'_days365GroupedMean_forGLMR125.csv']); %table with the mean for each day of the year
%% Integral Time Scale Calculation 
ts = dayTable.HoursProp;
its_cont = IntegralTimeScaleCalc(ts);

ts = meantab365.HoursProp;
its_nonCont = IntegralTimeScaleCalc(ts);