clearvars
close all

%% Parameters defined by user
filePrefix = 'PS'; % File name to match. 
siteabrev = 'PS1'; %abbreviation of site.
region = 'CCE1'; %region
sp = 'Pm'; % your species code
GDrive = 'I'; %Google Drive
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
DutyCy = 0; %if this data is duty cycled, make this equal to 1
%% load workspace
load([saveDir,'\',siteabrev,'_workspace125.mat']);

effortXls(1) = 'I'; %Correct GDrive for SWAL1
GDrive = 'I';
saveDir(1) = 'I';
tpwsPath(1) = 'I';
%% Set up duty cycled dates
% If only one or two deployments are duty cycled, adjust accordingly
if DutyCy == 1
    startTime = datetime(2007,07,02);
    endTime = datetime(2008,06,16);
else
end
%% Remove bin data with less than 5 clicks
binData(binData.Count < 5,:) = []; %identify any bins with less than 5 clicks and delete them
%% group data by 5min bins, hourly, days, weeks, and seasons 
%group data by 5 minute bins
binTable = synchronize(binData,binEffort);
binTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
binTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
binTable.maxPP = [];
binidx1 = (binTable.Count >= 5); %identify any bins with less than 5 clicks and delete them
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
%% accounting for the duty cycle and effort (hourly data)
[p,~]=size(hourlyTab);
hourlyTab.MaxEffort_Bin = ones(p,1)*(12); %total number of bins possible in one day
hourlyTab.MaxEffort_Sec = ones(p,1)*(3600); %seconds in one day

if DutyCy == 1
    DutyCycleIdxStart =  hourlyTab.tbin < startTime;
    startVal = find(DutyCycleIdxStart == 0,1);
    DutyCycleIdxEnd =  hourlyTab.tbin > endTime;
    endVal = find(DutyCycleIdxEnd == 1,1);
else
end

%dealing with duty cycled data
if strcmp(siteabrev,'CSM');
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * .2); %for Cross_01 and 02 only, 5 on 20 off (25 minute cycle)-- meaning you're recording 20% (0.2) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PG');
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * .333); %for Pagan_01 only, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CORC');
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * .5); %for CORC_01 and 02 only, 15 on 15 off (30 minute cycle)-- meaning you're recording 50% (0.5) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
    elseif strcmp(siteabrev,'CORC');
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * (5/35)); %for CORC_01 and 02 only, 15 on 15 off (30 minute cycle)-- meaning you're recording 50% (0.5) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CA');
    hourlyTab.Effort_Bin = hourlyTab.Effort_Bin * (4/12); %for GofCA10 and 11, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'QC');
    hourlyTab.Effort_Bin(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * (0.33); %for QC06 I evaluated the duty cycle using continous deployments and I should adjust by 33%
    hourlyTab.Effort_Sec(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
else
    hourlyTab.MaxEffort_Bin = ones(p,1)*(288);
end

%two ways to account for the difference in effort..
%proportion of hours with clicks
hourlyTab.Minutes = hourlyTab.Count_Bin * 5; %convert bins to minutes
hourlyTab.Hours = (hourlyTab.Count_Bin * 5) ./ 60; %convert the number of bins sperm whales were detected in to hours per day
hourlyTab.HoursProp = hourlyTab.Hours./(hourlyTab.Effort_Sec ./ (60 * 60)); %proportion of hours per day w/clicks

%normalizing bin count with duty cycle
hourlyTab.NormEffort_Bin = hourlyTab.Effort_Bin./hourlyTab.MaxEffort_Bin; %what proportion of the day was there effort
hourlyTab.NormEffort_Sec = hourlyTab.Effort_Sec./hourlyTab.MaxEffort_Sec; %what proportion of the day was there effort
hourlyTab.NormBin = hourlyTab.Count_Bin ./ hourlyTab.NormEffort_Bin; %what would the normalized bin count be given the amount of effort
hourlyTab.NormClick = hourlyTab.Count_Click ./ hourlyTab.NormEffort_Sec; %what would be the normalized click count given the amount of effort
hourlyTab.HoursNorm = (hourlyTab.NormBin * 5) ./ 60; %convert the number of 5-min bins per day to hours

%Winter starts on January (closest to the real thing, which is Dec. 21st)
hourlyTab.Season = zeros(p,1);
hourlyTab.month = month(hourlyTab.tbin);
summeridxD = (hourlyTab.month == 7  | hourlyTab.month == 8 | hourlyTab.month == 9);
fallidxD = (hourlyTab.month == 10  | hourlyTab.month == 11 | hourlyTab.month == 12);
winteridxD = (hourlyTab.month == 1  | hourlyTab.month == 2 | hourlyTab.month == 3);
springidxD = (hourlyTab.month == 4  | hourlyTab.month == 5 | hourlyTab.month == 6);

%adds the season according to the month the data was collected
hourlyTab.Season(summeridxD) = 1;
hourlyTab.Season(fallidxD) = 2;
hourlyTab.Season(winteridxD) = 3;
hourlyTab.Season(springidxD) = 4;

%add year and day to data
hourlyTab.Year = year(hourlyTab.tbin); 
hourlyTab.day = day(hourlyTab.tbin,'dayofyear');

NANidx = ismissing(hourlyTab(:,{'NormBin'}));
hourlyTab{:,{'NormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
hourlyTab{:,{'NormClick'}}(NANidx) = 0; %if there was effort, but no detections change the NormClick column to zero
hourlyTab{:,{'HoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero
hourlyTab{:,{'HoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero

%save day table to csv for R
writetable(timetable2table(hourlyTab),[saveDir,'\',siteabrev,'_binData_forGAMGEE.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
%% accounting for the duty cycle and effort (day data)
[p,~]=size(dayTable);
dayTable.MaxEffort_Bin = ones(p,1)*(288); %total number of bins possible in one day
dayTable.MaxEffort_Sec = ones(p,1)*(86400); %seconds in one day

if DutyCy == 1
    DutyCycleIdxStart =  dayTable.tbin < startTime;
    startVal = find(DutyCycleIdxStart == 0,1);
    DutyCycleIdxEnd =  dayTable.tbin > endTime;
    endVal = find(DutyCycleIdxEnd == 1,1);
else
end

%dealing with duty cycled data
if strcmp(siteabrev,'CSM');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .2); %for Cross_01 and 02 only, 5 on 20 off (25 minute cycle)-- meaning you're recording 20% (0.2) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PG');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .333); %for Pagan_01 only, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'HOKE');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * (5/35)); %for Pagan_01 only, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CORC');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .5); %for CORC_01 and 02 only, 15 on 15 off (30 minute cycle)-- meaning you're recording 50% (0.5) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CA');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .333); %for GofCA10 and 11, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'QC');
    dayTable.Effort_Bin(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * (0.33); %for QC06 I evaluated the duty cycle using continous deployments and I should adjust by 33%
    dayTable.Effort_Sec(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
else    
    dayTable.MaxEffort_Bin = ones(p,1)*(288);
end

%two ways to account for the difference in effort..
%proportion of hours with clicks
dayTable.Minutes = dayTable.Count_Bin * 5; %convert bins to minutes
dayTable.Hours = (dayTable.Count_Bin * 5) ./ 60; %convert the number of bins sperm whales were detected in to hours per day
dayTable.HoursProp = dayTable.Hours./(dayTable.Effort_Sec ./ (60 * 60)); %proportion of hours per day w/clicks

%normalizing bin count with duty cycle
dayTable.NormEffort_Bin = dayTable.Effort_Bin./dayTable.MaxEffort_Bin; %what proportion of the day was there effort
dayTable.NormEffort_Sec = dayTable.Effort_Sec./dayTable.MaxEffort_Sec; %what proportion of the day was there effort
dayTable.NormBin = dayTable.Count_Bin ./ dayTable.NormEffort_Bin; %what would the normalized bin count be given the amount of effort
dayTable.NormClick = dayTable.Count_Click ./ dayTable.NormEffort_Sec; %what would be the normalized click count given the amount of effort
dayTable.HoursNorm = (dayTable.NormBin * 5) ./ 60; %convert the number of 5-min bins per day to hours

%Winter starts on January (closest to the real thing, which is Dec. 21st)
dayTable.Season = zeros(p,1);
dayTable.month = month(dayTable.tbin);
summeridxD = (dayTable.month == 7  | dayTable.month == 8 | dayTable.month == 9);
fallidxD = (dayTable.month == 10  | dayTable.month == 11 | dayTable.month == 12);
winteridxD = (dayTable.month == 1  | dayTable.month == 2 | dayTable.month == 3);
springidxD = (dayTable.month == 4  | dayTable.month == 5 | dayTable.month == 6);

%adds the season according to the month the data was collected
dayTable.Season(summeridxD) = 1;
dayTable.Season(fallidxD) = 2;
dayTable.Season(winteridxD) = 3;
dayTable.Season(springidxD) = 4;

%add year and day to data
dayTable.Year = year(dayTable.tbin); 
dayTable.day = day(dayTable.tbin,'dayofyear');

NANidx = ismissing(dayTable(:,{'NormBin'}));
dayTable{:,{'NormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
dayTable{:,{'NormClick'}}(NANidx) = 0; %if there was effort, but no detections change the NormClick column to zero
dayTable{:,{'HoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero
dayTable{:,{'HoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero

%save day table to csv for R
dayBinTAB = timetable2table(dayTable);
writetable(dayBinTAB,[saveDir,'\',siteabrev,'_dayData_forGLMR125.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
%% day table with days grouped together (summed and averaged) ** USE THIS **
[MD,~] = findgroups(dayTable.day);

if length(MD) < 365
    meantab365 = table(dayTable.day(:), dayTable.HoursProp(:));
    meantab365.Properties.VariableNames = {'day' 'HoursProp'};
    sumtab365 = meantab365;
    meantab365.HoursProp(isnan(meantab365.HoursProp)) = 0;
    dayTable.HoursProp(isnan(dayTable.HoursProp)) = 0; 
else
dayTable.day = categorical(dayTable.day);
%mean
% [mean, sem, std, var, range] = grpstats(dayTable.HoursProp, dayTable.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
% meantable = array2table(mean);
% semtable = array2table(sem);
% stdtable = array2table(std);
% vartable = array2table(var);
% rangetable = array2table(range);
% newcol_mean = (1:length(mean))';
% meanarray365 = [newcol_mean mean sem std var range];
% meantab365 = array2table(meanarray365);
% meantab365.Properties.VariableNames = {'Day' 'HoursProp' 'SEM' 'Std' 'Var' 'Range'};
meantab365 = grpstats(timetable2table(dayTable),'day',{'mean','sem','std','var','range'},'DataVars',{'HoursProp','HoursNorm','NormBin','NormClick'}); %takes the mean of each day of the year
end

writetable(meantab365, [saveDir,'\',siteabrev,'_days365GroupedMean_forGLMR125.csv']); %table with the mean for each day of the year
%% Save workspace variable
save([saveDir,'\',siteabrev,'_workspaceStep2.mat']);