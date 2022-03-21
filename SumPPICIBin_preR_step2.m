clearvars
close all

%% Parameters defined by user
filePrefix = 'HOKE'; % File name to match. 
siteabrev = 'HOKE'; %abbreviation of site.
region = 'CCE'; %region
sp = 'Pm'; % your species code
itnum = '2'; % which iteration you are looking for
srate = 200; % sample rate
GDrive = 'H'; %Google Drive

effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx']; % specify excel file with effort times
% effortXls = ['I:\My Drive\',region,'_TPWS_metadataReduced\',siteabrev,'\SeasonalityAnalysis\Pm_Effort.xlsx']; % specify excel file with effort times

saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
% saveDir = ['I:\My Drive\',region,'_TPWS_metadataReduced\',siteabrev,'\SeasonalityAnalysis']; %specify directory to save files
%% load workspace
load([saveDir,'\',siteabrev,'_workspace125.mat']);
%effortXls = 'I:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\KS\Pm_Effort.xlsx'; % specify excel file with effort times
%saveDir = 'I:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\KS'; %specify directory to save files
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

writetable(timetable2table(hourlyTab),[saveDir,'\',siteabrev,'_binData_forGAMGEE.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R

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
%% accounting for the duty cycle and effort
[p,~]=size(dayTable);
dayTable.MaxEffort_Bin = ones(p,1)*(288); %total number of bins possible in one day
dayTable.MaxEffort_Sec = ones(p,1)*(86400); %seconds in one day

%dealing with duty cycled data
if strcmp(siteabrev,'CB');
    ge = dayTable.Effort_Bin(222:516); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(222:516) = ge * 240; %for CB02 10 on 2 off (12 minute cycle) -- meaning you're recording 0.8333 percent of the time
    dayTable.Effort_Sec(222:516) = dayTable.Effort_Bin(222:516) * 5 * 60; %convert from bins into efforts in seconds per day
    else
if strcmp(siteabrev,'BD');
        ge = dayTable.Effort_Bin(222:516); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(222:516) = ge * 144; %for ALEUT03BD ONLY 5 on 5 off (10 minute cycle) -- meaning you're recording 0.5 percent of the time
    dayTable.Effort_Sec(222:516) = dayTable.Effort_Bin(222:516) * 5 * 60; %convert from bins into efforts in seconds per day
else
    if strcmp(siteabrev,'CSM');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .2); %for Cross_01 and 02 only, 5 on 20 off (25 minute cycle)-- meaning you're recording 20% (0.2) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
else
    if strcmp(siteabrev,'PG');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .333); %for Pagan_01 only, 5 on 10 off (15 minute cycle)-- meaning you're recording 20% (0.2) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
else
    if strcmp(siteabrev,'CORC');
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .5); %for CORC_01 and 02 only, 15 on 15 off (30 minute cycle)-- meaning you're recording 20% (0.2) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
else
    if strcmp(siteabrev,'TIN');
    ge = dayTable.Effort_Bin(1:224); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(1:224) = ge * 72; %for TIN02 5 on 15 off (20 minute cycle) -- meaning you're recording 0.25 percent of the time
    dayTable.Effort_Sec(1:224) = dayTable.Effort_Bin(1:224) * 5 * 60; %convert from bins into efforts in seconds per day
    ge = dayTable.Effort_Bin(224:550); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(224:550) = ge * 240; %for TIN03 5 on 1 off (6 minute cycle) -- meaning you're recording 0.8333 percent of the time
    dayTable.Effort_Sec(224:550) = dayTable.Effort_Bin(224:550) * 5 * 60; %convert from bins into efforts in seconds per day
    ge = dayTable.Effort_Bin(551:1858); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(551:1858) = ge * 206; %for TIN04, 05, 06, 07, 08, and 09 5 on 2 off (7 minute cycle) -- meaning you're recording 0.714 percent of the time
    dayTable.Effort_Sec(551:1858) = dayTable.Effort_Bin(551:1858) * 5 * 60; %convert from bins into efforts in seconds per day
    else
if strcmp(siteabrev,'SAP');
    ge = dayTable.Effort_Bin(1:174); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(1:174) = ge * 36; %for TIN01 5 on 35 off (40 minute cycle) -- meaning you're recording 0.125 percent of the time
    dayTable.Effort_Sec(1:174) = dayTable.Effort_Bin(1:174) * 5 * 60; %convert from bins into efforts in seconds per day
    ge = dayTable.Effort_Bin(175:351); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(175:351) = ge * 72; %for TIN02 5 on 15 off (20 minute cycle) -- meaning you're recording 0.25 percent of the time
    dayTable.Effort_Sec(175:351) = dayTable.Effort_Bin(175:351) * 5 * 60; %convert from bins into efforts in seconds per day
    ge = dayTable.Effort_Bin(352:613); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(352:613) = ge * 240; %for TIN03 5 on 1 off (6 minute cycle) -- meaning you're recording 0.8333 percent of the time
    dayTable.Effort_Sec(352:613) = dayTable.Effort_Bin(352:613) * 5 * 60; %convert from bins into efforts in seconds per day
    ge = dayTable.Effort_Bin(614:2380); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(614:2380) = ge * 206; %for TIN04, 05, 06, 07, 08, and 09 5 on 2 off (7 minute cycle) -- meaning you're recording 0.714 percent of the time
    dayTable.Effort_Sec(614:2380) = dayTable.Effort_Bin(614:2380) * 5 * 60; %convert from bins into efforts in seconds per day
    else 
if strcmp(siteabrev,'QC');
    ge = dayTable.Effort_Bin(3:349); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(3:349) = ge * 41; %for QC06 5 on 30 off (35 minute cycle) -- meaning you're recording 0.143 percent of the time
    dayTable.Effort_Sec(3:349) = dayTable.Effort_Bin(3:349) * 5 * 60; %convert from bins into efforts in seconds per day
else
    if strcmp(siteabrev,'Wake');
    ge = dayTable.Effort_Bin(1:94); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(1:94) = ge * 144; %for Wake01 5 on 5 off (0 minute cycle) -- meaning you're recording 0.5 percent of the time
    dayTable.Effort_Sec(1:94) = dayTable.Effort_Bin(1:94) * 5 * 60; %convert from bins into efforts in seconds per day
    ge = dayTable.Effort_Bin(94:97); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(94:97) = ge * 240; %for Wake03 5 on 1 off (6 minute cycle) -- meaning you're recording 0.8333 percent of the time
    dayTable.Effort_Sec(94:97) = dayTable.Effort_Bin(94:97) * 5 * 60; %convert from bins into efforts in seconds per day
    ge = dayTable.Effort_Bin(97:1638); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(97:1638) = ge * 48; %for TIN04, 05, 06, 07, and 08 5 on 25 off (30 minute cycle) -- meaning you're recording 0.167 percent of the time
    dayTable.Effort_Sec(97:1638) = dayTable.Effort_Bin(97:1638) * 5 * 60; %convert from bins into efforts in seconds per day
    else
dayTable.MaxEffort_Bin = ones(p,1)*(288);
      end
end
end 
end
end
    end
end 
end
end

%two ways to account for the difference in effort..
%proportion of hours with clics
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

%save day table to csv for R
dayBinTAB = timetable2table(dayTable);
writetable(dayBinTAB,[saveDir,'\',siteabrev,'_dayData_forGLMR125.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
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

writetable(meantab365, [saveDir,'\',siteabrev,'_days365GroupedMean_forGLMR125.csv']); %table with the mean for each day of the year
%% Integral Time Scale Calculation 
dayTable.HoursProp(isnan(dayTable.HoursProp)) = 0; 
ts = dayTable.HoursProp;
its_cont = IntegralTimeScaleCalc(ts);

ts = meantab365.HoursProp;
its_nonCont = IntegralTimeScaleCalc(ts);
%% Save workspace variable
save([saveDir,'\',siteabrev,'_workspaceStep2.mat']);