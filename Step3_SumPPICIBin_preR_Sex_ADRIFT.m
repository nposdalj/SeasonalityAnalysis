clearvars
close all; clear all;clc;

% Step 3 - Similar to Step 2 but for sex specific groupings
% Groups data into 5-min bins and daily data accounting for effort(including duty cycle)...
% deals with duty cycle by either taking the proportion of presence based on effort or normalizes...
    % presence based on effort (usually these are pretty similar)
    
% Saves everything as workspaceStep3 to be used for subsequent steps

%IMPORTANT OUTPUTS:
% sexbinPresence - data binned in daily bins
% sexhourlyTab - data binned in 1-hr bins
% meantabFE365 - data averaged by Julian days for social groups
% meantabJU365 - data averaged by Julian days for mid-size
% meantabMA365 - data averaged by Julian days males
%.csv Outputs for future statistical analysis in R
    % bin Data - '*_binPresence.csv' %daily data
    % '_365GroupedMeanSocialGroup.csv']); %table with the mean for each day of the year
    % '*_365GroupedMeanMidSize.csv']); %table with the mean for each day of the year
    % '*_365GroupedMeanMale.csv']); %table with the mean for each day of the year
    % '*_binData_forGAMGEE_sexClasses.csv' %hourly bin data for modeling
%% Parameters defined by user
%Site names and data paths
filePrefix = 'ADRIFT_103'; % File name to match. 
genderFileName = 'ADRIFT_103'; %File name to match gender file
siteabrev = 'ADRIFT_103'; %abbreviation of site
sp = 'Pm'; % your species code

effortXls = 'J:\SeasonalityAnalysis\Pm_Effort.xlsx';% specify excel file with effort times
saveDir = ['J:\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
dayBinCSV= [saveDir,'\',siteabrev,'_dayData_forGLMR125.csv']; % specify csv document with general PM information
%% Load sex specific data from ICIgrams
filename = [saveDir,'\',genderFileName,'_',sp,'_gender_filled.mat'];
load(filename);
sexData = binData_mod; %not to get in the way of other bin data
%% load workspace 2
load([saveDir,'\',siteabrev,'_workspaceStep2.mat']);
p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin'); % get default parameters -- make sure these match for your species
%% Adjust bin data not to include less than 5 clicks or ICIs over 2000 ms
%This should have been done in previous steps but some of the old data
%needs adjusting still so I've kept this here
if any(isnan(sexData.SocialGroup)) %Identifies if this is 'new' data or 'old'
dataICIgram_noNaN = dataICIgram;
dataICIgram_noNaN(isnan(dataICIgram_noNaN.ICISel),:) = []; %remove NaNs
idx_ICItoohigh = find(dataICIgram_noNaN.ICISel > MaxICI);
dataICIgram_noNaN(idx_ICItoohigh, :) = [];

dataICIgramPP = varfun(@max,dataICIgram_noNaN,'GroupingVariable','tbin','InputVariable','PPall');
dataICIgramPP.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
dataICIgramPP.Properties.VariableNames{'max_PPall'} = 'maxPP';

dataICIgramPF = varfun(@mean,dataICIgram_noNaN,'GroupingVariable','tbin','InputVariable',{'PeakFrall','ICISel'});
dataICIgramPF.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
dataICIgramPF.Properties.VariableNames{'mean_PeakFrall'} = 'avgPeakFr';
dataICIgramPF.Properties.VariableNames{'mean_ICISel'} = 'meanICI';

dataICIgramFin = synchronize(dataICIgramPP,dataICIgramPF(:,2:3)); %merge two tables
dataICIgramFin(dataICIgramFin.Count < 5,:) = []; %identify any bins with less than 5 clicks and delete them

%merge original bin data and cleaned up dataICIgram table
sexData = synchronize(dataICIgramFin,sexData(:,4:8));
sexData(isnan(sexData.Count),:) = []; %remove NaNs
else
end
%% Extract effort from workspace 2 dayTable that already accounts for duty cycle, etc.
dayBinTAB = timetable2table(dayTable);
sexBinEffort = dayBinTAB(:,{'tbin','Effort_Bin','Effort_Sec','MaxEffort_Bin','MaxEffort_Sec'});
sexBinEffort = table2timetable(sexBinEffort);
%% group data by days and add effort
%Retime sex data daily
sexDEdata = sexData; %save just in case
sexData = retime(sexData,'daily','sum');
sexbinPresence = synchronize (sexData, sexBinEffort);

% Exclude columns that we don't care about
sexbinPresence.maxPP = [];
sexbinPresence.avgPeakFr = [];
sexbinPresence.mean_ICISel = [];
sexbinPresence.OtherA = [];
sexbinPresence.OtherB = [];
sexbinPresence.Count = [];
%% Two ways to account for effort
%PROPORTION OF HOURS WITH CLICKS
sexbinPresence.FeMinutes = sexbinPresence.SocialGroup *5;
sexbinPresence.JuMinutes = sexbinPresence.MidSize *5;
sexbinPresence.MaMinutes = sexbinPresence.Male *5;
sexbinPresence.FeHours = sexbinPresence.FeMinutes ./60;
sexbinPresence.JuHours = sexbinPresence.JuMinutes ./60;
sexbinPresence.MaHours = sexbinPresence.MaMinutes ./60;
sexbinPresence.FeHoursProp = sexbinPresence.FeHours ./(sexbinPresence.Effort_Sec ./ (60*60));
sexbinPresence.JuHoursProp = sexbinPresence.JuHours ./(sexbinPresence.Effort_Sec ./ (60*60));
sexbinPresence.MaHoursProp = sexbinPresence.MaHours ./(sexbinPresence.Effort_Sec ./ (60*60));

%NORMALIZING BIN COUNT ACCORDING TO EFFORT
sexbinPresence.NormEffort_Bin = sexbinPresence.Effort_Bin./sexbinPresence.MaxEffort_Bin; %what proportion of the day was there effort
sexbinPresence.NormEffort_Sec = sexbinPresence.Effort_Sec./sexbinPresence.MaxEffort_Sec; %what proportion of the day was there effort
sexbinPresence.SocialGroupNormBin = round(sexbinPresence.SocialGroup ./ sexbinPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for SocialGroups
sexbinPresence.MidSizeNormBin = round(sexbinPresence.MidSize ./ sexbinPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for MidSizes
sexbinPresence.MaleNormBin = round(sexbinPresence.Male ./ sexbinPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for Males
sexbinPresence.SocialGroupHoursNorm = round(sexbinPresence.SocialGroupNormBin ./ sexbinPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
sexbinPresence.MidSizeHoursNorm = round(sexbinPresence.MidSizeNormBin ./ sexbinPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
sexbinPresence.MaleHoursNorm = round(sexbinPresence.MaleNormBin ./ sexbinPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
%% Find month and season
[p,~]=size(sexbinPresence);
sexbinPresence.Season = zeros(p,1);
sexbinPresence.month = month(sexbinPresence.tbin);

%Winter starts on January (closest to the real thing, which is Dec. 21st)
summeridxD = (sexbinPresence.month == 7  | sexbinPresence.month == 8 | sexbinPresence.month == 9);
fallidxD = (sexbinPresence.month == 10  | sexbinPresence.month == 11 | sexbinPresence.month == 12);
winteridxD = (sexbinPresence.month == 1  | sexbinPresence.month == 2 | sexbinPresence.month == 3);
springidxD = (sexbinPresence.month == 4  | sexbinPresence.month == 5 | sexbinPresence.month == 6);

sexbinPresence.Season(summeridxD) = 1;
sexbinPresence.Season(fallidxD) = 2;
sexbinPresence.Season(winteridxD) = 3;
sexbinPresence.Season(springidxD) = 4;

%add year and day to data
sexbinPresence.Year = year(sexbinPresence.tbin);
sexbinPresence.day = day(sexbinPresence.tbin,'dayofyear');
%% Replace NANs where there was recording effort but no detections (this usually happens in the beginning or end of the recordings)
NANidx = ismissing(sexbinPresence(:,{'SocialGroupNormBin'}));
sexbinPresence{:,{'SocialGroupNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
sexbinPresence{:,{'FeHoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero
sexbinPresence{:,{'SocialGroupHoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero

NANidx = ismissing(sexbinPresence(:,{'MidSizeNormBin'}));
sexbinPresence{:,{'MidSizeNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
sexbinPresence{:,{'JuHoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero
sexbinPresence{:,{'MidSizeHoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero

NANidx = ismissing(sexbinPresence(:,{'MaleNormBin'}));
sexbinPresence{:,{'MaleNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
sexbinPresence{:,{'MaHoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero
sexbinPresence{:,{'MaleHoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero

writetable(timetable2table(sexbinPresence), [saveDir,'\', siteabrev, '_binPresence.csv']); %table with bin presence for each sex (timeseries)
%% Day table with days grouped together (summed and averaged)
[MD,~] = findgroups(sexbinPresence.day);

if length(MD) < 365
    meantab365 = table(sexbinPresence.day(:), sexbinPresence.FeHoursProp(:),sexbinPresence.JuHoursProp(:),sexbinPresence.MaHoursProp(:));
    meantab365.Properties.VariableNames = {'Day' 'HoursPropFE' 'HoursPropJU' 'HoursPropMA'};
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
    writetable(meantab365, [saveDir,'\',siteabrev,'_365GroupedMean.csv']); %table with the mean for each day of the year
else
    sexbinPresence.day = categorical(sexbinPresence.day);
    %mean for SocialGroups
    [meann, sem, std, var, range] = grpstats(sexbinPresence.FeHoursProp, sexbinPresence.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
    meantable = array2table(meann);
    semtable = array2table(sem);
    stdtable = array2table(std);
    vartable = array2table(var);
    rangetable = array2table(range);
    newcol_mean = (1:length(meann))';
    meanarrayFE365 = [newcol_mean meann sem std var range];
    meantabFE365 = array2table(meanarrayFE365);
    meantabFE365.Properties.VariableNames = {'Day' 'HoursPropFE' 'SEM' 'Std' 'Var' 'Range'};
    %mean for MidSizes
    [meann, sem, std, var, range] = grpstats(sexbinPresence.JuHoursProp, sexbinPresence.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
    meantable = array2table(meann);
    semtable = array2table(sem);
    stdtable = array2table(std);
    vartable = array2table(var);
    rangetable = array2table(range);
    newcol_mean = (1:length(meann))';
    meanarrayJU365 = [newcol_mean meann sem std var range];
    meantabJU365 = array2table(meanarrayJU365);
    meantabJU365.Properties.VariableNames = {'Day' 'HoursPropJU' 'SEM' 'Std' 'Var' 'Range'};
    %mean for males
    [meann, sem, std, var, range] = grpstats(sexbinPresence.MaHoursProp, sexbinPresence.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
    meantable = array2table(meann);
    semtable = array2table(sem);
    stdtable = array2table(std);
    vartable = array2table(var);
    rangetable = array2table(range);
    newcol_mean = (1:length(meann))';
    meanarrayMA365 = [newcol_mean meann sem std var range];
    meantabMA365 = array2table(meanarrayMA365);
    meantabMA365.Properties.VariableNames = {'Day' 'HoursPropMA' 'SEM' 'Std' 'Var' 'Range'};    

[pp,~]=size(meantabFE365);
meantabFE365.Season = zeros(pp,1);
meantabFE365.month = month(meantabFE365.Day);

[pp,~]=size(meantabJU365);
meantabJU365.Season = zeros(pp,1);
meantabJU365.month = month(meantabJU365.Day);

[pp,~]=size(meantabMA365);
meantabMA365.Season = zeros(pp,1);
meantabMA365.month = month(meantabMA365.Day);

%Winter starts on January 1st
summeridxD = (meantabFE365.month == 7  | meantabFE365.month == 8 | meantabFE365.month == 9);
fallidxD = (meantabFE365.month == 10  | meantabFE365.month == 11 | meantabFE365.month == 12);
winteridxD = (meantabFE365.month == 1  | meantabFE365.month == 2 | meantabFE365.month == 3);
springidxD = (meantabFE365.month == 4  | meantabFE365.month == 5 | meantabFE365.month == 6);

%adds the season according to the month the data was collected
meantabFE365.Season(summeridxD) = 1;
meantabFE365.Season(fallidxD) = 2;
meantabFE365.Season(winteridxD) = 3;
meantabFE365.Season(springidxD) = 4;

%Winter starts on January 1st
summeridxD = (meantabJU365.month == 7  | meantabJU365.month == 8 | meantabJU365.month == 9);
fallidxD = (meantabJU365.month == 10  | meantabJU365.month == 11 | meantabJU365.month == 12);
winteridxD = (meantabJU365.month == 1  | meantabJU365.month == 2 | meantabJU365.month == 3);
springidxD = (meantabJU365.month == 4  | meantabJU365.month == 5 | meantabJU365.month == 6);

%adds the season according to the month the data was collected
meantabJU365.Season(summeridxD) = 1;
meantabJU365.Season(fallidxD) = 2;
meantabJU365.Season(winteridxD) = 3;
meantabJU365.Season(springidxD) = 4;

%Winter starts on January 1st
summeridxD = (meantabMA365.month == 7  | meantabMA365.month == 8 | meantabMA365.month == 9);
fallidxD = (meantabMA365.month == 10  | meantabMA365.month == 11 | meantabMA365.month == 12);
winteridxD = (meantabMA365.month == 1  | meantabMA365.month == 2 | meantabMA365.month == 3);
springidxD = (meantabMA365.month == 4  | meantabMA365.month == 5 | meantabMA365.month == 6);

%adds the season according to the month the data was collected
meantabMA365.Season(summeridxD) = 1;
meantabMA365.Season(fallidxD) = 2;
meantabMA365.Season(winteridxD) = 3;
meantabMA365.Season(springidxD) = 4;

writetable(meantabFE365, [saveDir,'\',siteabrev,'_365GroupedMeanSocialGroup.csv']); %table with the mean for each day of the year
writetable(meantabJU365, [saveDir,'\',siteabrev,'_365GroupedMeanMidSize.csv']); %table with the mean for each day of the year
writetable(meantabMA365, [saveDir,'\',siteabrev,'_365GroupedMeanMale.csv']); %table with the mean for each day of the year
end
%% Hourly binary Data for GAMGEE Models
% Extract hourly effort from workspace 2 dayTable that already accounts for duty cycle, etc.
hourlytab = timetable2table(hourlyTab);
sexHourEffort = hourlytab(:,{'tbin','Effort_Bin','Effort_Sec','MaxEffort_Bin','MaxEffort_Sec'});
sexHourEffort = table2timetable(sexHourEffort);

%hourly
Click = retime(sexDEdata(:,5:7),'hourly','sum'); % #5-min bins per hour
sexhourlyTab = synchronize(Click,sexHourEffort);

%SocialGroups
binidx_hourlyF = (sexhourlyTab.SocialGroup >= 1);
[y,~]=size(sexhourlyTab);
sexhourlyTab.PreAbsF = zeros(y,1);
sexhourlyTab.PreAbsF(binidx_hourlyF) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

%MidSizes
binidx_hourlyJ = (sexhourlyTab.MidSize >= 1);
sexhourlyTab.PreAbsJ = zeros(y,1);
sexhourlyTab.PreAbsJ(binidx_hourlyJ) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

%Males
binidx_hourlyM = (sexhourlyTab.Male >= 1);
sexhourlyTab.PreAbsM = zeros(y,1);
sexhourlyTab.PreAbsM(binidx_hourlyM) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

writetable(timetable2table(sexhourlyTab),[saveDir,'\',siteabrev,'_binData_forGAMGEE_sexClasses.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
%% Save workspace variable
save([saveDir,'\',siteabrev,'_workspaceStep3.mat']);