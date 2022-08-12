clearvars
close all

%% Parameters defined by user
filePrefix = 'BC'; % File name to match. 
genderFileName = 'WAT_BC'; %File name to match gender file
siteabrev = 'BC'; %abbreviation of site
region = 'WAT';
sp = 'Pm'; % your species code
GDrive = 'I'; %Google Drive

effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx']; % specify excel file with effort times
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
dayBinCSV= [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\',siteabrev,'_dayData_forGLMR125.csv']; % specify csv document with general PM information
%% Load sex specific data from ICIgrams
filename = [saveDir,'\',genderFileName,'_',sp,'_gender.mat'];
load(filename);
sexData = binData; %change file name so it doesn't match general sperm whale data
%% load workspace 2
GDrive_correct = GDrive; % Preserve correct GDrive as it was entered above
load([saveDir,'\',siteabrev,'_workspaceStep2.mat']);

% Overwrite some path names
GDrive = GDrive_correct; %Correct GDrive if overwritten by loading workspace
effortXls(1) = GDrive;
saveDir(1) = GDrive;
tpwsPath(1) = GDrive;
%% Extract effort from workspace 2 dayTable that already accounts for duty cycle, etc.
sexBinEffort = dayBinTAB(:,{'tbin','Effort_Bin','Effort_Sec','MaxEffort_Bin','MaxEffort_Sec'});
sexBinEffort = table2timetable(sexBinEffort);
%% group data by days and add effort
%Retime sex data daily
sexDEdata = sexData;
sexData = retime(sexData,'daily','sum');
sexbinPresence = synchronize (sexData, sexBinEffort);
sexbinPresence.maxPP = [];
sexbinPresence.meanICI = [];
sexbinPresence.BigMale = [];
sexbinPresence.Other = [];
sexbinPresence.Count = [];
%% accounting for effort/dutycycle
%Proportion of hours
sexbinPresence.FeMinutes = sexbinPresence.Female *5;
sexbinPresence.JuMinutes = sexbinPresence.Juvenile *5;
sexbinPresence.MaMinutes = sexbinPresence.Male *5;
sexbinPresence.FeHours = sexbinPresence.FeMinutes ./60;
sexbinPresence.JuHours = sexbinPresence.JuMinutes ./60;
sexbinPresence.MaHours = sexbinPresence.MaMinutes ./60;
sexbinPresence.FeHoursProp = sexbinPresence.FeHours ./(sexbinPresence.Effort_Sec ./ (60*60));
sexbinPresence.JuHoursProp = sexbinPresence.JuHours ./(sexbinPresence.Effort_Sec ./ (60*60));
sexbinPresence.MaHoursProp = sexbinPresence.MaHours ./(sexbinPresence.Effort_Sec ./ (60*60));

%normalize bin counts for each sex based on bin effort
sexbinPresence.NormEffort_Bin = sexbinPresence.Effort_Bin./sexbinPresence.MaxEffort_Bin; %what proportion of the day was there effort
sexbinPresence.NormEffort_Sec = sexbinPresence.Effort_Sec./sexbinPresence.MaxEffort_Sec; %what proportion of the day was there effort
sexbinPresence.FemaleNormBin = round(sexbinPresence.Female ./ sexbinPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for Females
sexbinPresence.JuvenileNormBin = round(sexbinPresence.Juvenile ./ sexbinPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for Juveniles
sexbinPresence.MaleNormBin = round(sexbinPresence.Male ./ sexbinPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for Males
sexbinPresence.FemaleHoursNorm = round(sexbinPresence.FemaleNormBin ./ sexbinPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
sexbinPresence.JuvenileHoursNorm = round(sexbinPresence.JuvenileNormBin ./ sexbinPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
sexbinPresence.MaleHoursNorm = round(sexbinPresence.MaleNormBin ./ sexbinPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
%% find month and season
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

NANidx = ismissing(sexbinPresence(:,{'FemaleNormBin'}));
sexbinPresence{:,{'FemaleNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
sexbinPresence{:,{'FeHoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero
sexbinPresence{:,{'FemaleHoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero

NANidx = ismissing(sexbinPresence(:,{'JuvenileNormBin'}));
sexbinPresence{:,{'JuvenileNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
sexbinPresence{:,{'JuHoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero
sexbinPresence{:,{'JuvenileHoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero

NANidx = ismissing(sexbinPresence(:,{'MaleNormBin'}));
sexbinPresence{:,{'MaleNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
sexbinPresence{:,{'JuHoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero
sexbinPresence{:,{'JuvenileHoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero

writetable(timetable2table(sexbinPresence), [saveDir,'\', siteabrev, '_binPresence.csv']); %table with bin presence for each sex (timeseries)
%% day table with days grouped together (summed and averaged) ** USE THIS **
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
    %mean for females
    [mean, sem, std, var, range] = grpstats(sexbinPresence.FeHoursProp, sexbinPresence.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
    meantable = array2table(mean);
    semtable = array2table(sem);
    stdtable = array2table(std);
    vartable = array2table(var);
    rangetable = array2table(range);
    newcol_mean = (1:length(mean))';
    meanarrayFE365 = [newcol_mean mean sem std var range];
    meantabFE365 = array2table(meanarrayFE365);
    meantabFE365.Properties.VariableNames = {'Day' 'HoursPropFE' 'SEM' 'Std' 'Var' 'Range'};
    %mean for juveniles
    [mean, sem, std, var, range] = grpstats(sexbinPresence.JuHoursProp, sexbinPresence.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
    meantable = array2table(mean);
    semtable = array2table(sem);
    stdtable = array2table(std);
    vartable = array2table(var);
    rangetable = array2table(range);
    newcol_mean = (1:length(mean))';
    meanarrayJU365 = [newcol_mean mean sem std var range];
    meantabJU365 = array2table(meanarrayJU365);
    meantabJU365.Properties.VariableNames = {'Day' 'HoursPropJU' 'SEM' 'Std' 'Var' 'Range'};
    %mean for males
    [mean, sem, std, var, range] = grpstats(sexbinPresence.MaHoursProp, sexbinPresence.day, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
    meantable = array2table(mean);
    semtable = array2table(sem);
    stdtable = array2table(std);
    vartable = array2table(var);
    rangetable = array2table(range);
    newcol_mean = (1:length(mean))';
    meanarrayMA365 = [newcol_mean mean sem std var range];
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

writetable(meantabFE365, [saveDir,'\',siteabrev,'_365GroupedMeanFemale.csv']); %table with the mean for each day of the year
writetable(meantabJU365, [saveDir,'\',siteabrev,'_365GroupedMeanJuvenile.csv']); %table with the mean for each day of the year
writetable(meantabMA365, [saveDir,'\',siteabrev,'_365GroupedMeanMale.csv']); %table with the mean for each day of the year
end
%% Hourly binary Data for GAMGEE Models
% Extract hourly effort from workspace 2 dayTable that already accounts for duty cycle, etc.
hourlytab = timetable2table(hourlyTab);
sexHourEffort = hourlytab(:,{'tbin','Effort_Bin','Effort_Sec','MaxEffort_Bin','MaxEffort_Sec'});
sexHourEffort = table2timetable(sexHourEffort);

%hourly
Click = retime(sexDEdata(:,4:6),'hourly','sum'); % #5-min bins per hour
sexhourlyTab = synchronize(Click,sexHourEffort);

%Females
binidx_hourlyF = (sexhourlyTab.Female >= 1);
[y,~]=size(sexhourlyTab);
sexhourlyTab.PreAbsF = zeros(y,1);
sexhourlyTab.PreAbsF(binidx_hourlyF) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

%Juveniles
binidx_hourlyJ = (sexhourlyTab.Juvenile >= 1);
sexhourlyTab.PreAbsJ = zeros(y,1);
sexhourlyTab.PreAbsJ(binidx_hourlyJ) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

%Males
binidx_hourlyM = (sexhourlyTab.Male >= 1);
sexhourlyTab.PreAbsM = zeros(y,1);
sexhourlyTab.PreAbsM(binidx_hourlyM) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

writetable(timetable2table(sexhourlyTab),[saveDir,'\',siteabrev,'_binData_forGAMGEE_sexClasses.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
%% Save workspace variable
save([saveDir,'\',siteabrev,'_workspaceStep3.mat']);