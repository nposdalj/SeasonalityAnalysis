clearvars
close all
%% Parameters defined by user
filePrefix = 'GofAK_CB'; % File name to match
siteabrev = 'CB'; %abbreviation of site
sp = 'Pm'; % your species code
srate = 200; % sample rate
tpwsPath = 'E:\Project_Sites\CB\TPWS_125\TPWS2\'; %directory of TPWS files
saveDir = 'E:\Project_Sites\CB\Seasonality'; %specify directory to save files
effortXls = 'E:\Project_Sites\CB\Pm_Effort_CB.xlsx'; % specify excel file with effort times
dayBinCSV= 'E:\Project_Sites\CB\Seasonality\CB_dayData_forGLMR125.csv'; % specify csv document with general PM information
%% Get effort times matching prefix file
%when multiple sites in the effort table
allEfforts = readtable(effortXls); %read effort table
site = siteabrev; %abbreviation used in effort table

siteNUM = unique(allEfforts.Sites);
[sr,~] = size(siteNUM);

if sr > 1
    effTable = allEfforts(ismember(allEfforts.Sites,site),:); %effort is for multiple sites
else
    effTable = allEfforts; %effort is for one site only
end

% make Variable Names consistent
startVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'Start.*Effort'))>0,1,'first');
endVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'End.*Effort'))>0,1,'first');
effTable.Properties.VariableNames{startVar} = 'Start';
effTable.Properties.VariableNames{endVar} = 'End';

Start = datetime(x2mdate(effTable.Start),'ConvertFrom','datenum');
End = datetime(x2mdate(effTable.End),'ConvertFrom','datenum');

effort = table(Start,End);
%% get default parameters
p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin');
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
%% load data
filename = [tpwsPath,filePrefix,'_',sp,'_gender.mat'];
load(filename);
%% group data by days and add effort
binPresence = synchronize (binData, binEffort);
binPresence.Properties.VariableNames{'bin'} = 'Effort_Bin';
binPresence.Properties.VariableNames{'sec'} = 'Effort_Sec';
binPresence = retime(binPresence,'daily','sum');
binPresence.maxPP = [];
binPresence.meanICI = [];
binPresence.BigMale = [];
binPresence.Other = [];
binPresence.Count = [];
binPresence(~binPresence.Effort_Bin,:) = []; %removes days with no effort, NOT days with no presence
%% accounting for effort
[p,q]=size(binPresence);
binPresence.MaxEffort_Bin = ones(p,1)*(288);
binPresence.MaxEffort_Sec = ones(p,1) * (86400); %seconds in one day

%dealing with duty cycled data
if strcmp(siteabrev,'CB');
    %for CB02 ONLY - only .44 of each hour is recorded so effort of 5 min bins for each day is 127 bins
    ge = binPresence.Effort_Bin(222:507); %bin effort (excluding ships but not considering the duty cycle)
    ge = ge/288; %%proportion of data that was not 'ships' considering full recording effort
    binPresence.Effort_Bin(222:507) = ge * 127;
    binPresence.MaxEffort_Bin(222:507) = 127;
    binPresence.MaxEffort_Sec(222:507) = 127 * 5 * 60;
    binPresence.Effort_Sec(222:507) = binPresence.Effort_Bin(222:507) * 5 * 60;
    else
if strcmp(siteabrev,'BD');
    %for ALEUT03BD ONLY - only 0.33 of each hour is recorded so effort of 5 min bins for each day is 96
    ge = binPresence.Effort_Bin(274:end); %bin effort (excluding ships but not considering the duty cycle)
    ge = ge/288; %%proportion of data that was not 'ships' considering full recording effort
    binPresence.Effort_Bin(274:end) = ge * 96;
    binPresence.MaxEffort_Bin(274:end) = 96;
    binPresence.MaxEffort_Sec(274:end) = 96 * 5 * 60;
    binPresence.Effort_Sec(274:end) = binPresence.Effort_Bin(274:end) * 5 * 60;
    else
dayTable.MaxEffort_Bin = ones(p,1)*(288);
end
end

%Proportion of hours
binPresence.FeMinutes = binPresence.Female *5;
binPresence.JuMinutes = binPresence.Juvenile *5;
binPresence.MaMinutes = binPresence.Male *5;
binPresence.FeHours = binPresence.FeMinutes ./60;
binPresence.JuHours = binPresence.JuMinutes ./60;
binPresence.MaHours = binPresence.MaMinutes ./60;
binPresence.FeHoursProp = binPresence.FeHours ./(binPresence.Effort_Sec ./ (60*60));
binPresence.JuHoursProp = binPresence.JuHours ./(binPresence.Effort_Sec ./ (60*60));
binPresence.MaHoursProp = binPresence.MaHours ./(binPresence.Effort_Sec ./ (60*60));

%normalize bin counts for each sex based on bin effort
binPresence.NormEffort_Bin = binPresence.Effort_Bin./binPresence.MaxEffort_Bin; %what proportion of the day was there effort
binPresence.NormEffort_Sec = binPresence.Effort_Sec./binPresence.MaxEffort_Sec; %what proportion of the day was there effort
binPresence.FemaleNormBin = round(binPresence.Female ./ binPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for Females
binPresence.JuvenileNormBin = round(binPresence.Juvenile ./ binPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for Juveniles
binPresence.MaleNormBin = round(binPresence.Male ./ binPresence.NormEffort_Bin); %what would the normalized bin count be given the amount of effort for Males
binPresence.FemaleHoursNorm = round(binPresence.FemaleNormBin ./ binPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
binPresence.JuvenileHoursNorm = round(binPresence.JuvenileNormBin ./ binPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
binPresence.MaleHoursNorm = round(binPresence.MaleNormBin ./ binPresence.NormEffort_Bin); %convert the number of 5-min bins per day to hours
%% find month and season
binPresence.Season = zeros(p,1);
binPresence.month = month(binPresence.tbin);

%Winter starts on January (closest to the real thing, which is Dec. 21st)
summeridxD = (binPresence.month == 7  | binPresence.month == 8 | binPresence.month == 9);
fallidxD = (binPresence.month == 10  | binPresence.month == 11 | binPresence.month == 12);
winteridxD = (binPresence.month == 1  | binPresence.month == 2 | binPresence.month == 3);
springidxD = (binPresence.month == 4  | binPresence.month == 5 | binPresence.month == 6);

binPresence.Season(summeridxD) = 1;
binPresence.Season(fallidxD) = 2;
binPresence.Season(winteridxD) = 3;
binPresence.Season(springidxD) = 4;

%add year and day to data
binPresence.Year = year(binPresence.tbin);
binPresence.day = day(binPresence.tbin,'dayofyear');

NANidx = ismissing(binPresence(:,{'FemaleNormBin'}));
binPresence{:,{'FemaleNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero

NANidx = ismissing(binPresence(:,{'JuvenileNormBin'}));
binPresence{:,{'JuvenileNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero

NANidx = ismissing(binPresence(:,{'MaleNormBin'}));
binPresence{:,{'MaleNormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
%% day table with days grouped together (summed and averaged) ** USE THIS **
[MD,~] = findgroups(binPresence.day);

if length(MD) < 365
    meantab365 = table(binPresence.day(:), binPresence.FemaleNormBin(:),binPresence.JuvenileNormBin(:),binPresence.MaleNormBin(:));
    meantab365.Properties.VariableNames = {'Day' 'Female' 'Juvenile' 'Male'};
    sumtab365 = meantab365;
else
    binPresence.day = categorical(binPresence.day);
    %mean
    Fmeanarray = grpstats(binPresence.FemaleNormBin, binPresence.day, @mean); %takes the mean of each day of the year
    Jmeanarray = grpstats(binPresence.JuvenileNormBin, binPresence.day, @mean); %takes the mean of each day of the year
    Mmeanarray = grpstats(binPresence.MaleNormBin, binPresence.day, @mean); %takes the mean of each day of the year
    combinedMeanArray = [Fmeanarray Jmeanarray Mmeanarray];
    meantable = array2table(combinedMeanArray);
    newcol_mean = (1:length(Fmeanarray))';
    meanarray365 = [newcol_mean combinedMeanArray];
    meantab365 = array2table(meanarray365);
    meantab365.Properties.VariableNames = {'Day' 'Female' 'Juvenile' 'Male'};
end

%writetable(meantab365, [saveDir,'\',siteabrev,'_365GroupedMeanSexes.csv']); %table with the mean for each day of the year
%writetable(sumtab365, [saveDir,'\',siteabrev,'_365GroupedSumSexes.csv']); %table with the sum for each day of the year
%% Plotting

%Plot daily presence in 5-min bins for each class seperately
figure
subplot(3,1,1)
bar(binPresence.tbin,binPresence.FemaleNormBin,'FaceColor','y','BarWidth',3)
title(['Daily Presence of Social Units in the ',titleNAME])
subplot(3,1,2)
bar(binPresence.tbin,binPresence.JuvenileNormBin,'FaceColor','b','BarWidth',3)
ylabel('Daily Presence (5-min bins)')
ylim([0 50])
title(['Daily Presence of Mid-Size Animals in the ',titleNAME])
subplot(3,1,3)
bar(binPresence.tbin,binPresence.MaleNormBin,'FaceColor','c','BarWidth',3)
title(['Daily Presence of Males in the ',titleNAME])
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses_Subplots130.png']);

%Plot daily presence in 5-min bins for all classes in one plot
figure
bar(binPresence.tbin,binPresence.MaleNormBin,'FaceColor','c','BarWidth',1)
ylabel('Daily Presence (5-min bins)')
title(['Daily Presence of Each Size Class at ',titleNAME])
hold on
bar(binPresence.tbin,binPresence.JuvenileNormBin,'FaceColor','b','BarWidth',1)
bar(binPresence.tbin,binPresence.FemaleNormBin,'FaceColor','y','BarWidth',1)
legend('Males','Mid-Size','Social Units')
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses130.png']);

%Plotting sperm whale presence with presence of each size class over it
dayBinTAB = readtable(dayBinCSV); %read general PM table with presence information
figure
bar(dayBinTAB.tbin,dayBinTAB.NormBin,'k')
hold on
bar(binPresence.tbin,binPresence.MaleNormBin,'FaceColor','c','BarWidth',1)
bar(binPresence.tbin,binPresence.JuvenileNormBin,'FaceColor','b','BarWidth',1)
bar(binPresence.tbin,binPresence.FemaleNormBin,'FaceColor','y','BarWidth',1)
ylabel('Daily Presence (5-min bins)')
title(['Daily Presence with Each Size Class at ',titleNAME])
legend('All','Males','Mid-Size','Social Units')
saveas(gcf,[saveDir,'\',siteabrev,'AllDailyPresence_with_SizeClasses130.png']);

%% save binPresence table with M/F/J 
writetimetable(binPresence, [saveDir,'\',siteabrev,'_binPresence130.csv']); %table with the mean for each day of the year

