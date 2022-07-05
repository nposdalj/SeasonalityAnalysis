clearvars
close all

%% Parameters defined by user
filePrefix = 'BS'; % File name to match.
siteabrev = 'BS'; %abbreviation of site.
region = 'WAT'; %region
sp = 'Pm'; % your species code
itnum = '3'; % which iteration you are looking for
srate = 200; % sample rate
GDrive = 'H'; %Google Drive

effortXls = 'H:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis\NC\Pm_effort.xlsx'; % specify excel file with effort times

saveDir = 'H:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis'; % specify main directory to save files
exportDir = 'H:\My Drive\WAT_TPWS_metadataReduced\ComparingMethods';  % specify alternative directory to save files

%% load workspace with NP data
load([saveDir,'\',siteabrev,'\',siteabrev,'_workspace125.mat']);

%% group NP data by 5min bins
%group data by 5 minute bins
binTable_NP = synchronize(binData,binEffort);
binTable_NP.Properties.VariableNames{'bin'} = 'Effort_Bin';
binTable_NP.Properties.VariableNames{'sec'} = 'Effort_Sec';
binTable_NP.maxPP = [];
binidx1 = (binTable_NP.Count >= 1);
[y,~]=size(binTable_NP);
binTable_NP.PreAbs = zeros(y,1);
binTable_NP.PreAbs(binidx1) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 5min bin
%no effort bins are excluded 

%% load RC data
binTable_RC_allsites = readtable('H:\My Drive\WAT_TPWS_metadataReduced\ComparingMethods\SpermWhale_5minBin_RC.csv');
namesitedata_RC = [region,'_',siteabrev];
bins_RC = timetable(binTable_RC_allsites.Bin, binTable_RC_allsites.(namesitedata_RC));

%make edits to resolve anomalous issues
find(bins_RC.Var1 == Inf,1)
if strcmp(siteabrev, 'GS')
    bins_RC.Var1(227129) = NaN; %RC's data has an Inf here
end

%% combine data
bins_NP = timetable(binTable_NP.tbin,binTable_NP.Count);

timebins = synchronize(bins_NP,bins_RC);
timebins.Properties.VariableNames{'Var1_bins_NP'} = 'Count_NP';
timebins.Properties.VariableNames{'Var1_bins_RC'} = 'Count_RC';

figure(1) % Plot raw timebins data
stackedplot(timebins)

% Adjust timebins to only include period where both timeseries have data
start_RC = find(~isnan(timebins.Count_RC),1);
end_RC = find(~isnan(timebins.Count_RC),1, "last");
start_NP = find(~isnan(timebins.Count_NP),1);
end_NP = find(~isnan(timebins.Count_NP),1, "last");
startTime = max(start_RC, start_NP);
endTime = min(end_RC, end_NP);
timebins = timebins(startTime:endTime,:);

% Represent times with effort but no presence as 0, and represent times
% with no effort as NaN (according to gaps in RC's data)
timebins.Count_NP(isnan(timebins.Count_NP)) = 0;
timebins.Count_NP(isnan(timebins.Count_RC)) = NaN;
timebins = rmmissing(timebins); %remove periods without effort

figure(2) % Plot final timebins data
stackedplot(timebins)

%% OPTIONAL: Set NP bins w/ <50 clicks as 0 to compare w/ RC
rmNP50 = '1'; % Set 1 to set bins w/ <50 clicks as 0 in NP data; Set 0 to leave values intact

if rmNP50 == '1'
    timebins.Count_NP(timebins.Count_NP < 50) = 0;
end

%% Things to compare (compare plots and values)
NPdat(1) = length(timebins.Count_NP(timebins.Count_NP ~= 0)); % How many bins does RC have w/ presence vs us
RCdat(1) = length(timebins.Count_RC(timebins.Count_RC ~= 0)); % something I just thought of-- this will ignore nan's right???
NPdat(2) = mean(timebins.Count_NP); % Avg clicks at each site
RCdat(2) = mean(timebins.Count_RC);
daysduration = datenum(timebins.Time(end)) - datenum(timebins.Time(1)) + 1;
NPdat(3) = length(timebins.Count_NP(timebins.Count_NP ~= 0)) / daysduration / 7; % Avg number bins per week
RCdat(3) = length(timebins.Count_RC(timebins.Count_RC ~= 0)) / daysduration / 7;
NPdat(4) = sum(timebins.Count_NP) / daysduration / 7; % Avg number clicks per week
RCdat(4) = sum(timebins.Count_RC) / daysduration / 7;

summary_rows = {'Number of bins with presence','Average count of clicks per bin','Average number of bins per week','Average number clicks per week'};
expSummary = table(NPdat.',RCdat.','rowNames',summary_rows);
expSummary.Properties.VariableNames(1) = "NP";
expSummary.Properties.VariableNames(2) = "RC";
if rmNP50 == '1'
        writetable(expSummary,[exportDir,'/',siteabrev,'_MethodCompare_Summary_rmNP50.txt'],'Delimiter','\t','WriteRowNames',true)
else
    writetable(expSummary,[exportDir,'/',siteabrev,'_MethodCompare_Summary.txt'],'Delimiter','\t','WriteRowNames',true)
end

% Additional metrics (calculate below and log manually)
corrcoef(timebins.Count_NP,timebins.Count_RC, 'Rows', 'complete') % Correlation coef (R) of the two datasets at this site

%% Total number clicks per week [Data setup takes a while to run :)]
% I totally misunderstood "avg number clicks per week" but I spent so long on this so I refuse to throw it away :D
weeklydat = nan(4*52,3);
disp('Beginning compiling weekly totals, please be patient...')
try
for j = 2016:2019
    for i = 1:52
        k = 52*(j-2016) + i;
        if isempty(find(year(timebins.Time) == j & week(timebins.Time) == i, 1)) == 0
            weeklydat(k,1) = find(year(timebins.Time) == j & week(timebins.Time) == i, 1);
            weeklydat(k,2) = sum(timebins.Count_NP(weeklydat(k,1):(weeklydat(k,1)+2015)));
            weeklydat(k,3) = sum(timebins.Count_RC(weeklydat(k,1):(weeklydat(k,1)+2015)));
        else
        end
    end
    %disp(['Finished compiling weekly totals for ', char(string(j))])
    entertainme(1)
end
catch % If I had learned about catch three hours before I did I would've saved an incredible amount of time gah
    disp(['Finished compiling weekly totals for ', char(string(j)), '.'])
    disp('Final week is incomplete and will be ignored. Removing weeks that are incomplete / without data...')
    pause(2)
    weeklydat = rmmissing(weeklydat);
end
weeklysums = timetable(timebins.Time(weeklydat(:,1)),weeklydat(:,2),weeklydat(:,3));
weeklysums.Properties.VariableNames(1) = "NP";
weeklysums.Properties.VariableNames(2) = "RC";
disp('Done!')

maxy12 = max([weeklysums.NP;weeklysums.RC])*5/4; % set an appropriate y-limit for Figure 12

figure(12) % Plot total weekly click counts
set(gcf, 'Position', [100 100 1300 700])
subplot(4,1,1)                              % SUBPLOT 1: 2016
yearselect = year(weeklysums.Time)==2016;
plotweekly = bar(weeklysums.Time(yearselect), [weeklysums.NP(yearselect) weeklysums.RC(yearselect)], "BarWidth",1.5, "FaceAlpha",.75);
xlim([datetime('25-Dec-2015') datetime('31-Dec-2016')]); ylim([0 maxy12])
plotweekly(1).FaceColor = "#FF5964"; plotweekly(1).EdgeColor = "#FF5964";
plotweekly(2).FaceColor = "#35A7FF"; plotweekly(2).EdgeColor = "#35A7FF";
ylabel("Count")
legend("NP","RC")
title(['Weekly Total Click Count at ', siteabrev])

subplot(4,1,2)                              % SUBPLOT 2: 2017
yearselect = year(weeklysums.Time)==2017;
plotweekly = bar(weeklysums.Time(yearselect), [weeklysums.NP(yearselect) weeklysums.RC(yearselect)], "BarWidth",1.5, "FaceAlpha",.75);
xlim([datetime('25-Dec-2016') datetime('31-Dec-2017')]); ylim([0 maxy12])
plotweekly(1).FaceColor = "#FF5964"; plotweekly(1).EdgeColor = "#FF5964";
plotweekly(2).FaceColor = "#35A7FF"; plotweekly(2).EdgeColor = "#35A7FF";
ylabel("Count")

subplot(4,1,3)                              % SUBPLOT 3: 2018
yearselect = year(weeklysums.Time)==2018;
plotweekly = bar(weeklysums.Time(yearselect), [weeklysums.NP(yearselect) weeklysums.RC(yearselect)], "BarWidth",1.5, "FaceAlpha",.75);
xlim([datetime('25-Dec-2017') datetime('31-Dec-2018')]); ylim([0 maxy12])
plotweekly(1).FaceColor = "#FF5964"; plotweekly(1).EdgeColor = "#FF5964";
plotweekly(2).FaceColor = "#35A7FF"; plotweekly(2).EdgeColor = "#35A7FF";
ylabel("Count")

subplot(4,1,4)                              % SUBPLOT 4: 2019
yearselect = year(weeklysums.Time)==2019;
plotweekly = bar(weeklysums.Time(yearselect), [weeklysums.NP(yearselect) weeklysums.RC(yearselect)], "BarWidth",1.5, "FaceAlpha",.75);
xlim([datetime('25-Dec-2018') datetime('31-Dec-2019')]); ylim([0 maxy12])
plotweekly(1).FaceColor = "#FF5964"; plotweekly(1).EdgeColor = "#FF5964";
plotweekly(2).FaceColor = "#35A7FF"; plotweekly(2).EdgeColor = "#35A7FF";
ylabel("Count")
xlabel("Week")

if rmNP50 == '1'
    saveas(gcf,[exportDir,'\',siteabrev,'_TotalPresence_Weekly_rmNP50'],'png');
else
    saveas(gcf,[exportDir,'\',siteabrev,'_TotalPresence_Weekly'],'png');
end