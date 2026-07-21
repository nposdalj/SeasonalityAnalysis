 %% clear workspace
clear all;close all;clc;
% Needs to be run in 2018B or later
% AZORES_Plots_All - single site per run (edit siteName and rerun for the other AZORES site).
% Matches the flat F:\AZORES\SeasonalityAnalysis\<siteName>\ layout used by the AZORES Step1/Step2/Step3 scripts.
%
% Produces, for the selected site:
%   1) A CSV export of the 5-min bin size-class data (binData_mod), with a companion README
%   2) The bubble chart time series (general pattern + sex classes; both the
%      "general at top" and "difference/uncategorized" versions)
%   3) A stacked weekly timeseries (Social Groups / Mid-size / Adult Males in one panel,
%      as an alternative to the bubble chart)
%   4) Euler-style overlap diagrams comparing daily and hourly presence of the three classes
%   5) Overlapping ICI histograms comparing the ICI distribution of each size class
%
% Same color scheme throughout (matches OverlappingHistograms_Step1.m):
%   Social Groups = #66c2a5 (green), Mid-size = #fc8d62 (orange), Adult Males = #8da0cb (blue)
%% load data
siteName = 'AZORES_B_01'; % or 'AZORES_B_01' 'AZORES_A_04_DEEP'
siteNameDisp = strrep(siteName,'_','\_'); % escaped for use in titles - MATLAB's default tex interpreter
% treats a bare '_' as a subscript marker, which mangles siteName in plot titles
NumBub = 3;
DataDir = ['F:\AZORES\SeasonalityAnalysis\',siteName];
saveDirectory = DataDir; %plots/CSVs are saved alongside the site's data - no separate Plots tree for AZORES
Scale = 0; %(1) all scales are based on the lowest 'max' value (0) each size class has a respective scale

%Color scheme used throughout this script (matches the bubble chart / OverlappingHistograms convention)
colorSG = '#66c2a5'; %Social Groups
colorMS = '#fc8d62'; %Mid-size
colorMA = '#8da0cb'; %Adult Males
%% Load Workspace Step 2
load([DataDir,'\',siteName,'_workspaceStep2.mat']);
clear mean c

datesALL = unique(dayTable.Year)'; %unique dates for plot
dates = datesALL(1):datesALL(end); %fill in missing years

%Add first week of the year and last to account for x's
%Find the first and last day of the year
firstDay = datetime(dayTable.Year(1),1,1,0,0,0);
lastDay = datetime(dayTable.Year(end),12,31,0,0,0);
tbin = [firstDay; dayTable.tbin; lastDay];
ZeroCol = zeros(length(tbin), 1);
allDayTable = timetable(tbin,ZeroCol);
allDays = retime(allDayTable,'daily','fillwithmissing');
dayTableFull = synchronize(allDays,dayTable);
dayTableFull = removevars(dayTableFull, {'ZeroCol','Season','month','Year','day'});

%General
WeekData = retime(dayTableFull, 'weekly',@(x) mean(x, 'omitnan'));
WeekData.Noeffort = isnan(WeekData.Effort_Bin);
WeekData.Year = year(WeekData.tbin);
WeekData.wN = weeknum(WeekData.tbin,2,1); % Day of week begins on Monday and uses European standard
%deal with leap years, merge last day to week 52
leap = find(WeekData.wN == 53);
for l = 1:length(leap)
   WeekData(leap(l)-1,6:7) = varfun(@(x)mean(x,'omitnan'),WeekData(leap(l)-1:leap(l),6:7));
   WeekData(leap(l),:) = [];
end
% correct year for leap year
WeekData.diffYear = zeros(height(WeekData),1);
WeekData.diffYear(2:end) = diff(WeekData.Year);
adjYear = find(WeekData.diffYear == 1 & WeekData.wN == 52);
WeekData.Year(adjYear) = WeekData.Year(adjYear-1);
% convert bins to minutes
WeekData.NormBin = WeekData.NormBin *5;
%% Load Workspace Step 3
load([DataDir,'\',siteName,'_workspaceStep3.mat']);
clear mean
% NOTE: this load also brings binData_mod (the 5-min bin size-class table)
% into the workspace, since Step 3 loads it and never clears it before its
% own blanket "save workspaceStep3" - reused below for the CSV export and ICI histograms.

sexbinPresence.day = grp2idx(sexbinPresence.day);
binPresenceAll = synchronize(allDays,sexbinPresence);
binPresenceAll = removevars(binPresenceAll, {'ZeroCol'});
dataWeek = retime(binPresenceAll,'weekly',@(x) mean(x, 'omitnan'));
dataWeek.Noeffort = isnan(dataWeek.Effort_Bin);
dataWeek.Year = year(dataWeek.tbin);
dataWeek.wN = weeknum(dataWeek.tbin,2,1); % Day of week begins on Monday and uses European standard
%deal with leap years, merge last day to week 52
leap = find(dataWeek.wN == 53);
for l = 1:length(leap)
   dataWeek(leap(l)-1,6:7) = varfun(@nanmean,dataWeek(leap(l)-1:leap(l),6:7));
   dataWeek(leap(l),:) = [];
end
% correct year for leap year
dataWeek.diffYear = zeros(height(dataWeek),1);
dataWeek.diffYear(2:end) = diff(dataWeek.Year);
adjYear = find(dataWeek.diffYear == 1 & dataWeek.wN == 52);
dataWeek.Year(adjYear) = dataWeek.Year(adjYear-1);

% convert bins to minutes
dataWeek.SocialGroupNormBin = dataWeek.SocialGroupNormBin *5;
dataWeek.MidSizeNormBin = dataWeek.MidSizeNormBin *5;
dataWeek.MaleNormBin = dataWeek.MaleNormBin *5;
%% Delete first week for specific sites so they don't show up in the plots
dataWeek.year = year(dataWeek.tbin);
[DWyear,occurrenceDW] = unique(dataWeek.year);
if height(DWyear) > 1
if occurrenceDW(2) == 2
    dataWeek(1,:) = [];
end
dataWeek.Year = year(dataWeek.tbin);

WeekData.year = year(WeekData.tbin);
[WDyear,occurrenceWD] = unique(WeekData.year);
if occurrenceWD(2) == 2
    WeekData(1,:) = [];
end
end
WeekData.year = year(WeekData.tbin);
%% Checking to see how much was missed
CombinedWeek = dataWeek(:,19:21);
CombinedWeek.NormBin = WeekData.NormBin;
CombinedWeek.Added = CombinedWeek.SocialGroupNormBin + CombinedWeek.MidSizeNormBin + CombinedWeek.MaleNormBin;
CombinedWeek.Difference = CombinedWeek.NormBin - CombinedWeek.Added;
CombinedWeek.Difference( CombinedWeek.Difference <= 0 ) = 0;
%% Find Max Bubble Size and Distribution
%(1) Either use max and min values
Femax = max(dataWeek.SocialGroupNormBin);
Jumax = max(dataWeek.MidSizeNormBin);
Mamax = max(dataWeek.MaleNormBin);
Comax = max(CombinedWeek.Difference);
Gemax = max(CombinedWeek.NormBin);
MAXdiff = max([Femax, Jumax, Mamax, Comax]);
MINdiff = min([Femax, Jumax, Mamax, Comax]);
MAXgen = max([Femax, Jumax, Mamax, Gemax]);
MINgen = min([Femax, Jumax, Mamax, Gemax]);

% (2) Use the distribution of all of the values and ignore bins with very
% little data
TotalDataDiff = [dataWeek.SocialGroupNormBin; dataWeek.MidSizeNormBin; dataWeek.MaleNormBin; CombinedWeek.Difference]; %combine values
TotalDataGen = [dataWeek.SocialGroupNormBin; dataWeek.MidSizeNormBin; dataWeek.MaleNormBin; CombinedWeek.NormBin]; %combine values

%delete NAs and zeros for mean calculation
%Social groups
TotalDataSG = dataWeek.SocialGroupNormBin;
nanRows = isnan(TotalDataSG);
zeroRows = TotalDataSG==0;
badRows = nanRows | zeroRows;
TotalDataSG(badRows) = [];
%Mid-Size
TotalDataMS = dataWeek.MidSizeNormBin;
nanRows = isnan(TotalDataMS);
zeroRows = TotalDataMS==0;
badRows = nanRows | zeroRows;
TotalDataMS(badRows) = [];
%Males
TotalDataM = dataWeek.MaleNormBin;
nanRows = isnan(TotalDataM);
zeroRows = TotalDataM==0;
badRows = nanRows | zeroRows;
TotalDataM(badRows) = [];
%General
TotalDataG = CombinedWeek.NormBin;
nanRows = isnan(TotalDataG);
zeroRows = TotalDataG==0;
badRows = nanRows | zeroRows;
TotalDataG(badRows) = [];
%Difference
TotalDataD = CombinedWeek.Difference;
nanRows = isnan(TotalDataD);
zeroRows = TotalDataD==0;
badRows = nanRows | zeroRows;
TotalDataD(badRows) = [];
%All data with diff
nanRows = isnan(TotalDataDiff);
zeroRows = TotalDataDiff==0;
badRows = nanRows | zeroRows;
TotalDataDiff(badRows) = [];
%All data with gen
nanRows = isnan(TotalDataGen);
zeroRows = TotalDataGen==0;
badRows = nanRows | zeroRows;
TotalDataGen(badRows) = [];

%Find the histogram bins and set a cutoff
%Difference
figure, hist(TotalDataDiff,20);
[N,edges] = histcounts(TotalDataDiff,20);
lessThan = find(N<3); %find bins with less than 5
CutOffDiff = edges(lessThan(1)); %values over this value will be considered 'greater than' in the bubble sizes

%General
figure, hist(TotalDataGen,20);
[N,edges] = histcounts(TotalDataGen,20);
lessThan = find(N<3); %find bins with less than 5
CutOffGen = edges(lessThan(1)); %values over this value will be considered 'greater than' in the bubble sizes
%% Edit by AD: Manually mark periods of no effort for WAT sites
if strcmp(siteName, 'NC')
    WeekData.Noeffort(WeekData.Count_Click == 0) = 1;
    dataWeek.Noeffort(WeekData.Count_Click == 0) = 1;
elseif strcmp(siteName, 'BC')
    WeekData.Noeffort(WeekData.Count_Click == 0) = 1;
    dataWeek.Noeffort(WeekData.Count_Click == 0) = 1;
end
% AD: This adjustment may need to be made for other sites too, not sure
%% ===================================================================
%% Save 5-min bin size-class data as CSV, with a README
%% ===================================================================
binDataExport = timetable2table(binData_mod);
binDataExport = binDataExport(:,{'tbin','Count','mean_ICISel','SocialGroup','MidSize','Male'});

csvName = [siteName,'_5minBinData_SizeClass.csv'];
writetable(binDataExport,fullfile(saveDirectory,csvName));

readmeName = [siteName,'_5minBinData_SizeClass_README.txt'];
readmefile = fullfile(saveDirectory,readmeName);
fileid = fopen(readmefile, 'w');
fclose(fileid);
fileid = fopen(readmefile, 'at');
fprintf(fileid, [...
    '5-Minute Bin Size-Class Data - Column Descriptions\n' ...
    '====================================================\n' ...
    'File: ' csvName '\n' ...
    'Source: binData_mod (ICIgram_max classification, filled in by ICIgram_PostProcessing.m)\n\n' ...
    'tbin         - Start time of the 5-minute bin\n' ...
    'Count        - Number of validated clicks in this bin (ICI within the accepted threshold)\n' ...
    'mean_ICISel  - Mean inter-click interval (ms) of the validated clicks in this bin\n' ...
    'SocialGroup  - 1 if this bin was classified as Social Group presence, 0 otherwise\n' ...
    'MidSize      - 1 if this bin was classified as Mid-Size presence, 0 otherwise\n' ...
    'Male         - 1 if this bin was classified as Adult Male presence, 0 otherwise\n\n' ...
    'Note: maxPP and avgPeakFr (received level / peak frequency) are excluded from this\n' ...
    'export because AZORES TPWS files only contain click times (no MPP/MSP), so these\n' ...
    'columns are NaN for every bin.\n' ...
    'Note: OtherA/OtherB (undefined-class flags from ICIgram_max) are excluded since they\n' ...
    'were not used in classification for this dataset.\n']);
fclose(fileid);
%% Plot data with 'all clicks/general pattern' at the top of the subplot
%No Social Group Data
if nansum(dataWeek.SocialGroup) == 0
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
subplot(3,1,1)
years = unique(WeekData.Year);
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
black = [0,0,0];
absence = WeekData.NormBin == 0 & WeekData.Noeffort == 0;
for y = 1:length(years)
    hold on
    idxYear = WeekData.Year == years(y);
    idxNoeffort = WeekData.Noeffort == 1;
    scatter(WeekData.wN(idxYear & idxNoeffort),WeekData.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(WeekData.wN(idxYear & absence),WeekData.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',black)
    bubblechart(WeekData.wN(idxYear & keep),WeekData.Year(idxYear & keep),...
        round(WeekData.NormBin(idxYear & keep)),black)
end
bubblesize([2 15])
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffGen)]);
else
    bubblelim([1 round(max(WeekData.NormBin))]);
end
xlim([0,53])
ylabel('General')
set(gca,'xticklabel',[])
yticks(dates)

subplot(3,1,2)
years = unique(dataWeek.Year);
blue = colorMS;
keep = ~isnan(dataWeek.MidSizeNormBin) & dataWeek.MidSizeNormBin > 0;
absence = dataWeek.MidSizeNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.MidSizeNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffGen)]);
else
    bubblelim([1 round(max(dataWeek.MidSizeNormBin))]);
end
xlim([0,53])
ylabel('Mid-size')
set(gca,'xticklabel',[])
yticks(dates)

subplot(3,1,3)
blue = colorMA;
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
absence = dataWeek.MaleNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.MaleNormBin(idxYear & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffGen)]);
else
    bubblelim([1 round(max(dataWeek.MaleNormBin))]);
end
xlim([0,53])
xlabel('Week of the year')
ylabel('Adult Males')
yticks(dates)

else

fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
years = unique(dataWeek.Year);
% General Presence
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
subplot(4,1,1)
black = [0,0,0];
absence = WeekData.NormBin == 0 & WeekData.Noeffort == 0; % Change by AD
for y = 1:length(years)
    hold on
    idxYear = WeekData.Year == years(y);
    idxNoeffort = WeekData.Noeffort == 1;
    scatter(WeekData.wN(idxYear & idxNoeffort),WeekData.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(WeekData.wN(idxYear & absence),WeekData.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',black)
    bubblechart(WeekData.wN(idxYear & keep),WeekData.Year(idxYear & keep),...
        round(WeekData.NormBin(idxYear & keep)),black)
end
bubblesize([2 15])
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffGen)]);
else
    bubblelim([1 round(max(WeekData.NormBin))]);
end
xlim([0,53])
ylabel('General')
set(gca,'xticklabel',[])
yticks(dates)

% Social Group
subplot(4,1,2)
blue = colorSG;
keep = ~isnan(dataWeek.SocialGroupNormBin) & dataWeek.SocialGroupNormBin > 0;
absence = dataWeek.SocialGroupNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.SocialGroupNormBin(idxYear & keep)),blue)
end
bubblesize([2 15])
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffGen)]);
else
    bubblelim([1 round(max(dataWeek.SocialGroupNormBin))]);
end
xlim([0,53])
ylabel('Social Groups')
set(gca,'xticklabel',[])
yticks(dates)

subplot(4,1,3)
blue = colorMS;
keep = ~isnan(dataWeek.MidSizeNormBin) & dataWeek.MidSizeNormBin > 0;
absence = dataWeek.MidSizeNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.MidSizeNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffGen)]);
else
    bubblelim([1 round(max(dataWeek.MidSizeNormBin))]);
end
xlim([0,53])
ylabel('Mid-size')
set(gca,'xticklabel',[])
yticks(dates)

subplot(4,1,4)
blue = colorMA;
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
absence = dataWeek.MaleNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.MaleNormBin(idxYear & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffGen)]);
else
    bubblelim([1 round(max(dataWeek.MaleNormBin))]);
end
xlim([0,53])
xlabel('Week of the year')
ylabel('Adult Males')
yticks(dates)
end
%% save plot
set(gcf,'Position',[-1165         552         812         476])
if Scale == 1
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesScaled.png'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesScaled.pdf'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
else
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeries.png'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeries.pdf'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
end
%% Plotting the difference instead of the 'general pattern'
%No Social Group Data
if nansum(dataWeek.SocialGroup) == 0
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
subplot(3,1,1)
years = unique(dataWeek.year);
blue = colorMS;
keep = ~isnan(dataWeek.MidSizeNormBin) & dataWeek.MidSizeNormBin > 0;
absence = dataWeek.MidSizeNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.MidSizeNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffDiff)]);
else
    bubblelim([1 round(max(dataWeek.MidSizeNormBin))]);
end
xlim([0,53])
ylabel('Mid-size')
set(gca,'xticklabel',[])
yticks(dates)

subplot(3,1,2)
blue = colorMA;
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
absence = dataWeek.MaleNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.MaleNormBin(idxYear & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffDiff)]);
else
    bubblelim([1 round(max(dataWeek.MaleNormBin))]);
end
xlim([0,53])
ylabel('Adult Males')
yticks(dates)

subplot(3,1,3)
years = unique(WeekData.year);
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
black = '#C0C0C0'; %silver
absence = WeekData.NormBin == 0 & WeekData.Noeffort == 0; % Change by AD
for y = 1:length(years)
    hold on
    idxYear = WeekData.Year == years(y);
    idxNoeffort = WeekData.Noeffort == 1;
    scatter(WeekData.wN(idxYear & idxNoeffort),WeekData.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(WeekData.wN(idxYear & absence),WeekData.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',black)
    bubblechart(WeekData.wN(idxYear & keep),WeekData.Year(idxYear & keep),...
        round(WeekData.NormBin(idxYear & keep)),black)
end
bubblesize([2 8])
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffDiff)]);
else
    bubblelim([1 round(max(CombinedWeek.Difference))]);
end
xlim([0,53])
ylabel('Uncategorized')
xlabel('Week of the year')
set(gca,'xticklabel',[])
yticks(dates)

else

fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
years = unique(dataWeek.Year);
% Social Group
subplot(4,1,1)
blue = colorSG;
keep = ~isnan(dataWeek.SocialGroupNormBin) & dataWeek.SocialGroupNormBin > 0;
absence = dataWeek.SocialGroupNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.SocialGroupNormBin(idxYear & keep)),blue)
end
bubblesize([2 15])
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffDiff)]);
else
    bubblelim([1 round(max(dataWeek.SocialGroupNormBin))]);
end
xlim([0,53])
ylabel('Social Groups')
set(gca,'xticklabel',[])
yticks(dates)

subplot(4,1,2)
blue = colorMS;
keep = ~isnan(dataWeek.MidSizeNormBin) & dataWeek.MidSizeNormBin > 0;
absence = dataWeek.MidSizeNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.MidSizeNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffDiff)]);
else
    bubblelim([1 round(max(dataWeek.MidSizeNormBin))]);
end
xlim([0,53])
ylabel('Mid-size')
set(gca,'xticklabel',[])
yticks(dates)

subplot(4,1,3)
blue = colorMA;
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
absence = dataWeek.MaleNormBin == 0 & WeekData.Noeffort == 0; % Edit by AD
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.MaleNormBin(idxYear & keep)),blue)
end
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffDiff)]);
else
    bubblelim([1 round(max(dataWeek.MaleNormBin))]);
end
xlim([0,53])
ylabel('Adult Males')
yticks(dates)

% Difference
subplot(4,1,4)
black = '#C0C0C0'; %silver
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
absence = WeekData.NormBin == 0 & WeekData.Noeffort == 0; % Change by AD
for y = 1:length(years)
    hold on
    idxYear = WeekData.Year == years(y);
    idxNoeffort = WeekData.Noeffort == 1;
    scatter(WeekData.wN(idxYear & idxNoeffort),WeekData.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(WeekData.wN(idxYear & absence),WeekData.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',black)
    bubblechart(WeekData.wN(idxYear & keep),WeekData.Year(idxYear & keep),...
        round(WeekData.NormBin(idxYear & keep)),black)
end
bubblesize([2 15])
blgd= bubblelegend('Mean daily min');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
if Scale == 1
    bubblelim([1 round(CutOffDiff)]);
else
    bubblelim([1 round(max(CombinedWeek.Difference))]);
end
xlim([0,53])
ylabel('Uncategorized')
xlabel('Week of the year')
set(gca,'xticklabel',[])
yticks(dates)
end
%% save plot
set(gcf,'Position',[-1165         552         812         476])
if Scale == 1
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesDifferenceScaled.png'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesDifferenceScaled.pdf'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
else
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesDifference.png'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
    weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesDifference.pdf'];
    exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
end
%% ===================================================================
%% Stacked Weekly Timeseries (alternative to the bubble chart - all three
%% classes shown together in one panel instead of separate subplots)
%% Trimmed to the actual recording period - WeekData/dataWeek are padded
%% out to Jan 1 - Dec 31 of the first/last year for the bubble charts'
%% week-of-year axis, but that padding shouldn't show up here as empty bars.
%% ===================================================================
validWeeks = find(~dataWeek.Noeffort);
dataWeekTrim = dataWeek(validWeeks(1):validWeeks(end),:);

stackData = [dataWeekTrim.SocialGroupNormBin, dataWeekTrim.MidSizeNormBin, dataWeekTrim.MaleNormBin];
stackData(isnan(stackData)) = 0;

fStack = figure('Position',[296 417 900 380],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
hBar = bar(dataWeekTrim.tbin, stackData, 'stacked','EdgeColor','none','BarWidth',1);
hBar(1).FaceColor = colorSG;
hBar(2).FaceColor = colorMS;
hBar(3).FaceColor = colorMA;
legend({'Social Groups','Mid-size','Adult Males'},'Location','northeastoutside')
ylabel('Mean daily presence (min/week)')
xlabel('Date')
title(['Stacked Weekly Presence by Size Class - ',siteNameDisp])
xlim([dataWeekTrim.tbin(1), dataWeekTrim.tbin(end)])

stackfn = [saveDirectory,'\',siteName,'_StackedTimeSeries.png'];
exportgraphics(fStack,stackfn,'ContentType','vector','Resolution',300);
stackfn = [saveDirectory,'\',siteName,'_StackedTimeSeries.pdf'];
exportgraphics(fStack,stackfn,'ContentType','vector','Resolution',300);
%% ===================================================================
%% Same weekly presence data, but as three separate subplots (one per class)
%% instead of one stacked panel - easier to read each class's own trend/scale.
%% ===================================================================
fSubplots = figure('Position',[296 417 900 600],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
ax1 = subplot(3,1,1);
bar(dataWeekTrim.tbin, dataWeekTrim.SocialGroupNormBin,'FaceColor',colorSG,'EdgeColor','none','BarWidth',1)
ylabel('Social Groups')
title(['Weekly Presence by Size Class - ',siteNameDisp])
set(gca,'xticklabel',[])

ax2 = subplot(3,1,2);
bar(dataWeekTrim.tbin, dataWeekTrim.MidSizeNormBin,'FaceColor',colorMS,'EdgeColor','none','BarWidth',1)
ylabel('Mid-size (min/week)')
set(gca,'xticklabel',[])

ax3 = subplot(3,1,3);
bar(dataWeekTrim.tbin, dataWeekTrim.MaleNormBin,'FaceColor',colorMA,'EdgeColor','none','BarWidth',1)
ylabel('Adult Males')
xlabel('Date')

linkaxes([ax1,ax2,ax3],'x')
xlim(ax1,[dataWeekTrim.tbin(1), dataWeekTrim.tbin(end)])

% Show only every other month so the x-axis labels don't crowd each other
monthTicks = dateshift(dataWeekTrim.tbin(1),'start','month') : calmonths(2) : dataWeekTrim.tbin(end);
xticks(ax1,monthTicks); xticks(ax2,monthTicks); xticks(ax3,monthTicks);
xtickformat(ax3,'MMM yyyy');

subplotfn = [saveDirectory,'\',siteName,'_TimeSeries_Subplots.png'];
exportgraphics(fSubplots,subplotfn,'ContentType','vector','Resolution',300);
subplotfn = [saveDirectory,'\',siteName,'_TimeSeries_Subplots.pdf'];
exportgraphics(fSubplots,subplotfn,'ContentType','vector','Resolution',300);
%% ===================================================================
%% Area-Proportional Overlap Diagrams: Daily and Hourly presence overlap of the three classes
%% Circle sizes/positions are fit the same way vennX.m / PropVenn.m do it elsewhere in this
%% repo (pairwise area matching), but rendered as alpha-blended fills in the fixed class
%% colors (instead of vennX's colormap, which can't keep a fixed color per class and has
%% a region-aliasing quirk) with region counts computed directly and exactly (not approximated).
%% ===================================================================
classLabels = {'Social Groups','Mid-size','Adult Males'};
classColors = {colorSG, colorMS, colorMA};

% Daily presence (was the class present at all that day?)
dailyF = sexbinPresence.SocialGroup > 0;
dailyJ = sexbinPresence.MidSize > 0;
dailyM = sexbinPresence.Male > 0;
fVennDaily = plotClassVennProportional(dailyF,dailyJ,dailyM,classLabels,classColors, ...
    ['Daily Presence Overlap - ',siteNameDisp],'count');
vennfn = [saveDirectory,'\',siteName,'_VennDaily.png'];
exportgraphics(fVennDaily,vennfn,'ContentType','vector','Resolution',300);
vennfn = [saveDirectory,'\',siteName,'_VennDaily.pdf'];
exportgraphics(fVennDaily,vennfn,'ContentType','vector','Resolution',300);

fVennDailyPct = plotClassVennProportional(dailyF,dailyJ,dailyM,classLabels,classColors, ...
    ['Daily Presence Overlap (%) - ',siteNameDisp],'percent');
vennfn = [saveDirectory,'\',siteName,'_VennDaily_Percent.png'];
exportgraphics(fVennDailyPct,vennfn,'ContentType','vector','Resolution',300);
vennfn = [saveDirectory,'\',siteName,'_VennDaily_Percent.pdf'];
exportgraphics(fVennDailyPct,vennfn,'ContentType','vector','Resolution',300);

% Hourly presence (PreAbsF/PreAbsJ/PreAbsM were already computed in Step 3)
hourlyF = sexhourlyTab.PreAbsF > 0;
hourlyJ = sexhourlyTab.PreAbsJ > 0;
hourlyM = sexhourlyTab.PreAbsM > 0;
fVennHourly = plotClassVennProportional(hourlyF,hourlyJ,hourlyM,classLabels,classColors, ...
    ['Hourly Presence Overlap - ',siteNameDisp],'count');
vennfn = [saveDirectory,'\',siteName,'_VennHourly.png'];
exportgraphics(fVennHourly,vennfn,'ContentType','vector','Resolution',300);
vennfn = [saveDirectory,'\',siteName,'_VennHourly.pdf'];
exportgraphics(fVennHourly,vennfn,'ContentType','vector','Resolution',300);

fVennHourlyPct = plotClassVennProportional(hourlyF,hourlyJ,hourlyM,classLabels,classColors, ...
    ['Hourly Presence Overlap (%) - ',siteNameDisp],'percent');
vennfn = [saveDirectory,'\',siteName,'_VennHourly_Percent.png'];
exportgraphics(fVennHourlyPct,vennfn,'ContentType','vector','Resolution',300);
vennfn = [saveDirectory,'\',siteName,'_VennHourly_Percent.pdf'];
exportgraphics(fVennHourlyPct,vennfn,'ContentType','vector','Resolution',300);
%% ===================================================================
%% Overlapping ICI Histograms by Size Class (matches OverlappingHistograms_Step1.m convention)
%% ===================================================================
binDataT = timetable2table(binData_mod);
indSG = binDataT.SocialGroup == 1;
indMS = binDataT.MidSize == 1;
indMA = binDataT.Male == 1;

BinEdge = 300:50:2000;

fHist = figure('DefaultAxesFontSize',12,'DefaultTextFontName','Times');
histogram(binDataT.mean_ICISel(indSG),'BinEdges',BinEdge,'FaceColor',colorSG,'FaceAlpha',0.6,'EdgeColor','none')
hold on
histogram(binDataT.mean_ICISel(indMS),'BinEdges',BinEdge,'FaceColor',colorMS,'FaceAlpha',0.6,'EdgeColor','none')
histogram(binDataT.mean_ICISel(indMA),'BinEdges',BinEdge,'FaceColor',colorMA,'FaceAlpha',0.6,'EdgeColor','none')
legend({'Social Groups','Mid-size','Adult Males'})
xlabel('Interclick Interval (ms)')
ylabel('Count of 5-min bins')
title(['ICI Distribution by Size Class - ',siteNameDisp])
xlim([300 2000])
hold off

icifn = [saveDirectory,'\',siteName,'_ICIHistogram_BySizeClass.png'];
exportgraphics(fHist,icifn,'ContentType','vector','Resolution',300);
icifn = [saveDirectory,'\',siteName,'_ICIHistogram_BySizeClass.pdf'];
exportgraphics(fHist,icifn,'ContentType','vector','Resolution',300);

% Log-scale version - Adult Males tends to have far fewer 5-min bins than the
% other classes, so a linear y-axis can make its histogram nearly invisible.
fHistLog = figure('DefaultAxesFontSize',12,'DefaultTextFontName','Times');
histogram(binDataT.mean_ICISel(indSG),'BinEdges',BinEdge,'FaceColor',colorSG,'FaceAlpha',0.6,'EdgeColor','none')
hold on
histogram(binDataT.mean_ICISel(indMS),'BinEdges',BinEdge,'FaceColor',colorMS,'FaceAlpha',0.6,'EdgeColor','none')
histogram(binDataT.mean_ICISel(indMA),'BinEdges',BinEdge,'FaceColor',colorMA,'FaceAlpha',0.6,'EdgeColor','none')
legend({'Social Groups','Mid-size','Adult Males'})
xlabel('Interclick Interval (ms)')
ylabel('Count of 5-min bins (log scale)')
title(['ICI Distribution by Size Class (log scale) - ',siteNameDisp])
xlim([300 2000])
set(gca,'YScale','log')
hold off

icifn = [saveDirectory,'\',siteName,'_ICIHistogram_BySizeClass_Log.png'];
exportgraphics(fHistLog,icifn,'ContentType','vector','Resolution',300);
icifn = [saveDirectory,'\',siteName,'_ICIHistogram_BySizeClass_Log.pdf'];
exportgraphics(fHistLog,icifn,'ContentType','vector','Resolution',300);
%% save text file with max and mins for publication if needed
txtFileName = [siteName,'_','BubblePlotInfo.txt'];
paramfile = fullfile(saveDirectory,txtFileName);
fileid = fopen(paramfile, 'w');
fclose(fileid);
fileid = fopen(paramfile, 'at');
fprintf(fileid, ['Bubble Plot Time Series Information for Site ' siteName...
   '\n\nSocial Groups had a maximum mean daily presence (min) per week of \t' num2str(Femax)...
   '\n\nSocial Groups had a total mean daily presence (min) per week of \t' num2str(mean(dataWeek.SocialGroupNormBin,'omitnan')) '\n(Excludes NaNs)\t'...
   '\n\nWeeks when Social Groups were present had a total mean daily presence (min) per week of \t' num2str(mean(TotalDataSG)) '\n(Excludes NaNs and Zeros)\t'...
   '\n\nMid Size had a maximum mean daily presence (min) per week of \t' num2str(Jumax)...
   '\n\nMid Size had a total mean daily presence (min) per week of \t' num2str(mean(dataWeek.MidSizeNormBin,'omitnan')) '\n(Excludes NaNs)\t'...
   '\n\nWeeks when Mid Size were present had a total mean daily presence (min) per week of \t' num2str(mean(TotalDataMS)) '\n(Excludes NaNs and Zeros)\t'...
   '\n\nMales had a maximum mean daily presence (min) per week of \t' num2str(Mamax)...
   '\n\nMales had a total mean daily presence (min) per week of \t' num2str(mean(dataWeek.MaleNormBin,'omitnan')) '\n(Excludes NaNs)\t'...
   '\n\nWeeks when Males were present had a total mean daily presence (min) per week of \t' num2str(mean(TotalDataM)) '\n(Excludes NaNs and Zeros)\t'...
   '\n\nAll sperm whales had a maximum mean daily presence (min) per week of \t' num2str(Gemax)...
   '\n\nAll sperm whales had a total mean daily presence (min) per week of \t' num2str(mean(WeekData.NormBin,'omitnan')) '\n(Excludes NaNs)\t'...
   '\n\nWeeks when all sperm whales were present had a total mean daily presence (min) per week of \t' num2str(mean(TotalDataG)) '\n(Excludes NaNs and Zeros)\t'...
   '\n\nUncategorized had a maximum mean daily presence (min) per week of \t' num2str(Comax)...
   '\n\nUncategorized had a total mean daily presence (min) per week of \t' num2str(mean(CombinedWeek.Difference,'omitnan')) '\n(Excludes NaNs)\t'...
   '\n\nWeeks when uncategorized were present had a total mean daily presence (min) per week of \t' num2str(mean(TotalDataD)) '\n(Excludes NaNs and Zeros)\t']);
fclose(fileid);

%% ===================================================================
%% Local functions: area-proportional 3-set Euler diagram
%% Circle sizing/positioning reuses vennX.m's own pairwise area-matching
%% algorithm (its venn2 subroutine, adapted below as vennPairDist, and the
%% radius/center geometry from its plot_venn3), so the diagram is genuinely
%% area-proportional the same way PropVenn.m's vennX-based diagrams are
%% elsewhere in this repo. Rendering uses alpha-blended fills in fixed
%% per-class colors instead of vennX's colormap (which can't keep one fixed
%% color per class, and has a region-aliasing quirk in how it sums pixel
%% weights), and region counts are computed directly and exactly rather
%% than relying on vennX's own approximated region labels.
%% ===================================================================
function fig = plotClassVennProportional(A,B,C,labels,colors,titleStr,labelMode)
% labelMode: 'count' (default) shows raw region counts; 'percent' shows each
% region as a percentage of the total (the 7 regions are a mutually exclusive
% partition of "at least one class present", so the percentages sum to 100%).
if nargin < 7 || isempty(labelMode)
    labelMode = 'count';
end
A = logical(A); B = logical(B); C = logical(C);

% Region counts, in the "exclusive singles / inclusive-of-triple pairwise"
% convention vennX.m's fitting math expects (matches how PropVenn.m already
% feeds vennX.m elsewhere in this repo).
aOnly = sum(A & ~B & ~C);
abAll = sum(A & B);      % includes the triple overlap
bOnly = sum(B & ~A & ~C);
bcAll = sum(B & C);      % includes the triple overlap
cOnly = sum(C & ~A & ~B);
caAll = sum(C & A);      % includes the triple overlap
abc   = sum(A & B & C);

maxVal = max([aOnly,abAll,bOnly,bcAll,cOnly,caAll,abc,1]);
resolution = max(maxVal/300, 0.02); % ~300 steps across the largest region, per vennX.m's own guidance

distAB = vennPairDist( aOnly+caAll, abAll+abc, bOnly+bcAll, resolution );
distBC = vennPairDist( bOnly+abAll, bcAll+abc, cOnly+caAll, resolution );
distAC = vennPairDist( aOnly+abAll, caAll+abc, bcAll+cOnly, resolution );

r1 = sqrt( (aOnly+abAll+caAll+abc)/pi );
r2 = sqrt( (abAll+bOnly+bcAll+abc)/pi );
r3 = sqrt( (bcAll+cOnly+caAll+abc)/pi );

y = ( distAC^2 - distBC^2 + distAB^2 ) / 2 / distAB;
baseY = max(r1,r2);
centers = [ r1,            baseY; ...
            r1+distAB,     baseY; ...
            r1+y,          baseY+sqrt(distAC^2 - y^2) ];
radii = [r1, r2, r3];

fig = figure('DefaultAxesFontSize',12,'DefaultTextFontName','Times');
hold on
th = linspace(0,2*pi,200);
for k = 1:3
    c = colors{k};
    if ischar(c) || isstring(c) % fill() needs an RGB triplet, not a hex string like bar()/histogram() accept
        c = char(c);
        c = [hex2dec(c(2:3)) hex2dec(c(4:5)) hex2dec(c(6:7))]/255;
    end
    fill(centers(k,1)+radii(k)*cos(th), centers(k,2)+radii(k)*sin(th), c, 'FaceAlpha',0.5,'EdgeColor','none')
end

% Label each region at its true centroid, found on a coarse grid (same
% membership-mask technique vennX.m's plot_venn3 uses for its own labels).
gridRes = max(radii)/60;
xr = (min(centers(:,1)-radii)-0.1) : gridRes : (max(centers(:,1)+radii)+0.1);
yr = (min(centers(:,2)-radii)-0.1) : gridRes : (max(centers(:,2)+radii)+0.1);
[X,Y] = meshgrid(xr,yr);
inA = (X-centers(1,1)).^2 + (Y-centers(1,2)).^2 < r1^2;
inB = (X-centers(2,1)).^2 + (Y-centers(2,2)).^2 < r2^2;
inC = (X-centers(3,1)).^2 + (Y-centers(3,2)).^2 < r3^2;

regionVals = [aOnly, bOnly, cOnly, sum(A&B&~C), sum(B&C&~A), sum(A&C&~B), abc];
if strcmp(labelMode,'percent')
    total = sum(regionVals);
    if total > 0
        regionPct = round(100*regionVals/total);
        regionPct(regionVals > 0 & regionPct == 0) = 1; % don't let a nonzero region round down to 0%
        regionStrs = arrayfun(@(v) sprintf('%.0f%%',v), regionPct, 'UniformOutput', false);
    else
        regionStrs = repmat({'0%'},1,7);
    end
else
    regionStrs = arrayfun(@(v) num2str(v), regionVals, 'UniformOutput', false);
end

placeVennLabel(X,Y, inA & ~inB & ~inC, regionStrs{1}, false)
placeVennLabel(X,Y, inB & ~inA & ~inC, regionStrs{2}, false)
placeVennLabel(X,Y, inC & ~inA & ~inB, regionStrs{3}, false)
placeVennLabel(X,Y, inA & inB & ~inC, regionStrs{4}, false)
placeVennLabel(X,Y, inB & inC & ~inA, regionStrs{5}, false)
placeVennLabel(X,Y, inA & inC & ~inB, regionStrs{6}, false)
placeVennLabel(X,Y, inA & inB & inC, regionStrs{7}, true)

text(centers(1,1)-radii(1)*0.7, centers(1,2)+radii(1)*0.85, labels{1}, 'HorizontalAlignment','center','FontWeight','bold')
text(centers(2,1)+radii(2)*0.7, centers(2,2)+radii(2)*0.85, labels{2}, 'HorizontalAlignment','center','FontWeight','bold')
text(centers(3,1), centers(3,2)-radii(3)*0.85, labels{3}, 'HorizontalAlignment','center','FontWeight','bold')

axis equal off
title(titleStr)
hold off
end

function placeVennLabel(X,Y,mask,labelStr,bold)
% Places a region's (already-formatted) label at the centroid of its
% membership mask - skips regions that didn't materialize (e.g. zero-count
% overlaps at small resolution).
if ~any(mask(:))
    return
end
tx = mean(X(mask));
ty = mean(Y(mask));
if bold
    text(tx,ty,labelStr,'HorizontalAlignment','center','FontSize',12,'FontWeight','bold')
else
    text(tx,ty,labelStr,'HorizontalAlignment','center','FontSize',12)
end
end

function dist = vennPairDist(a,b,c,resolution)
% Adapted directly from vennX.m's venn2(): finds the center-to-center distance
% for two circles (areas a+b and b+c) whose overlap area best matches b, by
% shrinking the distance from fully-separated down to touching.
r1 = sqrt( (a+b)/pi );
r2 = sqrt( (b+c)/pi );

sizeY = max(2*r1, 2*r2);
[X,Y] = meshgrid(0:resolution:(2*r1+2*r2), 0:resolution:sizeY);
center1X = r1;
center1Y = sizeY/2;
center2Y = sizeY/2;

dist = r2; % fallback: fully overlapping (only reached if b is larger than either circle's own area)
for newCenter = (2*r1+r2):-resolution:r1
    img = (X-center1X).^2 + (Y-center1Y).^2 < r1^2 & ...
          (X-newCenter).^2 + (Y-center2Y).^2 < r2^2;
    if sum(img(:)) * resolution^2 > b
        dist = newCenter - center1X;
        break
    end
    dist = newCenter - center1X;
end
end
