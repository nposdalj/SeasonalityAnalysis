%% clear workspace
clear all;close all;clc;
% Needs to be run in 2018B or later
%% load data
siteName = 'KOA';
GDrive = 'I';
region = 'GofAK'; %region
NumBub = 3;
DataDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteName];
saveDirectory = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\',siteName];
Scale = 0; %(1) all scales are based on the lowest 'max' value (0) each size class has a respective scale
%% Load Workspace Step 2
GDrive_correctII = GDrive; % Store correct GDrive to overwrite some path names later
load([DataDir,'\',siteName,'_workspaceStep2.mat']);
clear mean

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

GDrive = GDrive_correctII; %Correct GDrive
effortXls(1) = GDrive;
saveDir(1) = GDrive;
tpwsPath(1) = GDrive;
dayBinCSV(1) = GDrive;
filename(1) = GDrive;

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
blue = '#fc8d62';
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
blue = '#8da0cb';
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
blue = '#66c2a5';%'#2e59a8';
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
blue = '#fc8d62';%'#349987';
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
blue = '#8da0cb';
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
blue = '#fc8d62';
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
blue = '#8da0cb';
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
blue = '#66c2a5';
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
blue = '#fc8d62';
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
blue = '#8da0cb';
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