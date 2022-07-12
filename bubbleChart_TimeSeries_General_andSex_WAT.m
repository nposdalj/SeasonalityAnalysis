%% clear workspace
clear all;close all;clc;

%% load data
siteName = 'WC';
NumBub = 3;
DataDir = 'H:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis';
saveDirectory = ['H:\My Drive\WAT_TPWS_metadataReduced\Plots\',siteName];
%% Retime data weekly
load([DataDir,'\',siteName,'\',siteName,'_workspaceStep2.mat']);
clear mean

%Add first week of the year and last to account for x's
%Find the first and last year
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
   WeekData(leap(l)-1,6:7) = varfun(@nanmean,WeekData(leap(l)-1:leap(l),6:7));
   WeekData(leap(l),:) = [];
end
% correct year for leap year
WeekData.diffYear = zeros(height(WeekData),1);
WeekData.diffYear(2:end) = diff(WeekData.Year);
adjYear = find(WeekData.diffYear == 1 & WeekData.wN == 52);
WeekData.Year(adjYear) = WeekData.Year(adjYear-1); 
% convert bins to minutes
WeekData.NormBin = WeekData.NormBin *5;

%Sex
load([DataDir,'\',siteName,'\',siteName,'_workspaceStep3.mat']);
clear mean
sexbinPresence.day = grp2idx(sexbinPresence.day);
sexbinPresenceAll = synchronize(allDays,sexbinPresence);
sexbinPresenceAll = removevars(sexbinPresenceAll, {'ZeroCol'});
dataWeek = retime(sexbinPresenceAll,'weekly',@(x) mean(x, 'omitnan'));
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
dataWeek.FemaleNormBin = dataWeek.FemaleNormBin *5;
dataWeek.JuvenileNormBin = dataWeek.JuvenileNormBin *5;
dataWeek.MaleNormBin = dataWeek.MaleNormBin *5;

%% Delete first week for specific sites
if (strcmp(siteName,'NC') | strcmp(siteName,'BC') | strcmp(siteName,'GS') | strcmp(siteName,'BP') | strcmp(siteName,'BS')...
        | strcmp(siteName,'WC') | strcmp(siteName,'OC') | strcmp(siteName,'HZ') | strcmp(siteName,'JAX'))
    WeekData(1,:) = [];
    dataWeek(1,:) = [];
end
dataWeek.year = year(dataWeek.tbin);
WeekData.year = year(WeekData.tbin); %All this is done so that you can generate the plots without an extra year on top :)
%% Checking to see how much was missed
if strcmp(siteName,'BD') % does this need to be changed depending on the site under investigation?
    CombinedWeek = dataWeek(:,20:22);
else
CombinedWeek = dataWeek(:,19:21);
end
CombinedWeek.NormBin = WeekData.NormBin;
CombinedWeek.Added = CombinedWeek.FemaleNormBin + CombinedWeek.JuvenileNormBin + CombinedWeek.MaleNormBin;
CombinedWeek.Difference = CombinedWeek.NormBin - CombinedWeek.Added;
CombinedWeek.Difference( CombinedWeek.Difference <= 0 ) = 0;
%% Find Max Bubble Size and Distribution
%Find max
Femax = max(dataWeek.FemaleNormBin);
Jumax = max(dataWeek.JuvenileNormBin);
Mamax = max(dataWeek.MaleNormBin);
Comax = max(CombinedWeek.Difference);
MAX = max([Femax, Jumax, Mamax, Comax]);
MIN = min([Femax, Jumax, Mamax, Comax]);

TotalData = [dataWeek.FemaleNormBin; dataWeek.JuvenileNormBin; dataWeek.MaleNormBin; CombinedWeek.Difference]; %combine values
%delete NAs and zeros
nanRows = isnan(TotalData);
zeroRows = TotalData==0;
badRows = nanRows | zeroRows;
TotalData(badRows) = [];

%Find the histogram bins and set a cutoff
figure, hist(TotalData,20)
[N,edges] = histcounts(TotalData,20);
lessThan = find(N<3); %find bins with less than 5
CutOff = edges(lessThan(1)); %values over this value will be considered 'greater than' in the bubble sizes

%Actual Cutoff (max bubble size)
%CutOff = MIN;

%ADDED BY AD: Manually mark periods of no effort for WAT sites
if siteName == "NC"
    WeekData.Noeffort(WeekData.Count_Click == 0) = 1;
    dataWeek.Noeffort(WeekData.Count_Click == 0) = 1;
elseif siteName == "BC"
    WeekData.Noeffort(WeekData.Count_Click == 0) = 1;
    dataWeek.Noeffort(WeekData.Count_Click == 0) = 1;
end
%% Plot data
%No Social Group Data
if strcmp(siteName, 'NOFEM') % For AD, don't run
    % Want to analyze social group presence, so skip to second version of
    % code (after "else")
    
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');

subplot(3,1,1)
years = unique(WeekData.Year);
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
black = [0,0,0];
absence = WeekData.NormBin == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
bubblelim([1 round(max(WeekData.NormBin))]);
xlim([0,53])
ylabel('General')
set(gca,'xticklabel',[])

subplot(3,1,2)
years = unique(dataWeek.Year);
blue = '#fc8d62';%'#349987';
keep = ~isnan(dataWeek.JuvenileNormBin) & dataWeek.JuvenileNormBin > 0;
absence = dataWeek.JuvenileNormBin == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.JuvenileNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.JuvenileNormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Mid-size')
set(gca,'xticklabel',[])

%fBubble3 = figure('Position',[411 1008 816 160]);
subplot(3,1,3)
blue = '#8da0cb';
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
absence = dataWeek.MaleNormBin == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.MaleNormBin))]);
xlim([0,53])
xlabel('Week of the year')
% ylim([2016.75,2017.25])
ylabel('Adult Males')

else
    
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
years = unique(dataWeek.Year);

% General Presence
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
subplot(4,1,1)
black = [0,0,0];
%absence = WeekData.NormBin == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
bubblelim([1 round(max(WeekData.NormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('General')
set(gca,'xticklabel',[])
yticks([2016 2017 2018 2019])

% Social Group
subplot(4,1,2)
blue = '#66c2a5';%'#2e59a8';
keep = ~isnan(dataWeek.FemaleNormBin) & dataWeek.FemaleNormBin > 0;
%absence = dataWeek.FemaleNormBin == 0;
absence = dataWeek.FemaleNormBin == 0 & WeekData.Noeffort == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.FemaleNormBin(idxYear & keep)),blue)
end
bubblesize([2 15])
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.FemaleNormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Social Groups')
set(gca,'xticklabel',[])
yticks([2016 2017 2018 2019])

subplot(4,1,3)
blue = '#fc8d62';%'#349987';
keep = ~isnan(dataWeek.JuvenileNormBin) & dataWeek.JuvenileNormBin > 0;
%absence = dataWeek.JuvenileNormBin == 0;
absence = dataWeek.JuvenileNormBin == 0 & WeekData.Noeffort == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.JuvenileNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.JuvenileNormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Mid-size')
set(gca,'xticklabel',[])
yticks([2016 2017 2018 2019])

%fBubble3 = figure('Position',[411 1008 816 160]);
subplot(4,1,4)
blue = '#8da0cb';
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
%absence = dataWeek.MaleNormBin == 0;
absence = dataWeek.MaleNormBin == 0 & WeekData.Noeffort == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.MaleNormBin))]);
xlim([0,53])
xlabel('Week of the year')
% ylim([2016.75,2017.25])
ylabel('Adult Males')
yticks([2016 2017 2018 2019])
end % don't run this line

%% save plot
set(gcf,'Position',[-1165         552         812         476])
%weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeries.pdf'];
%weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeries.png'];
weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeries_Style2.png'];
exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);
%% Plotting the difference instead of the 'general pattern'
%No Social Group Data
if strcmp(siteName, 'NOFEM')
    
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');

subplot(3,1,1)
years = unique(dataWeek.year);
blue = '#fc8d62';%'#349987';
keep = ~isnan(dataWeek.JuvenileNormBin) & dataWeek.JuvenileNormBin > 0;
absence = dataWeek.JuvenileNormBin == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.JuvenileNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(CombinedWeek.Difference))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Mid-size')
set(gca,'xticklabel',[])

%fBubble3 = figure('Position',[411 1008 816 160]);
subplot(3,1,2)
blue = '#8da0cb';
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
absence = dataWeek.MaleNormBin == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.MaleNormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Adult Males')

subplot(3,1,3)
years = unique(WeekData.year);
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
black = '#C0C0C0'; %silver
absence = WeekData.NormBin == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
bubblelim([1 round(max(CombinedWeek.Difference))]);
xlim([0,53])
ylabel('Difference')
xlabel('Week of the year')
set(gca,'xticklabel',[])

else
    
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
years = unique(dataWeek.Year);

% Social Group
subplot(4,1,1)
blue = '#66c2a5';%'#2e59a8';
keep = ~isnan(dataWeek.FemaleNormBin) & dataWeek.FemaleNormBin > 0;
absence = dataWeek.FemaleNormBin == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.FemaleNormBin(idxYear & keep)),blue)
end
bubblesize([2 15])
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
%bubblelim([1 round(max(dataWeek.FemaleNormBin))]);
bubblelim([1 round(CutOff)]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Social Groups')
set(gca,'xticklabel',[])

subplot(4,1,2)
blue = '#fc8d62';%'#349987';
keep = ~isnan(dataWeek.JuvenileNormBin) & dataWeek.JuvenileNormBin > 0;
absence = dataWeek.JuvenileNormBin == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.JuvenileNormBin(idxYear  & keep)),blue)
end
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
% bubblelim([1 round(max(dataWeek.JuvenileNormBin))]);
bubblelim([1 round(CutOff)]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Mid-size')
set(gca,'xticklabel',[])

%fBubble3 = figure('Position',[411 1008 816 160]);
subplot(4,1,3)
blue = '#8da0cb';
keep = ~isnan(dataWeek.MaleNormBin) & dataWeek.MaleNormBin > 0;
absence = dataWeek.MaleNormBin == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
bubblesize([2 15])
set(gca,'ydir','reverse')
%bubblelim([1 round(max(dataWeek.MaleNormBin))]);
bubblelim([1 round(CutOff)]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Adult Males')

% Difference %Q from AD - what's this
subplot(4,1,4)
black = '#C0C0C0'; %silver
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
absence = WeekData.NormBin == 0;
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
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = NumBub;
set(gca,'ydir','reverse')
%bubblelim([1 round(max(CombinedWeek.Difference))]);
bubblelim([1 round(CutOff)]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Difference')
xlabel('Week of the year')
set(gca,'xticklabel',[])
end

%% save plot
set(gcf,'Position',[-1165         552         812         476])
%weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesDifference.pdf'];
weeklyfn = [saveDirectory,'\',siteName,'_BubbleTimeSeriesDifference.png'];
exportgraphics(gcf,weeklyfn,'ContentType','vector','Resolution',300);

