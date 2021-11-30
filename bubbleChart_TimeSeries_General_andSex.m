%% clear workspace
clear all;close all;clc;

%% load data
siteName = 'KS';
DataDir = 'I:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis';
saveDirectory = 'I:\My Drive\GofAK_TPWS_metadataReduced\Plots';
%% Retime data weekly
load([DataDir,'\',siteName,'\',siteName,'_workspaceStep2.mat']);
clear mean
dayTable = removevars(dayTable, {'Season','month','Year','day'});
%General
WeekData = retime(dayTable, 'weekly',@(x) mean(x, 'omitnan'));
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
binPresence.day = grp2idx(binPresence.day);
dataWeek = retime(binPresence,'weekly',@(x) mean(x, 'omitnan'));
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
%% Plot data
%No Social Group Data
if strcmp(siteName, 'KOA') | strcmp(siteName, 'KS')
    
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
blgd.NumBubbles = 3;
set(gca,'ydir','reverse')
bubblelim([1 round(max(WeekData.NormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('General')

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
blgd.NumBubbles = 3;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.JuvenileNormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Mid-size')


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
blgd.NumBubbles = 3;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.MaleNormBin))]);
xlim([0,53])
xlabel('Week of the year')
% ylim([2016.75,2017.25])
ylabel('Adult Male')

else
    
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');
years = unique(dataWeek.Year);

% General Presence
keep = ~isnan(WeekData.NormBin) & WeekData.NormBin > 0;
subplot(4,1,1)
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
blgd.NumBubbles = 3;
set(gca,'ydir','reverse')
bubblelim([1 round(max(WeekData.NormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('General')

% Social Group
subplot(4,1,2)
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
blgd.NumBubbles = 3;
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.FemaleNormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Social Group')


subplot(4,1,3)
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
blgd.NumBubbles = 3;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.JuvenileNormBin))]);
xlim([0,53])
% ylim([2016.75,2017.25])
ylabel('Mid-size')


%fBubble3 = figure('Position',[411 1008 816 160]);
subplot(4,1,4)
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
blgd.NumBubbles = 3;
bubblesize([2 15])
set(gca,'ydir','reverse')
bubblelim([1 round(max(dataWeek.MaleNormBin))]);
xlim([0,53])
xlabel('Week of the year')
% ylim([2016.75,2017.25])
ylabel('Adult Male')
end

%% save plot
weeklyfn = [saveDirectory,'\',siteName,'\',siteName,'_BubbleTimeSeries.png'];
saveas(gcf,fullfile(weeklyfn),'png')
