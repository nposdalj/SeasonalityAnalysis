%% clear workspace
clear all;close all;clc;

%% load data
siteName = 'AB';
DataDir = 'I:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis';
load([DataDir,'\',siteName,'\',siteName,'_workspaceStep3.mat']);
saveDir = 'I:\My Drive\GofAK_TPWS_metadataReduced\Plots';

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

% Social Unit
fBubble = figure('Position',[296 417 766 378.5000],'DefaultAxesFontSize',12,'DefaultTextFontName','Times');

years = unique(dataWeek.Year);
subplot(3,1,1)
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
% bubblelim([9 1296]);
xlim([0,53])
ylabel('Social Unit')


subplot(3,1,2)
blue = '#fc8d62';%'#349987';
keep = ~isnan(dataWeek.Juvenile_adj) & dataWeek.Juvenile_adj > 0;
absence = dataWeek.Juvenile_adj == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear  & keep),dataWeek.Year(idxYear  & keep),...
        round(dataWeek.Juvenile_adj(idxYear  & keep)),blue)
end
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = 3;
bubblesize([2 15])
set(gca,'ydir','reverse')
% bubblelim([1 494]);
xlim([0,53])
ylabel('Mid-size')


%fBubble3 = figure('Position',[411 1008 816 160]);
subplot(3,1,3)
blue = '#8da0cb';
keep = ~isnan(dataWeek.Male_adj) & dataWeek.Male_adj > 0;
absence = dataWeek.Male_adj == 0;
for y = 1:length(years)
    hold on
    idxYear = dataWeek.Year == years(y);
    idxNoeffort = dataWeek.Noeffort == 1;
    scatter(dataWeek.wN(idxYear & idxNoeffort),dataWeek.Year(idxYear & idxNoeffort),7,...
        'x','MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8])
    scatter(dataWeek.wN(idxYear & absence),dataWeek.Year(idxYear & absence),3,...
        'o','MarkerEdgeColor',blue)
    bubblechart(dataWeek.wN(idxYear & keep),dataWeek.Year(idxYear & keep),...
        round(dataWeek.Male_adj(idxYear & keep)),blue)
end
blgd= bubblelegend('Avg. daily minutes');
blgd.Location = 'northeastoutside';
blgd.NumBubbles = 3;
bubblesize([2 15])
set(gca,'ydir','reverse')
% bubblelim([1 66]);
xlim([0,53])
xlabel('Week of the year')
ylabel('Adult Male')
