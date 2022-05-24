clearvars
close all

%% Parameters defined by user
filePrefix = 'WAT_BP'; % File name to match. 
siteabrev = 'BP'; %abbreviation of site.
SiteName_forPlots = 'Western Atlantic-BP';
saveDir = ['H:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis\BP']; %specify directory to save and load files
%% load workspace
load([saveDir,'\',siteabrev,'_workspace125.mat']);
saveDir = ['H:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis\BP']; %specify directory to save and load files
%% group data by 5min bins, days, weeks, and seasons
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

Click = retime(binData(:,1),'daily','sum'); % #click per day
Bin = retime(binData(:,1),'daily','count'); % #bin per day

%group data by day
dayData = synchronize(Click,Bin);
dayEffort = retime(binEffort,'daily','sum');
dayTab = synchronize(dayData,dayEffort);
dayTable = synchronize(dayData,dayEffort);
dayTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
dayTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
dayTableZeros = dayTable;
dayTable(~dayTable.Effort_Bin,:)=[]; %removes days with no effort, NOT days with no presence
%% statistical methods from Diogou, et al. 2019 - by daily 5 min bins ** USE THIS **
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
    dayTabe.Effort_Bin(222:516) = ge * 144; %for ALEUT03BD ONLY 5 on 5 off (10 minute cycle) -- meaning you're recording 0.5 percent of the time
    dayTable.Effort_Sec(222:516) = dayTable.Effort_Bin(222:516) * 5 * 60; %convert from bins into efforts in seconds per day
    else
dayTable.MaxEffort_Bin = ones(p,1)*(288);
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

%group data by week
weekTable = retime(dayTable,'weekly','mean');
%add month and year to week table
weekTable.year = year(weekTable.tbin);
weekTable.month = month(weekTable.tbin);

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
dayTable.year = year(dayTable.tbin); 
dayTable.day = day(dayTable.tbin,'dayofyear');

NANidx = ismissing(dayTable(:,{'NormBin'}));
dayTable{:,{'NormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
%%
weekTable.Percent = weekTable.Effort_Sec./86400;
tstart = weekTable.tbin(1)-7;
tend = weekTable.tbin(end)+7;
figure
yyaxis left
bar(weekTable.tbin,weekTable.HoursProp,'k')
xlim([tstart tend]);
title({'Average Weekly Presence of Sperm Whales',SiteName_forPlots})
ylabel('Proportion of Hours per Week with Clicks')
hold on
yyaxis right
plot(weekTable.tbin,weekTable.Percent,'.r')
ylabel('Percent Effort')
col = [0 0 0];
set(gcf,'defaultAxesColorOrder',[col;col])
weeklyfn = [siteabrev,'_','_weeklyPresence'];
saveas(gcf,fullfile(saveDir,weeklyfn),'png')
%% Average weekly click count (log scale)
if strcmp(siteabrev,'BD')
subplot(3,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2010),dayTable.HoursProp(dayTable.year==2010),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2010,01,01);
tend = datetime(2010,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(3,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2011),dayTable.HoursProp(dayTable.year==2011),'k')
tstart = datetime(2011,01,01);
tend = datetime(2011,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Average Daily # of 5-min Bins')
subplot(3,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2012),dayTable.HoursProp(dayTable.year==2012),'k')
tstart = datetime(2012,01,01);
tend = datetime(2012,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end

if strcmp(siteabrev,'CB')
figure
subplot(8,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2011),dayTable.HoursProp(dayTable.year==2011),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2011,01,01);
tend = datetime(2011,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(8,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2012),dayTable.HoursProp(dayTable.year==2012),'k')
tstart = datetime(2012,01,01);
tend = datetime(2012,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2013),dayTable.HoursProp(dayTable.year==2013),'k')
tstart = datetime(2013,01,01);
tend = datetime(2013,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,4) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2014),dayTable.HoursProp(dayTable.year==2014),'k')
tstart = datetime(2014,01,01);
tend = datetime(2014,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,5) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2015),dayTable.HoursProp(dayTable.year==2015),'k')
tstart = datetime(2015,01,01);
tend = datetime(2015,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(8,1,6) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2017),dayTable.HoursProp(dayTable.year==2017),'k')
tstart = datetime(2017,01,01);
tend = datetime(2017,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,7) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2018),dayTable.HoursProp(dayTable.year==2018),'k')
tstart = datetime(2018,01,01);
tend = datetime(2018,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,8) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2019),dayTable.HoursProp(dayTable.year==2019),'k')
tstart = datetime(2019,01,01);
tend = datetime(2019,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end

if strcmp(siteabrev,'QN')
figure
subplot(4,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2013),dayTable.HoursProp(dayTable.year==2013),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2013,01,01);
tend = datetime(2013,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(4,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2014),dayTable.HoursProp(dayTable.year==2014),'k')
tstart = datetime(2014,01,01);
tend = datetime(2014,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(4,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2015),dayTable.HoursProp(dayTable.year==2015),'k')
tstart = datetime(2015,01,01);
tend = datetime(2015,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(4,1,4) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2017),dayTable.HoursProp(dayTable.year==2017),'k')
tstart = datetime(2017,01,01);
tend = datetime(2017,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end


if strcmp(siteabrev,'GI')
figure
subplot(3,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2018),dayTable.HoursProp(dayTable.year==2018),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2018,01,01);
tend = datetime(2018,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(3,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2019),dayTable.HoursProp(dayTable.year==2019),'k')
tstart = datetime(2019,01,01);
tend = datetime(2019,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(3,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2020),dayTable.HoursProp(dayTable.year==2020),'k')
tstart = datetime(2020,01,01);
tend = datetime(2020,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end

if strcmp(siteabrev,'CORC')
figure
subplot(2,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2014),dayTable.HoursProp(dayTable.year==2014),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2014,01,01);
tend = datetime(2014,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(2,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2015),dayTable.HoursProp(dayTable.year==2015),'k')
tstart = datetime(2015,01,01);
tend = datetime(2015,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end


if strcmp(siteabrev,'TIN')
figure
subplot(9,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2011),dayTable.HoursProp(dayTable.year==2011),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2011,01,01);
tend = datetime(2011,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(9,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2012),dayTable.HoursProp(dayTable.year==2012),'k')
tstart = datetime(2012,01,01);
tend = datetime(2012,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(9,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2013),dayTable.HoursProp(dayTable.year==2013),'k')
tstart = datetime(2013,01,01);
tend = datetime(2013,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(9,1,4) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2014),dayTable.HoursProp(dayTable.year==2014),'k')
tstart = datetime(2014,01,01);
tend = datetime(2014,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(9,1,5) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2015),dayTable.HoursProp(dayTable.year==2015),'k')
tstart = datetime(2015,01,01);
tend = datetime(2015,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(9,1,6) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2016),dayTable.HoursProp(dayTable.year==2016),'k')
tstart = datetime(2016,01,01);
tend = datetime(2016,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(9,1,7) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2017),dayTable.HoursProp(dayTable.year==2017),'k')
tstart = datetime(2017,01,01);
tend = datetime(2017,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(9,1,8) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2018),dayTable.HoursProp(dayTable.year==2018),'k')
tstart = datetime(2018,01,01);
tend = datetime(2018,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(9,1,9) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2019),dayTable.HoursProp(dayTable.year==2019),'k')
tstart = datetime(2019,01,01);
tend = datetime(2019,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end


if strcmp(siteabrev,'SAP')
figure
subplot(10,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2010),dayTable.HoursProp(dayTable.year==2010),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2010,01,01);
tend = datetime(2010,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(10,1,2) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2011),dayTable.HoursProp(dayTable.year==2011),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2011,01,01);
tend = datetime(2011,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(10,1,3) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2012),dayTable.HoursProp(dayTable.year==2012),'k')
tstart = datetime(2012,01,01);
tend = datetime(2012,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(9,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2013),dayTable.HoursProp(dayTable.year==2013),'k')
tstart = datetime(2013,01,01);
tend = datetime(2013,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(10,1,5) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2014),dayTable.HoursProp(dayTable.year==2014),'k')
tstart = datetime(2014,01,01);
tend = datetime(2014,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(10,1,6) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2015),dayTable.HoursProp(dayTable.year==2015),'k')
tstart = datetime(2015,01,01);
tend = datetime(2015,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(10,1,7) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2016),dayTable.HoursProp(dayTable.year==2016),'k')
tstart = datetime(2016,01,01);
tend = datetime(2016,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(10,1,8) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2017),dayTable.HoursProp(dayTable.year==2017),'k')
tstart = datetime(2017,01,01);
tend = datetime(2017,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(10,1,9) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2018),dayTable.HoursProp(dayTable.year==2018),'k')
tstart = datetime(2018,01,01);
tend = datetime(2018,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(10,1,10) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2019),dayTable.HoursProp(dayTable.year==2019),'k')
tstart = datetime(2019,01,01);
tend = datetime(2019,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end

if strcmp(siteabrev,'QC')
figure
subplot(8,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2007),dayTable.HoursProp(dayTable.year==2007),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2007,01,01);
tend = datetime(2007,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(8,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2008),dayTable.HoursProp(dayTable.year==2008),'k')
tstart = datetime(2008,01,01);
tend = datetime(2008,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2009),dayTable.HoursProp(dayTable.year==2009),'k')
tstart = datetime(2009,01,01);
tend = datetime(2009,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,4) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2010),dayTable.HoursProp(dayTable.year==2010),'k')
tstart = datetime(2010,01,01);
tend = datetime(2010,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,5) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2011),dayTable.HoursProp(dayTable.year==2011),'k')
tstart = datetime(2011,01,01);
tend = datetime(2011,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(8,1,6) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2012),dayTable.HoursProp(dayTable.year==2012),'k')
tstart = datetime(2012,01,01);
tend = datetime(2012,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,7) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2013),dayTable.HoursProp(dayTable.year==2013),'k')
tstart = datetime(2013,01,01);
tend = datetime(2013,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,8) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2014),dayTable.HoursProp(dayTable.year==2014),'k')
tstart = datetime(2014,01,01);
tend = datetime(2014,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end

if strcmp(siteabrev,'Wake')
figure
subplot(8,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2010),dayTable.HoursProp(dayTable.year==2010),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2010,01,01);
tend = datetime(2010,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(8,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2011),dayTable.HoursProp(dayTable.year==2011),'k')
tstart = datetime(2011,01,01);
tend = datetime(2011,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2012),dayTable.HoursProp(dayTable.year==2012),'k')
tstart = datetime(2012,01,01);
tend = datetime(2012,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,4) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2013),dayTable.HoursProp(dayTable.year==2013),'k')
tstart = datetime(2013,01,01);
tend = datetime(2013,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,5) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2014),dayTable.HoursProp(dayTable.year==2014),'k')
tstart = datetime(2014,01,01);
tend = datetime(2014,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(8,1,6) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2015),dayTable.HoursProp(dayTable.year==2015),'k')
tstart = datetime(2015,01,01);
tend = datetime(2015,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,7) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2016),dayTable.HoursProp(dayTable.year==2016),'k')
tstart = datetime(2016,01,01);
tend = datetime(2016,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(8,1,8) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2017),dayTable.HoursProp(dayTable.year==2017),'k')
tstart = datetime(2017,01,01);
tend = datetime(2017,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end
if strcmp(siteabrev,'BS')
figure
subplot(4,1,1) %number of plots, column, row)
bar1 = bar(dayTable.tbin(dayTable.year==2016),dayTable.HoursProp(dayTable.year==2016),'k')
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots})
tstart = datetime(2016,01,01);
tend = datetime(2016,12,31);
xlim([tstart tend]);
ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);
ylim([ylim1 ylim2]);
subplot(4,1,2) %number of plots, column, row)
bar2 = bar(dayTable.tbin(dayTable.year==2017),dayTable.HoursProp(dayTable.year==2017),'k')
tstart = datetime(2017,01,01);
tend = datetime(2017,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
subplot(4,1,3) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2018),dayTable.HoursProp(dayTable.year==2018),'k')
tstart = datetime(2018,01,01);
tend = datetime(2018,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);
ylabel('Proportion of Hours/Day with Sperm Whales')
subplot(4,1,4) %number of plots, column, row)
bar3 = bar(dayTable.tbin(dayTable.year==2019),dayTable.HoursProp(dayTable.year==2019),'k')
tstart = datetime(2019,01,01);
tend = datetime(2019,12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);

xlabel('Time (months)')
% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')
end


