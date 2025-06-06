clearvars
close all

%% Parameters defined by user
GDrive = 'P';
Region = 'WAT';
filePrefix = 'WAT_BP'; % File name to match. 
siteabrev = 'BP'; %abbreviation of site.
SiteName_forPlots = 'Western Atlantic - BP';

inDir = [GDrive ':\My Drive\' Region '_TPWS_metadataReduced\SeasonalityAnalysis\' siteabrev]; % Directory to load files
saveDir = [GDrive ':\My Drive\' Region '_TPWS_metadataReduced\Plots\' siteabrev]; % Directory to save files

% Years to plot depending on site:
Years.CB   = [                     2011 2012 2013 2014 2015      2017 2018 2019      ];
Years.QN   = [                               2013 2014 2015      2017                ];
Years.BD   = [                2010 2011 2012                                         ];
Years.GI   = [                                                        2018 2019 2020 ];
Years.TIN  = [                     2011 2012 2013 2014 2015 2016 2017 2018 2019      ];
Years.SAP  = [                2010 2011 2012 2013 2014 2015 2016 2017 2018 2019      ];
Years.Wake = [                2010 2011 2012 2013 2014 2015 2016 2017                ];
Years.CORC = [                                    2014 2015                          ];
Years.QC   = [ 2007 2008 2009 2010 2011 2012 2013 2014                               ];
Years.HZ   = [                                         2015 2016 2017 2018 2019      ];
Years.NC   = [                                         2015 2016 2017 2018 2019      ];
Years.BC   = [                                              2016 2017 2018 2019      ];
Years.GS   = [                                              2016 2017 2018 2019      ];
Years.BP   = [                                              2016 2017 2018 2019      ];
Years.BS   = [                                              2016 2017 2018 2019      ];
Years.JAX  = [                                              2016 2017 2018 2019      ];
Years.OC   = [                                         2015 2016 2017 2018 2019      ];
Years.WC   = [                                              2016 2017 2018 2019      ];

%% load workspace
GDrive_correctII = GDrive; % store correct GDrive
saveDir_correctII = saveDir; % store correct saveDir

load([inDir,'\',siteabrev,'_workspace125.mat']);

GDrive = GDrive_correctII; % restore correct GDrive
saveDir = saveDir_correctII; % restore correct saveDir
effortXls(1) = GDrive;
tpwsPath(1) = GDrive;

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
if strcmp(siteabrev,'CB')
    ge = dayTable.Effort_Bin(222:516); %bin effort (excluding ships but not considering duty cycle)
    ge = ge/288; %proportion of data that was not ships if it were full recording effort
    dayTable.Effort_Bin(222:516) = ge * 240; %for CB02 10 on 2 off (12 minute cycle) -- meaning you're recording 0.8333 percent of the time
    dayTable.Effort_Sec(222:516) = dayTable.Effort_Bin(222:516) * 5 * 60; %convert from bins into efforts in seconds per day
    else
if strcmp(siteabrev,'BD')
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
plot(weekTable.tbin,weekTable.Percent,'.','Color',[.5 .5 .5])
ylabel('Percent Effort')
col = [0 0 0];
set(gcf,'defaultAxesColorOrder',[col;col])
weeklyfn = [siteabrev,'_','_weeklyPresence'];
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [.5 .5 .5];
saveas(gcf,fullfile(saveDir,weeklyfn),'png')

%% Average weekly click count (log scale)

numYears = length(Years.(siteabrev));
siteYears = Years.(siteabrev);

ylim1 = min(dayTable.HoursProp);
ylim2 = max(dayTable.HoursProp);

for i = 1:numYears
    
subplot(numYears,1,i) %number of plots, column, row)
bar_i = bar(dayTable.tbin(dayTable.year==siteYears(i)),dayTable.HoursProp(dayTable.year==siteYears(i)),'k');
tstart = datetime(siteYears(i),01,01);
tend = datetime(siteYears(i),12,31);
xlim([tstart tend]);
ylim([ylim1 ylim2]);

if i == 1  % Title, y-label, and x-label
title({'Average Daily Presence of Sperm Whales at ',SiteName_forPlots}) % Add title to top subplot
elseif i == round(numYears/2)
    ylabel('Average Daily # of 5-min Bins') % Add y-label to the middle subplot
elseif i == numYears
    xlabel('Time (months)') % Add x-label to bottom subplot
end

end

% Save plot
dailyfn = [siteabrev,'_','_subplots_daily_average_presence'];
saveas(gcf,fullfile(saveDir,dailyfn),'png')