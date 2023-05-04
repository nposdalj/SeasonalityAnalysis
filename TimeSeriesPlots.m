clearvars
close all
% This script creates time series plots for each site.
%% Parameters defined by user
filePrefix = 'HAT'; % File name to match. 
siteabrev = 'HAT'; %abbreviation of site.
GDrive = 'G'; %directory for Google Drive
region = 'WAT';
sp = 'Pm'; % your species code
titleNAME = 'Cape Hatteras';
dataDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory where workspaces are saved
%% load workspace
GDrive_corrected = GDrive; % Preserve correct GDrive as it was entered above
load([dataDir,'\',siteabrev,'_workspaceStep2.mat']);
load([dataDir,'\',siteabrev,'_workspaceStep3.mat']);

% Overwrite some path names
GDrive = GDrive_corrected; %Correct GDrive if overwritten by loading workspace
dayBinCSV(1) = GDrive;
effortXls(1) = GDrive;
filename(1) = GDrive;
tpwsPath(1) = GDrive;

saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\',siteabrev];
%% Fill in missing days
%day table
dayTable.day= double(dayTable.day);
dayTable = retime(dayTable,'daily','fillwithconstant');

%sex table
sexbinPresence.day = double(sexbinPresence.day);
sexbinPresence = retime(sexbinPresence,'daily','fillwithconstant');
%% Retime for weekly presence
%day table
weekTable = retime(dayTable,'weekly','sum');
weekTable.NormEffort_Bin = weekTable.Effort_Sec ./weekTable.MaxEffort_Sec;
weekTable.NormEffort_Bin(isnan(weekTable.NormEffort_Bin)) = 0;
weekTable.HoursProp = weekTable.Hours ./ (weekTable.Effort_Sec ./ (60*60));

%sex table
weekPresence = retime(sexbinPresence,'weekly','sum');
weekPresence.NormEffort_Bin = weekPresence.Effort_Sec ./weekPresence.MaxEffort_Sec;
weekPresence.NormEffort_Bin(isnan(weekPresence.NormEffort_Bin)) = 0;
weekPresence.FeHoursProp = weekPresence.FeHours ./(weekPresence.Effort_Sec ./ (60*60));
weekPresence.JuHoursProp = weekPresence.JuHours ./(weekPresence.Effort_Sec ./ (60*60));
weekPresence.MaHoursProp = weekPresence.MaHours ./(weekPresence.Effort_Sec ./ (60*60));
%% Plots
%Color scheme
gray = [.5 .5 .5];       % for effort
mint = [.4000 .7608 .6471];       % for social units
persimmon = [.9882 .5529 .3843];  % for mid-size animals
slate = [.5529 .6275 .7961];      % for males

%Plot proportion of hours per DAY with sperm whale presence
figure
yyaxis left
bar(dayTable.tbin, dayTable.HoursProp)
ylim([0 max(dayTable.HoursProp)]);
ylabel('Proportion of hours per day with sperm whale presence')
yyaxis right
plot(dayTable.tbin, (dayTable.Effort_Sec./dayTable.MaxEffort_Sec)*100, '.', 'Color',gray)
%ylim([-1 101])
xlim([dayTable.tbin(1) dayTable.tbin(end)]) %adjust x-axis to only show data range
ylabel('Percent effort')
ax = gca;
ax.YAxis(2).Color = gray;
title([{'Daily Presence of Sperm whales in the ',titleNAME}])
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence.png']);

%Plot proportion of hours per WEEK with sperm whale presence
figure
yyaxis left
bar(weekTable.tbin, weekTable.HoursNorm)
ylim([0 max(weekTable.HoursNorm)]);
xlim([weekTable.tbin(1),weekTable.tbin(end)])
ylabel('Normalized hours per week with sperm whale presence')
yyaxis right
plot(weekTable.tbin, weekTable.NormEffort_Bin*100,'.', 'Color',[.5 .5 .5])
ylim([-1 101])
ylabel('Percent effort')
ax = gca;
ax.YAxis(2).Color = gray;
title([{'Weekly Presence of Sperm whales in the ',titleNAME}])
saveas(gcf,[saveDir,'\',siteabrev,'WeeklyPresenceHoursNorm.png']);

%Plot proportion of hours per DAY with presence from each group
figure
set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[3 1 7.5 7])
subplot(3,1,1) % SUBPLOT 1: SOCIAL UNITS
yyaxis left
bar(sexbinPresence.tbin,sexbinPresence.FeHoursProp,'FaceColor',mint,'BarWidth',3)
xlim([sexbinPresence.tbin(1),sexbinPresence.tbin(end)])
ylim([0 max(sexbinPresence.FeHoursProp)])
yyaxis right
plot(sexbinPresence.tbin, sexbinPresence.NormEffort_Bin*100,'.','Color',gray)
ylim([-1 101])
title(['Social Units'])
ax = gca;
ax.YAxis(1).Color = mint;
ax.YAxis(2).Color = gray;
subplot(3,1,2) % SUBPLOT 2: MID-SIZE
yyaxis left
bar(sexbinPresence.tbin,sexbinPresence.JuHoursProp,'FaceColor',persimmon,'BarWidth',3)
xlim([sexbinPresence.tbin(1),sexbinPresence.tbin(end)])
ylim([0 max(sexbinPresence.JuHoursProp)])
ax = gca;
ax.YAxis(1).Color = persimmon;
label_y = ylabel('Proportion of hours per day with group presence','Color','k'); % Left y-axis label
label_y.Position(1) = -100;
yyaxis right
plot(sexbinPresence.tbin, sexbinPresence.NormEffort_Bin*100,'.','Color',gray)
ylim([-1 101])
title(['Mid-Size Animals']) % Right y-axis label
ax = gca;
ax.YAxis(2).Color = gray;
ylabel('Percent Effort','Color','k')
subplot(3,1,3) % SUBPLOT 3: MALES
yyaxis left
bar(sexbinPresence.tbin,sexbinPresence.MaHoursProp,'FaceColor',slate,'BarWidth',3)
xlim([sexbinPresence.tbin(1),sexbinPresence.tbin(end)])
ylim([0 max(sexbinPresence.MaHoursProp)])
title(['Males'])
yyaxis right
plot(sexbinPresence.tbin, sexbinPresence.NormEffort_Bin*100,'.','Color',gray)
ylim([-1 101])
ax = gca;
ax.YAxis(1).Color = slate;
ax.YAxis(2).Color = gray;
sgtitle(['Daily Presence of Sperm Whales in the ', titleNAME]) % Overarching title %
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses_Subplots.png']);

%Plot proportion of hours per WEEK with presence from each group
figure
set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[3 1 7.5 7])
subplot(3,1,1) % SUBPLOT 1: SOCIAL UNITS
yyaxis left
bar(weekPresence.tbin,weekPresence.FeHoursProp,'FaceColor',mint,'BarWidth',1)
xlim([weekPresence.tbin(1),weekPresence.tbin(end)])
ylim([0 max(weekPresence.FeHoursProp)])
yyaxis right
plot(weekPresence.tbin, weekPresence.NormEffort_Bin*100,'.','Color',gray) %changed col red -> gray
ylim([-1 101])
title(['Social Units'])
ax = gca;
ax.YAxis(1).Color = mint;
ax.YAxis(2).Color = gray;
subplot(3,1,2) % SUBPLOT 2: MID-SIZE
yyaxis left
bar(weekPresence.tbin,weekPresence.JuHoursProp,'FaceColor',persimmon,'BarWidth',1)
xlim([weekPresence.tbin(1),weekPresence.tbin(end)])
ylim([0 max(weekPresence.JuHoursProp)])
ylabel('Proportion of hours per week with group presence')
ax = gca;
ax.YAxis(1).Color = persimmon;
label_y = ylabel('Proportion of hours per day with group presence','Color','k'); % Left y-axis label
label_y.Position(1) = -100;
yyaxis right
plot(weekPresence.tbin, weekPresence.NormEffort_Bin*100,'.','Color',gray)
ylim([-1 101])
title(['Mid-Size Animals'])
ax = gca;
ax.YAxis(2).Color = gray;
ylabel('Percent Effort','Color','k')
subplot(3,1,3) % SUBPLOT 3: MALES
yyaxis left
bar(weekPresence.tbin,weekPresence.MaHoursProp,'FaceColor',slate,'BarWidth',1)
xlim([weekPresence.tbin(1),weekPresence.tbin(end)])
ylim([0 max(weekPresence.MaHoursProp)])
title(['Males'])
yyaxis right
plot(weekPresence.tbin, weekPresence.NormEffort_Bin*100,'.','Color',gray)
ylim([-1 101])
ax = gca;
ax.YAxis(1).Color = slate;
ax.YAxis(2).Color = gray;
sgtitle(['Weekly Presence of Sperm Whales in the ', titleNAME]) % Overarching title
saveas(gcf,[saveDir,'\',siteabrev,'WeeklyPresence_AllClasses_Subplots.png']);

%% Average yearly plots
%Average yearly presence of proportion of hours per DAY with sperm whale
%presence
figure
bar(meantab365.day, meantab365.mean_HoursProp)
%xlim([0 366])
xlabel('Day')
ylim([0 max(meantab365.mean_HoursProp)]);
ylabel('Average proportion of hours/day')
title([{'Average Daily Presence of Sperm Whales in the ',titleNAME}]);
saveas(gcf,[saveDir,'\',siteabrev,'AverageDailyPresence.png']);

%Average yearly presence of proportion of hours per WEEK with sperm whale
%presence
%retime average table
if length(MD) > 365
n=7; %average every 7th value to account for week
meantab365.day = double(string(meantab365.day));
meantab365ARRAY = table2array(meantab365);
mean_hoursNorm = meantab365ARRAY(:,8);
mean_SEM = meantab365ARRAY(:,9);

% HOURS NORM
s1 = size(mean_hoursNorm, 1);      % Find the next smaller multiple of n
m  = s1 - mod(s1, n);
y  = reshape(mean_hoursNorm(1:m,:), n, []);     % Reshape x to a [n, m/n] matrix
Avg = transpose(sum(y, 1) / n);  % Calculate the mean over the 1st dim

% SEM
y2  = reshape(mean_SEM(1:m,:), n, []);     % Reshape x to a [n, m/n] matrix
Avg2 = transpose(sum(y2, 1) / n);  % Calculate the mean over the 1st dim

meantab52 = (1:52)';
meantab365WEEK = array2table([meantab52 Avg Avg2]);
meantab365WEEK.Properties.VariableNames = {'Week' 'HoursProp' 'SEM'};

%Old code not ready to part with it yet
% meantab365WEEK = squeeze(mean(reshape(table2array(meantab365),[]),1));
% meantab365.datetime = datetime(meantab365.day, 'convertfrom','juliandate');
% meantab365.Week = week(meantab365.datetime);
% meantab365.Week = categorical(meantab365.Week);
% [mean, sem, std, var, range] = grpstats(meantab365.HoursProp, meantab365.Week, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
% meantable = array2table(mean);
% semtable = array2table(sem);
% stdtable = array2table(std);
% vartable = array2table(var);
% rangetable = array2table(range);
% newcol_mean = (1:length(mean))';
% meanarray365 = [newcol_mean mean sem std var range];
% WEEKmeantab365 = array2table(meanarray365);
% WEEKmeantab365.Properties.VariableNames = {'Week' 'HoursProp' 'SEM' 'Std' 'Var' 'Range'};

%make figure
figure
bar(meantab365WEEK.Week, meantab365WEEK.HoursProp,'FaceColor',[.6 .6 .6], 'EdgeColor','flat')
set(gcf,'position',[50 50 700 400])
xlim([0 52])
xlabel('Week')
ylabel('Average proportion of hours/week')
%title({'Average Weekly Presence of Sperm Whales in the ',titleNAME});
hold on
errorbar(meantab365WEEK.Week,meantab365WEEK.HoursProp, -(meantab365WEEK.SEM),meantab365WEEK.SEM, ...
    'Color','k','LineStyle','none')
ylim([0 inf])
saveas(gcf,[saveDir,'\',siteabrev,'AverageWeeklyPresence.png']);
%saveas(gcf,[saveDir,'\',siteabrev,'AverageWeeklyPresence_SansTitle.png']);
else
end

%Average yearly presence of proportion of hours per DAY with presence from
%each group

%Average yearly presence of proportion of hours per WEEK with presence from
%each group

figure
set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[3 1 7.5 7])
subplot(3,1,1)
bar(sexbinPresence.tbin,sexbinPresence.FeHoursProp,'FaceColor',mint,'BarWidth',3)
xlim([sexbinPresence.tbin(1),sexbinPresence.tbin(end)])
title(['Social Units'])
subplot(3,1,2)
bar(sexbinPresence.tbin,sexbinPresence.JuHoursProp,'FaceColor',persimmon,'BarWidth',3)
xlim([sexbinPresence.tbin(1),sexbinPresence.tbin(end)])
label_y = ylabel('Daily Presence (5-min bins)');
label_y.Position(1) = -100;
title(['Mid-Size Animals',])
subplot(3,1,3)
bar(sexbinPresence.tbin,sexbinPresence.MaHoursProp,'FaceColor',slate,'BarWidth',3)
xlim([sexbinPresence.tbin(1),sexbinPresence.tbin(end)])
title(['Males'])
sgtitle(['Daily Presence of Sperm Whales in the ', titleNAME]) % Overarching title
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses_Subplots130.png']);

%Plot daily presence in 5-min bins for all classes in one plot
figure
bar(sexbinPresence.tbin,sexbinPresence.MaleNormBin,'FaceColor',slate,'BarWidth',1)
ylabel('Daily Presence (5-min bins)')
title(['Daily Presence of Each Size Class at ',titleNAME])
hold on
bar(sexbinPresence.tbin,sexbinPresence.MidSizeNormBin,'FaceColor',persimmon,'BarWidth',1)
bar(sexbinPresence.tbin,sexbinPresence.SocialGroupNormBin,'FaceColor',mint,'BarWidth',1)
legend('Males','Mid-Size','Social Units')
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses130.png']);

%Plotting sperm whale presence with presence of each size class over it
dayBinTAB = readtable(dayBinCSV); %read general PM table with presence information
figure
bar(dayBinTAB.tbin,dayBinTAB.NormBin,'k')
hold on
bar(sexbinPresence.tbin,sexbinPresence.MaleNormBin,'FaceColor','c','BarWidth',1)
bar(sexbinPresence.tbin,sexbinPresence.MidSizeNormBin,'FaceColor','b','BarWidth',1)
bar(sexbinPresence.tbin,sexbinPresence.SocialGroupNormBin,'FaceColor','y','BarWidth',1)
ylabel('Daily Presence (5-min bins)')
title(['Daily Presence with Each Size Class at ',titleNAME])
legend('All','Males','Mid-Size','Social Units')
saveas(gcf,[saveDir,'\',siteabrev,'AllDailyPresence_with_SizeClasses130.png']);