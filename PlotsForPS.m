clearvars
close all
%% Parameters defined by user
filePrefix = 'PS'; % File name to match. 
siteabrev = 'PS'; %abbreviation of site.
sp = 'Pm'; % your species code
%saveDir = 'H:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis\BP'; %specify directory to save files
saveDir = 'I:\Shared drives\Pt. Sur\Analyzed data\Sperm whales\SeasonalityAnalysis';
titleNAME = 'Point Sur';
%% load workspace
load([saveDir,'\',siteabrev,'_workspaceStep2.mat']);
%% Retime for weekly presence
%day table
weekTable = retime(dayTable,'weekly','mean');
weekTable.NormEffort_Bin = weekTable.Effort_Sec ./weekTable.MaxEffort_Sec;
weekTable.NormEffort_Bin(isnan(weekTable.NormEffort_Bin)) = 0;
weekTable.HoursProp = weekTable.Hours ./ (weekTable.Effort_Sec ./ (60*60));

%% Plots
%Plot proportion of hours per DAY with sperm whale presence
figure
yyaxis left
bar(dayTable.tbin, dayTable.HoursProp)
ylim([0 max(dayTable.HoursProp)]);
ylabel('Proportion of hours per day with sperm whale presence')
yyaxis right
plot(dayTable.tbin, (dayTable.Effort_Sec./dayTable.MaxEffort_Sec)*100,'.r')
%ylim([-1 200])
ylabel('Percent effort')
title(['Daily Presence of Sperm whales in the ',titleNAME])
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence.png']);

%Plot proportion of hours per WEEK with sperm whale presence
figure
yyaxis left
bar(weekTable.tbin, weekTable.HoursProp)
ylim([0 max(weekTable.HoursProp)]);
ylabel('Proportion of hours per week with sperm whale presence')
yyaxis right
plot(weekTable.tbin, weekTable.NormEffort_Bin*100,'.r')
ylim([-1 101])
ylabel('Percent effort')
title([{'Weekly Presence of Sperm whales in the ',titleNAME}])
saveas(gcf,[saveDir,'\',siteabrev,'WeeklyPresence.png']);

%% Average yearly plots
%Average yearly presence of proportion of hours per DAY with sperm whale
%presence
figure
bar(meantab365.Day, meantab365.HoursProp)
xlim([0 366])
xlabel('Day')
ylim([0 max(meantab365.HoursProp)]);
ylabel('Average proportion of hours/day')
title([{'Average Daily Presence of Sperm Whales in the ',titleNAME}]);
saveas(gcf,[saveDir,'\',siteabrev,'AverageDailyPresence.png']);

%Average yearly presence of proportion of hours per WEEK with sperm whale
%presence
%retime average table
if length(MD) > 365
meantab365.datetime = datetime(meantab365.Day, 'convertfrom','juliandate');
meantab365.Week = week(meantab365.datetime);
meantab365.Week = categorical(meantab365.Week);
[mean, sem, std, var, range] = grpstats(meantab365.HoursProp, meantab365.Week, {'mean','sem','std','var','range'}); %takes the mean of each day of the year
meantable = array2table(mean);
semtable = array2table(sem);
stdtable = array2table(std);
vartable = array2table(var);
rangetable = array2table(range);
newcol_mean = (1:length(mean))';
meanarray365 = [newcol_mean mean sem std var range];
WEEKmeantab365 = array2table(meanarray365);
WEEKmeantab365.Properties.VariableNames = {'Week' 'HoursProp' 'SEM' 'Std' 'Var' 'Range'};
%make figure
figure
bar(WEEKmeantab365.Week, WEEKmeantab365.HoursProp)
xlim([0 52])
xlabel('Week')
ylabel('Average proportion of hours/week')
title({'Average Weekly Presence of Sperm Whales in the ',titleNAME});
hold on
errorbar(WEEKmeantab365.Week,WEEKmeantab365.HoursProp, -(WEEKmeantab365.SEM),WEEKmeantab365.SEM)
saveas(gcf,[saveDir,'\',siteabrev,'AverageWeeklyPresence.png']);
else
end
