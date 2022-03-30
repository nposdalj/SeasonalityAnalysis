clearvars
close all
%% Parameters defined by user
filePrefix = 'WAT_BS'; % File name to match. 
siteabrev = 'BS'; %abbreviation of site.
sp = 'Pm'; % your species code
saveDir = 'H:\My Drive\WAT_TPWS_metadataReduced\SeasonalityAnalysis\BS'; %specify directory to save files
titleNAME = 'Western Atlantic-Blake Spur';
%% load workspace
load([saveDir,'\',siteabrev,'_workspaceStep2.mat']);
load([saveDir,'\',siteabrev,'_workspaceStep3.mat']);
%% Fill in missing days
%day table
dayTable.day= double(dayTable.day);
dayTable = retime(dayTable,'daily','fillwithconstant');

%sex table
binPresence.day = double(binPresence.day);
binPresence = retime(binPresence,'daily','fillwithconstant');
%% Retime for weekly presence
%day table
weekTable = retime(dayTable,'weekly','sum');
weekTable.NormEffort_Bin = weekTable.Effort_Sec ./weekTable.MaxEffort_Sec;
weekTable.NormEffort_Bin(isnan(weekTable.NormEffort_Bin)) = 0;
weekTable.HoursProp = weekTable.Hours ./ (weekTable.Effort_Sec ./ (60*60));

%sex table
weekPresence = retime(binPresence,'weekly','sum');
weekPresence.NormEffort_Bin = weekPresence.Effort_Sec ./weekPresence.MaxEffort_Sec;
weekPresence.NormEffort_Bin(isnan(weekPresence.NormEffort_Bin)) = 0;
weekPresence.FeHoursProp = weekPresence.FeHours ./(weekPresence.Effort_Sec ./ (60*60));
weekPresence.JuHoursProp = weekPresence.JuHours ./(weekPresence.Effort_Sec ./ (60*60));
weekPresence.MaHoursProp = weekPresence.MaHours ./(weekPresence.Effort_Sec ./ (60*60));
%% Plots
%Plot proportion of hours per DAY with sperm whale presence
figure
yyaxis left
bar(dayTable.tbin, dayTable.HoursProp)
ylim([0 max(dayTable.HoursProp)]);
ylabel('Proportion of hours per day with sperm whale presence')
yyaxis right
plot(dayTable.tbin, (dayTable.Effort_Sec./dayTable.MaxEffort_Sec)*100,'.r')
ylim([-1 101])
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

%Plot proportion of hours per DAY with presence from each group
figure
subplot(3,1,1)
yyaxis left
bar(binPresence.tbin,binPresence.FeHoursProp,'FaceColor','y','BarWidth',3)
ylim([0 max(binPresence.FeHoursProp)])
yyaxis right
plot(binPresence.tbin, binPresence.NormEffort_Bin*100,'.r')
ylim([-1 101])
title(['Daily Presence of Social Units in the ',titleNAME])
subplot(3,1,2)
yyaxis left
bar(binPresence.tbin,binPresence.JuHoursProp,'FaceColor','b','BarWidth',3)
ylim([0 max(binPresence.JuHoursProp)])
ylabel('Proportion of hours per day with group presence')
yyaxis right
plot(binPresence.tbin, binPresence.NormEffort_Bin*100,'.r')
ylim([-1 101])
ylabel('Percent Effort')
title(['Daily Presence of Mid-Size Animals in the ',titleNAME])
subplot(3,1,3)
yyaxis left
bar(binPresence.tbin,binPresence.MaHoursProp,'FaceColor','c','BarWidth',3)
ylim([0 max(binPresence.MaHoursProp)])
title(['Daily Presence of Males in the ',titleNAME])
yyaxis right
plot(binPresence.tbin, binPresence.NormEffort_Bin*100,'.r')
ylim([-1 101])
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses_Subplots.png']);

%Plot proportion of hours per WEEK with presence from each group
figure
subplot(3,1,1)
yyaxis left
bar(weekPresence.tbin,weekPresence.FeHoursProp,'FaceColor','y','BarWidth',1)
ylim([0 max(weekPresence.FeHoursProp)])
yyaxis right
plot(weekPresence.tbin, weekPresence.NormEffort_Bin*100,'.r')
ylim([-1 101])
title(['Weekly Presence of Social Units in the ',titleNAME])
subplot(3,1,2)
yyaxis left
bar(weekPresence.tbin,weekPresence.JuHoursProp,'FaceColor','b','BarWidth',1)
ylim([0 max(weekPresence.JuHoursProp)])
ylabel('Proportion of hours per week with group presence')
yyaxis right
plot(weekPresence.tbin, weekPresence.NormEffort_Bin*100,'.r')
ylim([-1 101])
ylabel('Percent Effort')
title(['Weekly Presence of Mid-Size Animals in the ',titleNAME])
subplot(3,1,3)
yyaxis left
bar(weekPresence.tbin,weekPresence.MaHoursProp,'FaceColor','c','BarWidth',1)
ylim([0 max(weekPresence.MaHoursProp)])
title(['Weekly Presence of Males in the ',titleNAME])
yyaxis right
plot(weekPresence.tbin, weekPresence.NormEffort_Bin*100,'.r')
ylim([-1 101])
saveas(gcf,[saveDir,'\',siteabrev,'WeeklyPresence_AllClasses_Subplots.png']);

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

%Average yearly presence of proportion of hours per DAY with presence from
%each group

%Average yearly presence of proportion of hours per WEEK with presence from
%each group

figure
subplot(3,1,1)
bar(binPresence.tbin,binPresence.FeHoursProp,'FaceColor','y','BarWidth',3)
title(['Daily Presence of Social Units in the ',titleNAME])
subplot(3,1,2)
bar(binPresence.tbin,binPresence.JuHoursProp,'FaceColor','b','BarWidth',3)
ylabel('Daily Presence (5-min bins)')
title(['Daily Presence of Mid-Size Animals in the ',titleNAME])
subplot(3,1,3)
bar(binPresence.tbin,binPresence.MaHoursProp,'FaceColor','c','BarWidth',3)
title(['Daily Presence of Males in the ',titleNAME])
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses_Subplots130.png']);

%Plot daily presence in 5-min bins for all classes in one plot
figure
bar(binPresence.tbin,binPresence.MaleNormBin,'FaceColor','c','BarWidth',1)
ylabel('Daily Presence (5-min bins)')
title(['Daily Presence of Each Size Class at ',titleNAME])
hold on
bar(binPresence.tbin,binPresence.JuvenileNormBin,'FaceColor','b','BarWidth',1)
bar(binPresence.tbin,binPresence.FemaleNormBin,'FaceColor','y','BarWidth',1)
legend('Males','Mid-Size','Social Units')
saveas(gcf,[saveDir,'\',siteabrev,'DailyPresence_AllClasses130.png']);

%Plotting sperm whale presence with presence of each size class over it
dayBinTAB = readtable(dayBinCSV); %read general PM table with presence information
figure
bar(dayBinTAB.tbin,dayBinTAB.NormBin,'k')
hold on
bar(binPresence.tbin,binPresence.MaleNormBin,'FaceColor','c','BarWidth',1)
bar(binPresence.tbin,binPresence.JuvenileNormBin,'FaceColor','b','BarWidth',1)
bar(binPresence.tbin,binPresence.FemaleNormBin,'FaceColor','y','BarWidth',1)
ylabel('Daily Presence (5-min bins)')
title(['Daily Presence with Each Size Class at ',titleNAME])
legend('All','Males','Mid-Size','Social Units')
saveas(gcf,[saveDir,'\',siteabrev,'AllDailyPresence_with_SizeClasses130.png']);