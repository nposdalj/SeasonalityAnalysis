clearvars
close all
%% Parameters defined by user
filePrefix = 'Baja_GI'; % File name to match
siteabrev = 'GI'; %abbreviation of site
sp = 'Pm'; % your species code
saveDir = 'G:\Baja\Seasonality'; %specify directory to save files
titleNAME = 'Baja California - Guadalupe Island';
%% load workspace
load([saveDir,'\',siteabrev,'_workspaceStep2.mat']);
load([saveDir,'\',siteabrev,'_workspaceStep3.mat']);
%% Plots
%Plot proportion of hours per DAY with sperm whale presence
%Plot proportion of hours per WEEK with sperm whale presence

%Plot proportion of hours per DAY with presence from each group
%Plot proportion of hours per WEEK with presence from each group

%Average yearly presence of proportion of hours per DAY with sperm whale
%presence
%Average yearly presence of proportion of hours per WEEK with sperm whale
%presence

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