clearvars
close all
% This script creates time series plots for each site.
%% Parameters defined by user
filePrefixes = {'HZ' 'OC' 'NC' 'BC' 'WC' 'GS' 'BP' 'BS' 'JAX'}; % File name to match. 
siteabrevs = {'HZ' 'OC' 'NC' 'BC' 'WC' 'GS' 'BP' 'BS' 'JAX'}; %abbreviation of site.
GDrive = 'H'; %directory for Google Drive
region = 'WAT';
sp = 'Pm'; % your species code

%Pre-define meantab365WEEK summary arrays for general and size classes
mt365WEEK_reg = nan(52,length(siteabrevs));
mtFE365WEEK_reg = nan(52,length(siteabrevs));
mtJU365WEEK_reg = nan(52,length(siteabrevs));
mtMA365WEEK_reg = nan(52,length(siteabrevs));

for site = 1:length(siteabrevs)
    siteabrev = char(siteabrevs(site));

%% Load in site-specific data
%titleNAME = 'Oceanographers Canyon';
dataDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',char(siteabrevs(site))]; %specify directory where workspaces are saved

%% load workspace
load([dataDir,'\',siteabrev,'_workspaceStep2.mat']);
disp(['Workspace2 loaded for ' siteabrev])
load([dataDir,'\',siteabrev,'_workspaceStep3.mat']);
disp(['Workspace3 loaded for ' siteabrev]);

GDrive = 'H';       %Correct GDrive
dayBinCSV(1) = GDrive;
effortXls(1) = GDrive;
filename(1) = GDrive;
GDrive = GDrive;
tpwsPath(1) = GDrive;

saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\',siteabrev];

% %% Fill in missing days
% %day table
% dayTable.day= double(dayTable.day);
% dayTable = retime(dayTable,'daily','fillwithconstant');
% 
% %sex table
% sexbinPresence.day = double(sexbinPresence.day);
% sexbinPresence = retime(sexbinPresence,'daily','fillwithconstant');
% 
% %% Retime for weekly presence
% %day table
% weekTable = retime(dayTable,'weekly','sum');
% weekTable.NormEffort_Bin = weekTable.Effort_Sec ./weekTable.MaxEffort_Sec;
% weekTable.NormEffort_Bin(isnan(weekTable.NormEffort_Bin)) = 0;
% weekTable.HoursProp = weekTable.Hours ./ (weekTable.Effort_Sec ./ (60*60));
% 
% %sex table
% weekPresence = retime(sexbinPresence,'weekly','sum');
% weekPresence.NormEffort_Bin = weekPresence.Effort_Sec ./weekPresence.MaxEffort_Sec;
% weekPresence.NormEffort_Bin(isnan(weekPresence.NormEffort_Bin)) = 0;
% weekPresence.FeHoursProp = weekPresence.FeHours ./(weekPresence.Effort_Sec ./ (60*60));
% weekPresence.JuHoursProp = weekPresence.JuHours ./(weekPresence.Effort_Sec ./ (60*60));
% weekPresence.MaHoursProp = weekPresence.MaHours ./(weekPresence.Effort_Sec ./ (60*60));

% ADDED BY AD: Make a new version of meantab365 that is analogous to the class-specific tables
meantab365GEN = meantab365;

%Average yearly presence of proportion of hours per WEEK with sperm whale
%presence
%retime average tables - for overall AND for size classes
if length(MD) > 365
n=7; %average every 7th value to account for week

for i = 1:4
    
    if i == 1
        mt365_TYPE = meantabGEN365;
    elseif i == 2
        mt365_TYPE = meantabFE365;
    elseif i == 3
        mt365_TYPE = meantabJU365;
    elseif i == 4
        mt365_TYPE = meantabMA365;
    end

mt365_TYPE.day = double(string(mt365_TYPE.day));
mt365_TYPEARRAY = table2array(mt365_TYPE);
mean_hoursNorm = mt365_TYPEARRAY(:,8);
mean_SEM = mt365_TYPEARRAY(:,9);

% HOURS NORM
s1 = size(mean_hoursNorm, 1);      % Find the next smaller multiple of n % <- AHAHA THE CLASS-SPECIFIC TABLES DON'T HAVE THIS RIP AARON
m  = s1 - mod(s1, n);
y  = reshape(mean_hoursNorm(1:m,:), n, []);     % Reshape x to a [n, m/n] matrix
Avg = transpose(sum(y, 1) / n);  % Calculate the mean over the 1st dim

% SEM
y2  = reshape(mean_SEM(1:m,:), n, []);     % Reshape x to a [n, m/n] matrix
Avg2 = transpose(sum(y2, 1) / n);  % Calculate the mean over the 1st dim

meantab52 = (1:52)';
mt365_TYPE = array2table([meantab52 Avg Avg2]);
mt365_TYPE.Properties.VariableNames = {'Week' 'HoursProp' 'SEM'};

if i == 1
    mt365WEEK_reg(:,site) = meantab365WEEK.HoursProp; % Store the HoursProp data above in the array for all sites
elseif i == 2
    mtFE365WEEK_reg(:,site) = meantabFE365WEEK.HoursPropFE; % Store the HoursProp data above in the array for all sites
elseif i == 3
    mtJU365WEEK_reg(:,site) = meantabJU365WEEK.HoursPropJU; % Store the HoursProp data above in the array for all sites
elseif i == 4
    mtMA365WEEK_reg(:,site) = meantabMA365WEEK.HoursPropMA; % Store the HoursProp data above in the array for all sites
end
end

disp('Site Data stored summary arrays.')


% % Make figure
% figure(site)
% bar(meantab365WEEK.Week, meantab365WEEK.HoursProp,'FaceColor',[.6 .6 .6], 'EdgeColor','flat')
% set(gcf,'position',[50 50 700 400])
% xlim([0 52])
% xlabel('Week')
% ylabel('Average proportion of hours/week')
% %title({'Average Weekly Presence of Sperm Whales in the ',titleNAME});
% hold on
% errorbar(meantab365WEEK.Week,meantab365WEEK.HoursProp, -(meantab365WEEK.SEM),meantab365WEEK.SEM, ...
%     'Color','k','LineStyle','none')
% ylim([0 max(meantab365WEEK.HoursProp+meantab365WEEK.SEM)+1])
% %saveas(gcf,[saveDir,'\',siteabrev,'AverageWeeklyPresence.png']);
% %saveas(gcf,[saveDir,'\',siteabrev,'AverageWeeklyPresence_SansTitle.png']);
else
end

end

meantab365WEEK_regionalBACKUP = mt365WEEK_reg;
BACKUPmeantab365WEEK_regiona = mt365WEEK_reg;

meantab365WEEK_regPlot = nan(53,length(siteabrevs)+1);
meantab365WEEK_regPlot(1:52,1:length(siteabrevs)) = mt365WEEK_reg; % Make version that can be plotted, since pcolor plots vertices
    % (and will therefore be missing one row and one column of values)

allsites_cmap = figure(101);
colormap(parula)
allsites_cplot = pcolor(0:length(meantab365WEEK.Week), 0:length(siteabrevs), meantab365WEEK_regPlot.');
set(gca, 'YTick', 0.5:1:length(siteabrevs)-.5) % Make y-axis ticks line up with grid squares, rather than sit on gridlines
set(gca, 'YTickLabel', ['' siteabrevs]) % Assign site labels correctly to y-axis
set(gca, 'YDir', 'reverse')
set(gca, 'XTick', -0.5:5:50.5) % Make y-axis ticks line up with grid squares, rather than sit on gridlines
set(gca, 'XTickLabel', [0:5:50]) % Assign week numbers correctly to x-axis
set(allsites_cplot, 'EdgeColor', 'white')
colbar = colorbar;
title(['Average Weekly Presence of Sperm Whales in ' region ' Region'])
ylabel('Site')
xlabel('Week')
ylabel(colbar,'\fontsize{11} Average proportion of hours/week with presence')

%% Plot - BY CLASS
allsites_cmapF = figure(102);
colormap(parula)
allsites_cpotF = pcolor(0:52, 0:length(siteabrevs), weekPresenceF_regPlot.');