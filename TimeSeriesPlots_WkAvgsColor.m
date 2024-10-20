clearvars
close all

% This script plots weekly average presence of sperm whales at all sites,
% with the average number of hours with presence indicated by color.

%% Parameters defined by user
siteabrevs = {'HZ' 'OC' 'NC' 'BC' 'WC' 'GS' 'BP' 'BS' 'JAX'}; % List of abbreviations of sites in your region.
GDrive = 'P'; % Directory for Google Drive on this device
region = 'WAT';
sp = 'Pm'; % your species code

%% Pre-define meantab365WEEK summary arrays for general and size classes
mtGN365WEEK_reg = nan(52,length(siteabrevs));
mtFE365WEEK_reg = nan(52,length(siteabrevs));
mtJU365WEEK_reg = nan(52,length(siteabrevs));
mtMA365WEEK_reg = nan(52,length(siteabrevs));

%% Loop through sites

saveDir2 = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\All_Sites']; % Where tables and plots are saved

for site = 1:length(siteabrevs)
    siteabrev = char(siteabrevs(site));
    
    %% Load in site-specific workspaces
    dataDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',char(siteabrevs(site))];
        %specify directory where workspaces are saved
    GDrive_correctII = GDrive; % Temporarily store a protected version of GDrive, as the workspaces will overwrite it
        % named GDrive_correctII since workspaces also contain a GDrive_correct that will overwrite this one
    
    disp(['Loading Workspace2 for ' siteabrev '...'])
    load([dataDir,'\',siteabrev,'_workspaceStep2.mat']);
    disp(['Loading Workspace3 for ' siteabrev '...']);
    load([dataDir,'\',siteabrev,'_workspaceStep3.mat']);
    
    GDrive = GDrive_correctII;       % Restore the correct GDrive and apply to path names as needed
    dayBinCSV(1) = GDrive;
    effortXls(1) = GDrive;
    filename(1) = GDrive;
    tpwsPath(1) = GDrive;
    
    %% Make alternate version of sexbinPresence which, like dayTable, excludes days w/ no effort
    sexbinPresence_E = sexbinPresence;
    sexbinPresence_E(isnan(sexbinPresence_E.Effort_Bin),:) = [];
    
%     sexbinPresence_E.FemaleHoursNorm = (sexbinPresence_E.FemaleNormBin*5)./60; % CORRECT HOURSNORM!!!
%     sexbinPresence_E.JuvenileHoursNorm = (sexbinPresence_E.JuvenileNormBin*5)./60; % CORRECT HOURSNORM!!!
%     sexbinPresence_E.MaleHoursNorm = (sexbinPresence_E.MaleNormBin*5)./60; % CORRECT HOURSNORM!!!
    
    %% Make versions of meantabFE365,meantabJU365,meantabMA365 that are exactly analogous to meantab365
    % Just will exclude the stats for class-specific NormClick cuz sexbinPresence doesn't have that
    % Process of design: meantabFE365x is derived from sexbinPresence as
    % meantab365 is derived from dayTable.
    % How meantab365 is defined in Step2:
    % meantab365 = grpstats(timetable2table(dayTable),'day',{'mean','sem','std','var','range'},...
    %   'DataVars',{'HoursProp','HoursNorm','NormBin','NormClick'}); %takes the mean of each day of the year
    % Similar process in Step3 for class-specific, but just HoursProp:
    % [mean, sem, std, var, range] = grpstats(sexbinPresence.FeHoursProp, sexbinPresence.day,...
    %   {'mean','sem','std','var','range'}); %takes the mean of each day of the year
    
    meantabFE365x = grpstats(timetable2table(sexbinPresence_E),'day',{'mean','sem','std','var','range'},...
        'DataVars',{'FeHoursProp','FemaleHoursNorm','FemaleNormBin'});
    meantabJU365x = grpstats(timetable2table(sexbinPresence_E),'day',{'mean','sem','std','var','range'},...
        'DataVars',{'JuHoursProp','JuvenileHoursNorm','JuvenileNormBin'});
    meantabMA365x = grpstats(timetable2table(sexbinPresence_E),'day',{'mean','sem','std','var','range'},...
        'DataVars',{'MaHoursProp','MaleHoursNorm','MaleNormBin'});
    
    %% Get avg weekly presence data and store in summary arrays; modified from TimeSeriesPlots.m
    %Average yearly presence of proportion of hours per WEEK with sperm whale
    %presence
    %retime average tables - for overall AND for size classes
    if length(MD) > 365
        n=7; %average every 7th value to account for week
        
        for i = 1:4 % LOOP THROUGH CLASSES: 1 = GN (General), 2 = FE, 3 = JU, 4 = MA
            if i == 1   % Choose the class data for the current iteration
                mt365_TYPE = meantab365;
            elseif i == 2
                mt365_TYPE = meantabFE365x;
            elseif i == 3
                mt365_TYPE = meantabJU365x;
            elseif i == 4
                mt365_TYPE = meantabMA365x;
            end
            
            mt365_TYPE.day = double(string(mt365_TYPE.day));
            mt365_TYPEARRAY = table2array(mt365_TYPE);
            mean_hoursNorm = mt365_TYPEARRAY(:,8);
            mean_SEM = mt365_TYPEARRAY(:,9);
            
            % HOURS NORM
            s1 = size(mean_hoursNorm, 1);      % Find the next smaller multiple of n
            m  = s1 - mod(s1, n);
            y  = reshape(mean_hoursNorm(1:m,:), n, []);     % Reshape x to a [n, m/n] matrix
            Avg = transpose(sum(y, 1) / n);  % Calculate the mean over the 1st dim
            
            % SEM
            y2  = reshape(mean_SEM(1:m,:), n, []);     % Reshape x to a [n, m/n] matrix
            Avg2 = transpose(sum(y2, 1) / n);  % Calculate the mean over the 1st dim
            
            meantab52 = (1:52)';
            mt365_TYPEWEEK = array2table([meantab52 Avg Avg2]);
            mt365_TYPEWEEK.Properties.VariableNames = {'Week' 'HoursProp' 'SEM'};
            
            if i == 1
                mtGN365WEEK_reg(:,site) = mt365_TYPEWEEK.HoursProp; % Store the HoursProp data above in the array for all sites
            elseif i == 2
                mtFE365WEEK_reg(:,site) = mt365_TYPEWEEK.HoursProp; % Store the HoursProp data above in the array for all sites
            elseif i == 3
                mtJU365WEEK_reg(:,site) = mt365_TYPEWEEK.HoursProp; % Store the HoursProp data above in the array for all sites
            elseif i == 4
                mtMA365WEEK_reg(:,site) = mt365_TYPEWEEK.HoursProp; % Store the HoursProp data above in the array for all sites
            end
        end      
        disp(['Site Data for ' siteabrev ' stored in summary arrays.'])
        
    else
    end
    
end

%% Remake avg weekly presence summary arrays as tables, plot, and save tables and plots.

mtGN365WEEK_BACKUPreg = mtGN365WEEK_reg; % Make backups in case these variables get changed inadvertently
mtFE365WEEK_BACKUPreg = mtFE365WEEK_reg;
mtJU365WEEK_BACKUPreg = mtJU365WEEK_reg;
mtMA365WEEK_BACKUPreg = mtMA365WEEK_reg;

% FOR PLOTTING: Find absolute max across all plots so they can be set to the same color scale
presenceMax = max([max(mtGN365WEEK_reg(:)) max(mtFE365WEEK_reg(:)) max(mtJU365WEEK_reg(:)) max(mtMA365WEEK_reg(:))]);

for i = 1:4 % LOOP THROUGH CLASSES: 1 = GN (General), 2 = FE, 3 = JU, 4 = MA
    
    if i == 1       % Cycle through GN, FE, JU, and MA (generalized as XX) regional weekly meantabs
        mtXX365WEEK_reg = mtGN365WEEK_reg;
        classname = 'GN';
    elseif i == 2
        mtXX365WEEK_reg = mtFE365WEEK_reg;
        classname = 'FE';
    elseif i == 3
        mtXX365WEEK_reg = mtJU365WEEK_reg;
        classname = 'JU';
    elseif i == 4
        mtXX365WEEK_reg = mtMA365WEEK_reg;
        classname = 'MA';
    end
    
    %% Create and save tables
    
    mtXX365WEEK_regTab = [array2table((1:52).') array2table(mtGN365WEEK_reg)]; % Convert to table & include column for week numbers
    colnames = ['HoursProp' classname '_' siteabrevs(1:9)];
    mtXX365WEEK_regTab.Properties.VariableNames(1) = {'Week'};
    for vari = 1:length(siteabrevs) % Give table columns names specifying class and site
        mtXX365WEEK_regTab.Properties.VariableNames(vari+1) = {strcat('HoursProp', classname, '_', char(siteabrevs(vari)))};
    end
    writetable(mtXX365WEEK_regTab, [saveDir2 '\' region '_AverageWeeklyPresence_' classname '.xlsx'])
    
    if i == 1
        mtGN365WEEK_regTab = mtXX365WEEK_regTab;    % Generate the General table as a variable
    elseif i == 2
        mtFE365WEEK_regTab = mtXX365WEEK_regTab;    % Generate the Female table as a variable
    elseif i == 3
        mtJU365WEEK_regTab = mtXX365WEEK_regTab;    % Generate the Juvenile table as a variable
    elseif i == 4
        mtMA365WEEK_regTab = mtXX365WEEK_regTab;    % Generate the Male table as a variable
    end
    
    %% Create and save plots
    
    mtXX365WEEK_regPlot = nan(53,length(siteabrevs)+1);
    mtXX365WEEK_regPlot(1:52,1:length(siteabrevs)) = mtXX365WEEK_reg; % Make version that can be plotted, since pcolor
    % plots vertices (and will therefore be missing one row and one column of values without this adjustment)
    
    allsites_cmap = figure(100+i);
    colormap(parula)
    allsites_cplot = pcolor(0:52, 0:length(siteabrevs), mtXX365WEEK_regPlot.');
    set(gca, 'YTick', 0.5:1:length(siteabrevs)-.5) % Make y-axis ticks line up with grid squares, rather than sit on gridlines
    set(gca, 'YTickLabel', ['' siteabrevs]) % Assign site labels correctly to y-axis
    set(gca, 'YDir', 'reverse') % This is somewhat specific to WAT, where sites can be listed from north to south
    set(gca, 'XTick', -0.5:5:50.5) % Make y-axis ticks line up with grid squares, rather than sit on gridlines
    set(gca, 'XTickLabel', 0:5:50) % Assign week numbers correctly to x-axis
    set(gcf, 'Position', [100 100 750 350])
    set(allsites_cplot, 'EdgeColor', 'white')
    colbar = colorbar;
    %caxis([0 presenceMax])
    title(['Avg Weekly Presence of Sperm Whales in ' region ' - ' classname])
    ylabel('Site')
    xlabel('Week')
    ylabel(colbar,{'\fontsize{11} Average proportion of';'\fontsize{11} hours/week with presence'})
    saveas(allsites_cplot, [saveDir2, '\' 'WAT_AverageWeeklyPresence_' classname], 'png');
    
end