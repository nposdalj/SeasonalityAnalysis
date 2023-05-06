clearvars
close all

% AD's Manuscript Plots
%% Params defined by user

GDrive = 'L';

Sites = {'HZ', 'OC', 'NC', 'BC', 'WC', 'NFC', 'HAT', 'GS', 'BP', 'BS', 'JAX'};
Region = 'WAT';
sp = 'Pm'; % your species code

% Load in data for plotting parameter values (pDet, cue rate, etc.)
plotData = readtable([GDrive ':\My Drive\WAT_TPWS_metadataReduced\Plots\Manuscript_PlottingData.xlsx']);

% Plot HAT A and B in separate lanes? (Not yet configured)
HAT_split = 0; % 1 = Yes; 0 = No

% Color palette
palette.mint = [.4000 .7608 .6471];     palette.persimmon = [.9882 .5529 .3843];    palette.slate = [.5529 .6275 .7961];
palette.grey = [.7529 .7529 .7529];     palette.darkgrey = [.5000 .5000 .5000];     palette.magenta = [.8582 .3500 .6500];
palette.lightgrey = [.8700 .8700 .8700];

% Using GDrive shortcut?
using_shortcut = false;
shortcut = ':\.shortcut-targets-by-id\1FGSX39xqOmreo9qPfPoqhlhUNm1STQB9';

%% Prepare paths and load in data
if using_shortcut == false % Deal with GDrive shortcut if needed
    DrivePath = ':\My Drive';
elseif using_shortcut == true
    DrivePath = shortcut;
end

%% Probability of Detection - by Site and Sex

pDet_data = [plotData.pDetF plotData.pDetJ plotData.pDetM];
pDet_CVdata = [plotData.pDetF_CV plotData.pDetJ_CV plotData.pDetM_CV];

figure('name', 'Probability of Detection by Site and by Size Class')
pDetBar = bar(pDet_data, 'BarWidth', 1, 'EdgeColor', 'k');
ylim([0 0.6])
set(gca, 'xticklabel', plotData.Site)
pDetBar(1).FaceColor = palette.mint;
pDetBar(2).FaceColor = palette.persimmon;
pDetBar(3).FaceColor = palette.slate;
legend({'Social Groups' 'Mid-Size Animals' 'Males'}, 'FontSize', 10)
xlabel('Site')
ylabel('pDet')
set(gcf, 'Position', [100 100 1100 400])
hold on
F_CV_Bars = errorbar((1:height(plotData))-.225, plotData.pDetF, plotData.pDetF_CV, plotData.pDetF_CV, '.k', 'Marker', 'none', 'CapSize', 10, 'LineWidth', 1);
J_CV_Bars = errorbar((1:height(plotData)), plotData.pDetJ, plotData.pDetJ_CV, plotData.pDetJ_CV, '.k', 'Marker', 'none', 'CapSize', 10, 'LineWidth', 1);
M_CV_Bars = errorbar((1:height(plotData))+.225, plotData.pDetM, plotData.pDetM_CV, plotData.pDetM_CV, '.k', 'Marker', 'none', 'CapSize', 10, 'LineWidth', 1);
hold off

saveas(gcf, [GDrive DrivePath '\WAT_TPWS_metadataReduced\Manuscript\Figures\Timeseries\pDetBar.png'])

%% LOAD IN DATA FOR SITE PRESENCE PLOTS

GDrive_correctII = GDrive; % Preserve correct GDrive as it was entered above

for i = 1:length(Sites)
    
    siteabrev = char(Sites(i));
    dataDir = [GDrive,DrivePath,'\',Region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory where workspaces are saved
    
    disp(['Loading Site ' siteabrev '...']) % Load in workspaces
    load([dataDir,'\',siteabrev,'_workspaceStep2.mat']);
    load([dataDir,'\',siteabrev,'_workspaceStep3.mat']);
    
    % Overwrite some path names
    GDrive = GDrive_correctII; % Correct GDrive if overwritten by loading workspace
    dayBinCSV(1) = GDrive; effortXls(1) = GDrive; filename(1) = GDrive; tpwsPath(1) = GDrive;
    
    saveDir = [GDrive,DrivePath,'\',Region,'_TPWS_metadataReduced\Plots\',siteabrev];
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
    
    %% Save site weekPresence, clear weekPresence, and save site no-effort intervals
    SitesWeekPresence.(['Site', num2str(i)]) = weekPresence;
    clear weekTable weekPresence

    if strcmp(siteabrev, 'NFC') % Remove apparently erroneous effort time entries in NFC effort
        effort([1150 1151 1977 1978], :) = [];
    end
    
    SitesEffort.(['Site', num2str(i)]) = timetable2table(effort); % Preserve site-specific effort table for later use (for plotting large periods of no effort)
    noEffortIntervals = effort.Start(2:end) - effort.End(1:end-1); % Get durations of gaps between periods of effort
    largeEffortGapIdx.(['Site', num2str(i)]) = find(noEffortIntervals > days(7)); % Get index of the bins before large (>7 days) effort gaps

end
disp('Finished loading site workspaces and retaining needed variables')

%% Before plotting, remove some rows from NFC's effort table that seem to cause issues when plotting
oldEffort_NFC = SitesEffort.Site6; % Save old effort table - represents how the effort was originally logged
SitesEffort.Site6([1150, 1151, 1977, 1978], :) = []; % In table used for plotting, remove problematic rows
% Basically they have a start or end date that is way off, so a huge non-existent effort gap gets plotted

%% PLOTTING!
%% Stacked Sites - Acoustic Presence of Size Classes

HoursPropClasses = {'FeHoursProp', 'JuHoursProp', 'MaHoursProp'};
PlotColors = {palette.mint, palette.persimmon, palette.slate};
SaveNames = {'SocialGroup', 'MidSize', 'Male'};
PlotTitles = {'Social Groups', 'Mid-Size Animals', 'Adult Males'};

for j = 1:3 % Cycle through size classes
    
    HoursPropClass = char(HoursPropClasses(j)); % Select HoursProp, PlotColor, SaveName, and PlotTitle for this size class
    PlotColor = cell2mat(PlotColors(j));
    SaveName = char(SaveNames(j));
    PlotTitle = char(PlotTitles(j));
    
    figure('name', (HoursPropClass))
    set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[3 0 10 length(Sites)])
    
    x0 = dateshift(datetime('2015-04-26 0:00:00'), 'start', 'week'); % Earliest start time of all sites (OC)
    xf = dateshift(datetime('2019-06-18 14:17:09'), 'end', 'week'); % Latest end time of all sites (GS)
    for i = 1:length(Sites)
        subplot(length(Sites),1,i)
        yyaxis left
        bar(SitesWeekPresence.(['Site', num2str(i)]).tbin+days(3.5),SitesWeekPresence.(['Site', num2str(i)]).(HoursPropClass),'FaceColor',PlotColor,'BarWidth',1)
            % x-positions above shifted forward 3.5 days to center each bar
            % on its respective week, rather than centering it on the
            % first day of its week. Same shift applied to Percent Effort
        xlim([x0 xf])
        % ylim([0 max(SitesWeekPresence.(['Site', num2str(i)]).(HoursPropClass))])
        if j == 1 % Set social group yMax to 0.8
            ylim([0 .8])
        elseif j == 2 % Set mid-size yMax to 0.6
            ylim([0 .6])
        elseif j == 3 % Set male yMax to 0.11 and remove TickLabel at 0.05
            ylim([0 .11])
            ax_current = gca;
            ax_current.YAxis(1).TickLabels(2) = {' '};
        end
        ax = gca;
        ax.YAxis(1).Color = 'k'; % ax.YAxis(1).Color = PlotColor;
        if i == round(length(Sites) / 2)
            % ylabel('Proportion of hours per week with group presence', 'Color', 'k');
            labely = ylabel('Proportion of hours per week with size class presence', 'Color', 'k');
            labely.Position(1) = -430;
        end
        yyaxis right
        plot(SitesWeekPresence.(['Site', num2str(i)]).tbin+days(3.5), SitesWeekPresence.(['Site', num2str(i)]).NormEffort_Bin*100,'.','Color',palette.darkgrey) %changed col red -> gray
        ylim([-1 101])
        ax_current = gca;
        ax_current.YAxis(2).TickLabels(2) = {' '};
        if i == 1
            % title(PlotTitle) % Currently just not including title
        end
        if i ~= length(Sites)
            set(gca, 'XTicklabel', [])
        end
        ax.YAxis(2).Color = 'k'; % ax.YAxis(2).Color = palette.darkgrey;
        if i == round(length(Sites) / 2)
            ylabel('Percent Effort (%) ','Color','k'); % Left y-axis label
        end
        
        hold on
        
        % Shade regions before and after effort begins
        fill([x0, dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(1), 'end', 'week', 'previous'),...
            dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(1), 'end', 'week', 'previous'), x0],...
            [100,100,0,0], palette.lightgrey, 'LineStyle','none') % Shade before effort begins
        fill([xf, dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(end), 'start', 'week', 'next'),...
            dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(end), 'start', 'week', 'next'), xf],...
            [100,100,0,0], palette.lightgrey, 'LineStyle','none') % Shade after effort begins
        
        % Conditionally shade special deployment gaps with no effort
        largeEffortGaps = largeEffortGapIdx.(['Site', num2str(i)]);
        SiteEffort = SitesEffort.(['Site', num2str(i)]);
        for k = 1:(length(largeEffortGaps)) % Plot all no-effort periods with a length of more than a week (those which are indexed in largeEffortGapIdx
            largeEffortGapStart = datetime(SiteEffort.End(largeEffortGaps(k))); % Start time of large effort gap
            largeEffortGapEnd = datetime(SiteEffort.Start(largeEffortGaps(k)+1));
            fill([dateshift(largeEffortGapStart, 'start', 'week', 'next'), dateshift(largeEffortGapStart, 'start', 'week', 'next'),...
                dateshift(largeEffortGapEnd, 'end', 'week', 'previous'), dateshift(largeEffortGapEnd, 'end', 'week', 'previous')],...
                [100,0,0,100], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
        end
        
% Alternative code for shading gaps with no effort -- keeping this because it actually looked different in a few ways        
%         if i == 1 % Site 1, HZ
%             fill([dateshift(datetime(2016,03,25,2,21,21), 'end', 'week'), dateshift(datetime(2016,04,22,18,00,00), 'start', 'week'),...
%                 dateshift(datetime(2016,04,22,18,00,00), 'start', 'week'), dateshift(datetime(2016,03,25,2,21,21), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%             fill([dateshift(datetime(2017,06,19,07,05,06), 'end', 'week'), dateshift(datetime(2017,07,09,00,00,00), 'start', 'week'),...
%                 dateshift(datetime(2017,07,09,00,00,00), 'start', 'week'), dateshift(datetime(2017,06,19,07,05,06), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%             fill([dateshift(datetime(2018,01,13,15,25,06), 'end', 'week'), dateshift(datetime(2018,06,11,17,59,59), 'start', 'week'),...
%                 dateshift(datetime(2018,06,11,17,59,59), 'start', 'week'), dateshift(datetime(2018,01,13,15,25,06), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%         elseif i == 2 % Site 2, OC
%             fill([dateshift(datetime('2016-02-09 6:16:15'), 'end', 'week'), dateshift(datetime('2016-04-24 5:59:59'), 'start', 'week'),...
%                 dateshift(datetime('2016-04-24 5:59:59'), 'start', 'week'), dateshift(datetime('2016-02-09 6:16:15'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%             fill([dateshift(datetime('2017-05-18 6:37:25'), 'end', 'week'), dateshift(datetime('2017-07-06 23:59:59'), 'start', 'week'),...
%                 dateshift(datetime('2017-07-06 23:59:59'), 'start', 'week'), dateshift(datetime('2017-05-18 6:37:25'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%             fill([dateshift(datetime('2018-04-16 6:56:18'), 'end', 'week'), dateshift(datetime('2018-06-10 6:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2018-06-10 6:00:00'), 'start', 'week'), dateshift(datetime('2018-04-16 6:56:18'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%         elseif i == 3 % Site 3, NC
%             fill([dateshift(datetime('2015-09-18 16:28:51'), 'end', 'week'), dateshift(datetime('2016-04-21 18:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2016-04-21 18:00:00'), 'start', 'week'), dateshift(datetime('2015-09-18 16:28:51'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%             fill([dateshift(datetime('2017-5-24 14:53:51'), 'end', 'week'), dateshift(datetime('2017-07-16 18:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2017-07-16 18:00:00'), 'start', 'week'), dateshift(datetime('2017-5-24 14:53:51'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%         elseif i == 4 % Site 4, BC
%             fill([dateshift(datetime('2017-06-10 23:04:05'), 'end', 'week'), dateshift(datetime('2017-06-30 12:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2017-06-30 12:00:00'), 'start', 'week'), dateshift(datetime('2017-06-10 23:04:05'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%         elseif i == 9 % Site 9, JAX
%             fill([dateshift(datetime('2017-10-28 17:27:48'), 'end', 'week'), dateshift(datetime('2018-06-27 0:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2018-06-27 0:00:00'), 'start', 'week'), dateshift(datetime('2017-10-28 17:27:48'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%         end
        
        set(gca, 'Layer', 'top') % Move tick marks to top so they aren't covered by shaded regions
        text(datetime(2015,5,15),70,Sites(i), 'FontSize', 10, 'FontWeight', 'bold') % Site Label
        
        hold off
    end   
    saveas(gcf, [GDrive DrivePath '\' Region '_TPWS_metadataReduced\Manuscript\Figures\Timeseries\AllSites_' SaveName 'ClickPresence.png'])
    
end

%% Stacked Sites - Group Presence of Size Classes

HoursPropClasses = {'FeHoursProp', 'JuHoursProp', 'MaHoursProp'};
PlotColors = {palette.mint, palette.persimmon, palette.slate};
SaveNames = {'SocialGroup', 'MidSize', 'Male'};
PlotTitles = {'Social Groups', 'Mid-Size Animals', 'Adult Males'};

for j = 1:3 % Cycle through size classes
    
    HoursPropClass = char(HoursPropClasses(j)); % Select HoursProp, PlotColor, SaveName, and PlotTitle for this size class
    PlotColor = cell2mat(PlotColors(j));
    SaveName = char(SaveNames(j));
    PlotTitle = char(PlotTitles(j));
    
    figure('name', (HoursPropClass))
    set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[3 0 10 length(Sites)])
    
    x0 = dateshift(datetime('2015-04-26 0:00:00'), 'start', 'week'); % Earliest start time of all sites (OC)
    xf = dateshift(datetime('2019-06-18 14:17:09'), 'end', 'week'); % Latest end time of all sites (GS)
    for i = 1:length(Sites)
        subplot(length(Sites),1,i)
        yyaxis left
        errorbar(SitesWeekPresence.(['Site', num2str(i)]).tbin+days(3.5),SitesWeekPresence.(['Site', num2str(i)]).(HoursPropClass),...
            .1*ones(height(SitesWeekPresence.(['Site', num2str(i)])), 1))
        plot(SitesWeekPresence.(['Site', num2str(i)]).tbin+days(3.5),SitesWeekPresence.(['Site', num2str(i)]).(HoursPropClass),...
            'o', 'MarkerEdgeColor',PlotColor, 'MarkerSize', 4, 'MarkerFaceColor', 'w')
            % x-positions above shifted forward 3.5 days to center each bar
            % on its respective week, rather than centering it on the
            % first day of its week. Same shift applied to Percent Effort
        xlim([x0 xf])
        % ylim([0 max(SitesWeekPresence.(['Site', num2str(i)]).(HoursPropClass))])
        if j == 1 % Set social group yMax to 0.8
            ylim([0 .8])
        elseif j == 2 % Set mid-size yMax to 0.6
            ylim([0 .6])
        elseif j == 3 % Set male yMax to 0.11 and remove TickLabel at 0.05
            ylim([0 .11])
            ax_current = gca;
            ax_current.YAxis(1).TickLabels(2) = {' '};
        end
        ax = gca;
        ax.YAxis(1).Color = 'k'; % ax.YAxis(1).Color = PlotColor;
        if i == round(length(Sites) / 2)
            % ylabel('Proportion of hours per week with group presence', 'Color', 'k');
            labely = ylabel('Proportion of hours per week with size class presence', 'Color', 'k');
            labely.Position(1) = -430;
        end
        yyaxis right
        plot(SitesWeekPresence.(['Site', num2str(i)]).tbin+days(3.5), SitesWeekPresence.(['Site', num2str(i)]).NormEffort_Bin*100,'.','Color',palette.darkgrey) %changed col red -> gray
        ylim([-1 101])
        ax_current = gca;
        ax_current.YAxis(2).TickLabels(2) = {' '};
        if i == 1
            % title(PlotTitle) % Currently just not including title
        end
        if i ~= length(Sites)
            set(gca, 'XTicklabel', [])
        end
        ax.YAxis(2).Color = 'k'; % ax.YAxis(2).Color = palette.darkgrey;
        if i == round(length(Sites) / 2)
            ylabel('Percent Effort (%) ','Color','k'); % Left y-axis label
        end
        
        hold on
        
        % Shade regions before and after effort begins
        fill([x0, dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(1), 'end', 'week', 'previous'),...
            dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(1), 'end', 'week', 'previous'), x0],...
            [100,100,0,0], palette.lightgrey, 'LineStyle','none') % Shade before effort begins
        fill([xf, dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(end), 'start', 'week', 'next'),...
            dateshift(SitesWeekPresence.(['Site', num2str(i)]).tbin(end), 'start', 'week', 'next'), xf],...
            [100,100,0,0], palette.lightgrey, 'LineStyle','none') % Shade after effort begins
        
        % Conditionally shade special deployment gaps with no effort
        largeEffortGaps = largeEffortGapIdx.(['Site', num2str(i)]);
        SiteEffort = SitesEffort.(['Site', num2str(i)]);
        for k = 1:(length(largeEffortGaps)) % Plot all no-effort periods with a length of more than a week (those which are indexed in largeEffortGapIdx
            largeEffortGapStart = datetime(SiteEffort.End(largeEffortGaps(k))); % Start time of large effort gap
            largeEffortGapEnd = datetime(SiteEffort.Start(largeEffortGaps(k)+1));
            fill([dateshift(largeEffortGapStart, 'start', 'week', 'next'), dateshift(largeEffortGapStart, 'start', 'week', 'next'),...
                dateshift(largeEffortGapEnd, 'end', 'week', 'previous'), dateshift(largeEffortGapEnd, 'end', 'week', 'previous')],...
                [100,0,0,100], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
        end
        
% Alternative code for shading gaps with no effort -- keeping this because it actually looked different in a few ways        
%         if i == 1 % Site 1, HZ
%             fill([dateshift(datetime(2016,03,25,2,21,21), 'end', 'week'), dateshift(datetime(2016,04,22,18,00,00), 'start', 'week'),...
%                 dateshift(datetime(2016,04,22,18,00,00), 'start', 'week'), dateshift(datetime(2016,03,25,2,21,21), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%             fill([dateshift(datetime(2017,06,19,07,05,06), 'end', 'week'), dateshift(datetime(2017,07,09,00,00,00), 'start', 'week'),...
%                 dateshift(datetime(2017,07,09,00,00,00), 'start', 'week'), dateshift(datetime(2017,06,19,07,05,06), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%             fill([dateshift(datetime(2018,01,13,15,25,06), 'end', 'week'), dateshift(datetime(2018,06,11,17,59,59), 'start', 'week'),...
%                 dateshift(datetime(2018,06,11,17,59,59), 'start', 'week'), dateshift(datetime(2018,01,13,15,25,06), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%         elseif i == 2 % Site 2, OC
%             fill([dateshift(datetime('2016-02-09 6:16:15'), 'end', 'week'), dateshift(datetime('2016-04-24 5:59:59'), 'start', 'week'),...
%                 dateshift(datetime('2016-04-24 5:59:59'), 'start', 'week'), dateshift(datetime('2016-02-09 6:16:15'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%             fill([dateshift(datetime('2017-05-18 6:37:25'), 'end', 'week'), dateshift(datetime('2017-07-06 23:59:59'), 'start', 'week'),...
%                 dateshift(datetime('2017-07-06 23:59:59'), 'start', 'week'), dateshift(datetime('2017-05-18 6:37:25'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%             fill([dateshift(datetime('2018-04-16 6:56:18'), 'end', 'week'), dateshift(datetime('2018-06-10 6:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2018-06-10 6:00:00'), 'start', 'week'), dateshift(datetime('2018-04-16 6:56:18'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%         elseif i == 3 % Site 3, NC
%             fill([dateshift(datetime('2015-09-18 16:28:51'), 'end', 'week'), dateshift(datetime('2016-04-21 18:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2016-04-21 18:00:00'), 'start', 'week'), dateshift(datetime('2015-09-18 16:28:51'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%             fill([dateshift(datetime('2017-5-24 14:53:51'), 'end', 'week'), dateshift(datetime('2017-07-16 18:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2017-07-16 18:00:00'), 'start', 'week'), dateshift(datetime('2017-5-24 14:53:51'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none', 'Marker', 'none')
%         elseif i == 4 % Site 4, BC
%             fill([dateshift(datetime('2017-06-10 23:04:05'), 'end', 'week'), dateshift(datetime('2017-06-30 12:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2017-06-30 12:00:00'), 'start', 'week'), dateshift(datetime('2017-06-10 23:04:05'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%         elseif i == 9 % Site 9, JAX
%             fill([dateshift(datetime('2017-10-28 17:27:48'), 'end', 'week'), dateshift(datetime('2018-06-27 0:00:00'), 'start', 'week'),...
%                 dateshift(datetime('2018-06-27 0:00:00'), 'start', 'week'), dateshift(datetime('2017-10-28 17:27:48'), 'end', 'week')],...
%                 [100,100,0,0], palette.lightgrey, 'LineStyle','none')
%         end
        
        set(gca, 'Layer', 'top') % Move tick marks to top so they aren't covered by shaded regions
        text(datetime(2015,5,15),70,Sites(i), 'FontSize', 10, 'FontWeight', 'bold') % Site Label
        
        hold off
    end   
    saveas(gcf, [GDrive DrivePath '\' Region '_TPWS_metadataReduced\Manuscript\Figures\Timeseries\AllSites_' SaveName 'GroupPresence.png'])
    
end
