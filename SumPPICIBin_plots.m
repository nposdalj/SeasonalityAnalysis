clearvars
close all;clear all;clc;

% Uses output from Step 2 to create generic plots (weekly presence, ICI, peak to peak recieved level, etc.)
% Accounts for effort but doesn't normalize bin or click count
%% Parameters defined by user
%Site names and paths
filePrefix = 'OC'; % File name to match. 
siteabrev = 'OC'; %abbreviation of site.
region = 'WAT'; %region
sp = 'Pm'; % your species code
GDrive = 'I'; %Google Drive
PlotSiteName =  'Oceanographers Canyon';
saveDirr = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\',siteabrev,'\']; %specify directory to save files
%% load workspace
GDrive_correct = GDrive; % Preserve correct GDrive as it was entered above
load([saveDirr,'\',siteabrev,'_workspaceStep2.mat']);
GDrive = 'I'; %Correct GDrive for SWAL1
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\',siteabrev,'\']; %specify directory to save files
p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin'); % get default parameters -- make sure these match for your species

% Overwrite some path names
GDrive = GDrive_correct; %Correct GDrive if overwritten by loading workspace
effortXls(1) = GDrive;
saveDirr(1) = GDrive;
tpwsPath(1) = GDrive;
%% Group data by week for plotting
dayTable(:,{'Season','month','Year','day'}) = [];
weekData = retime(dayTable,'weekly','mean');
weekEffort = retime(dayTable,'weekly','sum');
weekTable = retime(dayTable,'weekly','sum');

%Percent effort
weekTable.Percent = (weekTable.Effort_Sec./weekTable.MaxEffort_Sec)*100;
idx1 = weekTable.MaxEffort_Bin == 0;
idx2 = isnan(weekTable.Percent);
idx = idx1 + idx2;
weekTable.Percent(idx == 2) = 0;
%% Calculate average peak to peak per 5 min bins
binData_MEAN = varfun(@mean,data,'GroupingVariable','tbin','InputVariable','PPall');
binData_MEAN(binData_MEAN.GroupCount < 5,:) = []; %identify any bins with less than 5 clicks and delete them
%% Create plots, binlog and pplog files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weekly group and click count with percent effort reported
figure(6); set(6, 'name','Weekly presence')
set(gca,'defaultAxesColorOrder',[0 0 0;0 0 0]);
subplot(2,1,1);
yyaxis left
bar(weekTable.tbin,weekTable.Count_Click,'k')
xlim([weekTable.tbin(1) weekTable.tbin(end)])
title({'Weekly Presence of Sperm Whales', ['at ' PlotSiteName]})
ylabel('# of Clicks');
hold on
yyaxis right
plot(weekTable.tbin,weekTable.Percent,'.r')
ylim([-1 101]);
ylabel('Percent Effort')
if length(weekTable.tbin) > 53 && length(weekTable.tbin) <= 104 % 2 years
    step = calmonths(1);
elseif length(weekTable.tbin) > 104 && length(weekTable.tbin) <= 209 % 4 years
    step = calmonths(2);
elseif length(weekTable.tbin) > 209 && length(weekTable.tbin) <= 313 % 6 years
    step = calmonths(3);
elseif length(weekTable.tbin) > 313 && length(weekTable.tbin) <= 417 % 8 years
    step = calmonths(4);
elseif length(weekTable.tbin) >= 417 % more than 8 years
    step = calyears(1);
end
% define tick steps only if more than 1 year of data
if length(weekTable.tbin) > 53
    xtickformat('MMMyy')
    xticks(weekTable.tbin(1):step:weekTable.tbin(end))
    xtickangle(45)
end
hold off

subplot(2,1,2);
yyaxis left
bar(weekTable.tbin,weekTable.Count_Bin,'k')
xlim([weekTable.tbin(1) weekTable.tbin(end)])
ylabel('# of Bins');
hold on
yyaxis right
plot(weekTable.tbin,weekTable.Percent,'.r')
ylim([-1 101]);
ylabel('Percent Effort')
col = [0 0 0];
set(gcf,'defaultAxesColorOrder',[col;col]);

if length(weekTable.tbin) > 53 && length(weekTable.tbin) <= 104 % 2 years
    step = calmonths(1);
elseif length(weekTable.tbin) > 104 && length(weekTable.tbin) <= 209 % 4 years
    step = calmonths(2);
elseif length(weekTable.tbin) > 209 && length(weekTable.tbin) <= 313 % 6 years
    step = calmonths(3);
elseif length(weekTable.tbin) > 313 && length(weekTable.tbin) <= 417 % 8 years
    step = calmonths(4);
elseif length(weekTable.tbin) >= 417 % more than 8 years
    step = calyears(1);
end
% define tick steps only if more than 1 year of data
if length(weekTable.tbin) > 53
    xtickformat('MMMyy')
    xticks(weekTable.tbin(1):step:weekTable.tbin(end))
    xtickangle(45)
end
hold off
 
% Save plot
weeklyfn = [filePrefix,'_',p.speName,'_weeklypresence_withPercentEffort'];
saveas(gcf,fullfile(saveDir,weeklyfn),'png')
%% Plot Inter-Click Interval
iciIdx = find(ICIall > p.iciRange(1) & ICIall < p.iciRange(2));
% statistics
miciSel = mean(ICIall(iciIdx));
sdiciSel = std(ICIall(iciIdx));
moiciSel = mode(ICIall(iciIdx));
meiciSel = median(ICIall(iciIdx));

figure(1); set(1,'name','Inter-Click Interval')
h1 = gca;
centerIci = p.iciRange(1):1:p.iciRange(2);
[nhist] = histc(ICIall(iciIdx),centerIci);
bar(h1,centerIci,nhist, 'barwidth', 1, 'basevalue', 1)
xlim(h1,p.iciRange);
title(h1,sprintf('%s N=%d',filePrefix,length(ICIall)), 'Interpreter', 'none')
xlabel(h1,'Inter-Click Interval (ms)')
ylabel(h1,'Counts')
% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', miciSel);
stdlabel = sprintf('Std = %0.2f', sdiciSel);
melabel = sprintf('Median = %0.2f', meiciSel);
molabel = sprintf('Mode = %0.2f', moiciSel);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% save ici data and figure
icifn = [filePrefix,'_',p.speName,'_ici'];
saveas(h1,fullfile(saveDir,icifn),'png')
%% Plot Peak-to-peak per click
% statistics
mpp = mean(PPall);
sdpp = std(PPall);
mopp = mode(PPall);
mepp = median(PPall);

% Plot histogram
figure(2); set(2,'name','Received Levels')
h2 = gca;
center = 120:1:p.p1Hi;
[nhist] = histc(PPall,center);
bar(h2,center, nhist, 'barwidth', 1, 'basevalue', 1)
title(h2,sprintf('%s N=%d',filePrefix,length(PPall)), 'Interpreter', 'none')
xlabel(h2,'Peak-Peak Amplitude (dB)')
ylabel(h2,'Click Counts')
xlim(h2,[min(PPall), p.p1Hi])
% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mpp);
stdlabel = sprintf('Std = %0.2f', sdpp);
melabel = sprintf('Median = %0.2f', mepp);
molabel = sprintf('Mode = %0.2f', mopp);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
ppfn = [filePrefix,'_',p.speName,'_pp'];
saveas(h2,fullfile(saveDir,ppfn),'png')
%% Plot Peak-to-peak per bin
% statistics
mppBin = mean(binData.maxPP);
sdppBin = std(binData.maxPP);
moppBin = mode(binData.maxPP);
meppBin = median(binData.maxPP);

% Plot histogram
figure(3); set(3,'name','Received Levels per Bin')
h3 = gca;
centerBin = 120:1:p.p1Hi;
[nhistBin] = histc(binData.maxPP,centerBin);
bar(h3,centerBin, nhistBin, 'barwidth', 1, 'basevalue', 1)
title(h3,sprintf('%s N=%d',filePrefix,length(binData.maxPP)), 'Interpreter', 'none')
xlabel(h3,'Peak-Peak Amplitude (dB)')
ylabel(h3,[num2str(p.binDur),' min Bin Counts'])
xlim(h3,[120, p.p1Hi])
% create labels and textbox
mnlabelBin = sprintf('Mean = %0.2f', mppBin);
stdlabelBin = sprintf('Std = %0.2f', sdppBin);
melabelBin = sprintf('Median = %0.2f', meppBin);
molabelBin = sprintf('Mode = %0.2f', moppBin);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabelBin,stdlabelBin,...
    melabelBin,molabelBin});

% Save plot
binfn = [filePrefix,'_',p.speName,'_bin'];
saveas(h3,fullfile(saveDir,binfn),'png')
%% Plot weekly mean of detections 
figure(5); set(5,'name','Weekly presence','DefaultAxesColor',[.8 .8 .8]) 
set(gca,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
h5(1) = subplot(2,1,1);
h5(2) = subplot(2,1,2);
hold(h5(1), 'on')
bar(h5(1),weekEffort.tbin,weekEffort.Effort_Sec,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
box on;
bar(h5(1),weekData.tbin,weekData.Count_Click,'FaceColor',[0,0,0],'BarWidth', 1)
hold(h5(1), 'off')
hold(h5(2), 'on')
bar(h5(2),weekEffort.tbin,weekEffort.Effort_Bin,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
bar(h5(2),weekData.tbin,weekData.Count_Bin,'FaceColor',[0,0,0],'BarWidth', 1)
box off
hold(h5(2), 'off')
%set(h5(1),'xticklabel', '');
ylabel(h5(1),'# of Clicks')
ylabel(h5(2),'# of Bins')
title(h5(1),'Weekly Click Counting')
title(h5(2),'Weekly Group Counting');
axis (h5(1),'tight')
axis (h5(2),'tight')
clickMax = max(weekData.Count_Click);
binMax = max(weekData.Count_Bin);
ylim(h5(1),[0 clickMax])
ylim(h5(2),[0 binMax])

% Save plot
weeklyfn2 = [filePrefix,'_',p.speName,'_weeklypresence'];
saveas(gcf,fullfile(saveDir,weeklyfn2),'png')
%% Plot time series peak-peak of all deployments
figure(7)
plot(data.tbin,data.PPall,'.')
title(sprintf('%s N=%d',filePrefix,length(PPall)), 'Interpreter', 'none')
xlabel('Time')
ylabel('Peak-Peak Amplitude (dB)')

% Save plot
timeseries_fn1 = [filePrefix,'_',p.speName,'_TimeSeries_PeakPeak'];
saveas(gcf,fullfile(saveDir,timeseries_fn1),'png');
%% Plot average peak frequency in bins vs time of all deployments
figure(8)
plot(binData_MEAN.tbin,binData_MEAN.mean_PPall,'.')
[pp,q] = size(binData_MEAN);
title(sprintf('%s N=%d',filePrefix,pp), 'Interpreter', 'none')
xlabel('Time')
ylabel('Averag Peak-Peak Amplitude per Bin (dB)')

% Save plot
timeseries_fn2 = [filePrefix,'_',p.speName,'_TimeSeries_AveragePPbin'];
saveas(gcf,fullfile(saveDir,timeseries_fn2),'png');