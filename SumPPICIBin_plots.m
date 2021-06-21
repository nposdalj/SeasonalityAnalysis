clearvars
close all
%% Parameters defined by user
filePrefix = 'Baja_GI'; % File name to match. 
siteabrev = 'GI'; %abbreviation of site.
key = 'Baja_GI_'; %after what identifying marker is the deployment number
PlotSiteName = 'Baja Guadualpe Island';
sp = 'Pm'; % your species code
itnum = '3'; % which iteration you are looking for
srate = 200; % sample rate
tpwsPath = 'G:\My Drive\CCE_TPWS_metadataReduced\Baja_GI\TPWS_130'; %directory of TPWS files
effortXls = 'G:\My Drive\CCE_TPWS_metadataReduced\Baja_GI\Pm_Effort.xls'; % specify excel file with effort times
saveDir = 'G:\My Drive\CCE_TPWS_metadataReduced\Baja_GI\Seasonality'; %specify directory to save files
%RC_data = 1; %If you're using RCs data, make this equal to 1, otherwise make it equal to 0.
%% define subfolder that fit specified iteration
if itnum > 1
   for id = 2: str2num(itnum) % iterate id times according to itnum
       subfolder = ['TPWS',num2str(id)];
       tpwsPath = (fullfile(tpwsPath,subfolder));
   end
end
%% Find all TPWS files that fit your specifications (does not look in subdirectories)
% Concatenate parts of file name
% if isempty(sp)
%     detfn = [filePrefix,'.*','TPWS',itnum,'.mat'];
% else
%     detfn = [filePrefix,'.*',sp,'.*TPWS',itnum,'.mat'];
% end

if exist('RC_data','var')
    detfn = [filePrefix,'.*','Delphin_','.*TPWS',itnum,'.*','.mat'];
% else
% if isempty(p.speName)
%     detfn = [p.filePrefix,'.*','TPWS',p.iterationNum,'.mat'];
else
    detfn = [filePrefix,'.*',sp,'.*TPWS',itnum,'.mat'];
end

% Get a list of all the files in the start directory
fileList = cellstr(ls(tpwsPath));

% Find the file name that matches the filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
if isempty(fileMatchIdx)
    % if no matches, throw error
    error('No files matching filePrefix found!')
end
%% Get effort times matching prefix file
%when multiple sites in the effort table
allEfforts = readtable(effortXls); %read effort table
site = siteabrev; %abbreviation used in effort table

siteNUM = unique(allEfforts.Sites);
[sr,~] = size(siteNUM);

if sr > 1
    effTable = allEfforts(ismember(allEfforts.Sites,site),:); %effort is for multiple sites
else
    effTable = allEfforts; %effort is for one site only
end

% make Variable Names consistent
startVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'Start.*Effort'))>0,1,'first');
endVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'End.*Effort'))>0,1,'first');
effTable.Properties.VariableNames{startVar} = 'Start';
effTable.Properties.VariableNames{endVar} = 'End';

Start = datetime(x2mdate(effTable.Start),'ConvertFrom','datenum');
End = datetime(x2mdate(effTable.End),'ConvertFrom','datenum');

effort = timetable(Start,End);
%% Concatenate all detections from the same site
concatFiles = fileList(fileMatchIdx);
%% get default parameters
p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin');
%% Concatenate variables
PPall = []; TTall = []; ICIall = []; % initialize matrices
for idsk = 1 : length(concatFiles)
    % Load file
    fprintf('Loading %d/%d file %s\n',idsk,length(concatFiles),fullfile(tpwsPath,concatFiles{idsk}))
    D = load(fullfile(tpwsPath,concatFiles{idsk}));
    
    % find times outside effort (sometimes there are detections
    % which are from the audio test at the beggining of the wav file)
    within = cell2mat(arrayfun(@(x)sum(isbetween(x,datenum(effort.Start),datenum(effort.End))),D.MTT,'uni',false));
    goodIdx = find(within ~= 0);
    MTT = D.MTT;
    MPP = D.MPP;
    MTT = D.MTT(goodIdx); % only keep the good detections
    MPP = D.MPP(goodIdx);
    MPP(:,2) = idsk;
    
    % concatenate
    TTall = [TTall; MTT];   % group start times
    PPall = [PPall; MPP];   % group peak-to-peak
    
    % Inter-Click Interval
    ici = diff(MTT)*24*60*60*1000; % in ms
    ICIall = [ICIall;[ici; nan]];  % group inter-click interval
end
%% After parfor data may not be sorted. Sort all the variables.
[~,sorted] = sort(TTall);
TTall = TTall(sorted);
PPall_deplo = PPall;
PPall = PPall(sorted);
ICIall = ICIall(sorted);
%% Create timetable per click
tbin = datetime(TTall,'ConvertFrom','datenum');
clickData = timetable(tbin,PPall,ICIall);
clear tbin
%% Convert times to bin vector times
vTT = datevec(TTall);
tbin = datetime([vTT(:,1:4), floor(vTT(:,5)/p.binDur)*p.binDur, ...
    zeros(length(vTT),1)]);
%% create table and get click counts and max pp per bin
data = timetable(tbin,TTall,PPall);
binData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
binData.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binData.Properties.VariableNames{'max_PPall'} = 'maxPP';
binDataIDX = binData.Count >= 5; %only keep bins that have more than 5 clicks
binData = binData(binDataIDX,:); %only keep bins that have more than 5 clicks

positiveCounts = sum(binData.Count);
positiveBins = length(binData.Count);

%% create table and get click counts and mean pp per bin
binData_MEAN = varfun(@mean,data,'GroupingVariable','tbin','InputVariable','PPall');
binData_MEAN.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binData_MEAN.Properties.VariableNames{'mean_PPall'} = 'meanPP';
binData_MEAN = binData_MEAN(binDataIDX,:);%only keep bins that have more than 5 clicks
%% group effort in bins
effort.diffSec = seconds(effort.End-effort.Start);
effort.bins = effort.diffSec/(60*p.binDur);
effort.roundbin = round(effort.diffSec/(60*p.binDur));

secMonitEffort = sum(effort.diffSec);
binMonitEffort = sum(effort.roundbin);

[er,~] = size(effort.Start);

if er > 1
    binEffort = intervalToBinTimetable(effort.Start,effort.End,p); % convert intervals in bins when there is multiple lines of effort
    binEffort.sec = binEffort.bin*(p.binDur*60);
else
    binEffort = intervalToBinTimetable_Only1RowEffort(effort.Start,effort.End,p); % convert intervals in bins when there is only one line of effort
    binEffort.sec = binEffort.bin*(p.binDur*60);
end
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

%group data by week
weekData = retime(dayData,'weekly','mean');

weekEffort = retime(binEffort,'weekly','sum');
weekTable = retime(dayTable,'weekly','sum');

%group data by month
monthData = retime(dayData, 'monthly', 'mean');
monthEffort = retime(binEffort,'monthly','sum');
monthTable = retime(dayTable, 'monthly','sum');
monthTable(~monthTable.Effort_Bin,:)=[]; %removes days with no effort, NOT days with no presence
%% Create plots, binlog and pplog files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weekly group and click count with percent effort reported
weekTable.Percent = (weekTable.Effort_Sec./604800)*100;
figure(6); set(6, 'name','Weekly presence')
set(gca,'defaultAxesColorOrder',[0 0 0;0 0 0]);
subplot(2,1,1);
yyaxis left
bar(weekTable.tbin,weekTable.Count_Click,'k')
xlim([weekTable.tbin(1) weekTable.tbin(end)])
title({'Average Weekly Presence of Sperm Whales', ['at ' PlotSiteName]})
ylabel({'Weekly Mean'; 'of clicks per day'});
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
ylabel({'Weekly Mean'; 'of 5-min bins per day'});
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

%% Plot Peak-to-peak of each deployment over the other
% detfn2 = [filePrefix,'.*','pp2','.mat'];
% fileMatchIdx2 = find(~cellfun(@isempty,regexp(fileList,detfn2))>0);
% concatFiles2 = fileList(fileMatchIdx2);
% 
% MPP_deplo = [];
% for idsk = 1:length(concatFiles2)
%     fullFileName = fullfile(tpwsPath,concatFiles2{idsk});
%     MPP_temp = load(fullFileName);
%     MPP_deplo = [MPP_deplo; MPP_temp.MPP_save]; 
% end
% 
% Value = regexp(MPP_deplo.disk,'\d*','Match');
% MPP_deplo.deployment = str2double(Value);
% 
% G = findgroups(MPP_deplo.deployment);
% [g, gN] = grp2idx(G);
% nhist2 = histcounts(MPP_deplo.MPP,G);
% 
% G = findgroups(MPP_deplo.deployment);
% func = @(x) histcounts(x);
% [G,L] = findgroups(MPP_deplo);
% [MPP_bydeplo,n] = splitapply(func,MPP_deplo.MPP,G);
% L.nhist = MPP_bydeplo;
% 
% nhist2 = splitapply(func,MPP_deplo.MPP,G);
% 
% for kk = 1:length(unique(MPP_deplo.deployment))
%     [nhist2{kk}] = histc(MPP_deplo.MPP(G)
% 
% end
% nhist2 = histc(MPP_deplo.MPP,center);
% figure
% plotHandles = zeros(1,length(concatFiles2));
% plotLabels = cell(1,length(concatFiles2));
% for idsk = 1: length(concatFiles2)
%     plotHandles(idsk) = plot(center,nhist2{idsk});
%     plotLabels{idsk} = [siteabrev,'\_',num2str(idsk)];
%     hold on
% end
% set(gca,'YScale','log')
% xlabel('Peak-Peak Amplitidue (dB)')
% ylabel('Click Counts (log)')
% legend(plotHandles,plotLabels);
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
bar(h5(1),weekEffort.tbin,weekEffort.sec,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
box on;
bar(h5(1),weekData.tbin,weekData.Count_Click,'FaceColor',[0,0,0],'BarWidth', 1)
hold(h5(1), 'off')
hold(h5(2), 'on')
bar(h5(2),weekEffort.tbin,weekEffort.bin,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
bar(h5(2),weekData.tbin,weekData.Count_Bin,'FaceColor',[0,0,0],'BarWidth', 1)
box off
hold(h5(2), 'off')
%set(h5(1),'xticklabel', '');
ylabel(h5(1),{'Weekly mean';'of clicks per day'})
ylabel(h5(2),{'Weekly mean';['of ', num2str(p.binDur), ' min bins per day']})
title(h5(1),'Click Counting')
title(h5(2),'Group Counting');
axis (h5(1),'tight')
axis (h5(2),'tight')
clickMax = max(weekData.Count_Click);
binMax = max(weekData.Count_Bin);
ylim(h5(1),[0 clickMax])
ylim(h5(2),[0 binMax])


% define step according to number of weeks
if length(weekData.tbin) > 53 && length(weekData.tbin) <= 104 % 2 years
    step = calmonths(1);
elseif length(weekData.tbin) > 104 && length(weekData.tbin) <= 209 % 4 years
    step = calmonths(2);
elseif length(weekData.tbin) > 209 && length(weekData.tbin) <= 313 % 6 years
    step = calmonths(3);
elseif length(weekData.tbin) > 313 && length(weekData.tbin) <= 417 % 8 years
    step = calmonths(4);
elseif length(weekData.tbin) >= 417 % more than 8 years
    step = calyears(1);
end
% define tick steps only if more than 1 year of data
if length(weekData.tbin) > 53
    xtickformat('MMMyy')
    xticks(weekData.tbin(1):step:weekData.tbin(end))
    xtickangle(45)
end

% Save plot
weeklyfn2 = [filePrefix,'_',p.speName,'_weeklypresence'];
saveas(gcf,fullfile(saveDir,weeklyfn2),'png')

%% Plot time series peak-peak of all deployments
figure(7)
plot(data.tbin,data.PPall,'.')
title(sprintf('%s N=%d',filePrefix,length(PPall)), 'Interpreter', 'none')
xlabel('Time')
ylabel('Peak-Peak Amplitude (dB)')
timeseries_fn1 = [filePrefix,'_',p.speName,'_TimeSeries_PeakPeak'];
saveas(gcf,fullfile(saveDir,timeseries_fn1),'png');

%% Plot average peak frequency in bins vs. time of all deployments
figure(8)
plot(binData_MEAN.tbin,binData_MEAN.meanPP,'.')
[pp,q] = size(binData_MEAN);
title(sprintf('%s N=%d',filePrefix,pp), 'Interpreter', 'none')
xlabel('Time')
ylabel('Averag Peak-Peak Amplitude per Bin (dB)')
timeseries_fn2 = [filePrefix,'_',p.speName,'_TimeSeries_AveragePPbin'];
saveas(gcf,fullfile(saveDir,timeseries_fn2),'png');

%% Plots to delete later - these were for the lab presentation on 10/19/2020
% figure; bar(binData.tbin,binData.Count)
% ylabel('# of Clicks in a 5-min Bin')
% title('Number of clicks in each 5-minute bin')
% figure;bar(dayData.tbin,dayData.Count_Click)
% ylabel('Daily # of Clicks')
% title('Number of clicks each day')
% figure;bar(dayData.tbin,dayData.Count_Bin)
% ylabel('Daily # of 5-min bins')
% title('Number of 5-min bins each day')
% 
% %plotting on top of one another
% figure
% bar(binData.tbin,binData.Count)
% hold on
% bar(binData_RC.tbin,binData_RC.Count)
% ylabel('# of Clicks in a 5-min Bin')
% title('Number of clicks in each 5-minute bin')
% legend
% 
% figure
% bar(dayData.tbin,dayData.Count_Click)
% hold on
% bar(dayData_RC.tbin,dayData_RC.Count_Click,'r')
% ylabel('Daily # of Clicks')
% title('Number of clicks each day')
% legend('Pm detector','De detector')
% 
% figure
% bar(dayData.tbin,dayData.Count_Bin)
% hold on
% bar(dayData_RC.tbin,dayData_RC.Count_Bin,'r')
% ylabel('Daily # of 5-min bins')
% title('Number of 5-min bins each day')
% legend('Pm detector','De detector')
% 
% %Plotting difference
% dayDiff_Click = dayData.Count_Click - dayData_RC.Count_Click;
% dayDiff_Bin = dayData.Count_Bin - dayData_RC.Count_Bin;
% figure;bar(dayData.tbin,dayDiff_Click)
% ylabel('Difference in Daily # of Clicks')
% title('Difference in the # of clicks each day')
% figure;bar(dayData.tbin,dayDiff_Bin)
% ylabel('Difference in Daily # of 5-min bins')
% title('Difference in the # of 5-min bins each day')