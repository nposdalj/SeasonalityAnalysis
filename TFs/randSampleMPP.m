% randSampleMPP

% Note: IDAll file has all click times, and all pp amplitudes for one
% deployment, could just be [MTT, MPP], merged cumulatively across all
% disks.

% This code was provided by KEF and modified by NP 08312022
close all;clear all; clc;
%% User Definied Parameters
siteabrevv = {'QN'}; %abbreviation of site.
GDrive = 'I'; %directory for Google Drive
region = 'GofAK';
sp = 'Pm'; % your species code
itnum = '3'; % which iteration you are looking for (which TPWS folder)
NN = 100000; %random clicks to choose from each deployment
RLmin = 125; %Recieved level threshold of your data
RLmax = 160; %Max recieved level you're intersted in plotting
ReRun = 1; %If workspace exists, load that instead of re-running (0) or create new workspace regardless if one already exists (1)
%% Calculate variable sizes
[~,gg] = size(siteabrevv);
rlVec = RLmin:1:RLmax;
%% Load HARP data
saveDIR = [GDrive,':\My Drive\TestTFs']; %directory where to save output
HARPsum = [saveDIR,'\HARPdataSummary.xlsx']; %HARP data summary sheet
dtable = readtable(HARPsum);
stxt = size(dtable); 
tfnum = [];
ifoundx = 0;
%% Load workspaces
for ss = 1:gg
siteabrev = siteabrevv{ss}; %choose site
tpwsDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\TPWS_125\',siteabrev]; %directory of TPWS files
dataDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory where workspaces are saved
effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx'];% specify excel file with effort times
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\',siteabrev]; % add desired saveDir back in here
file = dir([dataDir,['\',siteabrev,'_simulatedRL.mat']]); %if workspace exists, it would be saved here
%% Create empty variables
counterI = 1;
%% define subfolder that fit specified iteration
if itnum > 1
   for id = 2: str2num(itnum) % iterate id times according to itnum
       subfolder = ['TPWS',num2str(id)];
       tpwsPath = (fullfile(tpwsDir,subfolder));
   end
end
%% Find all TPWS files that fit your specifications (does not look in subdirectories)
% Concatenate parts of file name
detfn = [region,'_',siteabrev,'.*',sp,'.*TPWS',itnum,'.mat'];

% Get a list of all the files in the start directory
fileList = cellstr(ls(tpwsPath));

% Find the file name that matches the filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);

if isempty(fileMatchIdx)
    itnum = '2'; % which iteration you are looking for (which TPWS folder)
    % define subfolder that fit specified iteration
    if itnum > 1
        for id = 2: str2num(itnum) % iterate id times according to itnum
            subfolder = ['TPWS',num2str(id)];
            tpwsPath = (fullfile(tpwsDir,subfolder));
        end
    % Concatenate parts of file name
    detfn = ['.*',sp,'.*TPWS',itnum,'.mat'];

    % Get a list of all the files in the start directory
    fileList = cellstr(ls(tpwsPath));

    % Find the file name that matches the filePrefix
    fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
end
end
%% Get effort times matching prefix file
allEfforts = readtable(effortXls); %read effort table
effTable = allEfforts; %effort is for one site only

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
siteDiskList = cell(length(concatFiles),1);

for s = 1: length(concatFiles)
    a = cell2mat(strfind(concatFiles(s),'_disk'));
    siteDiskList{s,1} = concatFiles{s}(1:a-1);
end

siteDisk = unique(siteDiskList);
[qq,~] = size(siteDisk);
SiteMPP = cell(1,qq);
SiteMTT = cell(1,qq);
SiteICI = cell(1,qq);
histSubset = cell(1,qq);
spMean = cell(1,qq);
spStd = cell(1,qq);

if ReRun == 1 || ReRun == 0 && isempty(file)
    disp(['Loading data from site ',siteabrev])
for i = 1:length(siteDisk)
    disp(['loading times from: ', siteDisk{i}]);
    index = strfind(siteDiskList, siteDisk{i});
    siteDiskIdx = find(not(cellfun('isempty', index)));
%% Concatenate variables
PPall = []; TTall = []; ICIall = []; PeakFrall = [];% initialize matrices
for idsk = 1 : length(siteDiskIdx)
    % Load file
    disk = siteDiskIdx(idsk);
    fprintf('Loading %d/%d file %s\n',idsk,length(siteDiskIdx),fullfile(tpwsPath,concatFiles{disk}))
    D = load(fullfile(tpwsPath,concatFiles{disk}));
    
    % find times outside effort (sometimes there are detections
    % which are from the audio test at the beggining of the wav file)
    within = cell2mat(arrayfun(@(x)sum(isbetween(x,datetime(effort.Start),datetime(effort.End))),datetime(D.MTT,'ConvertFrom','datenum'),'uni',false));
    goodIdx = find(within ~= 0);
    MTT = D.MTT(goodIdx); % only keep the good detections
    MPP = D.MPP(goodIdx);
    
    % concatenate
    TTall = [TTall; MTT];   % group start times
    PPall = [PPall; MPP];   % group peak-to-peak
    
    % Inter-Click Interval
    ici = diff(MTT)*24*60*60*1000; % in ms
    ICIall = [ICIall;[ici; nan]];  % group inter-click interval
    %% After parfor data may not be sorted. Sort all the variables.
    [~,sorted] = sort(TTall);
    TTall = TTall(sorted);
    PPall = PPall(sorted);
    ICIall = ICIall(sorted);
    %% Save variables
    SiteMPP{i} = PPall;
    SiteMTT{i} = TTall;
    SiteICI{i} = ICIall;
end
end
N = round(min(cellfun(@(c) size(c,1), SiteMPP))/1.3); %come up with a sample size that's 1.3 the size of the min number of clicks

for i = 1:length(siteDisk)
%Find min length of SiteMPP to adjust N accordingly
for itr = 1:50
    p = randsample(length(SiteMPP{i}),N); %choose N random click samples
    histSubset{i} = [histSubset{i};hist(SiteMPP{i}(p),rlVec)];
end
end

figure
if length(siteDisk) > 12
    for i = 1:length(siteDisk)
    subplot(5,5,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDisk{i}])
    counterI = counterI+1;
    end
elseif length(siteDisk) > 9 && length(siteDisk) < 12
for i = 1:length(siteDisk)
    subplot(3,4,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDisk{i}])
    counterI = counterI+1;
end
elseif length(siteDisk) > 4 && length(siteDisk) < 10
    for i = 1:length(siteDisk)
    subplot(3,3,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDisk{i}])
    counterI = counterI+1;
    end
else
for i = 1:length(siteDisk)
    subplot(2,2,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDisk{i}])
    counterI = counterI+1;
end
end
plotName = [saveDir,'\',region,'_',siteabrev,'_SubPlotRecievedLevel'];
saveas(gcf,[plotName,'.png'])
saveas(gcf,[plotName,'.fig'])

%Find matching TF to include in the plot
for itab = 1:stxt(1)
    ifound = strfind(dtable.Data_ID(itab),siteabrevv);
    if cell2mat(ifound) >0
        ifoundx = ifoundx + 1;
        tfnum(ifoundx) = str2double(dtable.PreAmp(itab));
    end
end

c = distinguishable_colors(length(siteDisk)); %more colors

%legend table
legmat = [];
for iR = 1:length(siteDisk)
    legmat{iR} = [cell2mat(siteDisk(iR)),' ',num2str(tfnum(iR)),' n = ',num2str(length(SiteMPP{iR}))];
end

figure
for iR = 1:length(siteDisk)
    spMean{iR} = mean(histSubset{iR});
    spStd{iR} = std(histSubset{iR});
    y = log10(spMean{iR});
    ystd1 = log10(spMean{iR}+spStd{iR});
    ystd2 = log10(spMean{iR}-spStd{iR});
    badY = isinf(y);
    y(badY)=[];
    ystd1(badY)=[];
    ystd2(badY)=[];
    h = plot(rlVec(~badY),y);
    hold on
    x2 = [rlVec(~badY),fliplr(rlVec(~badY))];
    inBetween = ([ystd1,fliplr(ystd2)]);
    fill(x2, inBetween,c(iR,:),'edgecolor','none','facealpha',.5);
    delete(h)
end
    title(['Recieved Level at ',region,' ',siteabrev,' N = ',num2str(N)])
    xlim([125,150])
    grid on
    legend(legmat,'Location','Best')
    xlabel('Recieved Level')
    ylabel('Log(Counts)')
    plotName = [saveDir,'\',region,'_',siteabrev,'_StackedRecievedLevel'];
    saveas(gcf,[plotName,'.png'])
    saveas(gcf,[plotName,'.fig'])
    close all
%% save .mat files for future need
save([dataDir,['\',siteabrev,'_simulatedRL.mat']],'spStd','histSubset','spMean','SiteMPP','SiteICI','SiteMTT','rlVec','N','NN','-v7.3')
elseif ReRun == 0 && ~isempty(file)
    disp(['Loading workspace from site ',siteabrev])
    load([dataDir,['\',siteabrev,'_simulatedRL.mat']])
end
%% Do it again but remove the deployments with less data 
if strcmp('PS1',siteabrevv)
    %Remove PS01 and PS13
    SiteMPPx = SiteMPP;
    SiteMPPx(:,6) = [];
    SiteMPPx(:,1) = [];

    siteDiskk = siteDisk;
    siteDiskk(6,:) = [];
    siteDiskk(1,:) = [];
    
    N = round(min(cellfun(@(c) size(c,1), SiteMPPx))/1.3); %come up with a sample size that's 1.3 the size of the min number of clicks
    
    [qq,~] = size(siteDiskk);
    SiteMPP = cell(1,qq);
    SiteMTT = cell(1,qq);
    SiteICI = cell(1,qq);
    histSubset = cell(1,qq);
    spMean = cell(1,qq);
    spStd = cell(1,qq);

for i = 1:length(siteDiskk)
%Find min length of SiteMPP to adjust N accordingly
for itr = 1:50
    p = randsample(length(SiteMPPx{i}),N); %choose N random click samples
    histSubset{i} = [histSubset{i};hist(SiteMPPx{i}(p),rlVec)];
end
end

counterI = 1;
figure
if length(siteDiskk) > 12
    for i = 1:length(siteDiskk)
    subplot(5,5,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDiskk{i}])
    counterI = counterI+1;
    end
elseif length(siteDiskk) > 9 && length(siteDiskk) < 12
for i = 1:length(siteDiskk)
    subplot(3,4,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDiskk{i}])
    counterI = counterI+1;
end
elseif length(siteDiskk) > 4 && length(siteDiskk) < 10
    for i = 1:length(siteDiskk)
    subplot(3,3,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDiskk{i}])
    counterI = counterI+1;
    end
else
for i = 1:length(siteDiskk)
    subplot(2,2,counterI)
    spMean{i} = mean(histSubset{i});
    spStd{i} = std(histSubset{i});
    semilogy(rlVec,histSubset{i})
    hold on
    plot(rlVec,spMean{i},'-k','Linewidth',2)
    plot(rlVec,spMean{i}+spStd{i},'--k','Linewidth',2)
    plot(rlVec,spMean{i}-spStd{i},'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')
    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    plot([125,125],[1,10000],'--r')
    ylim([1,40000])
    xlim([125 150])
    xlim([RLmin,RLmax])
    title(['Recieved Level for ',siteDiskk{i}])
    counterI = counterI+1;
end
end
plotName = [saveDir,'\',region,'_',siteabrev,'_SubPlotRecievedLevel_higherDet'];
saveas(gcf,[plotName,'.png'])
saveas(gcf,[plotName,'.fig'])

%Find matching TF to include in the plot
tfnum = [];
ifoundx = 0;
for itab = 1:stxt(1)
    ifound = strfind(dtable.Data_ID(itab),siteabrevv);
    if cell2mat(ifound) >0
        ifoundx = ifoundx + 1;
        tfnum(ifoundx) = str2double(dtable.PreAmp(itab));
    end
end

%delete PS01 and PS13
tfnumx = tfnum;
tfnumx(:,6:7) = [];
tfnumx(:,1) = [];

c = distinguishable_colors(length(siteDiskk)); %more colors

%legend table
legmat = [];
for iR = 1:length(siteDiskk)
    legmat{iR} = [cell2mat(siteDisk(iR)),' ',num2str(tfnum(iR)),' n = ',num2str(length(SiteMPPx{iR}))];
end

figure
for iR = 1:length(siteDiskk)
    spMean{iR} = mean(histSubset{iR});
    spStd{iR} = std(histSubset{iR});
    y = log10(spMean{iR});
    ystd1 = log10(spMean{iR}+spStd{iR});
    ystd2 = log10(spMean{iR}-spStd{iR});
    badY = isinf(y);
    y(badY)=[];
    ystd1(badY)=[];
    ystd2(badY)=[];
    h = plot(rlVec(~badY),y);
    hold on
    x2 = [rlVec(~badY),fliplr(rlVec(~badY))];
    inBetween = ([ystd1,fliplr(ystd2)]);
    fill(x2, inBetween,c(iR,:),'edgecolor','none','facealpha',.5);
    delete(h)
end
    title(['Recieved Level at ',region,' ',siteabrev,' N = ',num2str(N)])
    xlim([125,150])
    grid on
    legend(legmat,'Location','Best')
    xlabel('Recieved Level')
    ylabel('Log(Counts)')
    plotName = [saveDir,'\',region,'_',siteabrev,'_StackedRecievedLevel_higherDet'];
    saveas(gcf,[plotName,'.png'])
    saveas(gcf,[plotName,'.fig'])
    close all
%% save .mat files for future need
save([dataDir,['\',siteabrev,'_simulatedRL_higherDet.mat']],'spStd','histSubset','spMean','SiteMPPx','siteDiskk','SiteICI','SiteMTT','rlVec','N','NN','-v7.3')
end
    
end