% randSampleMPP

% Note: IDAll file has all click times, and all pp amplitudes for one
% deployment, could just be [MTT, MPP], merged cumulatively across all
% disks.

% This code was provided by KEF and modified by NP 08312022
close all;clear all clc;
%% User Definied Parameters
siteabrevv = {'NC','WC','OC','BC','HZ','BS','BP','GS','JAX'}; %abbreviation of site.
GDrive = 'I'; %directory for Google Drive
region = 'WAT';
sp = 'Pm'; % your species code
itnum = '3'; % which iteration you are looking for (which TPWS folder)
RLmin = 125; %Recieved level threshold of your data
RLmax = 170; %Max recieved level you're intersted in plotting
%% Calculate variable sizes
[~,gg] = size(siteabrevv);
rlVec = RLmin:1:RLmax;
%% Load workspaces
for ss = 1:gg
siteabrev = siteabrevv{ss}; %choose site
disp(['Loading data from site ',siteabrev])
tpwsPath = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\TPWS_125\',siteabrev]; %directory of TPWS files
dataDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory where workspaces are saved
effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx'];% specify excel file with effort times
%load([dataDir,'\',siteabrev,'_workspaceStep2.mat']);
%load([dataDir,'\',siteabrev,'_workspace125.mat']); %load this in last to properly load the variable p
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\Plots\',siteabrev]; % add desired saveDir back in here

%% define subfolder that fit specified iteration
if itnum > 1
   for id = 2: str2num(itnum) % iterate id times according to itnum
       subfolder = ['TPWS',num2str(id)];
       tpwsPath = (fullfile(tpwsPath,subfolder));
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
    % if no matches, throw error
    error('No files matching filePrefix found!')
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
for i = 1:length(siteDisk)
    disp(['loading times from: ', siteDisk{i}]);
    index = strfind(siteDiskList, siteDisk{i});
    siteDiskIdx = find(not(cellfun('isempty', index)));
%% Concatenate variables
PPall = []; TTall = []; ICIall = []; PeakFrall = [];% initialize matrices
for idsk = 1 : length(siteDiskIdx)
    % Load file
    fprintf('Loading %d/%d file %s\n',idsk,length(siteDiskIdx),fullfile(tpwsPath,concatFiles{idsk}))
    D = load(fullfile(tpwsPath,concatFiles{idsk}));
    
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
end
end
%% After parfor data may not be sorted. Sort all the variables.
[~,sorted] = sort(TTall);
TTall = TTall(sorted);
PPall = PPall(sorted);
ICIall = ICIall(sorted);

%% Create empty variables
histSubset = {};
counterI = 1;
spMean = [];
spStd = [];
nameStr = {};
figure(9);clf

for itr = 1:50
     try
        isempty(histSubset{iD,iS})
     catch
         histSubset{iD,iS} = [];           
     end
     p = randsample(length(thisSet),100000);
     histSubset{iD,iS} = [histSubset{iD,iS};hist(IDAll(thisSet(p),3),rlVec)];
end

end



for iD = 1:length(IDList)
    [~,dirParts] = fileparts(IDList(iD).folder);
    load(fullfile(IDList(iD).folder,IDList(iD).name))
    if isempty(IDAll)
        continue
    end
    for iS = 9% 1:length(mySpID)
        
        
        thisSet = find(IDAll(:,2) >=9);
        for itr = 1:50
            try
                isempty(histSubset{iD,iS})
            catch
                
                histSubset{iD,iS} = [];
                
            end
            p = randsample(length(thisSet),100000);
            histSubset{iD,iS} = [histSubset{iD,iS};hist(IDAll(thisSet(p),3),rlVec)];
            
        end
    end
    
    subplot(4,4,counterI)
    spMean(counterI,:) = (mean(histSubset{iD,iS}));
    spStd(counterI,:) = (std(histSubset{iD,iS}));
    semilogy(rlVec,histSubset{iD,iS})

    hold on
   
    plot(rlVec,spMean(counterI,:),'-k','Linewidth',2)
    plot(rlVec,spMean(counterI,:)+spStd(counterI,:),'--k','Linewidth',2)
    plot(rlVec,spMean(counterI,:)-spStd(counterI,:),'--k','Linewidth',2)
    grid on
    ylabel('log10(counts)')

    xlabel('RL(dB_P_P)')
    plot([125,125],[1,10000],'--r')
    ylim([1,3000])
    xlim([120,140])
    nameTry = strfind(dirParts,'GOM');
    if isempty(nameTry)
        nameTry = strfind(dirParts,'GofMX');
        nameStr{counterI} = dirParts(7:10);
    else
        nameStr{counterI} = strrep(dirParts(5:9),'_','-');
    end
    title(nameStr{counterI})   
    counterI = counterI+1;

end

figure(19);clf
for iR = 1:size(spMean,1)
    y = log10(spMean(iR,:));
    ystd1 = log10(spMean(iR,:)+spStd(iR,:));
    ystd2 = log10(spMean(iR,:)-spStd(iR,:));
    badY = isinf(y);
    y(badY)=[];
    ystd1(badY)=[];
    
    ystd2(badY)=[];
    h = plot(rlVec(~badY),y);
    %ylim([1,1000])
    hold on
    x2 = [rlVec(~badY),fliplr(rlVec(~badY))];
    inBetween = ([ystd1,fliplr(ystd2)]);
    
    fill(x2, inBetween,h.Color,'edgecolor','none','facealpha',.5);
    delete(h)
end


xlim([120,140])
grid on
plot(rlVec,log10(mean(spMean)),'-k','linewidth',2)
ylim([0,3.5])
[~,cIdx] = min(abs(rlVec-125));
meanCross = log10(mean(spMean(:,cIdx)));
plot([125,125],[1,10000],'--r')
plot(rlVec,ones(size(rlVec))*meanCross,'--r')

[C,I] = min(abs(spMean - meanCross),[],2);
legend(cellstr(nameStr))