close all;clear all;clc;
%% select directory where ship files are located
siteabrev = 'GS';
siteNameMatch = 'GS';
shipDataType = 2; % 1 - old ship data, %2 - new ship data

shipDir = ['I:\My Drive\WAT_TPWS_metadataReduced\metadata_reduced\',siteabrev];
shipTimesDir = ['I:\My Drive\WAT_TPWS_metadataReduced\shipTimes\',siteabrev,]; % directory where to save ship times .mat files

IDDir = ['I:\My Drive\WAT_TPWS_metadataReduced\TPWS_125\',siteabrev];
IDTimesDir = ['I:\My Drive\WAT_TPWS_metadataReduced\IDTimes\',siteabrev]; % directory where to save ID times .mat files
maxDetEdit = 2; % number of TPWS folders (i.e. TPWS4 is 4)
saveTable = ['I:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx']; % directory where to save excel file with effort times
%% write ship and ID file times
%Check to see if ship files are in the old format or new format
if shipDataType == 1
    writeSFilesTimes(shipDir,shipTimesDir); %for old ship files
else
    writeSMatTimes(shipDir,shipTimesDir); %for new ship files
end
%run join_IDs to group no effort times
writeIDTimes(IDDir,IDTimesDir,maxDetEdit);
%% Get a list of all the files in the start directory
shipList = cellstr(ls(shipTimesDir));
shipfiles = shipList(3:end); % exclude dots

IDList = cellstr(ls(IDTimesDir));
IDfiles = IDList(3:end); % exclude dots
%% get start end dates of disks
[edgeffort,latLongs, depl, site] = NP_dates;
strdepl = num2str(depl,'%02d');
locations = strcat(site,strdepl);
%% Extract effort times
effTable = table();
for n = 1:length(site)
    if strcmp(site(n),siteNameMatch)
        loc = site{n};
        strdepl = num2str(depl(n),'%02d');
        deplCompare = [loc strdepl];
    % load ship times
    iS = contains(shipfiles,deplCompare);
    iID = contains(IDfiles,deplCompare);
     if any([iS;iID]) % only if there is ship and ID files
        if any(iS)
            %if length(iS) == 1
            ship = load(fullfile(shipTimesDir,shipfiles{iS}));
            %shiptimes = ship.times;
            %else
                %for jj = 1:length(iS)
                    %if jj == 1
                    %ship = load(fullfile(shipTimesDir,shipfiles{jj}));
                    %shiptimes = ship.times;
                    %else
                       %shipTEMP = load(fullfile(shipTimesDir,shipfiles{jj})); 
                       %shipTEMPtimes = shipTEMP.times;
                       %shiptimes = [shiptimes;shipTEMPtimes];
                    %end
        end
        else
            ship.times = [];
    end
end
        
        if any(iID)
%             if length(iID) == 1
            ID = load(fullfile(IDTimesDir,IDfiles{iID}));
%             IDtimes = ID.times;
%             else
%                 for jj = 1:length(iID)
%                     if jj == 1
%                         ID = load(fullfile(IDTimesDir,IDfiles{jj}));
%                         IDtimes = ID.times;
%                     else
%                         IDTEMPtimes = load(fullfile(IDTimesDir,IDfiles{jj}));
%                         IDTEMPtimes = IDTEMPtimes.times;
%                         IDtimes = [IDtimes;IDTEMPtimes];
%                     end
%                 end
            end
        else 
            ID.times = []; 
        end
      times = [ship.times; ID.times];
      times = groupoverlaps(times);
   
      % get times between no effort intervals
      effort = [];
      effort(:,1) = [edgeffort(n,1); times(:,2) + datenum(0,0,0,0,0,1/1000)]; % add a ms
      effort(:,2) = [times(:,1) - datenum(0,0,0,0,0,1/1000); edgeffort(n,2)]; % extract a ms
        
      S = cell(length(effort),1);
      LOC = cell(length(effort),1);
      lat = cell(length(effort),1);
      log = cell(length(effort),1);
        
      S(:) = {loc};
      LOC(:) = {strdepl};
      lat(:) = {latLongs(n,1)};
      log(:) = {latLongs(n,2)};
        
      effTable = [effTable; table(S,LOC,lat,log,m2xdate(effort(:,1)),m2xdate(effort(:,2)))];
    end
    else
    end
end

effTable.Properties.VariableNames = {'Sites','Deployments','Latitude','Longitude','StartEffort','EndEffort'};

% save effort times in excel file
writetable(effTable,saveTable)

disp('Effort times saved')
