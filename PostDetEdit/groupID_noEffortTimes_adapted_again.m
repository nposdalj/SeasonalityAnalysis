%% Clear all variables
close all;clear all;clc;
%% select directory where ship files are located
siteabrev = 'Palmyra';
siteNameMatch = 'Palmyra';
region = 'CentralPac'; %all of the WAT data has a space between the site and the deployment #
shipDataType = 1; % 1 - old ship data, %2 - new ship data
ShipIDReDo = 0; % If you want to re-run ship and ID times, change this to 1
GDrive = 'I';

%shipDir = 'I:\Shared drives\Pt. Sur\Analyzed data\Sperm whales\metadata_reduced';
%shipTimesDir = 'I:\Shared drives\Pt. Sur\Analyzed data\Sperm whales\ShipTimes';
shipDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\metadata_reduced\',siteabrev];
shipTimesDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\ShipTimes\',siteabrev]; % directory where to save ship times .mat files
%shipDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\',siteabrev,'\metadata_reduced'];
%shipTimesDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\',siteabrev,'\ShipTimes']; % directory where to save ship times .mat files
%shipDir = 'E:\WC\metadata_reduced';
%shipTimesDir = 'E:\WC\ShipTimes';

% IDDir = 'I:\Shared drives\Pt. Sur\Analyzed data\Sperm whales\TPWS_125';
% IDTimesDir = 'I:\Shared drives\Pt. Sur\Analyzed data\Sperm whales\IDTimes';
IDDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\TPWS_125\',siteabrev];
%IDDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\',siteabrev,'\TPWS_125'];
%IDDir = 'E:\WC\TPWS_125';
IDTimesDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\IDTimes\',siteabrev]; % directory where to save ID times .mat files
% IDTimesDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\',siteabrev,'\IDTimes']; % directory where to save ID times .mat files
%IDTimesDir = 'E:\WC\IDTimes';

maxDetEdit = 2; % number of TPWS folders (i.e. TPWS4 is 4)
%saveTable = 'E:\WC\SeasonalityAnalysis\Pm_Effort.xlsx';
%saveTable = 'I:\Shared drives\Pt. Sur\Analyzed data\Sperm whales\SeasonalityAnalysis\Pm_Effort.xlsx';
saveTable = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx'];
%saveTable = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\',siteabrev,'\Seasonality\Pm_Effort.xlsx'];
%% write ship and ID file times
%Check to see if ship files are in the old format or new format
if ShipIDReDo == 1
if shipDataType == 1
    writeSFilesTimes(shipDir,shipTimesDir); %for old ship files
else
    writeSMatTimes(shipDir,shipTimesDir); %for new ship files
end
%run join_IDs to group no effort times
writeIDTimes(IDDir,IDTimesDir,maxDetEdit);
else
end
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
    if strcmp(site(n),siteNameMatch);
        loc = site{n};
        strdepl = num2str(depl(n),'%02d');
        if strcmp(siteNameMatch,'ALEUT01KS');
            deplCompare = siteNameMatch;
        elseif strcmp(region,'WAT')
            deplCompare = [loc,'_',strdepl];
        elseif strcmp(region,'JAX')
            deplCompare = [loc,'_D_',strdepl];
        elseif strcmp(siteNameMatch,'CORC');
            deplCompare = ['OTSG_',siteNameMatch,'4_',strdepl];
        elseif strcmp(siteNameMatch,'DCPP01C');
            deplCompare = siteNameMatch;
        elseif strcmp(siteNameMatch,'HOKE');
            deplCompare = siteNameMatch;
        else
            deplCompare = [loc strdepl];
        end
    
        
    % load ship and ID times
    iS = contains(shipfiles,deplCompare);
    iID = contains(IDfiles,deplCompare);
    
     if any([iS;iID]) % only if there is ship and ID files
        if any(iS)
            ship = load(fullfile(shipTimesDir,shipfiles{iS}));
        else
            ship.times = [];
        end

        if any(iID)
            ID = load(fullfile(IDTimesDir,IDfiles{iID}));
        else 
            ID.times = []; 
        end
     else
         ship.times = [];
         ID.times = [];
     end
     
      times = [ship.times; ID.times];
      times = groupoverlaps(times);
   
      % get times between no effort intervals
      effort = [];
      effort(:,1) = [edgeffort(n,1); times(:,2) + datenum(0,0,0,0,0,1/1000)]; % add a ms
      effort(:,2) = [times(:,1) - datenum(0,0,0,0,0,1/1000); edgeffort(n,2)]; % extract a ms
      
      [p,q] = size(effort);
      
      S = cell(p,1);
      LOC = cell(p,1);
      lat = cell(p,1);
      log = cell(p,1);
        
      S(:) = {loc};
      LOC(:) = {strdepl};
      lat(:) = {latLongs(n,1)};
      log(:) = {latLongs(n,2)};
        
      effTable = [effTable; table(S,LOC,lat,log,m2xdate(effort(:,1)),m2xdate(effort(:,2)))];
    else
    end
end

effTable.Properties.VariableNames = {'Sites','Deployments','Latitude','Longitude','StartEffort','EndEffort'};

% save effort times in excel file
writetable(effTable,saveTable)

disp('Effort times saved')
