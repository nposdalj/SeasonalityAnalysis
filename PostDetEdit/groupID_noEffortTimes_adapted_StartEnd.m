close all;clear all;clc;
%% select directory where ship files are located
site = 'WC';
region = 'WAT'; %all of the WAT data has a space between the site and the deployment
saveTable = ['E:\',site,'\SeasonalityAnalysis\Effort.xlsx'];
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
