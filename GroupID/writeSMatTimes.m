function writeSMatTimes(cDir,saveDir)

% get deployment names
folders = dir(cDir) ;
folders = ({folders.name})' ;
folders = folders(3:end);

siteDiskList = cell(length(folders),1);
for s = 1: length(folders)
    a = cell2mat(strfind(folders(s),'_disk'));
    siteDiskList{s,1} = folders{s}(1:a-1);
end

siteDisk = unique(siteDiskList);

for i = 1:length(siteDisk)
    disp(['loading times from: ', siteDisk{i}]);
    index = strfind(siteDiskList, siteDisk{i});
    siteDiskIdx = find(not(cellfun('isempty', index)));
        % return detection times from the same site
    times = [];
    allfolders = [];
    for j = 1:length(siteDiskIdx)
        f = siteDiskIdx(j);
        fold = fullfile(cDir,folders{f});
        disp(['   ', folders{f}]);
        SearchFileMaskShip = {'*.mat'};
        SearchPathMaskShip = {[cDir,'\',siteDisk{i},'*']};
        SearchPathMask = {fold};
        SearchRecursiv = 1;
[PathFileListShip, FileListShip, PathListShip] = ...
    utFindFiles(SearchFileMaskShip, SearchPathMaskShip, SearchRecursiv);

if length(FileListShip) == length(siteDiskIdx)
    
    MatchIDX = strfind(PathListShip, siteDiskIdx(j));
    
if ~isempty(MatchIDX)
    MatchIndex = find(not(cellfun('isempty', index)));
    PathFileListShipMatch = PathFileListShip;
        
if length(siteDiskIdx) == 1
   load(PathFileListShipMatch{1});
   shipTAB = [table(shipLabels) table(num2cell(shipTimes))];
   shipArray = table2array(shipTAB);
   shipArray(strcmp(shipArray(:,1),'ambient'),:) = [];
   shipArray(:,1) = [];
else
    if j == 1
    shipTAB = [];
    load(PathFileListShipMatch{1});
    shipTAB = [table(shipLabels) table(num2cell(shipTimes))];
        else
    load(PathFileListShipMatch{j});
    ships = [table(shipLabels) table(num2cell(shipTimes))];
    shipTAB = [shipTAB ; ships];  
    end
end
else
end
else
    [PathFileListShipa, FileListShipa, PathListShipa] = ...
    utFindFiles(SearchFileMaskShip, SearchPathMask, SearchRecursiv);
if ~isempty(PathFileListShipa)

if length(siteDiskIdx) == 1
   load(char(PathFileListShipa));
   shipTAB = [table(shipLabels) table(num2cell(shipTimes))];
   shipArray = table2array(shipTAB);
   shipArray(strcmp(shipArray(:,1),'ambient'),:) = [];
   shipArray(:,1) = [];
else
    if j == 1
    shipTAB = [];
    load(char(PathFileListShipa));
    shipTAB = [table(shipLabels) table(num2cell(shipTimes))];
%     elseif i == 9 && j == 2
%     shipTAB = [];
%     load(char(PathFileListShipa));
%     shipTAB = [table(shipLabels) table(num2cell(shipTimes))];
    else            
    load(char(PathFileListShipa));
    ships = [table(shipLabels) table(num2cell(shipTimes))];
    shipTAB = [shipTAB ; ships];  
    end
end
else
end

end
end

    if exist('shipTAB','var')
shipArray = table2array(shipTAB);
shipArray(strcmp(shipArray(:,1),'ambient'),:) = [];
shipArray(:,1) = [];

shipStart = cell2mat(shipArray(:,1));
shipEnd = cell2mat(shipArray(:,2));
times = [shipStart shipEnd];
% write times in mat file
save([saveDir,'\',char(siteDisk{i}),'.mat'],'times');
clear shipTAB
    else
    end

end

disp('Done writing ship file times')
