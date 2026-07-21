clearvars
close all
%% Parameters defined by user
filePaths = {'F:\AZORES\SeasonalityAnalysis\AZORES_B_01\', ...
             'F:\AZORES\SeasonalityAnalysis\AZORES_A_04_DEEP\'}; %specify directories to search for output CSVs - add more roots here as needed
saveDir = filePaths{1}; %combined output is saved alongside the first listed directory
%% Find all files that fit your specifications, across all listed directories
files = [];
for iPath = 1:length(filePaths)
    theseFiles = dir([filePaths{iPath},'**\*_dayData_forGLMR125.csv']);
    files = [files; theseFiles];
end
n = length(files);
x = cell(1, numel(files));
%load all of the tables
for i=1:n
    fn = fullfile({files(i).folder},{files(i).name});
    fn_char = char(fn);
    x{i}=readtable(fn_char);
end

%add a new column for each table with the site name
for i=1:n
siteName = {files(i).name};
newSiteName = extractBefore(siteName,'_dayData_forGLMR125.csv'); % strip the known suffix instead of splitting on the
% first underscore, since site names like AZORES_B_01 / AZORES_A_04_DEEP contain underscores themselves
gg = height(x{1,i}); %length of table
x{1,i}.Site = repmat(newSiteName,gg,1);
end

%combine all the tables into one
table = vertcat(x{:});
writetable(table,[saveDir,'All_Data.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
