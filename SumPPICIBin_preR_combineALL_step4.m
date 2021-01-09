clearvars
close all
%% Parameters defined by user
filePath = 'G:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\'; %specify directory to save files
%% Find all files that fit your specifications
files = dir([filePath,'**\*_binPresence.csv']);
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
newSiteName = extractBefore(siteName,'_');
gg = height(x{1,i}); %length of table
x{1,i}.Site = repmat(newSiteName,gg,1);
end

%combine all the tables into one
table = vertcat(x{:});
writetable(table,[filePath,'All_Data.csv']); %save table to .csv to continue stats in R F:\Seasonality\Kruskal_RankSumSTATS.R
