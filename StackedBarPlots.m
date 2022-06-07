clearvars
close all

%this code should be run in MATLAB 2019b or later
%making a change here
%% Parameters defined by user
filePath = 'G:\My Drive\CCE_TPWS_metadataReduced\Baja_GI\Seasonality\'; %specify directory to save files
%% Find all files that fit your specifications for sites with less than a year
files = dir([filePath,'**\*days365GroupedMean_forGLMR125.csv']);
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
table_short = vertcat(x{:});
%% Find all files that fit your specifications for sites with more than a year
%Females
files = dir([filePath,'**\*_365GroupedMeanFemale.csv']);
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
table_females = vertcat(x{:});

%Juveniles
files = dir([filePath,'**\*_365GroupedMeanJuvenile.csv']);
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
table_juveniles = vertcat(x{:});

%Males
files = dir([filePath,'**\*_365GroupedMeanMale.csv']);
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
table_males = vertcat(x{:});
%% define stacked bar plot colors
blue = [0.2081    0.1663    0.5292];
cyan = [0.0383, 0.6742, 0.7435];
yellow = [0.9763    0.9831    0.0538];
grey = [0.6 0.6 0.6 ];
tilecolor = [blue; cyan; yellow; grey];
%% Stacked bar plot for all sites
%find proportion for sites with less than a year of data 
table_short.Fem = table_short.HoursPropFE > 0;
table_short.Juv = table_short.HoursPropJU > 0;
table_short.Mal = table_short.HoursPropMA > 0;
table_short.NA = table_short.HoursPropFE == 0 & table_short.HoursPropJU == 0 & table_short.HoursPropMA == 0;
table_short.Site = string(table_short.Site);
%table_short.Males = sum(table_short.Juv > 0 & table_short.Mal > 0);
xKOA = [sum(table_short.Juv(table_short.Site=='KOA'))  sum(table_short.Mal(table_short.Site=='KOA')) sum(table_short.Fem(table_short.Site=='KOA')) sum(table_short.NA(table_short.Site=='KOA'))];
xKS = [sum(table_short.Juv(table_short.Site=='KS'))  sum(table_short.Mal(table_short.Site=='KS')) sum(table_short.Fem(table_short.Site=='KS')) sum(table_short.NA(table_short.Site=='KS'))];
xAB = [sum(table_short.Juv(table_short.Site=='AB'))  sum(table_short.Mal(table_short.Site=='AB')) sum(table_short.Fem(table_short.Site=='AB')) sum(table_short.NA(table_short.Site=='AB'))];

%find proportion for sites with more than a year of data
columnIndicesToDelete = [1 3 4 5 6 7 8 9];
table_females(:,columnIndicesToDelete) = [];
table_juveniles(:,columnIndicesToDelete) = [];
columnIndicesToDelete = [3 4 5 6 7 8];
table_males(:,columnIndicesToDelete) = [];
masterTAB = [table_males table_females table_juveniles];
masterTAB = [masterTAB(:,1) masterTAB(:,3:5) masterTAB(:,2)];
masterTAB.Fem = masterTAB.HoursPropFE > 0;
masterTAB.Juv = masterTAB.HoursPropJU > 0;
masterTAB.Mal = masterTAB.HoursPropMA > 0;
masterTAB.NA = masterTAB.HoursPropFE == 0 & masterTAB.HoursPropJU == 0 & masterTAB.HoursPropMA == 0;
masterTAB.Site = string(masterTAB.Site);
xCB = [sum(masterTAB.Juv(masterTAB.Site=='CB'))  sum(masterTAB.Mal(masterTAB.Site=='CB')) sum(masterTAB.Fem(masterTAB.Site=='CB')) sum(masterTAB.NA(masterTAB.Site=='CB'))];
xPT = [sum(masterTAB.Juv(masterTAB.Site=='PT'))  sum(masterTAB.Mal(masterTAB.Site=='PT')) sum(masterTAB.Fem(masterTAB.Site=='PT')) sum(masterTAB.NA(masterTAB.Site=='PT'))];
xQN = [sum(masterTAB.Juv(masterTAB.Site=='QN'))  sum(masterTAB.Mal(masterTAB.Site=='QN')) sum(masterTAB.Fem(masterTAB.Site=='QN')) sum(masterTAB.NA(masterTAB.Site=='QN'))];
xBD = [sum(masterTAB.Juv(masterTAB.Site=='BD'))  sum(masterTAB.Mal(masterTAB.Site=='BD')) sum(masterTAB.Fem(masterTAB.Site=='BD')) sum(masterTAB.NA(masterTAB.Site=='BD'))];
xMASTER = [xCB; xBD; xQN; xPT; xKOA; xAB; xKS];

b = bar(xMASTER, 'stacked','FaceColor','flat');
b(1).CData = cyan;
b(2).CData = blue;
b(3).CData = yellow;
b(4).CData = grey;
xlabel('Site')
ylabel('Number of Recording Days')
title('Proprtion of Days with Presence of Each Group')
legend('Mid-Size','Males','Social Units','No Animals')
ax = gca;
ax.FontSize = 16;

%% Stacked bar plot for all sites - edited for Ally for Navy presentation
xMASTER2 = [xCB; xQN; xPT; xKOA; xAB];

b = bar(xMASTER2, 'stacked','FaceColor','flat');
b(1).CData = cyan;
b(2).CData = blue;
b(3).CData = yellow;
b(4).CData = grey;
xlabel('Site')
ylabel('Number of Recording Days')
title({'Proprtion of Days with Presence of', 'Sex Class in The Gulf of Alaska'})
legend('Mid-Size','Males','Social Units','No Animals')
ax = gca;
ax.FontSize = 16;
set(gca, 'XTickLabel',{'CB','QN','PT','KOA','AB'});
