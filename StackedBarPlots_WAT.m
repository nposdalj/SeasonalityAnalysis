clearvars
close all

%this code should be run in MATLAB 2019b or later
%making a change here
%% Parameters defined by user
sitename = 'BP'; %specify site; must rerun thru line 96 for each site. After each run comment the site (lines 128-132) to avoid replacing
filePath = ['G:\.shortcut-targets-by-id\1FGSX39xqOmreo9qPfPoqhlhUNm1STQB9\WAT_TPWS_metadataReduced\SeasonalityAnalysis\', sitename]; %specify directory to save files
%% Find all files that fit your specifications for sites with less than a year
files = dir([filePath,'\',sitename,'_days365GroupedMean_forGLMR125.csv']);
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
files = dir([filePath,'\',sitename,'_365GroupedMeanFemale.csv']);
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
files = dir([filePath,'\',sitename,'_365GroupedMeanJuvenile.csv']);
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
files = dir([filePath,'\',sitename,'_365GroupedMeanMale.csv']);
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
gray = [0.6 0.6 0.6];
mint = [0.4000 0.7608 0.6471];       % for social units
persimmon = [0.9882 0.5529 0.3843];  % for mid-size animals
slate = [0.5529 0.6275 0.7961];      % for males
tilecolor = [slate; persimmon; mint; gray];
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
%xNC = [sum(masterTAB.Juv(masterTAB.Site=='NC'))  sum(masterTAB.Mal(masterTAB.Site=='NC')) sum(masterTAB.Fem(masterTAB.Site=='NC')) sum(masterTAB.NA(masterTAB.Site=='NC'))];
%xBC = [sum(masterTAB.Juv(masterTAB.Site=='BC'))  sum(masterTAB.Mal(masterTAB.Site=='BC')) sum(masterTAB.Fem(masterTAB.Site=='BC')) sum(masterTAB.NA(masterTAB.Site=='BC'))];
%xGS = [sum(masterTAB.Juv(masterTAB.Site=='GS'))  sum(masterTAB.Mal(masterTAB.Site=='GS')) sum(masterTAB.Fem(masterTAB.Site=='GS')) sum(masterTAB.NA(masterTAB.Site=='GS'))];
%xBP = [sum(masterTAB.Juv(masterTAB.Site=='BP'))  sum(masterTAB.Mal(masterTAB.Site=='BP')) sum(masterTAB.Fem(masterTAB.Site=='BP')) sum(masterTAB.NA(masterTAB.Site=='BP'))];
%xBS = [sum(masterTAB.Juv(masterTAB.Site=='BS'))  sum(masterTAB.Mal(masterTAB.Site=='BS')) sum(masterTAB.Fem(masterTAB.Site=='BS')) sum(masterTAB.NA(masterTAB.Site=='BS'))];
%xMASTER = [xNC; xBC; xGS; xBP; xBS];

b = bar(xMASTER, 'stacked','FaceColor','flat');
b(1).FaceColor = persimmon;  %mid-size
b(2).FaceColor = slate;      %males
b(3).FaceColor = mint;    %social units
b(4).FaceColor = gray;      %no animals
xlabel('Site')
ylabel('Number of Recording Days')
title('Proportion of Days with Presence of Each Group')
legend('Mid-Size','Males','Social Units','No Animals')
ax = gca;
ax.FontSize = 11;
set(gca, 'XTickLabel',{'NC','BC','GS','BP','BS'});

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
