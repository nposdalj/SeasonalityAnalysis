close all;clear all;clc;

%%specify file path
%% Parameters defined by user
filePath = 'G:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\'; %specify directory to save files
%% load lat and longs for each site
CB_latLongs = [58.645683,-148.07; %01 in decimals
58.671333,-148.02; %02
58.673483,-148.00; %03
58.671867,-148.02; %04
58.671,-148.02; %05
58.670817,-148.02; %06
58.65525,-148.09; %07
58.6695,-148.03; %08
58.6695,-148.0338889]; %09

QN_latLongs = [56.339017,-145.18; %qn01
56.339383,-145.18; %qn02
56.3413,-145.18; %qn04
56.340683,-145.18; %qn05
56.339833,-145.19]; %qn06

PT_latLongs = [56.24345, -142.75; %pt01
56.243917, -142.75; %pt02
56.242917, -142.75; %pt03
56.243333, -142.75]; %pt04

OCNMSQC_latLongs = [47.466,-125.16; %6
47.50005,-125.36; %12
47.500433,-125.36; %14
47.500533,-125.36; %15
47.500633,-125.65]; %16

ALEUTBD_latLongs = [52.633333,-185.63; %BD02
52.076,-185.64]; %BD03
%% load lat and long with only 1 deployment
AB_mean = [57.513667,-146.50];
ABtext = repmat({'AB'},size(AB_mean,1),1);
AB = [ABtext num2cell(AB_mean)];

ALEUT01KS_mean = [52.316783,-188.52]; 
KStext = repmat({'KS'},size(ALEUT01KS_mean,1),1);
KS = [KStext num2cell(ALEUT01KS_mean)];

KOA_mean = [57.224,-150.53];
KOAtext = repmat({'KOA'},size(KOA_mean,1),1);
KOA = [KOAtext num2cell(KOA_mean)];

%find means of sites with multiple deployments
CB_mean = mean(CB_latLongs);
CBtext = repmat({'CB'},size(CB_mean,1),1);
CB = [CBtext num2cell(CB_mean)];

QN_mean = mean(QN_latLongs);
QNtext = repmat({'QN'},size(QN_mean,1),1);
QN = [QNtext num2cell(QN_mean)];

PT_mean = mean(PT_latLongs);
PTtext = repmat({'PT'},size(PT_mean,1),1);
PT = [PTtext num2cell(PT_mean)];

ALEUTBD_mean = mean(ALEUTBD_latLongs);
BDtext = repmat({'BD'},size(ALEUTBD_mean,1),1);
BD = [BDtext num2cell(ALEUTBD_mean)];

%create one table with all lat and longs
LL = [AB; BD; CB; KS; PT; QN; KOA];
LatLong = cell2mat(LL(:,2:3));

LatLongTAB = array2table(LatLong);
LatLongTAB.Properties.VariableNames = {'Latitude' 'Longitude'};
LatLongTAB{:,'Site'} = {'AB'; 'BD'; 'CB'; 'KS'; 'PT'; 'QN'; 'KOA'};
%LatLongTAB.Longitude(2) = -(LatLongTAB.Longitude(2) + 10);
%% load bin data
% Find all files that fit your specifications for sites with less than a year
files = dir([filePath,'**\*_365GroupedMean.csv']);
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

%% Make a master table with average proportion of hours on that day of the year
columnIndicesToDelete = [1 3 4 5 6 7 8 9];
table_females(:,columnIndicesToDelete) = [];
table_juveniles(:,columnIndicesToDelete) = [];
columnIndicesToDelete = [3 4 5 6];
table_males(:,columnIndicesToDelete) = [];
masterTAB = [table_males table_females table_juveniles];
masterTAB = [masterTAB(:,1:2) masterTAB(:,6:7) masterTAB(:,3:5)];
%add table_short to masterTAB
table_short = [table_short(:,1) table_short(:,4) table_short(:,2:3) table_short(:,5:7)];
masterTAB = [masterTAB;table_short];
%% Take the mean proportion of hours for each season
masterTABmean = varfun(@mean,masterTAB,'InputVariables',{'HoursPropMA',...
    'HoursPropFE','HoursPropJU'},'GroupingVariables',{'Season','Site'});
%% Fill in zeros
masterTABmean.Site = categorical(masterTABmean.Site);
masterTABmean.Season(23) = 2;
masterTABmean.Site(23) = 'AB';
masterTABmean.Season(24) = 2;
masterTABmean.Site(24) = 'KOA';
masterTABmean.Season(25) = 2;
masterTABmean.Site(25) = 'KS';
masterTABmean.Season(26) = 3;
masterTABmean.Site(26) = 'AB';
masterTABmean.Season(27) = 3;
masterTABmean.Site(27) = 'KOA';
masterTABmean.Season(28) = 3;
masterTABmean.Site(28) = 'KS';

%% sort tables alphabetically
LatLongTAB = sortrows(LatLongTAB,'Site');
masterTABmean = sortrows(masterTABmean,'Site');

%% assign bathymetry type to data
masterTABmean.Site = string(masterTABmean.Site);
[q,~]=size(masterTABmean);
masterTABmean.bathy = zeros(q,1);
masterTABmean.bathy(masterTABmean.Site == 'CB') = 1;
masterTABmean.bathy(masterTABmean.Site == 'KOA') = 1;
masterTABmean.bathy(masterTABmean.Site == 'AB') = 2;
masterTABmean.bathy(masterTABmean.Site == 'PT') = 3;
masterTABmean.bathy(masterTABmean.Site == 'QN') = 3;
masterTABmean.bathy(masterTABmean.Site == 'BD') = 4;
masterTABmean.bathy(masterTABmean.Site == 'KS') = 4;

LatLongTAB.Site = string(LatLongTAB.Site);
[q,~]=size(LatLongTAB);
LatLongTAB.bathy = zeros(q,1);
LatLongTAB.bathy(LatLongTAB.Site == 'CB') = 1;
LatLongTAB.bathy(LatLongTAB.Site == 'KOA') = 1;
LatLongTAB.bathy(LatLongTAB.Site == 'AB') = 2;
LatLongTAB.bathy(LatLongTAB.Site == 'PT') = 3;
LatLongTAB.bathy(LatLongTAB.Site == 'QN') = 3;
LatLongTAB.bathy(LatLongTAB.Site == 'BD') = 4;
LatLongTAB.bathy(LatLongTAB.Site == 'KS') = 4;
%% Geobubble with sizes for Female presence - subplot
masterTABmean.bathy = categorical(masterTABmean.bathy);
colors = [1 1 0; 1 0 0; 0.2 0 0; 0 1 0];
% colors = distinguishable_colors(7);
figure(1)
% set(gcf,'DefaultAxesColorOrder',colors)
subplot(2,2,1)
wi = geobubble(LatLongTAB.Latitude,LatLongTAB.Longitude,masterTABmean.mean_HoursPropFE(masterTABmean.Season ==3),masterTABmean.bathy(masterTABmean.Season ==3),'SizeLimits',[0 0.053],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
wi.BubbleColorList = colors;
title 'Winter (Jan-Mar)';
subplot(2,2,2)
sp = geobubble(LatLongTAB.Latitude,LatLongTAB.Longitude,masterTABmean.mean_HoursPropFE(masterTABmean.Season ==4),masterTABmean.bathy(masterTABmean.Season ==4),'SizeLimits',[0 0.053],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
sp.BubbleColorList = colors;
title 'Spring (May-Jun)';
% sp.SizeLegendTitle = 'Proportion of Hours per Day';
subplot(2,2,3)
su = geobubble(LatLongTAB.Latitude,LatLongTAB.Longitude,masterTABmean.mean_HoursPropFE(masterTABmean.Season ==1),masterTABmean.bathy(masterTABmean.Season ==1),'SizeLimits',[0 0.053],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
su.BubbleColorList = colors;
title 'Summer (Jul-Sep)';
subplot(2,2,4)
fa = geobubble(LatLongTAB.Latitude,LatLongTAB.Longitude,masterTABmean.mean_HoursPropFE(masterTABmean.Season ==2),masterTABmean.bathy(masterTABmean.Season ==2),'SizeLimits',[0 0.053],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
fa.BubbleColorList = colors;
title 'Fall (Oct-Dec)';
set(gcf,'Color','w');
save('FemalePresence.png');
export_fig FemalePresenceHQ.png
%% Geobubble with sizes for Female presence - subplot (GOA AND BSAI ONLY)
LatLongTABF_NL = LatLongTABF;
LatLongTABF_NL([4,5,6,7,9],:) = []; 
LatLongTABF_NL.Site = categorical(LatLongTABF_NL.Site);
lat_lims = [50 62];
long_lims = [-180 -120];

%subplots
NLcolors = distinguishable_colors(6);
figure(4)
subplot(2,2,1)
site = LatLongTABF_NL.Site;
wiNL = geobubble(LatLongTABF_NL.Latitude,LatLongTABF_NL.Longitude,LatLongTABF_NL.WinterBin,site,'SizeLimits',[0 50],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
wiNL.BubbleColorList = NLcolors;
title 'Female Social Unit Presence in the Winter (Dec-Feb)';
wiNL.SizeLegendTitle = 'Daily Presence (min)';
subplot(2,2,2)
sp = geobubble(LatLongTABF_NL.Latitude,LatLongTABF_NL.Longitude,LatLongTABF_NL.SpringBin,site,'SizeLimits',[0 50],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
sp.BubbleColorList = NLcolors;
title 'Female Social Unit Presence in the Spring (Mar-May)';
sp.SizeLegendTitle = 'Daily Presence (min)';
subplot(2,2,3)
su = geobubble(LatLongTABF_NL.Latitude,LatLongTABF_NL.Longitude,LatLongTABF_NL.SummerBin,site,'SizeLimits',[0 50],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
su.BubbleColorList = NLcolors;
title 'Female Social Unit Presence in the Summer (Jun-Aug)';
su.SizeLegendTitle = 'Daily Presence (min)';
subplot(2,2,4)
fa = geobubble(LatLongTABF_NL.Latitude,LatLongTABF_NL.Longitude,LatLongTABF_NL.FallBin,LatLongTABF_NL.Site,'SizeLimits',[0 50],'BubbleWidthRange',[1 20]);
% geolimits(lat_lims,long_lims);
fa.BubbleColorList = NLcolors;
title 'Female Social Unit Presence in the Fall (Sep-Nov)';
fa.SizeLegendTitle = 'Daily Presence (min)';
set(gcf,'Color','w');
save('FemalePresenceNL.png');
export_fig FemalePresenceNLHQ.png
%% Geobubble with sizes for Juvenile presence - subplot
LatLongTABJ.Site = categorical(LatLongTABJ.Site);
lat_lims = [28 62];
long_lims = [-180 -119];

%subplots
colors = distinguishable_colors(11);
figure(2)
set(gcf,'DefaultAxesColorOrder',colors)
subplot(2,2,1)
wi = geobubble(LatLongTABJ.Latitude,LatLongTABJ.Longitude,LatLongTABJ.WinterBin,LatLongTABJ.Site,'SizeLimits',[0 60],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
wi.BubbleColorList = colors;
title 'Mid-Size Sperm Whale Presence in the Winter (Dec-Feb)';
wi.SizeLegendTitle = 'Daily Presence (min)';
subplot(2,2,2)
sp = geobubble(LatLongTABJ.Latitude,LatLongTABJ.Longitude,LatLongTABJ.SpringBin,LatLongTABJ.Site,'SizeLimits',[0 60],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
sp.BubbleColorList = colors;
title 'Mid-Size Sperm Whale Presence in the Spring (Mar-May)';
sp.SizeLegendTitle = 'Daily Presence (min)';
subplot(2,2,3)
su = geobubble(LatLongTABJ.Latitude,LatLongTABJ.Longitude,LatLongTABJ.SummerBin,LatLongTABJ.Site,'SizeLimits',[0 60],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
su.BubbleColorList = colors;
title 'Mid-Size Sperm Whale Presence in the Summer (Jun-Aug)';
su.SizeLegendTitle = 'Daily Presence (min)';
subplot(2,2,4)
fa = geobubble(LatLongTABJ.Latitude,LatLongTABJ.Longitude,LatLongTABJ.FallBin,LatLongTABJ.Site,'SizeLimits',[0 60],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
fa.BubbleColorList = colors;
title 'Mid-Size Sperm Whale Presence in the Fall (Sep-Nov)';
fa.SizeLegendTitle = 'Daily Presence (min)';
set(gcf,'Color','w');
save('JuvenilePresence.png');
export_fig JuvenilePresenceHQ.png
%% Geobubble with sizes for presence
LatLongTABM.Site = categorical(LatLongTABM.Site);
lat_lims = [28 62];
long_lims = [-180 -119];

%subplots
colors = distinguishable_colors(11);
figure(3)
subplot(2,2,1)
wi = geobubble(LatLongTABM.Latitude,LatLongTABM.Longitude,LatLongTABM.WinterBin,LatLongTABM.Site,'SizeLimits',[0 70],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
wi.BubbleColorList = colors;
title 'Male Presence in the Winter (Dec-Feb)';
wi.SizeLegendTitle = 'Daily Presence (min)';
subplot(2,2,2)
sp = geobubble(LatLongTABM.Latitude,LatLongTABM.Longitude,LatLongTABM.SpringBin,LatLongTABM.Site,'SizeLimits',[0 70],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
title 'Male Presence in the Spring (Mar-May)';
sp.SizeLegendTitle = 'Daily Presence (min)';
sp.BubbleColorList = colors;
subplot(2,2,3)
su = geobubble(LatLongTABM.Latitude,LatLongTABM.Longitude,LatLongTABM.SummerBin,LatLongTABM.Site,'SizeLimits',[0 70],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
title 'Male Presence in the Summer (Jun-Aug)';
su.SizeLegendTitle = 'Daily Presence (min)';
su.BubbleColorList = colors;
subplot(2,2,4)
fa = geobubble(LatLongTABM.Latitude,LatLongTABM.Longitude,LatLongTABM.FallBin,LatLongTABM.Site,'SizeLimits',[0 70],'BubbleWidthRange',[1 20]);
geolimits(lat_lims,long_lims);
title 'Male Presence in the Fall (Sep-Nov)';
fa.SizeLegendTitle = 'Daily Presence (min)';
fa.BubbleColorList = colors;
set(gcf,'Color','w');
save('MalePresence.png');
export_fig MalePresenceHQ.png
