
close all;clear all;clc;

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
[CBlat,CBlong] = meanm(CB_latLongs(:,1),CB_latLongs(:,2));
CB_mean = [CBlat,CBlong];
CBtext = repmat({'CB'},size(CB_mean,1),1);
CB = [CBtext num2cell(CB_mean)];

[QNlat,QNlong] = meanm(QN_latLongs(:,1),QN_latLongs(:,2));
QN_mean = [QNlat, QNlong];
QNtext = repmat({'QN'},size(QN_mean,1),1);
QN = [QNtext num2cell(QN_mean)];

[PTlat,PTlong] = meanm(PT_latLongs(:,1),PT_latLongs(:,2));
PT_mean = [PTlat, PTlong];
PTtext = repmat({'PT'},size(PT_mean,1),1);
PT = [PTtext num2cell(PT_mean)];

[BDlat, BDlong] = meanm(ALEUTBD_latLongs(:,1),ALEUTBD_latLongs(:,2));
ALEUTBD_mean = [BDlat, BDlong];
BDtext = repmat({'BD'},size(ALEUTBD_mean,1),1);
BD = [BDtext num2cell(ALEUTBD_mean)];

%create one table with all lat and longs
LL = [AB; BD; CB; KS; PT; QN; KOA];
LatLong = cell2mat(LL(:,2:3));

LatLongTAB = array2table(LatLong);
LatLongTAB.Properties.VariableNames = {'Latitude' 'Longitude'};
LatLongTAB{:,'Site'} = {'AB'; 'BD'; 'CB'; 'KS'; 'PT'; 'QN'; 'KOA'};
LatLongTAB.Longitude(2) = -(LatLongTAB.Longitude(2) + 10);

lat_lims = [45 65];
long_lims = [-200 -120];
%% grey site map with no color distinction
figure(1)
LatLongTAB.Site = categorical(LatLongTAB.Site);
A = 500;
latitude = LatLongTAB.Latitude;
longitude = LatLongTAB.Longitude;
sites = LatLongTAB.Site;
gm = geoscatter(latitude,longitude,A,sites,'.');

geolimits([45 65],[-220 -120]);
text(latitude(1),longitude(1)-0.5,'AB','HorizontalAlignment','right','FontSize',10);
text(latitude(2)+0.75,longitude(2),'BD','HorizontalAlignment','right','FontSize',10);
text(latitude(3),longitude(3)-0.5,'CB','HorizontalAlignment','right','FontSize',10);
text(latitude(4)+.4,longitude(4)-0.5,'KS','HorizontalAlignment','right','FontSize',10);
text(latitude(5),longitude(5)+2.5,'PT','HorizontalAlignment','right','FontSize',10);
text(latitude(6),longitude(6)-0.75,'QN','HorizontalAlignment','right','FontSize',10);
text(latitude(7),longitude(7)-0.75,'KOA','HorizontalAlignment','right','FontSize',10);
geolimits([45 65],[-200 -120]);

%% grey site map
figure(2)
LatLongTAB.Site = categorical(LatLongTAB.Site);
A = 50;
latitude = LatLongTAB.Latitude;
longitude = LatLongTAB.Longitude;

gm = geoscatter(latitude,longitude,A,'.','k');  
geolimits([-3 74],[-180 -110]);