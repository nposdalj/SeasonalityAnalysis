
close all;clear all;clc;

%% load lat and longs for each site
NC_latLongs = [39.83248333,-69.98; %01 in decimals
39.83238333, -69.98; %nc02
39.83258333, -69.98; %nc03
39.83295, -69.98]; %nc04

BC_latLongs = [39.19105, -72.23; %bc01
39.1905, -72.23; %bc02
39.19191667, -72.23]; %bc03

GS_latLongs = [33.66563333, -76.00; %gs01
33.66701667, -76.00; %gs02
33.66991667, -76.00]; %gs03
%gs04???

BP_latLongs = [32.10603333, -77.09; %bp01
32.10695, -77.09; %bp02
32.10526667, -77.09]; %bp03
%bp04???

BS_latLongs = [30.58378333, -77.39; %bs01
30.58303333, -77.39; %bs02
30.58295, -77.39]; %bs03

WC_latLongs = [38.37415, -73.37; %wc01
38.37385, -73.37; %wc02
38.37336667, -73.37]; %wc03

OC_latLongs = [40.2633, -67.99; %oc01
40.26331667, -67.99; %oc02
40.26333333, -67.99; %oc03
40.23, -67.98]; %oc04

HZ_latLongs = [41.06191667, -66.35; %hz01
41.06183333, -66.35; %hz02
41.06165, -66.35; %hz03
41.06165, -66.35]; %hz04

JAX_latLongs = [30.15183333, -79.77; %JAX_D_13
30.15268333, -79.77; %JAX_D_14
30.15225, -79.77]; %JAX_D_15
%Not sure if other JAX deployments will be included but this is it for now,
%update if needed.

%% load lat and long with only 1 deployment
%AB_mean = [57.513667,-146.50];
%ABtext = repmat({'AB'},size(AB_mean,1),1);
%AB = [ABtext num2cell(AB_mean)];

%ALEUT01KS_mean = [52.316783,-188.52];

%KStext = repmat({'KS'},size(ALEUT01KS_mean,1),1);
%KS = [KStext num2cell(ALEUT01KS_mean)];

%KOA_mean = [57.224,-150.53];
%KOAtext = repmat({'KOA'},size(KOA_mean,1),1);
%KOA = [KOAtext num2cell(KOA_mean)];

%find means of sites with multiple deployments
[NClat,NClong] = meanm(NC_latLongs(:,1),NC_latLongs(:,2));
NC_mean = [NClat,NClong];
NCtext = repmat({'NC'},size(NC_mean,1),1);
NC = [NCtext num2cell(NC_mean)];

[BClat,BClong] = meanm(BC_latLongs(:,1),BC_latLongs(:,2));
BC_mean = [BClat, BClong];
BCtext = repmat({'BC'},size(BC_mean,1),1);
BC = [BCtext num2cell(BC_mean)];

[GSlat,GSlong] = meanm(GS_latLongs(:,1),GS_latLongs(:,2));
GS_mean = [GSlat, GSlong];
GStext = repmat({'GS'},size(GS_mean,1),1);
GS = [GStext num2cell(GS_mean)];

[BPlat, BPlong] = meanm(BP_latLongs(:,1),BP_latLongs(:,2));
BP_mean = [BPlat, BPlong];
BPtext = repmat({'BP'},size(BP_mean,1),1);
BP = [BPtext num2cell(BP_mean)];

[BSlat, BSlong] = meanm(BS_latLongs(:,1),BS_latLongs(:,2));
BS_mean = [BSlat, BSlong];
BStext = repmat({'BS'},size(BS_mean,1),1);
BS = [BStext num2cell(BS_mean)];

[WClat, WClong] = meanm(WC_latLongs(:,1),WC_latLongs(:,2));
WC_mean = [WClat, WClong];
WCtext = repmat({'WC'},size(WC_mean,1),1);
WC = [WCtext num2cell(WC_mean)];

[OClat, OClong] = meanm(OC_latLongs(:,1),OC_latLongs(:,2));
OC_mean = [OClat, OClong];
OCtext = repmat({'OC'},size(OC_mean,1),1);
OC = [OCtext num2cell(OC_mean)];

[HZlat, HZlong] = meanm(HZ_latLongs(:,1),HZ_latLongs(:,2));
HZ_mean = [HZlat, HZlong];
HZtext = repmat({'HZ'},size(HZ_mean,1),1);
HZ = [HZtext num2cell(HZ_mean)];

[JAXlat, JAXlong] = meanm(JAX_latLongs(:,1),JAX_latLongs(:,2));
JAX_mean = [JAXlat, JAXlong];
JAXtext = repmat({'JAX'},size(JAX_mean,1),1);
JAX = [JAXtext num2cell(JAX_mean)];

%create one table with all lat and longs
LL = [NC; BC; GS; BP; BS; WC; OC; HZ; JAX];
LatLong = cell2mat(LL(:,2:3));

LatLongTAB = array2table(LatLong);
LatLongTAB.Properties.VariableNames = {'Latitude' 'Longitude'};
LatLongTAB{:,'Site'} = {'NC'; 'BC'; 'GS'; 'BP'; 'BS'; 'WC'; 'OC'; 'HZ'; 'JAX'};
%LatLongTAB.Longitude(2) = -(LatLongTAB.Longitude(2) + 10);

lat_lims = [28.5 42];
long_lims = [-80 -65];
%% grey site map with no color distinction
figure(1)
LatLongTAB.Site = categorical(LatLongTAB.Site);
A = 500;
latitude = LatLongTAB.Latitude;
longitude = LatLongTAB.Longitude;
sites = LatLongTAB.Site;
gm = geoscatter(latitude,longitude,A,sites,'.');

geolimits([28.5 42],[-80 -65]); %geolimits([45 65],[-220 -120]);
text(latitude(1),longitude(1)-0.5,'NC','HorizontalAlignment','right','FontSize',10);
text(latitude(2),longitude(2)-0.5,'BC','HorizontalAlignment','right','FontSize',10);
text(latitude(3),longitude(3)-0.5,'GS','HorizontalAlignment','right','FontSize',10);
text(latitude(4),longitude(4)-0.5,'BP','HorizontalAlignment','right','FontSize',10);
text(latitude(5),longitude(5)-0.5,'BS','HorizontalAlignment','right','FontSize',10);
text(latitude(6),longitude(6)-0.5,'WC','HorizontalAlignment','right','FontSize',10);
text(latitude(7),longitude(7)-0.5,'OC','HorizontalAlignment','right','FontSize',10);
text(latitude(8),longitude(8)-0.5,'HZ','HorizontalAlignment','right','FontSize',10);
text(latitude(9),longitude(9)-0.5,'JAX','HorizontalAlignment','right','FontSize',10);
geolimits([28.5 42],[-80 -65]);

%% grey site map
figure(2)
LatLongTAB.Site = categorical(LatLongTAB.Site);
A = 50;
latitude = LatLongTAB.Latitude;
longitude = LatLongTAB.Longitude;

gm = geoscatter(latitude,longitude,A,'.','k');  
geolimits([28.5 42],[-80 -65]); %geolimits([-3 74],[-180 -110]);

%% Added by AD: grey site map with no color distinction and with labels
figure(3)
LatLongTAB.Site = categorical(LatLongTAB.Site);
A = 500;
latitude = LatLongTAB.Latitude;
longitude = LatLongTAB.Longitude;
sites = LatLongTAB.Site;
gm = geoscatter(latitude,longitude,A,sites,'.','k');

geolimits([28.5 42],[-80 -65]); %geolimits([45 65],[-220 -120]);
text(latitude(1),longitude(1)-0.5,'NC','HorizontalAlignment','right','FontSize',10);
text(latitude(2),longitude(2)-0.5,'BC','HorizontalAlignment','right','FontSize',10);
text(latitude(3),longitude(3)-0.5,'GS','HorizontalAlignment','right','FontSize',10);
text(latitude(4),longitude(4)-0.5,'BP','HorizontalAlignment','right','FontSize',10);
text(latitude(5),longitude(5)+1.5,'BS','HorizontalAlignment','right','FontSize',10);
text(latitude(6),longitude(6)-0.5,'WC','HorizontalAlignment','right','FontSize',10);
text(latitude(7),longitude(7)-0.5,'OC','HorizontalAlignment','right','FontSize',10);
text(latitude(8),longitude(8)-0.5,'HZ','HorizontalAlignment','right','FontSize',10);
text(latitude(9),longitude(9)+1.5,'JAX','HorizontalAlignment','right','FontSize',10);
geolimits([28.5 42],[-80 -65]);