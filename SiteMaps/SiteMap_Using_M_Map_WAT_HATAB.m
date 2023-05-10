%% This code was modified by NP on 4/4/2022 from the M_Map website to plot bathymetry in the GOA/BSAI
 % Designed to work in R2016b
%source - https://www.eoas.ubc.ca/~rich/map.html#examples

% To get m_map functions, on manatee14, the code currently references:
% C:\Users\HARP\Documents\AD_Working\m_map1.4\m_map

close all;clear all;clc;
%% Load Site data
% load lat and longs for each site
SaveDir = 'G:\My Drive\WAT_TPWS_metadataReduced\Plots';
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

HAT_latLongs = [35.34218333, -74.86;
35.30183333, -74.86;
35.58413333, -74.75;
35.58351667, -74.75;
35.58413333, -74.75;
35.58351667, -74.75;
35.58976667, -74.75;
35.5893, -74.75];

NFC_latLongs = [37.16651667, -74.47;
37.1674, - 74.47;
37.16451667, -74.47];

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

[NFClat, NFClong] = meanm(NFC_latLongs(:,1),NFC_latLongs(:,2));
NFC_mean = [NFClat, NFClong];
NFCtext = repmat({'NFC'},size(NFC_mean,1),1);
NFC = [NFCtext num2cell(NFC_mean)];

[HATlat, HATlong] = meanm(HAT_latLongs(:,1),HAT_latLongs(:,2));
HAT_mean = [HATlat, HATlong];
HATtext = repmat({'HAT'},size(HAT_mean,1),1);
HAT = [HATtext num2cell(HAT_mean)];

%create one table with all lat and longs
LL = [NC; BC; GS; BP; BS; WC; OC; HZ; JAX; NFC; HAT];
LatLong = cell2mat(LL(:,2:3));

SiteData = array2table(LatLong);
SiteData.Properties.VariableNames = {'Latitude' 'Longitude'};
SiteData{:,'Site'} = {'NC'; 'BC'; 'GS'; 'BP'; 'BS'; 'WC'; 'OC'; 'HZ'; 'JAX'; 'NFC'; 'HAT'};
%% Create map 
%The projections that successfully work: UTM, Transverse mercator (or this), Mercator (probably the best),...
%Miller Cylindrical, Albers Equal-Area Conic, Lambert Conformal Conic, Hammer-Aitoff, Mollweide, Robinson
m_proj('Robinson','long',[-82 -62],'lat',[28.5 42]); %identify ranges of the map
[CS,CH]=m_etopo2('contourf',[-7000:1000:-1000 -500 -200 0 ],'edgecolor','none'); %load bathymetry from ETOPO 1
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_line(SiteData.Longitude(1),SiteData.Latitude(1),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(1)+0.4,SiteData.Latitude(1),SiteData.Site(1),'FontWeight','Bold','FontSize',12);
% m_line(SiteData.Longitude(1),SiteData.Latitude(1),'marker','s',...   % Previous point formatting style
%           'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(1)+0.1,SiteData.Latitude(1),SiteData.Site(1),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(2),SiteData.Latitude(2),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(2)+0.4,SiteData.Latitude(2),SiteData.Site(2),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(3),SiteData.Latitude(3),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(3)+0.4,SiteData.Latitude(3),SiteData.Site(3),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(4),SiteData.Latitude(4),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(4)+0.4,SiteData.Latitude(4),SiteData.Site(4),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(5),SiteData.Latitude(5),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(5)+0.4,SiteData.Latitude(5),SiteData.Site(5),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(6),SiteData.Latitude(6),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(6)+0.4,SiteData.Latitude(6),SiteData.Site(6),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(7),SiteData.Latitude(7),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(7)+0.4,SiteData.Latitude(7),SiteData.Site(7),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(8),SiteData.Latitude(8),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(8)+0.4,SiteData.Latitude(8),SiteData.Site(8),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(9),SiteData.Latitude(9),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(9)+0.4,SiteData.Latitude(9),SiteData.Site(9),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(10),SiteData.Latitude(10),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(10)+0.4,SiteData.Latitude(10),SiteData.Site(10),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(11),SiteData.Latitude(11),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(11)+0.4,SiteData.Latitude(11),SiteData.Site(11),'FontWeight','Bold','FontSize',12);
m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
colormap(m_colmap('blues'));  
caxis([-7000 000]);
%[ax,h]=m_contfbar([.68 .95],.35,CS,CH,'endpiece','no','axfrac',.05,'YColor','white');
[ax,h]=m_contfbar(.18,[.45 .85],CS,CH,'endpiece','no','axfrac',.03,'YColor','black');
title(ax,'meters','FontWeight','Bold')
ax.FontSize = 11;
%ax.XColor = [1, 1, 1];
%ax.YColor = [1, 1, 1];
% ax.FontWeight = 'Bold';
set(gcf,'color','w');  % otherwise 'print' turns lakes black
%% Save plot
x0=10;
y0=10;
width=550*1.5;
height=400*1.5;
set(gcf,'position',[x0,y0,width,height])

%png
filename = [SaveDir,'\SiteMap_M_Map_WAT_HATAB.png'];
saveas(gcf,filename, 'png')

%pdf
filename = [SaveDir,'\SiteMap_M_Map_WAT_HATAB.pdf'];
saveas(gcf,filename, 'pdf')