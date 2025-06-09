%% This code was modified by NP on 4/4/2022 from the M_Map website to plot bathymetry in the GOA/BSAI
 % Designed to work in R2016b
%source - https://www.eoas.ubc.ca/~rich/map.html#examples

% To get m_map functions, on manatee14, the code currently references:
% C:\Users\HARP\Documents\AD_Working\m_map1.4\m_map

close all;clear all;clc;
%% Load Site data
% load lat and longs for each site
SaveDir = 'H:\WAT_TPWS_metadataReduced\Plots\Site Maps';

HATA_latLongs = [35.34218333, -74.86;
35.30183333, -74.86;
35.58413333, -74.75;
35.58351667, -74.75];

HATB_latLongs = [35.58413333, -74.75;
35.58351667, -74.75;
35.58976667, -74.75;
35.5893, -74.75];

[HATAlat, HATAlong] = meanm(HATA_latLongs(:,1),HATA_latLongs(:,2));
HATA_mean = [HATAlat, HATAlong];
HATAtext = repmat({'HAT_A'},size(HATA_mean,1),1);
HATA = [HATAtext num2cell(HATA_mean)];

[HATBlat, HATBlong] = meanm(HATB_latLongs(:,1),HATB_latLongs(:,2));
HATB_mean = [HATBlat, HATBlong];
HATBtext = repmat({'HAT_B'},size(HATB_mean,1),1);
HATB = [HATBtext num2cell(HATB_mean)];

%create one table with all lat and longs
LL = [HATA; HATB];
LatLong = cell2mat(LL(:,2:3));

SiteData = array2table(LatLong);
SiteData.Properties.VariableNames = {'Latitude' 'Longitude'};
SiteData{:,'Site'} = {'HAT A'; 'HAT B'};
%% Create map 
%The projections that successfully work: UTM, Transverse mercator (or this), Mercator (probably the best),...
%Miller Cylindrical, Albers Equal-Area Conic, Lambert Conformal Conic, Hammer-Aitoff, Mollweide, Robinson
m_proj('Robinson','long',[-74.85 -74.65],'lat',[35.35 35.65]); %identify ranges of the map
[CS,CH]=m_etopo2('contourf',[-6000:500:-1000 -700 -500 -200 0 ],'edgecolor','none'); %load bathymetry from ETOPO 1
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_line(SiteData.Longitude(1),SiteData.Latitude(1),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(1)+0.02,SiteData.Latitude(1),SiteData.Site(1),'FontWeight','Bold','FontSize',16);
m_line(SiteData.Longitude(2),SiteData.Latitude(2),'marker','.','markersize',15,...
          'linest','none','color','k','clip','point');
m_text(SiteData.Longitude(2)+0.02,SiteData.Latitude(2),SiteData.Site(2),'FontWeight','Bold','FontSize',16);
m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
colormap(m_colmap('blues'));  
caxis([-6000 000]);
m_ruler(.7,[.05 .4],3,'fontsize',16)
set(gcf,'color','w');  % otherwise 'print' turns lakes black
%% Save plot
x0=10;
y0=10;
width=550*1.5;
height=400*1.5;
set(gcf,'position',[x0,y0,width,height])

%png
filename = [SaveDir,'\SiteMap_M_Map_WAT_HAT_inset.png'];
saveas(gcf,filename, 'png')

%pdf
filename = [SaveDir,'\SiteMap_M_Map_WAT_HAT_inset.pdf'];
saveas(gcf,filename, 'pdf')