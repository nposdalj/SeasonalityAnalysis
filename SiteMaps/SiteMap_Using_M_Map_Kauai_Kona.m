close all;clear all;clc;
GDrive = 'L';
SaveDir = [GDrive,':\My Drive\FourLoko\SiteMaps']
%% load lat and longs for each site
HAWAIIK_latLongs = [19.5815,-156.02; %01
19.58058333, -156.02; %02
19.58156667, -156.01; %03
19.57761667, -156.01; %05
19.58261667, -156.02; %06
19.58196667, -156.02; %07
19.58146667, -156.02; %08
19.58148333, -156.02; %09
19.58231667, -156.02; %10
19.58241667, -156.02; %11
19.58275, -156.02; %13
19.58293333, -156.02; %14
19.58293333, -156.02; %15
19.58293333, -156.02; %16
19.583083335, -153.02; %17
19.58308333, -156.02; %18
19.58323333, -156.02; %19
19.58323333, -156.02; %20
%missing 22
19.58306667, -156.02; %23_01
19.0097175, -156.02; %25
%missing 26
19.58298333, -156.00]; %27

Kauai_latLongs = [21.952,-159.89;
    21.953,-159.89;
    21.9492,-159.89];
%% find means of sites with multiple deployments
[HAWAIIKlat, HAWAIIKlong] = meanm(HAWAIIK_latLongs(:,1),HAWAIIK_latLongs(:,2));
HAWAIIK_mean = [HAWAIIKlat, HAWAIIKlong];
HAWAIIKtext = repmat({'Hawaii'},size(HAWAIIK_mean,1),1);
HAWAIIK = [HAWAIIKtext num2cell(HAWAIIK_mean)];

[Kauailat, Kauailong] = meanm(Kauai_latLongs(:,1),Kauai_latLongs(:,2));
Kauai_mean = [Kauailat, Kauailong];
Kauaitext = repmat({'Kauai'},size(Kauai_mean,1),1);
Kauai = [Kauaitext num2cell(Kauai_mean)];
%% create one table with all lat and longs
LL = [HAWAIIK; Kauai];
LatLong = cell2mat(LL(:,2:3));

LatLongTAB = array2table(LatLong);
LatLongTAB.Properties.VariableNames = {'Latitude' 'Longitude'};

LatLongTAB{:,'Site'} = {'Kona'; 'Kauai'};
%% grey site map
m_proj('Robinson','long',[-122.5 -122.2],'lat',[36.25 36.45]); %identify ranges of the map
[CS,CH]=m_etopo2('contourf',[-2500:100:-1000 -500 -200 0 ],'edgecolor','none'); %load bathymetry from ETOPO 1
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_line(SiteData.Longitude(1),SiteData.Latitude(1),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(1)+0.003,SiteData.Latitude(1),SiteData.Deployment(1),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(2),SiteData.Latitude(2),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(2),SiteData.Latitude(2),SiteData.Deployment(2),'FontWeight','Bold','FontSize',12);
m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
colormap(m_colmap('blues'));  
caxis([-2500 000]);
[ax,h]=m_contfbar([.6 .8],.2,CS,CH,'endpiece','no','axfrac',.05,'YColor','white');
title(ax,'meters','FontWeight','Bold')
ax.FontSize = 12;
ax.FontWeight = 'Bold';
set(gcf,'color','w');  % otherwise 'print' turns lakes black
%% Save plot
x0=10;
y0=10;
width=550*1.5;
height=400*1.5;
set(gcf,'position',[x0,y0,width,height])

%png
filename = [SaveDir,'\SiteMap_M_Map_Kona_Kauai.png'];
saveas(gcf,filename, 'png')

%pdf
filename = [SaveDir,'\SiteMap_M_Map_Kona_Kauai.pdf'];
saveas(gcf,filename, 'pdf')