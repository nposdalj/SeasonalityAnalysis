%% This code was modified by NP on 4/4/2022 from the M_Map website to plot bathymetry in the GOA/BSAI
 % Designed to work in R2016b
%source - https://www.eoas.ubc.ca/~rich/map.html#examples
%To use this code, you need to have m_map added to the Matlab path located
%on Natalie's Google drive I:\My Drive\m_map
close all;clear all;clc;
%% Load Site data
% load lat and longs for each site
SaveDir = 'I:\My Drive\CCE_TPWS_metadataReduced\Plots\PS';
PS_latLongs = [36.2991, -122.3938;
36.29873333, -122.39305;
36.29868333, -122.3933;
36.3146, -122.3910833;
36.39045, -122.30585;
36.38886667, -122.3066;
36.38893333, -122.3068167;
36.38945, -122.3076167;
36.39131667, -122.3075;
36.3912, -122.3069833;
36.29908333, -122.3938833;
36.29911667, -122.3938333;
36.3703, -122.314767;
36.37023333, -122.3144;
36.36963333, -122.31495;
36.37066667, -122.3146167;];

SiteData = array2table(PS_latLongs);
SiteData.Properties.VariableNames = {'Latitude' 'Longitude'};
SiteData{:,'Deployment'} = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '11'; '12'; '13'; '14'; '15'; '16'; '17'};
%% Create map 
%The projections that successfully work: UTM, Transverse mercator (or this), Mercator (probably the best),...
%Miller Cylindrical, Albers Equal-Area Conic, Lambert Conformal Conic, Hammer-Aitoff, Mollweide, Robinson
m_proj('Robinson','long',[-122.5 -122.2],'lat',[36.25 36.45]); %identify ranges of the map
[CS,CH]=m_etopo2('contourf',[-2500:100:-1000 -500 -200 0 ],'edgecolor','none'); %load bathymetry from ETOPO 1
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_line(SiteData.Longitude(1),SiteData.Latitude(1),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(1)+0.003,SiteData.Latitude(1),SiteData.Deployment(1),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(2),SiteData.Latitude(2),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(2),SiteData.Latitude(2),SiteData.Deployment(2),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(3),SiteData.Latitude(3),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(3),SiteData.Latitude(3),SiteData.Deployment(3),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(4),SiteData.Latitude(4),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(4)+0.003,SiteData.Latitude(4),SiteData.Deployment(4),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(5),SiteData.Latitude(5),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(5),SiteData.Latitude(5),SiteData.Deployment(5),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(6),SiteData.Latitude(6),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(6),SiteData.Latitude(6),SiteData.Deployment(6),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(7),SiteData.Latitude(7),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(7),SiteData.Latitude(7),SiteData.Deployment(7),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(8),SiteData.Latitude(8),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(8),SiteData.Latitude(8),SiteData.Deployment(8),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(9),SiteData.Latitude(9),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(9),SiteData.Latitude(9),'5-9,11','FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(10),SiteData.Latitude(10),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(10),SiteData.Latitude(10),SiteData.Deployment(10),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(11),SiteData.Latitude(11),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(11),SiteData.Latitude(11),SiteData.Deployment(11),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(12),SiteData.Latitude(12),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(12),SiteData.Latitude(12),'1-3,12-13','FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(13),SiteData.Latitude(13),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(13),SiteData.Latitude(13),SiteData.Deployment(13),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(14),SiteData.Latitude(14),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(14),SiteData.Latitude(14),SiteData.Deployment(14),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(15),SiteData.Latitude(15),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
% m_text(SiteData.Longitude(15),SiteData.Latitude(15),SiteData.Deployment(15),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(16),SiteData.Latitude(16),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(16),SiteData.Latitude(16),'14-17','FontWeight','Bold','FontSize',12);
m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
colormap(m_colmap('blues'));  
caxis([-2500 000]);
%[ax,h]=m_contfbar([.68 .95],.35,CS,CH,'endpiece','no','axfrac',.05,'YColor','white');
[ax,h]=m_contfbar([.6 .8],.2,CS,CH,'endpiece','no','axfrac',.05,'YColor','white');
title(ax,'meters','FontWeight','Bold')
ax.FontSize = 12;
%ax.XColor = [1, 1, 1];
%ax.YColor = [1, 1, 1];
ax.FontWeight = 'Bold';
set(gcf,'color','w');  % otherwise 'print' turns lakes black
%% Save plot
x0=10;
y0=10;
width=550*1.5;
height=400*1.5;
set(gcf,'position',[x0,y0,width,height])

%png
filename = [SaveDir,'\SiteMap_M_Map_PS.png'];
saveas(gcf,filename, 'png')

%pdf
filename = [SaveDir,'\SiteMap_M_Map_PS.pdf'];
saveas(gcf,filename, 'pdf')