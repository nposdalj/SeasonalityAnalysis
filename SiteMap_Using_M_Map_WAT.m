%% This code was modified by NP on 4/4/2022 from the M_Map website to plot bathymetry in the GOA/BSAI
%source - https://www.eoas.ubc.ca/~rich/map.html#examples
close all;clear all;clc;
%% Load Site data
% load lat and longs for each site
SaveDir = 'I:\My Drive\Manuscripts\GOA\Figures';
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

%create one table with all lat and longs
LL = [NC; BC; GS; BP; BS];
LatLong = cell2mat(LL(:,2:3));

SiteData = array2table(LatLong);
SiteData.Properties.VariableNames = {'Latitude' 'Longitude'};
SiteData{:,'Site'} = {'NC'; 'BC'; 'GS'; 'BP'; 'BS'};
%% Create map 
%The projections that successfully work: UTM, Transverse mercator (or this), Mercator (probably the best),...
%Miller Cylindrical, Albers Equal-Area Conic, Lambert Conformal Conic, Hammer-Aitoff, Mollweide, Robinson
m_proj('Robinson','long',[-80 -65],'lat',[28.5 42]); %identify ranges of the map
[CS,CH]=m_etopo2('contourf',[-7000:1000:-1000 -500 -200 0 ],'edgecolor','none'); %load bathymetry from ETOPO 1
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_line(SiteData.Longitude(1),SiteData.Latitude(1),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(1)+0.1,SiteData.Latitude(1),SiteData.Site(1),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(2),SiteData.Latitude(2),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(2)+0.1,SiteData.Latitude(2),SiteData.Site(2),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(3),SiteData.Latitude(3),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(3)+0.1,SiteData.Latitude(3),SiteData.Site(3),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(4),SiteData.Latitude(4),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(4)+0.1,SiteData.Latitude(4),SiteData.Site(4),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(5),SiteData.Latitude(5),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(5)+0.1,SiteData.Latitude(5),SiteData.Site(5),'FontWeight','Bold','FontSize',12);
m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
colormap(m_colmap('blues'));  
caxis([-7000 000]);
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
filename = [SaveDir,'\SiteMap_M_Map_WAT.png'];
saveas(gcf,filename, 'png')

%pdf
filename = [SaveDir,'\SiteMap_M_Map_WAT.pdf'];
saveas(gcf,filename, 'pdf')