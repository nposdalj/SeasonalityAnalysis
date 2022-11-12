%% This code was modified by NP on 4/4/2022 from the M_Map website to plot bathymetry in the GOA/BSAI
 % Designed to work in R2016b
%source - https://www.eoas.ubc.ca/~rich/map.html#examples
close all; clc;

%% Parameters defined by user + Load Site data

% Parameters defined by user
GDrive = 'G';
Region = 'WAT';
Sites = {'NC' 'BC' 'GS' 'BP' 'BS' 'WC' 'OC' 'HZ' 'JAX'};

% No need to change the following
numsites = length(Sites);

% load lat and longs for each site
SaveDir = [GDrive ':\My Drive\' Region '_TPWS_metadataReduced\Plots'];
Site_latLongs.NC = [39.83248333,-69.98; %01 in decimals
39.83238333, -69.98; %nc02
39.83258333, -69.98; %nc03
39.83295, -69.98]; %nc04

Site_latLongs.BC = [39.19105, -72.23; %bc01
39.1905, -72.23; %bc02
39.19191667, -72.23]; %bc03

Site_latLongs.GS = [33.66563333, -76.00; %gs01
33.66701667, -76.00; %gs02
33.66991667, -76.00]; %gs03
%gs04???

Site_latLongs.BP = [32.10603333, -77.09; %bp01
32.10695, -77.09; %bp02
32.10526667, -77.09]; %bp03
%bp04???

Site_latLongs.BS = [30.58378333, -77.39; %bs01
30.58303333, -77.39; %bs02
30.58295, -77.39]; %bs03

Site_latLongs.WC = [38.37415, -73.37; %wc01
38.37385, -73.37; %wc02
38.37336667, -73.37]; %wc03

Site_latLongs.OC = [40.2633, -67.99; %oc01
40.26331667, -67.99; %oc02
40.26333333, -67.99; %oc03
40.23, -67.98]; %oc04

Site_latLongs.HZ = [41.06191667, -66.35; %hz01
41.06183333, -66.35; %hz02
41.06165, -66.35; %hz03
41.06165, -66.35]; %hz04

Site_latLongs.JAX = [30.15183333, -79.77; %JAX_D_13
30.15268333, -79.77; %JAX_D_14
30.15225, -79.77]; %JAX_D_15
%Not sure if other JAX deployments will be included but this is it for now,
%update if needed.

%find means of sites with multiple deployments
[NClat,NClong] = meanm(Site_latLongs.NC(:,1),Site_latLongs.NC(:,2));
NC_mean = [NClat,NClong];
NCtext = repmat({'NC'},size(NC_mean,1),1);
NC = [NCtext num2cell(NC_mean)];

[BClat,BClong] = meanm(Site_latLongs.BC(:,1),Site_latLongs.BC(:,2));
BC_mean = [BClat, BClong];
BCtext = repmat({'BC'},size(BC_mean,1),1);
BC = [BCtext num2cell(BC_mean)];

[GSlat,GSlong] = meanm(Site_latLongs.GS(:,1),Site_latLongs.GS(:,2));
GS_mean = [GSlat, GSlong];
GStext = repmat({'GS'},size(GS_mean,1),1);
GS = [GStext num2cell(GS_mean)];

[BPlat, BPlong] = meanm(Site_latLongs.BP(:,1),Site_latLongs.BP(:,2));
BP_mean = [BPlat, BPlong];
BPtext = repmat({'BP'},size(BP_mean,1),1);
BP = [BPtext num2cell(BP_mean)];

[BSlat, BSlong] = meanm(Site_latLongs.BS(:,1),Site_latLongs.BS(:,2));
BS_mean = [BSlat, BSlong];
BStext = repmat({'BS'},size(BS_mean,1),1);
BS = [BStext num2cell(BS_mean)];

[WClat, WClong] = meanm(Site_latLongs.WC(:,1),Site_latLongs.WC(:,2));
WC_mean = [WClat, WClong];
WCtext = repmat({'WC'},size(WC_mean,1),1);
WC = [WCtext num2cell(WC_mean)];

[OClat, OClong] = meanm(Site_latLongs.OC(:,1),Site_latLongs.OC(:,2));
OC_mean = [OClat, OClong];
OCtext = repmat({'OC'},size(OC_mean,1),1);
OC = [OCtext num2cell(OC_mean)];

[HZlat, HZlong] = meanm(Site_latLongs.HZ(:,1),Site_latLongs.HZ(:,2));
HZ_mean = [HZlat, HZlong];
HZtext = repmat({'HZ'},size(HZ_mean,1),1);
HZ = [HZtext num2cell(HZ_mean)];

[JAXlat, JAXlong] = meanm(Site_latLongs.JAX(:,1),Site_latLongs.JAX(:,2));
JAX_mean = [JAXlat, JAXlong];
JAXtext = repmat({'JAX'},size(JAX_mean,1),1);
JAX = [JAXtext num2cell(JAX_mean)];

%create one table with all lat and longs
LL = [NC; BC; GS; BP; BS; WC; OC; HZ; JAX];
LatLong = cell2mat(LL(:,2:3));

SiteData = array2table(LatLong);
SiteData.Properties.VariableNames = {'Latitude' 'Longitude'};
SiteData{:,'Site'} = {'NC'; 'BC'; 'GS'; 'BP'; 'BS'; 'WC'; 'OC'; 'HZ'; 'JAX'};

%% Loop through sites, create map of deployments
for i = 1:numsites
    site_i = Sites(i);
    latLongs_i = Site_latLongs.(char(site_i));
    latLongs_i = [latLongs_i (1:size(latLongs_i, 1))'];
    latLongs_i = array2table(latLongs_i);
    latLongs_i.Properties.VariableNames = {'Latitude', 'Longitude', 'DeplNum'};
    
    plotLatLims = [min(latLongs_i.Latitude)-.01 max(latLongs_i.Latitude+.01)];
    plotLongLims = [min(latLongs_i.Longitude)-.02 max(latLongs_i.Longitude+.02)];
    
    figure(i)
    
    m_proj('Robinson','long',plotLongLims,'lat',plotLatLims); %identify ranges of the map
    [CS,CH]=m_etopo2('contourf',[-1200:10:-700],'edgecolor','none'); %load bathymetry from ETOPO 1
    m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
    
    numDepls = size(latLongs_i, 1);
    for j = 1:numDepls
        m_line(latLongs_i.Longitude(j),latLongs_i.Latitude(j),'marker','s',...
            'linest','none','markerfacecolor','w','clip','point');
        m_text(latLongs_i.Longitude(j)+0.002,latLongs_i.Latitude(j),[char(site_i) '-0' num2str(latLongs_i.DeplNum(j))],'FontWeight','Bold','FontSize',12);
        
    end
    
    
    m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
    colormap(m_colmap('blues'));
    caxis([-1200 -700]);
    [ax,h]=m_contfbar([.6 .8],.2,CS,CH,'endpiece','no','axfrac',.05,'YColor','white');
    title(ax,'meters','FontWeight','Bold')
    ax.FontSize = 12;
    ax.FontWeight = 'Bold';
    set(gcf,'color','w');  % otherwise 'print' turns lakes black
    
    m_ruler(.8,[.5 .9],3,'fontsize',12)
    
    saveas(gcf, ['G:\My Drive\WAT_TPWS_metadataReduced\Plots\Site Maps\DeplMap_' char(site_i) '.png'])
end

%% Create map 
%The projections that successfully work: UTM, Transverse mercator (or this), Mercator (probably the best),...
%Miller Cylindrical, Albers Equal-Area Conic, Lambert Conformal Conic, Hammer-Aitoff, Mollweide, Robinson
m_proj('Robinson','long',[-82 -62],'lat',[28.5 42]); %identify ranges of the map
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
m_line(SiteData.Longitude(6),SiteData.Latitude(6),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(6)+0.1,SiteData.Latitude(6),SiteData.Site(6),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(7),SiteData.Latitude(7),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(7)+0.1,SiteData.Latitude(7),SiteData.Site(7),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(8),SiteData.Latitude(8),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(8)+0.1,SiteData.Latitude(8),SiteData.Site(8),'FontWeight','Bold','FontSize',12);
m_line(SiteData.Longitude(9),SiteData.Latitude(9),'marker','s',...
          'linest','none','markerfacecolor','w','clip','point');
m_text(SiteData.Longitude(9)+0.1,SiteData.Latitude(9),SiteData.Site(9),'FontWeight','Bold','FontSize',12);
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