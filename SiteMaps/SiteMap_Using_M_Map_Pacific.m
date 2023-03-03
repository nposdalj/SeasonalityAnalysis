
close all;clear all;clc;
%% load lat and longs for each site
CB_latLongs = [58.645683,-148.07; %01
58.671333,-148.02; %02
58.673483,-148.00; %03
58.671867,-148.02; %04
58.671,-148.02; %05
58.670817,-148.02; %06
58.65525,-148.09; %07
58.6695,-148.03; %08
58.6695,-148.0338889; %09
58.66961667,-148.03]; %10

QN_latLongs = [56.339017,-145.18; %01
56.339383,-145.18; %02
56.3413,-145.18; %04
56.340683,-145.18; %05
56.339833,-145.19]; %06

PT_latLongs = [56.24345, -142.75; %01
56.243917, -142.75; %02
56.242917, -142.75; %03
56.243333, -142.75]; %04

OCNMSQC_latLongs = [47.466,-125.16; %6
47.50005,-125.36; %12
47.500433,-125.36; %14
47.500533,-125.36; %15
47.500633,-125.65]; %16

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

PAL_latLongs = [5.8641, -162.17; %2
5.8647, -162.16; %3
5.86545, -162.16; %4
5.8651, -162.16; %5
5.86295, -162.16
5.9042, -162.04; %7
5.895316667, -162.04; %8
5.894833333, -162.04; %9
5.895083333, -162.04]; %10

PHR_latLongs = [27.72528333, -175.64; %1
27.727, -175.63; %2
27.72531667, -175.64; %4
27.72535, -175.64; %5
27.72515, -175.64; %6
27.74103333, -175.56; %7
% missing 8
% missing 9
27.74098333, -175.56]; %10

SAIPAN_latLongs = [15.31663333, 145.46; %1
15.3171, 145.46; %2
15.31778333, 145.46; %3
15.32125, 145.46; %4
15.32125, 145.46; %5
15.31743333, 145.46]; %6
%missing 7]; %7

TIN_latLongs = [15.03906667, 145.75; %2
15.0398, 145.75; %3
15.03735, 145.75; %4
15.03735, 145.75; %5
15.04003333, 145.75]; %6
%missing 7]; %7

Wake_latLongs = [19.22, -166.69; %1
19.2209, -166.69; %3
19.22156667, -166.69; %4
19.22216667, -166.69; %5
19.22346667, -166.69; %6
19.37206667, -166.69]; %7

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

ALEUTBD_latLongs = [52.633333,-174.37; %BD02
52.076,-174.36]; %BD03

GofCA_latLongs = [29.028,-113.38;
    29.027, -113.38];

GI_latLongs = [29.14,-118.26;
    29.14,-118.26];    

Kauai_latLongs = [21.952,-159.89;
    21.953,-159.89;
    21.9492,-159.89];
%% find means of sites with multiple deployments
[CBlat,CBlong] = meanm(CB_latLongs(:,1),CB_latLongs(:,2));
CB_mean = [CBlat,CBlong];
CBtext = repmat({'CB'},size(CB_mean,1),1);
CB = [CBtext num2cell(CB_mean)];

[QNlat,QNlong] = meanm(QN_latLongs(:,1),QN_latLongs(:,2));
QN_mean = [QNlat, QNlong];
QNtext = repmat({'QN'},size(QN_mean,1),1);
QN = [QNtext num2cell(QN_mean)];

[BDlat, BDlong] = meanm(ALEUTBD_latLongs(:,1),ALEUTBD_latLongs(:,2));
ALEUTBD_mean = [BDlat, BDlong];
BDtext = repmat({'BD'},size(ALEUTBD_mean,1),1);
BD = [BDtext num2cell(ALEUTBD_mean)];

[PTlat,PTlong] = meanm(PT_latLongs(:,1),PT_latLongs(:,2));
PT_mean = [PTlat, PTlong];
PTtext = repmat({'PT'},size(PT_mean,1),1);
PT = [PTtext num2cell(PT_mean)];

[QClat,QClong] = meanm(OCNMSQC_latLongs(:,1),OCNMSQC_latLongs(:,2));
OCNMSQC_mean = [QClat, QClong];
OCNMStext = repmat({'OCNMS'},size(OCNMSQC_mean,1),1);
OCNMS = [OCNMStext num2cell(OCNMSQC_mean)];

[HAWAIIKlat, HAWAIIKlong] = meanm(HAWAIIK_latLongs(:,1),HAWAIIK_latLongs(:,2));
HAWAIIK_mean = [HAWAIIKlat, HAWAIIKlong];
HAWAIIKtext = repmat({'Hawaii'},size(HAWAIIK_mean,1),1);
HAWAIIK = [HAWAIIKtext num2cell(HAWAIIK_mean)];

[PALlat, PALlong] = meanm(PAL_latLongs(:,1),PAL_latLongs(:,2));
PAL_mean = [PALlat, PALlong];
PALtext = repmat({'PAL'},size(PAL_mean,1),1);
PAL = [PALtext num2cell(PAL_mean)];

[PHRlat, PHRlong] = meanm(PHR_latLongs(:,1),PHR_latLongs(:,2));
PHR_mean = [PHRlat, PHRlong];
PHRtext = repmat({'PHR'},size(PHR_mean,1),1);
PHR = [PHRtext num2cell(PHR_mean)];

[SAIPANlat, SAIPANlong] = meanm(SAIPAN_latLongs(:,1),SAIPAN_latLongs(:,2));
SAIPAN_mean = [SAIPANlat, SAIPANlong];
SAIPANtext = repmat({'Saipan'},size(SAIPAN_mean,1),1);
SAIPAN = [SAIPANtext num2cell(SAIPAN_mean)];

[TINlat, TINlong] = meanm(TIN_latLongs(:,1),TIN_latLongs(:,2));
TIN_mean = [TINlat, TINlong];
TINtext = repmat({'Tinian'},size(TIN_mean,1),1);
TIN = [TINtext num2cell(TIN_mean)];

[Wakelat, Wakelong] = meanm(Wake_latLongs(:,1),Wake_latLongs(:,2));
Wake_mean = [Wakelat, Wakelong];
Waketext = repmat({'Wake'},size(Wake_mean,1),1);
Wake = [Waketext num2cell(Wake_mean)];

[PSlat, PSlong] = meanm(PS_latLongs(:,1),PS_latLongs(:,2));
PS_mean = [PSlat, PSlong];
PStext = repmat({'PS'},size(PS_mean,1),1);
PS = [PStext num2cell(PS_mean)];

[GofCAlat, GofCAlong] = meanm(GofCA_latLongs(:,1),GofCA_latLongs(:,2));
GofCA_mean = [GofCAlat, GofCAlong];
GofCAtext = repmat({'GofCA'},size(GofCA_mean,1),1);
GofCA = [GofCAtext num2cell(GofCA_mean)];

[GIlat, GIlong] = meanm(GI_latLongs(:,1),GI_latLongs(:,2));
GI_mean = [GIlat, GIlong];
GItext = repmat({'GI'},size(GI_mean,1),1);
GI = [GItext num2cell(GI_mean)];

[Kauailat, Kauailong] = meanm(Kauai_latLongs(:,1),Kauai_latLongs(:,2));
Kauai_mean = [Kauailat, Kauailong];
Kauaitext = repmat({'Kauai'},size(Kauai_mean,1),1);
Kauai = [Kauaitext num2cell(Kauai_mean)];
%% Sites with only one deployment
AB_mean = [57.513667,-146.50];
ABtext = repmat({'AB'},size(AB_mean,1),1);
AB = [ABtext num2cell(AB_mean)];

ALEUT01KS_mean = [52.316783,-188.52];
KStext = repmat({'KS'},size(ALEUT01KS_mean,1),1);
KS = [KStext num2cell(ALEUT01KS_mean)];

KOA_mean = [57.224,-150.53];
KOAtext = repmat({'KOA'},size(KOA_mean,1),1);
KOA = [KOAtext num2cell(KOA_mean)];

DCPP_mean = [35.4,-122.4];
DCPPtext = repmat({'DCPP'},size(DCPP_mean,1),1);
DCPP = [DCPPtext num2cell(DCPP_mean)];

CCE_mean = [33.48,-121.42];
CCEtext = repmat({'CCE'},size(CCE_mean,1),1);
CCE = [CCEtext num2cell(CCE_mean)];

HOKE_mean = [32.1,-126.91];
HOKEtext = repmat({'HOKE'},size(HOKE_mean,1),1);
HOKE = [HOKEtext num2cell(HOKE_mean)];

CORC_mean = [31.747,-121.38];
CORCtext = repmat({'CORC'},size(CORC_mean,1),1);
CORC = [CORCtext num2cell(CORC_mean)];

CSM_mean = [18.722,-158.25];
CSMtext = repmat({'CSM'},size(CSM_mean,1),1);
CSM = [CSMtext num2cell(CSM_mean)];

LSM_mean = [28.62745,-176.73];
LSMtext = repmat({'LSM'},size(LSM_mean,1),1);
LSM = [LSMtext num2cell(LSM_mean)];

EQ_mean = [0.44345,-164];
EQtext = repmat({'EQ'},size(EQ_mean,1),1);
EQ = [EQtext num2cell(EQ_mean)];

KR_mean = [6.365,-162.2923];
KRtext = repmat({'KR'},size(KR_mean,1),1);
KR = [KRtext num2cell(KR_mean)];

Pagan_mean = [17.963,-145.48];
Pagantext = repmat({'Pagan'},size(Pagan_mean,1),1);
Pagan = [Pagantext num2cell(Pagan_mean)];
%% create one table with all lat and longs
LL = [CB; HAWAIIK; OCNMS; PAL; PHR; PT; QN; SAIPAN; TIN; Wake; PS; DCPP;...
    CCE; HOKE; CORC; GI; GofCA; Kauai; LSM; KR; EQ; Pagan; AB; BD; KOA; KS; CSM];
LatLong = cell2mat(LL(:,2:3));

LatLongTAB = array2table(LatLong);
LatLongTAB.Properties.VariableNames = {'Latitude' 'Longitude'};

LatLongTAB{:,'Site'} = {'CB'; 'Kona'; 'QC'; 'Palmyra'; 'PHR'; 'PT'; 'QN';...
    'Saipan'; 'Tinian'; 'Wake'; 'PS'; 'DCPP'; 'CCE'; 'Hoke'; 'CORC'; 'GI';...
    'GofCA'; 'Kauai'; 'LSM'; 'KR'; 'EQ'; 'Pagan'; 'AB'; 'BD'; 'KOA'; 'KS'; 'CSM'};

LatLongTAB.Longitude(8) = -2.154600000000000e+02;
LatLongTAB.Longitude(9) = -2.157500000000000e+02;
%% grey site map
figure(1)
LatLongTAB.Site = categorical(LatLongTAB.Site);
A = 200;
latitude = LatLongTAB.Latitude;
longitude = LatLongTAB.Longitude;
gm = geoscatter(latitude,longitude,A,'.','k');  
% text(latitude(1),longitude(1)-0.5,'CB','HorizontalAlignment','right','FontSize',16);
% text(latitude(2),longitude(2)+14,'Kona','HorizontalAlignment','right','FontSize',16);
% text(latitude(3),longitude(3)-0.5,'QC','HorizontalAlignment','right','FontSize',16);
% text(latitude(5),longitude(5)-0.5,'PHR','HorizontalAlignment','right','FontSize',16);
% text(latitude(4),longitude(4)-0.5,'Palmyra','HorizontalAlignment','right','FontSize',16);
% text(latitude(6),longitude(6),'PT','HorizontalAlignment','right','FontSize',16);
% text(latitude(7),longitude(7)+7.5,'QN','HorizontalAlignment','right','FontSize',16);
% text(latitude(8),longitude(8)-1.5,'Saipan','HorizontalAlignment','right','FontSize',16);
% text(latitude(9)+2,longitude(9)-0.5,'Tinian','HorizontalAlignment','right','FontSize',16);
% text(latitude(10),longitude(10)+10,'Wake','HorizontalAlignment','right','FontSize',16);
% text(latitude(11),longitude(11),'PS','HorizontalAlignment','right','FontSize',16);
% text(latitude(12),longitude(12)-0.5,'DCPP','HorizontalAlignment','right','FontSize',16);
% text(latitude(13),longitude(13)-0.5,'CCE','HorizontalAlignment','right','FontSize',16);
% text(latitude(14),longitude(14)-0.5,'Hoke','HorizontalAlignment','right','FontSize',16);
% text(latitude(15),longitude(15)-0.5,'CORC','HorizontalAlignment','right','FontSize',16);
% text(latitude(16),longitude(16)-0.5,'GI','HorizontalAlignment','right','FontSize',16);
% text(latitude(17),longitude(17)-0.5,'GofCA','HorizontalAlignment','right','FontSize',16);
% text(latitude(18),longitude(18)-0.5,'Kauai','HorizontalAlignment','right','FontSize',16);
% text(latitude(19),longitude(19)-0.5,'LSM','HorizontalAlignment','right','FontSize',16);
% text(latitude(20),longitude(20)-0.5,'KR','HorizontalAlignment','right','FontSize',16);
% text(latitude(21),longitude(21)-0.5,'EQ','HorizontalAlignment','right','FontSize',16);
% text(latitude(22),longitude(22)-0.5,'Pagan','HorizontalAlignment','right','FontSize',16);
% text(latitude(23),longitude(23)-0.5,'AB','HorizontalAlignment','right','FontSize',16);
% text(latitude(24),longitude(24)-0.5,'BD','HorizontalAlignment','right','FontSize',16);
% text(latitude(25),longitude(25)-0.5,'KOA','HorizontalAlignment','right','FontSize',16);
% text(latitude(26),longitude(26)-0.5,'KS','HorizontalAlignment','right','FontSize',16);
geolimits([-3 65],[-220 -120]);
set(gcf,'Color','w');
save('Site_map_Chapter4.png');
export_fig Site_map_Chapter4.png