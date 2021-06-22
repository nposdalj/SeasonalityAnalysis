
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

ALEUTBD_latLongs = [52.633333,-175.63; %02
52.076,-175.64]; %03

GCAT_latLongs = [28.01012222, -112.01; %GofCA07-T
    28.01010722, -112.01]; %GofCA08-T

GCAP_latLongs = [23.8305, -109.63; %1
23.82986667, -109.63; %2
23.8304, -109.63; %3
23.92633333, -109.68; %4
23.82741667, -109.63; %4
23.82413333, -109.63]; %6  

GCACB_latLongs = [29.02805, -113.38; %10
29.02753333, -113.38]; %11

MC_latLongs = [28+(50.746/60),-(88+(27.927/60)); %01
    28+(50.771/60),-(88+(27.907/60)); %02
    28+(50.775/60),-(88+(27.909/60)); %03
    28+(50.775/60),-(88+(27.946/60)); %04_disk01
    28+(50.775/60),-(88+(27.946/60)); %04_disk03-08
    28+(50.797/60),-(88+(27.991/60)); %05
    28+(50.853/60),-(88+(28.041/60)); %06
    28+(50.781/60),-(88+(28.059/60)); %07_disk01-03
    28+(50.781/60),-(88+(28.059/60)); %07_disk05-14
    28+(58.850/60),-(88+(28.101/60)); %09
    28+(58.732/60),-(88+(28.082/60)); %10
    28+(58.724/60),-(88+(28.079/60)); %11
    28+(58.727/60),-(88+(28.084/60)); %12
    28+(50.832/60),-(88+(27.973/60))]; %13

GC_latLongs = [27+(33.470/60),-(91+(10.010/60)); %01
    27+(33.466/60),	-(91+(10.014/60)); %02
    27+(33.424/60),	-(91+(10.073/60)); %03
    27+(33.426/60),	-(91+(10.060/60)); %04
    27+(33.440/60),	-(91+(10.562/60)); %05
    27+(33.347/60),	-(91+(10.092/60)); %06_disk01-12
    27+(33.347/60),	-(91+(10.092/60)); %06_disk14-15
    27+(33.347/60),	-(91+(10.092/60)); %07
    27+(33.366/60),	-(91+(10.073/60)); %08
    27+(33.364/60),	-(91+(10.096/60)); %09
    27+(33.367/60),	-(91+(10.083/60)); %10
    27+(33.413/40),-(91+(10.339/60))]; %11

DT_latLongs = [25+(31.911/60),-(84+(38.251/60));  %01
25+(31.911/60),-(84+(38.251/60)); %02
25+(31.859/60),-(84+(38.262/60)); %03
25+(31.867/60),-(84+(38.265/60)); %04
25+(31.938/60),-(84+(38.041/60)); %05
25+(31.941/60),-(84+(38.046/60)); %06
25+(32.219/60),-(84+(38.150/60)); %07
25+(32.320/60),-(84+(37.873/60)); %08
25+(32.316/60),-(84+(37.878/60)); %09
25+(32.360/60),-(84+(37.743/60)); %10
25+(32.355/60),-84+(37.733/60)]; %11

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

CSM_latLongs = [18.72208333, -158.23; %1-A
18.72238333, -158.221; %2-A
18.72223333, -158.224]; %3-A

KAU_latLongs = [21.95273333, -159.89; %1-A
21.95373333, -159.89; %2-A
21.9492, -159.89]; %5-A

LSM_latLongs = [28.61335, -176.70]; %1

PAL_WT_latLongs = [5.8641, -162.17; %2
5.8647, -162.16; %3
5.86545, -162.16; %4
5.8651, -162.16; %5
5.86295, -162.16]; %6

PAL_NS_latLongs = [5.9042, -162.04; %7
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

HZ_latLongs = [41.06191667, -66.35; %1
41.06183333, -66.35; %2
41.06165, -66.35; %3
41.06165, -66.35]; %4

NC_latLongs = [39.83248333, -69.98; %1
39.83238333, -69.98; %2
39.83258333, -69.98; %3
39.83295, -69.98]; %4

BC_latLongs = [39.19105, -72.23; %1
39.1905, -72.23; %2
39.19191667, -72.23]; %3

GS_latLongs = [33.66563333, -76.00; %1
33.66701667, -76.00; %2
33.66991667, -76.00]; %3

BP_latLongs = [32.10603333, -77.09; %1
32.10695, -77.09; %2
32.10526667, -77.09]; %3

BS_latLongs = [30.58378333, -77.39; %1
30.58303333, -77.39; %2
30.58295, -77.39]; %3

NFC_latLongs = [37.16623333, -74.47; %1
37.16651667, -74.47; %2
37.1674, -74.47; %3
37.16451667, -74.47]; %4

HATa_latLongs = [35.34053333, -74.86; %1
35.3406, -74.86; %2
35.34445, -74.86; %3
35.34676667, -74.86; %4
35.34218333, -74.86; %5
35.30183333, -74.86]; %6

HATb_latLongs = [35.58413333, -74.75; %1
35.58351667, -74.74; %3
35.58976667, -74.75; %4
35.5893, -74.75; %5
35.58441667, -74.75]; %6

JAXA_latLongs = [30.2771, -80.22; %1
30.28051667 -80.22; %2
30.27533333 -80.22; %3
30.2682 -80.22; %5
30.27818333 -80.22; %6
30.28501667 -80.22]; %8

JAXB_latLongs = [30.2582, -80.43; %1
30.25916667, -80.43; %4
30.25708333, -80.43; %5
30.25768333, -80.43]; %6

JAXC_latLongs = [30.33286667, -80.20; %9
30.32643333, -80.20]; %10

JAXD_latLongs = [30.1506, -79.77; %11
30.14893333, -79.77; %12
30.15183333, -79.77; %13
30.15268333, -79.77; %14
30.15225, -79.77]; %15

OC_latLongs = [40.2633, -67.99; %1
40.26331667, -67.99; %2
40.26333333, -67.99; %3
40.23, -67.98]; %4

WC_latLongs = [38.37415, -73.37; %1
38.37385, -73.37; %2
38.37336667, -73.37]; %3

PI_latLongs = [72.72433333, -76.23; %1
72.72496667, -76.23; %2
72.72483333, -76.23; %3
72.72283333, -76.23]; %4
%% load lat and long with only 1 deployment
EQ_mean = [0.44345, -164.00]; %1-A
EQtext = repmat({'EQ'},size(EQ_mean,1),1);
EQ = [EQtext num2cell(EQ_mean)];

GI_mean = [29.14103333, -118.26]; %1
GItext = repmat({'GI'},size(GI_mean,1),1);
GI = [GItext num2cell(GI_mean)];

GCABLA_mean = [29.02698333, 113.38]; %9
GCABLAtext = repmat({'GCABLA'},size(GCABLA_mean,1),1);
GCABLA = [GCABLAtext num2cell(GCABLA_mean)];

AB_mean = [57.513667,-150.53];
ABtext = repmat({'AB'},size(AB_mean,1),1);
AB = [ABtext num2cell(AB_mean)];

KOA_mean = [57.224,-146.50];
KOAtext = repmat({'KOA'},size(KOA_mean,1),1);
KOA = [KOAtext num2cell(KOA_mean)];

CORC_mean = [31.747,-121.38]; 
CORCtext = repmat({'CORC'},size(CORC_mean,1),1);
CORC = [CORCtext num2cell(CORC_mean)];

HOKE_mean = [32.1,-126.91];
HOKEtext = repmat({'HOKE'},size(HOKE_mean,1),1);
HOKE = [HOKEtext num2cell(HOKE_mean)];

CCE_mean = [33.482833,-123.42]; 
CCEtext = repmat({'CCE'},size(CCE_mean,1),1);
CCE = [CCEtext num2cell(CCE_mean)];

DCPP01C_mean = [35.4,-122.44];
DCPP01Ctext = repmat({'DCPP01C'},size(DCPP01C_mean,1),1);
DCPP01C = [DCPP01Ctext num2cell(DCPP01C_mean)];

ALEUT01KS_mean = [52.316783,-178.52]; 
KStext = repmat({'KS'},size(ALEUT01KS_mean,1),1);
KS = [KStext num2cell(ALEUT01KS_mean)];

KR_mean = [6.365133333,-162.2923167];
KRtext = repmat({'KR'},size(KR_mean,1),1);
KR = [KRtext num2cell(KR_mean)];

%find means of sites with multiple deployments
[GCAPlat,GCAPlong] = meanm(GCAP_latLongs(:,1),GCAP_latLongs(:,2));
GCAP_mean = [GCAPlat,GCAPlong];
GCAPtext = repmat({'GCAP'},size(GCAP_mean,1),1);
GCAP = [GCAPtext num2cell(GCAP_mean)];

[OClat,OClong] = meanm(OC_latLongs(:,1),OC_latLongs(:,2));
OC_mean = [OClat,OClong];
OCtext = repmat({'OC'},size(OC_mean,1),1);
OC = [OCtext num2cell(OC_mean)];

[WClat,WClong] = meanm(WC_latLongs(:,1),WC_latLongs(:,2));
WC_mean = [WClat,WClong];
WCtext = repmat({'WC'},size(WC_mean,1),1);
WC = [WCtext num2cell(WC_mean)];

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

[QClat,QClong] = meanm(OCNMSQC_latLongs(:,1),OCNMSQC_latLongs(:,2));
OCNMSQC_mean = [QClat, QClong];
OCNMStext = repmat({'OCNMS'},size(OCNMSQC_mean,1),1);
OCNMS = [OCNMStext num2cell(OCNMSQC_mean)];

[BDlat, BDlong] = meanm(ALEUTBD_latLongs(:,1),ALEUTBD_latLongs(:,2));
ALEUTBD_mean = [BDlat, BDlong];
BDtext = repmat({'BD'},size(ALEUTBD_mean,1),1);
BD = [BDtext num2cell(ALEUTBD_mean)];

[BClat, BClong] = meanm(BC_latLongs(:,1),BC_latLongs(:,2));
BC_mean = [BClat, BClong];
BCtext = repmat({'BC'},size(BC_mean,1),1);
BC = [BCtext num2cell(BC_mean)];

[BPlat, BPlong] = meanm(BP_latLongs(:,1),BP_latLongs(:,2));
BP_mean = [BPlat, BPlong];
BPtext = repmat({'BP'},size(BP_mean,1),1);
BP = [BPtext num2cell(BP_mean)];

[BSlat, BSlong] = meanm(BS_latLongs(:,1),BS_latLongs(:,2));
BS_mean = [BSlat, BSlong];
BStext = repmat({'BS'},size(BS_mean,1),1);
BS = [BStext num2cell(BS_mean)];

[CSMlat, CSMlong] = meanm(CSM_latLongs(:,1),CSM_latLongs(:,2));
CSM_mean = [CSMlat, CSMlong];
CSMtext = repmat({'CSM-A'},size(CSM_mean,1),1);
CSM = [BStext num2cell(CSM_mean)];

[DTlat, DTlong] = meanm(DT_latLongs(:,1),DT_latLongs(:,2));
DT_mean = [DTlat, DTlong];
DTtext = repmat({'DT'},size(DT_mean,1),1);
DT = [DTtext num2cell(DT_mean)];

[GCATlat, GCATlong] = meanm(GCAT_latLongs(:,1),GCAT_latLongs(:,2));
GCAT_mean = [GCATlat, GCATlong];
GCATtext = repmat({'GofCA-T'},size(GCAT_mean,1),1);
GCAT = [GCATtext num2cell(GCAT_mean)];

[GSlat, GSlong] = meanm(GS_latLongs(:,1),GS_latLongs(:,2));
GS_mean = [GSlat, GSlong];
GStext = repmat({'GS'},size(GS_mean,1),1);
GS = [GStext num2cell(GS_mean)];

[GClat, GClong] = meanm(GC_latLongs(:,1),GC_latLongs(:,2));
GC_mean = [GClat, GClong];
GCtext = repmat({'GC'},size(GC_mean,1),1);
GC = [GCtext num2cell(GC_mean)];

[HATalat, HATalong] = meanm(HATa_latLongs(:,1),HATa_latLongs(:,2));
HATa_mean = [HATalat, HATalong];
HATatext = repmat({'HAT-A'},size(HATa_mean,1),1);
HATa = [HATatext num2cell(HATa_mean)];

[HATblat, HATblong] = meanm(HATb_latLongs(:,1),HATb_latLongs(:,2));
HATb_mean = [HATblat, HATblong];
HATbtext = repmat({'HAT-B'},size(HATb_mean,1),1);
HATb = [HATbtext num2cell(HATb_mean)];

[HAWAIIKlat, HAWAIIKlong] = meanm(HAWAIIK_latLongs(:,1),HAWAIIK_latLongs(:,2));
HAWAIIK_mean = [HAWAIIKlat, HAWAIIKlong];
HAWAIIKtext = repmat({'Hawaii'},size(HAWAIIK_mean,1),1);
HAWAIIK = [HAWAIIKtext num2cell(HAWAIIK_mean)];

[HZlat, HZlong] = meanm(HZ_latLongs(:,1),HZ_latLongs(:,2));
HZ_mean = [HZlat, HZlong];
HZtext = repmat({'HZ'},size(HZ_mean,1),1);
HZ = [HZtext num2cell(HZ_mean)];

[JAXAlat, JAXAlong] = meanm(JAXA_latLongs(:,1),JAXA_latLongs(:,2));
JAXA_mean = [JAXAlat, JAXAlong];
JAXAtext = repmat({'JAX-A'},size(JAXA_mean,1),1);
JAXA = [JAXAtext num2cell(JAXA_mean)];

[JAXBlat, JAXBlong] = meanm(JAXB_latLongs(:,1),JAXB_latLongs(:,2));
JAXB_mean = [JAXBlat, JAXBlong];
JAXBtext = repmat({'JAX-B'},size(JAXB_mean,1),1);
JAXB = [JAXBtext num2cell(JAXB_mean)];

[JAXClat, JAXClong] = meanm(JAXC_latLongs(:,1),JAXC_latLongs(:,2));
JAXC_mean = [JAXClat, JAXClong];
JAXCtext = repmat({'JAX-C'},size(JAXC_mean,1),1);
JAXC = [JAXCtext num2cell(JAXC_mean)];

[JAXDlat, JAXDlong] = meanm(JAXD_latLongs(:,1),JAXD_latLongs(:,2));
JAXD_mean = [JAXDlat, JAXDlong];
JAXDtext = repmat({'JAX-D'},size(JAXD_mean,1),1);
JAXD = [JAXDtext num2cell(JAXD_mean)];

[KAUlat, KAUlong] = meanm(KAU_latLongs(:,1),KAU_latLongs(:,2));
KAU_mean = [KAUlat, KAUlong];
KAUtext = repmat({'Kauii'},size(KAU_mean,1),1);
KAU = [KAUtext num2cell(KAU_mean)];

[LSMlat, LSMlong] = meanm(LSM_latLongs(:,1),LSM_latLongs(:,2));
LSM_mean = [LSMlat, LSMlong];
LSMtext = repmat({'LSM'},size(LSM_mean,1),1);
LSM = [LSMtext num2cell(LSM_mean)];

[MClat, MClong] = meanm(MC_latLongs(:,1),MC_latLongs(:,2));
MC_mean = [MClat, MClong];
MCtext = repmat({'MC'},size(MC_mean,1),1);
MC = [MCtext num2cell(MC_mean)];

[NClat, NClong] = meanm(NC_latLongs(:,1),NC_latLongs(:,2));
NC_mean = [NClat, NClong];
NCtext = repmat({'NC'},size(NC_mean,1),1);
NC = [NCtext num2cell(NC_mean)];

[NFClat, NFClong] = meanm(NFC_latLongs(:,1),NFC_latLongs(:,2));
NFC_mean = [NFClat, NFClong];
NFCtext = repmat({'NFC'},size(NFC_mean,1),1);
NFC = [NFCtext num2cell(NFC_mean)];

[PALNSlat, PALNSlong] = meanm(PAL_NS_latLongs(:,1),PAL_NS_latLongs(:,2));
PALNS_mean = [PALNSlat, PALNSlong];
PALNStext = repmat({'PAL-NS'},size(PALNS_mean,1),1);
PALNS = [PALNStext num2cell(PALNS_mean)];

[PALWTlat, PALWTlong] = meanm(PAL_WT_latLongs(:,1),PAL_WT_latLongs(:,2));
PALWT_mean = [PALWTlat, PALWTlong];
PALWTtext = repmat({'PAL-WT'},size(PALWT_mean,1),1);
PALWT = [PALWTtext num2cell(PALWT_mean)];

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

[QClat,QClong] = meanm(OCNMSQC_latLongs(:,1),OCNMSQC_latLongs(:,2));
OCNMSQC_mean = [QClat, QClong];
OCNMStext = repmat({'OCNMS'},size(OCNMSQC_mean,1),1);
OCNMS = [OCNMStext num2cell(OCNMSQC_mean)];

[PIlat,PIlong] = meanm(PI_latLongs(:,1),PI_latLongs(:,2));
PI_mean = [PIlat, PIlong];
PItext = repmat({'PI'},size(PI_mean,1),1);
PI = [PItext num2cell(PI_mean)];
%% create one table with all lat and longs
LL = [AB; BC; BD; BP; BS; CB; CCE; CORC; CSM; DCPP01C; DT; EQ; GC; GCABLA; GCAP; GCAT; GS; HATa; HATb; HAWAIIK;...
    HOKE; HZ; JAXA; JAXB; JAXC; JAXD; KAU; KOA; KS; LSM; MC; NC; NFC; OC; OCNMS; PALNS; PALWT; PHR; PT; QN; SAIPAN; TIN; Wake; WC; KR; PI; GI];
LatLong = cell2mat(LL(:,2:3));

LatLongTAB = array2table(LatLong);
LatLongTAB.Properties.VariableNames = {'Latitude' 'Longitude'};

LatLongTAB{:,'Site'} = {'AB'; 'BC'; 'BD'; 'BP'; 'BS'; 'CB'; 'CCE'; 'CORC'; 'CSM'; 'DCPP01C'; 'DT'; 'EQ'; 'GC';...
    'GCABLA';'GCAP'; 'GCAT'; 'GS'; 'HATa'; 'HATb'; 'HAWAIIK'; 'HOKE'; 'HZ'; 'JAXA'; 'JAXB'; 'JAXC'; 'JAXD'; 'KAU';...
    'KOA'; 'KS'; 'LSM'; 'MC'; 'NC'; 'NFC'; 'OC'; 'OCNMS'; 'PALNS'; 'PALWT'; 'PHR'; 'PT'; 'QN'; 'SAIPAN'; 'TIN'; 'Wake'; 'WC'; 'KR'; 'PI'; 'GI'};

LatLongTAB.Longitude(41) = -2.154600000000000e+02; %Saipan
LatLongTAB.Longitude(42) = -2.157500000000000e+02; %Tinian
%% grey site map
figure(3)
LatLongTAB.Site = categorical(LatLongTAB.Site);
A = 200;
latitude = LatLongTAB.Latitude;
longitude = LatLongTAB.Longitude;
gm = geoscatter(latitude,longitude,A,'.','k');  
text(latitude(1),longitude(1)-0.5,'AB','HorizontalAlignment','right','FontSize',16);
text(latitude(2)-1,longitude(2)+4,'BC','HorizontalAlignment','right','FontSize',16);
text(latitude(3),longitude(3)+5,'BD','HorizontalAlignment','right','FontSize',16);
text(latitude(4),longitude(4)-0.5,'BP','HorizontalAlignment','right','FontSize',16);
text(latitude(5),longitude(5)+5,'BS','HorizontalAlignment','right','FontSize',16);
text(latitude(6),longitude(6)-0.5,'CB','HorizontalAlignment','right','FontSize',16);
text(latitude(7),longitude(7)-0.5,'CCE','HorizontalAlignment','right','FontSize',16);
text(latitude(8),longitude(8)+10,'CORC','HorizontalAlignment','right','FontSize',16);
text(latitude(9),longitude(9)-0.5,'CSM','HorizontalAlignment','right','FontSize',16);
text(latitude(10),longitude(10)-0.5,'DCPP','HorizontalAlignment','right','FontSize',16);
text(latitude(11),longitude(11)-0.5,'DT','HorizontalAlignment','right','FontSize',16);
text(latitude(12),longitude(12)-0.5,'EQ','HorizontalAlignment','right','FontSize',16);
text(latitude(13),longitude(13)-0.5,'GC','HorizontalAlignment','right','FontSize',16);
%text(latitude(14),longitude(14)-0.5,'GCAB','HorizontalAlignment','right','FontSize',16);
text(latitude(15),longitude(15)-0.5,'GCAPP','HorizontalAlignment','right','FontSize',16);
text(latitude(16)-1,longitude(16)-0.5,'GCAT','HorizontalAlignment','right','FontSize',16);
text(latitude(17),longitude(17)-0.5,'GS','HorizontalAlignment','right','FontSize',16);
text(latitude(18),longitude(18)-0.5,'HAT','HorizontalAlignment','right','FontSize',16);
text(latitude(20),longitude(20)+8,'HAW','HorizontalAlignment','right','FontSize',16);
text(latitude(21)-1.5,longitude(21)-0.5,'HOKE','HorizontalAlignment','right','FontSize',16);
text(latitude(22)+0.5,longitude(22)-0.5,'HZ','HorizontalAlignment','right','FontSize',16);
text(latitude(23),longitude(23)-0.5,'JAX','HorizontalAlignment','right','FontSize',16);
text(latitude(27),longitude(27)-0.5,'KAU','HorizontalAlignment','right','FontSize',16);
text(latitude(28)+1,longitude(28)+6.5,'KOA','HorizontalAlignment','right','FontSize',16);
text(latitude(29),longitude(29)-0.5,'KS','HorizontalAlignment','right','FontSize',16);
text(latitude(30),longitude(30)-0.5,'LSM','HorizontalAlignment','right','FontSize',16);
text(latitude(31)+0.5,longitude(31)-0.5,'MC','HorizontalAlignment','right','FontSize',16);
text(latitude(32)+0.5,longitude(32)-0.5,'NC','HorizontalAlignment','right','FontSize',16);
text(latitude(33),longitude(33)-0.5,'NFC','HorizontalAlignment','right','FontSize',16);
text(latitude(34)-1,longitude(34)+4.5,'OC','HorizontalAlignment','right','FontSize',16);
text(latitude(35),longitude(35)-0.5,'QC','HorizontalAlignment','right','FontSize',16);
text(latitude(36),longitude(36)-0.5,'PAL','HorizontalAlignment','right','FontSize',16);
text(latitude(38),longitude(38)+7,'PHR','HorizontalAlignment','right','FontSize',16);
text(latitude(39),longitude(39)+4.5,'PT','HorizontalAlignment','right','FontSize',16);
text(latitude(40),longitude(40)-0.5,'QN','HorizontalAlignment','right','FontSize',16);
text(latitude(41)+2,longitude(41)-0.5,'SAI','HorizontalAlignment','right','FontSize',16);
text(latitude(42),longitude(42)+6,'TIN','HorizontalAlignment','right','FontSize',16);
text(latitude(43),longitude(43)-0.5,'WAK','HorizontalAlignment','right','FontSize',16);
text(latitude(44)+0.5,longitude(44)-1,'WC','HorizontalAlignment','right','FontSize',16);
text(latitude(45)+0.5,longitude(45)+5,'KR','HorizontalAlignment','right','FontSize',16);
text(latitude(46),longitude(46)-0.5,'PI','HorizontalAlignment','right','FontSize',16);
text(latitude(47),longitude(47)-0.5,'GI','HorizontalAlignment','right','FontSize',16);
geolimits([-3 74],[-180 -110]);
set(gcf,'Color','w');
save('Site_map.png');
export_fig Site_mapHQ.png