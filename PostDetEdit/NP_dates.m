function [edgeffort,latLongs, depl, site] = NP_dates
edgeffort = [];
site = {};
%% CB site info
CB_dates = [datenum([2011 07 13]), datenum([2012 02 19]); %cb01
datenum([2012 05 03]), datenum([2013 02 21]); %cb02
datenum([2013 06 06]), datenum([2013 09 05]); %cb03
datenum([2013 09 05]), datenum([2014 04 28]); %cb04
datenum([2014 04 29]), datenum([2014 09 09]); %cb05
datenum([2014 09 09]), datenum([2015 05 02]); %cb06
datenum([2015 05 01]), datenum([2015 09 06]); %cb07
datenum([2017 04 30 12 00 00]), datenum([2017 09 12 12 06 21]); %cb08
datenum([2017 09 14 12 00 00]), datenum([2018 06 16 14 53 48]);%cb09
datenum([2019 04 25 08 00 00]), datenum([2019 09 27 03 27 30])]; %cb10

CB_depl = [1,2,3,4,5,6,7,8,9,10];

CB_site(1:length(CB_depl),1) = {'CB'};

CB_latLongs = [58.645683,-148.07; %01 in decimals
58.671333,-148.02; %02
58.673483,-148.00; %03
58.671867,-148.02; %04
58.671,-148.02; %05
58.670817,-148.02; %06
58.65525,-148.09; %07
58.6695,-148.03; %08
58.6695,-148.0338889; %09
58.66961667,-148.03];%10
%% QN site info

QN_dates = [datenum([2013 09 11 21 00 00]),datenum([2014 04 16 23 50 00]);%qn02 %datenum([2013 06 11]),datenum([2013 09 10]);%qn01
datenum([2014 09 10 18 00 00]),datenum([2015 05 02 18 41 45]);%qn04
datenum([2015 05 02 06 00 00]),datenum([2015 08 18 18 16 15]);%qn05
datenum([2017 04 30 00 00 00]),datenum([2017 09 14 20 38 51])];%qn06

QN_depl = [2,4,5,6];

QN_site(1:length(QN_depl),1) = {'QN'};

QN_latLongs = [56.339383,-145.18; %qn02 56.339017,-145.18; %qn01
56.3413,-145.18; %qn04
56.340683,-145.18; %qn05
56.339833,-145.19]; %qn06
%% AB site info
AB_dates = [datenum([2017 04 29 12 00 00]), datenum([2017 09 14 18 37 36])];%ab0101hp

AB_depl = 1;

AB_site(1:length(AB_depl),1) = {'AB'};

AB_latLongs = [57.513667,-146.50]; %01
%% KOA site info
KOA_dates = [datenum([2019 4 24 19 00 00]), datenum([2019 9 27 19 06 28])];%KOA

KOA_depl = 1;

KOA_site(1:length(KOA_depl),1) = {'KOA'};

KOA_latLongs = [57.224,-150.53]; %01
%% CCE site info
CCE_dates = [datenum([2016 10 10 00 00 00]), datenum([2017 11 08 22 26 21])];%CCE
CCE_depl = 1;
CCE_site(1:length(CCE_depl),1) = {'CCE'};
CCE_latLongs = [57.224,-150.53]; %01
%% PT site info

PT_dates = [datenum([2012 09 09]), datenum([2013 06 10]); %pt01
    datenum([2013 06 11]), datenum([2013 08 20]); %pt01
    datenum([2013 09 03]), datenum([2014 03 21]); %pt01
    datenum([2014 04 30]), datenum([2014 09 10])]; %pt01

PT_depl = [1,2,3,4];
PT_site(1:length(PT_depl),1) = {'PT'};

PT_latLongs = [56.24345, -142.75; %pt01
56.243917, -142.75; %pt02
56.242917, -142.75; %pt03
56.243333, -142.75]; %pt04
%% Kauai site info

Kauai_dates = [datenum([2009 10 08]), datenum([2010 05 13 01 00 00]); %03
    datenum([2010 06 04]), datenum([2010 08 20 01 43 00]); %02
    datenum([2016 07 09]), datenum([2017 08 09 00 58 15])]; %05

Kauai_depl = [1,2,5];
Kauai_site(1:length(Kauai_depl),1) = {'Kauai'};

Kauai_latLongs = [56.24345, -142.75; %pt01
56.243917, -142.75; %pt02
56.242917, -142.75]; %pt03
%% CORC

CORC_dates = [datenum([2014 04 30 23 59 59]), datenum([2015 12 04 12 00 00])];
CORC_depl = 2;
CORC_site(1:length(CORC_depl),1) = {'CORC'};
CORC_latLongs = [31.747,-121.38]; 
%% HOKE
HOKE_dates = [datenum([2008 09 15]), datenum([2009 06 06 14 13 13])];
HOKE_depl = 1;
HOKE_site(1:length(HOKE_depl),1) = {'HOKE'};
HOKE_latLongs = [32.1,-126.91];
%% Equator
Equator_dates = [datenum([2012 03 06 00 00 00]), datenum([2012 06 17 05 35 00])];
Equator_depl = 1;
Equator_site(1:length(Equator_depl),1) = {'Equator'};
Equator_latLongs = [32.1,-126.91];
%% LSM
LSM_dates = [datenum([2009 05 18 00 00 00]), datenum([2009 08 15 10 33 45])];
LSM_depl = 1;
LSM_site(1:length(LSM_depl),1) = {'LSM'};
LSM_latLongs = [32.1,-126.91];
%% DCPP01C
DCPP01C_dates = [datenum([2012 11 07 02 00 00]), datenum([2013 03 19 22 48 45])];
DCPP01C_depl = 1;
DCPP01C_site(1:length(DCPP01C_depl),1) = {'DCPP01C'};
DCPP01C_latLongs = [35.4,-122.44];
%% ALEUTBD

ALEUT_dates = [datenum([2010 08 27]), datenum([2011 05 26 08 05 00]); %BD02
    datenum([2011 05 31]), datenum([2012 08 11])]; %BD03

ALEUT_depl = [2,3];
ALEUT_site(1:length(ALEUT_depl),1) = {'ALEUT'};

ALEUT_latLongs = [52.633333,-175.63; %BD02
52.076,-175.64]; %BD03
%% ALEUT01KS

ALEUT01KS_dates = [datenum([2010 06 03]), datenum([2010 08 25])];

ALEUT01KS_depl = 1;

ALEUT01KS_site(1:length(ALEUT01KS_depl),1) = {'ALEUT01KS'};

ALEUT01KS_latLongs = [52.316783,-178.52]; 
%% OCNMS_QC

OCNMSQC_dates = [datenum([2007 07 03]), datenum([2008 06 15]); %6
datenum([2011 01 27]), datenum([2011 10 07]); %12
datenum([2011 12 07]), datenum([2012 07 11]); %14
datenum([2012 09 14]), datenum([2013 06 30]); %15
datenum([2013 07 17]), datenum([2014 05 02])]; %16

OCNMSQC_depl = [6,12,14,15,16];
OCNMSQC_site(1:length(OCNMSQC_depl),1) = {'OCNMS'};

OCNMSQC_latLongs = [47.466,-125.16; %6
47.50005,-125.36; %12
47.500433,-125.36; %14
47.500533,-125.36; %15
47.500633,-125.65]; %16

%% Baja_GI_01
GI_dates = [datenum([2018,11,19,03,00,00]), datenum([2019,10,22,15,58,45]); %1
    datenum([2019,10,23,18,00,00]), datenum([2020,10,03,13,51,15])]; %2
GI_depl = [1,2];
GI_site(1:length(GI_depl),1) = {'GI'};
GI_latLongs = [29.14103333, -118.26; %1
   29.14243333, - 118.26]; %1
%% Antarc01El
EI_dates = [datenum([2018,11,19]), datenum([2019,10,22])];
EI_depl = [1];
EI_site(1:length(EI_depl),1) = {'Antarc01EI'};
EI_latLongs = [-60.8869, -57.238]; %1

%% Saipan
Saipan_dates = [datenum([2010,03,05]), datenum([2010,08,25]); %1
    datenum([2011,04,27]), datenum([2011,10,20]); %2
    datenum([2012,06,20]), datenum([2013,03,08]); %3
    datenum([2013,07,23]), datenum([2014,01,17]); %4
    datenum([2014,06,18]), datenum([2015,04,17]); %5
    datenum([2015,05,13]), datenum([2016,05,02]); %6
    datenum([2016,05,30]), datenum([2017,05,17]); %7
    datenum([2017,05,29]), datenum([2018,06,02]); %8
    datenum([2018,07,12]), datenum([2019,02,01])]; %9
Saipan_depl = [1,2,3,4,5,6,7,8,9];
Saipan_site(1:length(Saipan_depl),1) = {'Saipan'};
Saipan_latLongs = [15.31663333, -145.46; %1
15.3171, -145.46; %2
15.31778333, -145.46; %3
15.32125, -145.46; %4
15.32125, -145.46; %5
15.31743333, -145.46; %6
15.31671667, -145.46; %7
15.31701667, -145.46; %8
15.31615, -145.46]; %9

%% Tinian
Tinian_dates = [datenum([2011,04,13]), datenum([2011,11,22]); %2
    datenum([2012,06,23]), datenum([2013,05,14]); %3
    datenum([2013,07,23]), datenum([2014,06,15]); %4
    datenum([2014,06,27]), datenum([2014,11,11]); %5
    datenum([2015,05,13]), datenum([2016,05,23]); %6
    datenum([2016,05,30]), datenum([2016,11,05]); %7
    datenum([2018,07,12]), datenum([2019,05,12])]; %9
Tinian_depl = [2,3,4,5,6,7,9];
Tinian_site(1:length(Tinian_depl),1) = {'Tinian'};
Tinian_latLongs = [15.03906667, -145.76; %2
15.0398, -145.76; %3
15.04003333, -145.76; %4
15.04003333, -145.76; %5
15.000674, -145.76; %6
15.04003333, -145.76; %7
15.04003333, -145.76]; %9

%% Wake
Wake_dates = [datenum([2010,01,31]), datenum([2010,05,04]); %1
    datenum([2011,03,25]), datenum([2011,05,27]); %3
    datenum([2012,02,25]), datenum([2013,01,03]); %4
    datenum([2014,06,20]), datenum([2015,05,08]); %5
    datenum([2015,05,05]), datenum([2016,05,24]); %6
    datenum([2016,04,12]), datenum([2016,12,16]); %7
    datenum([2017,05,02]), datenum([2017,10,28])]; %8
Wake_depl = [1,3,4,5,6,7,8];
Wake_site(1:length(Wake_depl),1) = {'Wake'};
Wake_latLongs = [19.22, -166.68; %1
19.2209, -166.69; %3
19.22156667, -166.69; %4
19.22156667, -166.69; %5
19.22346667, -166.69; %6
19.37206667, -166.69; %7
19.37206667, -167.16]; %8

%% OC
OC_dates = [datenum([2015,04,26,00,00,00]),datenum([2016,02,09,06,16,15]);
    datenum([2016,04,26,05,59,59]), datenum([2017,05,18,06,37,35]);
    datenum([2017,07,09,00,00,00]),datenum([2018,04,16,05,56,18]);
    datenum([2018,06,10,06,00,00]),datenum([2019,05,19,04,33,45])];
    OC_depl = [1,2,3,4];
OC_site(1:length(OC_depl),1) = {'OC'};
OC_latLongs = [40.26331667, 67.99;
    40.26331667, 67.99;
    40.26331667, 67.99
    40.26331667, 67.99]; 

%% GS
GS_dates = [datenum([2016,04,29, 00, 00, 00]), datenum([2017, 06, 27, 18, 35, 06]);
    datenum([2017,06,28 00, 00, 00]), datenum([2018,06,26,11,31,21]);
    datenum([2018,06,28,23,59,59]), datenum([2019,06,18,14,17,09])];
GS_depl = [1,2,3];
GS_site(1:length(GS_depl),1) = {'GS'};
GS_latLongs = [39.938, - 76;
40.021, -76;
40.195, - 76];
%% WC
WC_dates = [datenum([2016,04,20, 06, 00, 00]), datenum([2017, 06, 29, 20,57,36]);
    datenum([2017,06,30 00, 00, 00]), datenum([2018,06,02,20,42,36]);
    datenum([2018,06,02,22,00,00]), datenum([2019,05,19,08,32,30])];
WC_depl = [1,2,3];
WC_site(1:length(WC_depl),1) = {'WC'};
WC_latLongs = [39.938, - 76;
40.021, -76;
40.195, - 76];
%% BP
BP_dates = [datenum([2016,04,20, 18, 00, 00]), datenum([2017, 06, 10, 23,04,05]);
    datenum([2017,06,30 12, 00, 00]), datenum([2018,06,03,11,31,21]);
    datenum([2018,06,03,12,00,00]), datenum([2019,05,19,19,30,00])];
BP_depl = [1,2,3];
BP_site(1:length(BP_depl),1) = {'BP'};
BP_latLongs = [39.938, - 76;
40.021, -76;
40.195, - 76];
%% BS
BS_dates = [datenum([2016,04,27, 18, 00, 00]), datenum([2017, 06, 26, 15,22,05]);
    datenum([2017,06,26 18, 00, 00]), datenum([2018,06,23,07,32,33]);
    datenum([2018,06,28,00,00,00]), datenum([2019,06,16,20,13,45])];
BS_depl = [1,2,3];
BS_site(1:length(BS_depl),1) = {'BS'};
BS_latLongs = [30.583, 77.39;
30.583, 77.39;
30.582, 77.39];
%% JAX
JAX_dates = [datenum([2016,04,26, 18, 00, 00]), datenum([2017, 06, 25,19,23,35]);
    datenum([2017,06,25 18,03,57]), datenum([2017,10,28,17,27,48]);
    datenum([2018,06,27,00,00,00]), datenum([2019,06,15,11,03,45])];
JAX_depl = [13,14,15];
JAX_site(1:length(JAX_depl),1) = {'JAX'};
JAX_latLongs = [39.938, - 76;
40.021, -76;
40.195, - 76];
%% BC
BC_dates = [datenum([2016,04,20, 18, 00, 00]), datenum([2017, 06, 10, 23, 04, 05]);
    datenum([2017,06,30 12, 00, 00]), datenum([2018,06,03,11,31,21]);
    datenum([2018,06,03,12,00,00]), datenum([2019,05,19,19,30,00])];
BC_depl = [1,2,3];
BC_site(1:length(BC_depl),1) = {'BC'};
BC_latLongs = [39.938, - 76;
40.021, -76;
40.195, - 76];
%% HZ
HZ_dates = [datenum([2015,6,27,18,00,00]), datenum([2016,03,25,02,21,21]);
    datenum([2016,04,22, 18, 00, 00]), datenum([2017, 06, 19, 07,05,06]);
    datenum([2017,07,09 00, 00, 00]), datenum([2018,01,13,15,25,06]);
    datenum([2018,06,11,17,59,59]), datenum([2019,05,10,06,33,44])];
HZ_depl = [1,2,3,4];
HZ_site(1:length(HZ_depl),1) = {'HZ'};
HZ_latLongs = [39.938, - 76;
40.021, -76;
40.195, - 76;
40.195, - 76];

%% NC
NC_dates = [datenum([2015,4,27,0,00,00]), datenum([2015,09,18,16,28,51]);
    datenum([2016,04,21, 18, 00, 00]), datenum([2017, 05, 24, 14,53,51]);
    datenum([2017,07,16,18, 00, 00]), datenum([2018,06,09,13,02,06]);
    datenum([2018,06,10,00,00,00]), datenum([2019,06,03,04,43,45])];
NC_depl = [1,2,3,4];
NC_site(1:length(NC_depl),1) = {'NC'};
NC_latLongs = [39.938, - 76;
40.021, -76;
40.195, - 76;
40.195, - 76];
%% GofCA
CA_dates = [datenum([2009,4,26,20,00,00]), datenum([2009,09,12,20,06,05]);
    datenum([2009,12,07, 00, 00, 00]), datenum([2010, 05, 18, 18,47,00])];
CA_depl = [10,11];
CA_site(1:length(CA_depl),1) = {'GofCA'};
CA_latLongs = [39.938, - 76;
40.195, - 76];
%% PS
PS_dates = [datenum([2006 10 03]), datenum([2007 01 16]); %PS01
datenum([2007 01 25]), datenum([2007 06 03]); %PS02
datenum([2007 07 19]), datenum([2007 10 29]); %PS03
datenum([2008 01 25]), datenum([2008 07 08]); %PS04
datenum([2008 08 04]), datenum([2009 01 06]); %PS05
datenum([2009 02 01]), datenum([2009 04 30]); %PS06
datenum([2009 05 01]), datenum([2009 09 22]); %PS07
datenum([2009 09 23]), datenum([2010 01 06]); %PS08
datenum([2010 02 26]), datenum([2010 11 02]); %PS09
%datenum([2011 06 21]), datenum([2012 04 07]); %PS11
datenum([2011 06 21]), datenum([2011 11 29]); %PS11
datenum([2011 11 30]), datenum([2012 06 24]); %PS12
datenum([2012 07 03]), datenum([2012 08 26]); %PS13
datenum([2018 11 14]), datenum([2019 06 10]); %PS14
datenum([2019 06 11]), datenum([2020 01 25])]; %PS15

PS_depl = [1,2,3,4,5,6,7,8,9,11,12,13,14,15];

PS_site(1:length(PS_depl),1) = {'PS'};

%DONT FIX THIS ONE
PS_latLongs = [58.645683,-148.07; %01 in decimals
58.671333,-148.02; %02
58.673483,-148.00; %03
58.671867,-148.02; %04
58.671,-148.02; %05
58.670817,-148.02; %06
58.65525,-148.09; %07
58.6695,-148.03; %08
58.6695,-148.0338889; %09
58.66961667,-148.03;
58.66961667,-148.03;
58.66961667,-148.03;
58.66961667,-148.03;
58.66961667,-148.03];%10

%% Palmyra
PAL_dates = [datenum([2006 10 19 04 00 00]), datenum([2007 03 23 14 00 00]); %PAL02
    datenum([2007 04 09 00 00 00]), datenum([2007 09 20 21 25 00]); %PAL03
    datenum([2007 09 22 00 00 00]), datenum([2008 02 13 16 00 21]); %PAL04
    datenum([2008 05 27 12 00 00]), datenum([2008 09 17 09 26 13]); %PAL05
    datenum([2008 10 22 00 00 00]), datenum([2009 04 02 01 23 13]); %PAL06
    datenum([2009 06 02 00 00 00]), datenum([2009 09 27 21 00 00]); %PAL07
    datenum([2009 10 05 00 00 00]), datenum([2009 11 12 08 23 46]); %PAL08
    datenum([2010 06 12 00 00 00]), datenum([2010 08 25 02 08 46]); %PAL09
    datenum([2010 08 27 00 00 00]), datenum([2010 12 09 16 06 15])]; %PAL10

PAL_depl = [2,3,4,5,6,7,8,9,10];

PAL_site(1:length(PAL_depl),1) = {'Palmyra'};

%DONT FIX THIS ONE
PAL_latLongs = [58.645683,-148.07; %02
58.671333,-148.02; %03
58.673483,-148.00; %05
58.671867,-148.02; %05
58.671,-148.02; %06
58.670817,-148.02; %07
58.65525,-148.09; %08
58.6695,-148.03; %09
58.6695,-148.0338889]; %10

%% PHR
PHR_dates = [datenum([2009 10 20 00 00 00]), datenum([2010 05 24 00 42 30]); %PHR01
    datenum([2010 06 01 00 00 00]), datenum([2010 09 17 12 10 00]); %PHR02
    datenum([2011 04 12 00 00 00]), datenum([2011 07 29 10 50 00]); %PHR04
    datenum([2011 08 15 12 00 00]), datenum([2012 01 07 10 03 45]); %PHR05
    datenum([2014 09 12 00 00 00]), datenum([2015 07 16 15 20 00]); %PHR08
    datenum([2015 10 15 00 00 00]), datenum([2016 08 14 02 33 45]); %PHR09
    datenum([2016 08 20 00 00 00]), datenum([2017 03 14 15 03 45]); %PHR10
    datenum([2017 08 20 00 00 00]), datenum([2018 03 05 15 01 20]); %PHR11
    datenum([2018 10 15 00 00 00]), datenum([2019 06 10 08 32 30])]; %PHR12

PHR_depl = [1,2,4,5,8,9,10,11,12];

PHR_site(1:length(PHR_depl),1) = {'PHR'};

%DONT FIX THIS ONE
PHR_latLongs = [58.645683,-148.07; %02
58.671333,-148.02; %03
58.673483,-148.00; %05
58.671867,-148.02; %05
58.671,-148.02; %06
58.670817,-148.02; %07
58.65525,-148.09; %08
58.6695,-148.03; %09
58.6695,-148.0338889]; %10
%% Hawaii
Hawaii_dates = [datenum([2009 02 10]), datenum([2009 04 01 17 19 16]); %05
datenum([2009 04 23]), datenum([2009 08 18 17 48 35]); %06
datenum([2009 10 25]), datenum([2009 12 15 22 16 20]); %07
datenum([2009 12 20]), datenum([2010 03 05 16 03 47]); %08
datenum([2010 05 01]), datenum([2010 06 16 16 56 41]); %09
datenum([2010 09 30]), datenum([2011 03 12 11 55 45]); %10
datenum([2011 05 12]), datenum([2011 10 22 08 18 30]); %11
datenum([2013 10 23]), datenum([2014 03 25 00 00 00]); %16
datenum([2014 03 25]), datenum([2014 07 14 05 57 34]);%17
datenum([2014 07 28]), datenum([2014 10 12 16 41 46]); %18
datenum([2014 12 06]), datenum([2015 03 06 14 34 52]); %19
datenum([2015 04 25]), datenum([2015 08 18 00 38 13]); %20
datenum([2015 11 07]), datenum([2016 02 27 00 00 00]); %22
datenum([2016 07 04]), datenum([2016 09 14 02 45 06]); %23
datenum([2017 07 12]), datenum([2017 10 25 01 08 51]); %26
datenum([2017 10 26]), datenum([2018 04 25 02 51 21]); %27
datenum([2018 04 29]), datenum([2018 11 19 00 00 00]); %28
datenum([2018 11 23]), datenum([2019 03 31 00 00 00]); %29
datenum([2019 04 04]), datenum([2019 09 29 00 00 00])]; %30

Hawaii_depl = [5,6,7,8,9,10,11,16,17,18,19,20,22,23,26,27,28,29,30];

Hawaii_site(1:length(Hawaii_depl),1) = {'Hawaii'};

Hawaii_latLongs = [58.645683,-148.07; %01 in decimals
58.671333,-148.02; %02
58.673483,-148.00; %03
58.671867,-148.02; %04
58.671,-148.02; %05
58.670817,-148.02; %06
58.65525,-148.09; %07
58.6695,-148.03; %08
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.6695,-148.0338889; %09
58.66961667,-148.03];%10
%% concatenate all information
edgeffort = [CB_dates; QN_dates; AB_dates; PT_dates; CORC_dates;...
    HOKE_dates; ALEUT_dates; DCPP01C_dates; ALEUT01KS_dates; OCNMSQC_dates; KOA_dates; GI_dates;...
    EI_dates; Saipan_dates; Tinian_dates; Wake_dates; OC_dates; GS_dates; NC_dates; JAX_dates; HZ_dates; BC_dates;...
    BP_dates; BS_dates; WC_dates; CA_dates; CCE_dates; PS_dates; PAL_dates; Equator_dates; LSM_dates; PHR_dates; Hawaii_dates;...
    Kauai_dates];
latLongs = [CB_latLongs; QN_latLongs; AB_latLongs; PT_latLongs;...
    CORC_latLongs; HOKE_latLongs; ALEUT_latLongs; DCPP01C_latLongs; ALEUT01KS_latLongs; OCNMSQC_latLongs;...
    KOA_latLongs; GI_latLongs; EI_latLongs; Saipan_latLongs; Tinian_latLongs; Wake_latLongs; OC_latLongs; GS_latLongs;...
    NC_latLongs; JAX_latLongs; HZ_latLongs; BC_latLongs; BP_latLongs; BS_latLongs; WC_latLongs; CA_latLongs; CCE_latLongs;...
    PS_latLongs; PAL_latLongs; Equator_latLongs; LSM_latLongs; PHR_latLongs; Hawaii_latLongs; Kauai_latLongs];
depl = [CB_depl'; QN_depl'; AB_depl'; PT_depl'; CORC_depl'; HOKE_depl';...
    ALEUT_depl'; DCPP01C_depl'; ALEUT01KS_depl';OCNMSQC_depl'; KOA_depl';GI_depl'; EI_depl'; Saipan_depl';...
    Tinian_depl'; Wake_depl'; OC_depl'; GS_depl'; NC_depl'; JAX_depl'; HZ_depl'; BC_depl'; BP_depl'; BS_depl'; WC_depl';...
    CA_depl'; CCE_depl'; PS_depl'; PAL_depl';Equator_depl'; LSM_depl'; PHR_depl'; Hawaii_depl'; Kauai_depl'];
site = [CB_site; QN_site; AB_site; PT_site; CORC_site; HOKE_site;...
    ALEUT_site; DCPP01C_site; ALEUT01KS_site; OCNMSQC_site; KOA_site; GI_site; EI_site; Saipan_site;...
    Tinian_site; Wake_site; OC_site; GS_site; NC_site; JAX_site; HZ_site; BC_site; BP_site; BS_site; WC_site;...
    CA_site; CCE_site; PS_site; PAL_site; Equator_site; LSM_site; PHR_site; Hawaii_site; Kauai_site];