function [edgeffort, depl, site, duty] = DutyCycleInfo
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

CB__duty = [0,0;10,12;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;];
%% Kauai site info

Kauai_dates = [datenum([2009 10 08]), datenum([2010 05 13 01 00 00]); %03
    datenum([2010 06 04]), datenum([2010 08 20 01 43 00]); %02
    datenum([2016 07 09]), datenum([2017 08 09 00 58 15])]; %05

Kauai_depl = [1,2,5];
Kauai_site(1:length(Kauai_depl),1) = {'Kauai'};

Kauai_latLongs = [56.24345, -142.75; %pt01
56.243917, -142.75; %pt02
56.242917, -142.75]; %pt03

Kauai_duty = [5,20;0,0;5,7];
%% CORC

CORC_dates = [datenum([2014 04 30 23 59 59]), datenum([2015 12 04 12 00 00])];
CORC_depl = 2;
CORC_site(1:length(CORC_depl),1) = {'CORC'};
CORC_latLongs = [31.747,-121.38]; 
CORC_duty = [15,30];
%% HOKE
HOKE_dates = [datenum([2008 09 15]), datenum([2009 06 06 14 13 13])];
HOKE_depl = 1;
HOKE_site(1:length(HOKE_depl),1) = {'HOKE'};
HOKE_latLongs = [32.1,-126.91];
HOKE_duty = [5,35];
%% LSM
LSM_dates = [datenum([2009 05 18 00 00 00]), datenum([2009 08 15 10 33 45])];
LSM_depl = 1;
LSM_site(1:length(LSM_depl),1) = {'LSM'};
LSM_latLongs = [32.1,-126.91];
LSM_duty = [5,10];
%% ALEUTBD

ALEUT_dates = [datenum([2010 08 27]), datenum([2011 05 26 08 05 00]); %BD02
    datenum([2011 05 31]), datenum([2012 08 11])]; %BD03

ALEUT_depl = [2,3];
ALEUT_site(1:length(ALEUT_depl),1) = {'BD'};
ALEUT_duty = [0,0;5,10];
%% OCNMS_QC

OCNMSQC_dates = [datenum([2007 07 03]), datenum([2008 06 15]); %6
datenum([2011 01 27]), datenum([2011 10 07]); %12
datenum([2011 12 07]), datenum([2012 07 11]); %14
datenum([2012 09 14]), datenum([2013 06 30]); %15
datenum([2013 07 17]), datenum([2014 05 02])]; %16

OCNMSQC_depl = [6,12,14,15,16];
OCNMSQC_site(1:length(OCNMSQC_depl),1) = {'QC'};
OCNMSQC_duty = [5,3;0,0;0,0;0,0;0,0;];
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
Saipan_duty = [5,40;5,20;5,6;5,7;5,7;5,7;5,7;5,7;5,7;];
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
Tinian_duty = [5,20;5,6;5,7;5,7;5,7;5,7;5,7;];
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
Wake_duty = [5,10;5,6;5,30;5,30;5,30;5,30;5,30;];
%% GofCA
CA_dates = [datenum([2009,4,26,20,00,00]), datenum([2009,09,12,20,06,05]);
    datenum([2009,12,07, 00, 00, 00]), datenum([2010, 05, 18, 18,47,00])];
CA_depl = [10,11];
CA_site(1:length(CA_depl),1) = {'CA'};
CA_duty = [5,15;5,15];
%% PS1
PS1_dates = [datenum([2006 10 03]), datenum([2007 01 16]); %PS01
datenum([2007 01 25]), datenum([2007 06 03]); %PS02
datenum([2007 07 19]), datenum([2007 10 29]); %PS03
datenum([2008 01 25]), datenum([2008 07 08]); %PS04
datenum([2011 11 30]), datenum([2012 06 24]); %PS12
datenum([2012 07 03]), datenum([2012 08 26])]; %PS13

PS1_depl = [1,2,3,4,12,13];

PS1_site(1:length(PS1_depl),1) = {'PS1'};

PS1_duty = [5,15;5,15;5,20;5,20;0,0;5,10];
%% PS2
PS2_dates = [datenum([2008 08 04]), datenum([2009 01 06]); %PS05
datenum([2009 02 01]), datenum([2009 04 30]); %PS06
datenum([2009 05 01]), datenum([2009 09 22]); %PS07
datenum([2009 09 23]), datenum([2010 01 06]); %PS08
datenum([2010 02 26]), datenum([2010 11 02]); %PS09
datenum([2011 06 21]), datenum([2011 11 29]); %PS11
datenum([2018 11 14]), datenum([2019 06 10]); %PS14
datenum([2019 06 11]), datenum([2020 01 25])]; %PS15

PS2_depl = [5,6,7,8,9,11,14,15];

PS2_site(1:length(PS2_depl),1) = {'PS2'};

PS2_duty = [5,15;5,10;5,15;5,10;5,25;5,10;0,0;0,0];
%% Palmyra
PAL_dates = [datenum([2006 10 19 04 00 00]), datenum([2007 03 23 14 00 00]); %PAL02
    datenum([2007 04 09 00 00 00]), datenum([2007 09 20 21 25 00]); %PAL03
    datenum([2008 05 27 12 00 00]), datenum([2008 09 17 09 26 13]); %PAL05
    datenum([2008 10 22 00 00 00]), datenum([2009 04 02 01 23 13]); %PAL06
    datenum([2009 06 02 00 00 00]), datenum([2009 09 27 21 00 00]); %PAL07
    datenum([2009 10 05 00 00 00]), datenum([2009 11 12 08 23 46])]; %PAL08

PAL_depl = [2,3,5,6,7,8];

PAL_site(1:length(PAL_depl),1) = {'Palmyra'};

PAL_duty = [5,20;5,20;5,20;5,20;5,20;5,20];
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

PHR_duty = [5,20;0,0;0,0;5,8;5,20;5,30;5,30;5,30;5,30;];
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

Hawaii_site(1:length(Hawaii_depl),1) = {'Kona'};

Hawaii_duty = [0,0;5,15;0,0;5,12;5,25;5,8;5,8;5,15;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;];
%% concatenate all information
edgeffort = [CB_dates; CORC_dates; HOKE_dates; ALEUT_dates; OCNMSQC_dates;...
    Saipan_dates; Tinian_dates; Wake_dates; CA_dates; PS1_dates; PS2_dates;...
    PAL_dates; Equator_dates; LSM_dates; PHR_dates; Hawaii_dates; Kauai_dates];
depl = [CB_depl';CORC_depl'; HOKE_depl'; ALEUT_depl'; OCNMSQC_depl'; Saipan_depl';...
    Tinian_depl'; Wake_depl'; CA_depl'; PS1_depl'; PS2_depl'; PAL_depl';...
    LSM_depl'; PHR_depl'; Hawaii_depl'; Kauai_depl'];
site = [CB_site; CORC_site; HOKE_site; ALEUT_site; OCNMSQC_site; Saipan_site;...
    Tinian_site; Wake_site; CA_site; PS1_site; PS2_site; PAL_site;...
    LSM_site; PHR_site; Hawaii_site; Kauai_site];
duty = [CB_duty; CORC_duty; HOKE_duty; ALEUT_duty; OCNMSQC_duty; Saipan_duty;...
    Tinian_duty; Wake_duty; CA_duty; PS1_duty; PS2_duty; PAL_duty;...
    LSM_duty; PHR_duty; Hawaii_duty; Kauai_duty];