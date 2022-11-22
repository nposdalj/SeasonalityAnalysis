clearvars
close all; clear all;clc;

% Step 2 - Groups data into 5-min bins and daily data accounting for effort(including duty cycle)...
% deals with duty cycle by either taking the proportion of presence based on effort or normalizes...
    % presence based on effort (usually these are pretty similar)
    
% Saves everything as workspaceStep2 to be used for subsequent steps

%IMPORTANT OUTPUTS:
% hourlyTab - data binned in 1-hr bins
% dayTable - data binned in daily bins
% meanTab365 - data averaged by Julian days
%.csv Outputs for future statistical analysis in R
    % bin Data - '*_binData_forGAMGEE.csv'
    % day Data = '*_dayData_forGLMR125.csv'
    % mean of julian day - '*_days365GroupedMean_forGLMR125.csv'
%% Parameters defined by user
%Site names and data paths
filePrefix = 'CB'; % File name to match. 
siteabrev = 'CB'; %abbreviation of site.
region = 'GofAK'; %region
sp = 'Pm'; % your species code
GDrive = 'I'; %Google Drive
saveDir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
DutyCy = 1; %if this data only has 1 deployment that is duty cycled make it equal to 1 otherwise, make it equal to the number...
% of deployments that have different duty cycles that must be accounted for; if this data is NOT duty cycled,...
% or if the entire deployment is duty cycled, make it equal to 0
%% load workspace
GDrive_correct = GDrive; % Preserve correct GDrive as it was entered above
load([saveDir,'\',siteabrev,'_workspace125.mat']);
GDrive = 'I'; %Correct GDrive for SWAL1

% Overwrite some path names
GDrive = GDrive_correct; %Correct GDrive if overwritten by loading workspace
effortXls(1) = GDrive;
saveDir(1) = GDrive;
tpwsPath(1) = GDrive;
%% Set up duty cycled dates
% If only one or two deployments are duty cycled, adjust accordingly
if DutyCy == 1
    if strcmp(siteabrev,'QC')
        startTime = datetime(2007,07,05); %QC
        endTime = datetime(2008,06,15); %QC
    elseif strcmp(siteabrev,'CSM')
        startTime = datetime(2005,04,26); %CSM
        endTime = datetime(2006,05,11); %CSM
    elseif strcmp(siteabrev,'Palmyra')
        startTime = datetime(2006,10,19); %Palmyra
        endTime = datetime(2009,11,12); %Palmyra
    elseif strcmp(siteabrev,'BD')
        startTime = datetime(2011,05,31); 
        endTime = datetime(2012,08,11);
    elseif strcmp(siteabrev,'CB')
        startTime = datetime(2012,05,03);
        endTime = datetime(2013,02,21);
    end
elseif DutyCy == 3
    if strcmp(siteabrev,'PS1')
        startTime1 = datetime(2006,10,03); %PS1
        endTime1 = datetime(2007,06,03); %PS1
        startTime2 = datetime(2007,07,19); %PS1
        endTime2 = datetime(2008,07,08); %PS1
        startTime3 = datetime(2012,07,03); %PS1
        endTime3 = datetime(2012,08,26); %PS1
    elseif strcmp(siteabrev,'Tinian')
        startTime1 = datetime(2011,04,13); %Tinian
        endTime1 = datetime(2011,11,22); %Tinian
        startTime2 = datetime(2012,06,23); %Tinian
        endTime2 = datetime(2013,05,14); %Tinian
        startTime3 = datetime(2013,07,23); %Tinian
        endTime3 = datetime(2019,05,12); %Tinian
    elseif  strcmp(siteabrev,'Wake')
        startTime1 = datetime(2010,01,31); %Wake
        endTime1 = datetime(2010,04,25); %Wake
        startTime2 = datetime(2011,03,25); %Wake
        endTime2 = datetime(2011,05,27); %Wake
        startTime3 = datetime(2012,02,25); %Wake
        endTime3 = datetime(2017,10,28); %Wake
    end
elseif DutyCy == 2 %Kauai
    startTime1 = datetime(2009,10,08); %Kauai
    endTime1 = datetime(2010,05,13); %Kauai
    startTime2 = datetime(2016,07,09); %Kauai
    endTime2 = datetime(2017,08,09); %Kauai
elseif DutyCy == 4 %PHR
    if strcmp(siteabrev,'PHR')
        startTime1 = datetime(2009,10,20);%PHR
        endTime1 = datetime(2010,05,23);%PHR
        startTime2 = datetime(2011,08,15);%PHR
        endTime2 = datetime(2012,01,07);%PHR
        startTime3 = datetime(2014,09,12);%PHR
        endTime3 = datetime(2015,07,16);%PHR
        startTime4 = datetime(2015,10,15);%PHR
        endTime4 = datetime(2019,06,10);%PHR
    elseif strcmp(siteabrev,'Saipan')
        startTime1 = datetime(2010,03,05);%Saipan
        endTime1 = datetime(2010,08,25);%Saipan
        startTime2 = datetime(2011,04,27);%Saipan
        endTime2 = datetime(2011,10,20);%Saipan
        startTime3 = datetime(2012,06,20);%Saipan
        endTime3 = datetime(2013,03,08);%Saipan
        startTime4 = datetime(2013,07,23);%Saipan
        endTime4 = datetime(2019,02,01);%Saipan
    end
elseif DutyCy == 6 %PS2
    startTime1 = datetime(2008,08,04);
    endTime1 = datetime(2009,01,06);
    startTime2 = datetime(2009,02,01);
    endTime2 = datetime(2009,04,30);
    startTime3 = datetime(2009,05,01);
    endTime3 = datetime(2009,09,22);
    startTime4 = datetime(2009,09,23);
    endTime4 = datetime(2010,01,06);
    startTime5 = datetime(2010,02,26);
    endTime5 = datetime(2010,11,02);
    startTime6 = datetime(2011,06,21);
    endTime6 = datetime(2012,04,07);
elseif DutyCy == 5 %Kona
    startTime1 = datetime(2009,04,23);
    endTime1 = datetime(2009,08,18);
    startTime2 = datetime(2009,12,20);
    endTime2 = datetime(2010,03,05);
    startTime3 = datetime(2010,05,01);
    endTime3 = datetime(2010,06,16);
    startTime4 = datetime(2010,09,30);
    endTime4 = datetime(2011,10,22);
    startTime5 = datetime(2013,10,22);
    endTime5 = datetime(22014,03,25);
end
%% Remove bin data with less than 5 clicks -- this should be happening in Step 1, but just in case
binData(binData.Count < 5,:) = []; %identify any bins with less than 5 clicks and delete them
%% group data by 5min bins, hourly, days, weeks, and seasons 
% STILL DOESN'T INCLUDE EFFORT FOR DUTY CYCLE
% BINS (5-min)
binTable = synchronize(binData,binEffort);
binTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
binTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
binTable.maxPP = [];
binidx1 = (binTable.Count >= 5); %identify any bins with less than 5 clicks and delete them (another safety)
[y,~]=size(binTable);
%Add binary data
binTable.PreAbs = zeros(y,1);
binTable.PreAbs(binidx1) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 5min bin

% HOURLY
Click = retime(binData(:,1),'hourly','sum'); % #click per day
Bin = retime(binData(:,1),'hourly','count'); % #bin per day
hourData = synchronize(Click,Bin);
hourlyEffort = retime(binEffort,'hourly','sum');
hourlyTab = synchronize(hourData,hourlyEffort);
hourlyTab.Properties.VariableNames{'bin'} = 'Effort_Bin';
hourlyTab.Properties.VariableNames{'sec'} = 'Effort_Sec';
hourlyTab(~hourlyTab.Effort_Bin,:)=[]; %removes days with no effort, NOT days with no presence (where there are big gaps in data between deployments
%Add binary data
binidx_hourly = (hourlyTab.Count_Bin >= 1);
[y,~]=size(hourlyTab);
hourlyTab.PreAbs = zeros(y,1);
hourlyTab.PreAbs(binidx_hourly) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 1 hour bin

% DAILY
Click = retime(binData(:,1),'daily','sum'); % #click per day
Bin = retime(binData(:,1),'daily','count'); % #bin per day
dayData = synchronize(Click,Bin);
dayEffort = retime(binEffort,'daily','sum');
dayTab = synchronize(dayData,dayEffort);
dayTable = synchronize(dayData,dayEffort);
dayTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
dayTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
dayTableZeros = dayTable;
dayTable(~dayTable.Effort_Bin,:)=[]; %removes days with no effort, NOT days with no presence
%% Accounting for the duty cycle and effort (hourly data)
[p,~]=size(hourlyTab);
hourlyTab.MaxEffort_Bin = ones(p,1)*(12); %total number of 5-min bins possible in one day
hourlyTab.MaxEffort_Sec = ones(p,1)*(3600); %seconds in an hour

if DutyCy == 1
    DutyCycleIdxStart =  hourlyTab.tbin < startTime;
    startVal = find(DutyCycleIdxStart == 0,1);
    DutyCycleIdxEnd =  hourlyTab.tbin > endTime;
    endVal = find(DutyCycleIdxEnd == 1,1);
elseif DutyCy == 2
    DutyCycleIdxStart1 =  hourlyTab.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  hourlyTab.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  hourlyTab.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  hourlyTab.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
elseif DutyCy == 3
    DutyCycleIdxStart1 =  hourlyTab.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  hourlyTab.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  hourlyTab.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  hourlyTab.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  hourlyTab.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  hourlyTab.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
elseif DutyCy == 4
    DutyCycleIdxStart1 =  hourlyTab.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  hourlyTab.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  hourlyTab.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  hourlyTab.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  hourlyTab.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  hourlyTab.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
    DutyCycleIdxStart4 =  hourlyTab.tbin < startTime4;
    startVal4 = find(DutyCycleIdxStart4 == 0,1);
    DutyCycleIdxEnd4 =  hourlyTab.tbin > endTime4;
    endVal4 = find(DutyCycleIdxEnd4 == 1,1);
elseif DutyCy == 5
    DutyCycleIdxStart1 =  hourlyTab.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  hourlyTab.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  hourlyTab.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  hourlyTab.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  hourlyTab.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  hourlyTab.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
    DutyCycleIdxStart4 =  hourlyTab.tbin < startTime4;
    startVal4 = find(DutyCycleIdxStart4 == 0,1);
    DutyCycleIdxEnd4 =  hourlyTab.tbin > endTime4;
    endVal4 = find(DutyCycleIdxEnd4 == 1,1);
    DutyCycleIdxStart5 =  hourlyTab.tbin < startTime5;
    startVal5 = find(DutyCycleIdxStart5 == 0,1);
    DutyCycleIdxEnd5 =  hourlyTab.tbin > endTime5;
    endVal5 = find(DutyCycleIdxEnd5 == 1,1);
elseif DutyCy == 6
    DutyCycleIdxStart1 =  hourlyTab.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  hourlyTab.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  hourlyTab.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  hourlyTab.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  hourlyTab.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  hourlyTab.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
    DutyCycleIdxStart4 =  hourlyTab.tbin < startTime4;
    startVal4 = find(DutyCycleIdxStart4 == 0,1);
    DutyCycleIdxEnd4 =  hourlyTab.tbin > endTime4;
    endVal4 = find(DutyCycleIdxEnd4 == 1,1);
    DutyCycleIdxStart5 =  hourlyTab.tbin < startTime5;
    startVal5 = find(DutyCycleIdxStart5 == 0,1);
    DutyCycleIdxEnd5 =  hourlyTab.tbin > endTime5;
    endVal5 = find(DutyCycleIdxEnd5 == 1,1);
    DutyCycleIdxStart6 =  hourlyTab.tbin < startTime6;
    startVal6 = find(DutyCycleIdxStart6 == 0,1);
    DutyCycleIdxEnd6 =  hourlyTab.tbin > endTime6;
    endVal6 = find(DutyCycleIdxEnd6 == 1,1);
end

% Adjusting effort for each site
if strcmp(siteabrev,'CSM'); 
    hourlyTab.Effort_Bin(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * (0.2); %for Cross_01 and 02 only, 5 on 20 off (25 minute cycle)-- meaning you're recording 20% (0.2) of the time
    hourlyTab.Effort_Sec(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'LSM'); %all duty cycled
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * .5); %for LSM only, 5 on 5 off (10 minute cycle)-- meaning you're recording 50% (0.5) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Pagan'); %all duty cycled
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * .33); %for Pagan_01 only, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
     hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CORC'); %all duty cycled
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * .5); %for CORC_01 and 02 only, 15 on 15 off (30 minute cycle)-- meaning you're recording 50% (0.5) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'HOKE'); %all duty cycled
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * (5/35)); %for CORC_01 and 02 only, 15 on 15 off (30 minute cycle)-- meaning you're recording 50% (0.5) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CA'); %all duty cycled
    hourlyTab.Effort_Bin = hourlyTab.Effort_Bin * (4/12); %for GofCA10 and 11, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'QC');
    hourlyTab.Effort_Bin(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * (0.33); %for QC06 I evaluated the duty cycle using continous deployments and I should adjust by 33%
    hourlyTab.Effort_Sec(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'BD');
    hourlyTab.Effort_Bin(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * (0.49); %for BD02,  I evaluated the duty cycle using continous deployments and I should adjust by 49%
    hourlyTab.Effort_Sec(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CB');
    hourlyTab.Effort_Bin(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * (0.83); %for CB02,  I evaluated the duty cycle using continous deployments and I should adjust by %
    hourlyTab.Effort_Sec(startVal:endVal) = hourlyTab.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Palmyra'); %all duty cycled
    hourlyTab.Effort_Bin = floor(hourlyTab.Effort_Bin * 0.25); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec = hourlyTab.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Kauai');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.26); %for Kauai01 I evaluated the duty cycle using Kauai02 (continous deployment) and I should adjust by 26%
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.71); %for Kaua05 I evaluated the duty cycle using Kauai02 (continous deployment) and I should adjust by 71%
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PS1');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.33); %for PS1_01 and 02 I evaluated the duty cycle using PS 12 (continous deployment) and I should adjust by 33%
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.25); %for PS1_03 and 04 I evaluated the duty cycle using PS 12 (continous deployment) and I should adjust by 25%
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * (0.5); %for PS1_13 I evaluated the duty cycle using PS 12 (continous deployment) and I should adjust by 50%
    hourlyTab.Effort_Sec(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Tinian');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.25); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.83); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * (0.71); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Wake');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.5); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.83); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * (0.17); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PS2');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.33); %for PS2_05 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 33%
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.5); %for PS2_06 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 25%
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * (0.33); %for PS2_07 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 50%
    hourlyTab.Effort_Sec(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * (0.5); %for PS2_08 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 50%
    hourlyTab.Effort_Sec(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal5:endVal5) = hourlyTab.Effort_Bin(startVal5:endVal5) * (0.2); %for PS2_09 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 50%
    hourlyTab.Effort_Sec(startVal5:endVal5) = hourlyTab.Effort_Bin(startVal5:endVal5) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal6:endVal6) = hourlyTab.Effort_Bin(startVal6:endVal6) * (0.5); %for PS2_11 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 33%
    hourlyTab.Effort_Sec(startVal6:endVal6) = hourlyTab.Effort_Bin(startVal6:endVal6) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Kona');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.33); %for Hawaii05 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 33%
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.42); %for Hawaii08 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 42%
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * (0.20); %for Hawaii09 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 20%
    hourlyTab.Effort_Sec(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * (0.63); %for Hawaii10-11 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 63%
    hourlyTab.Effort_Sec(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal5:endVal5) = hourlyTab.Effort_Bin(startVal5:endVal5) * (0.26); %for Hawaii16 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 26%
    hourlyTab.Effort_Sec(startVal5:endVal5) = hourlyTab.Effort_Bin(startVal5:endVal5) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PHR');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.25); %for PHR01 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 25%
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.62); %%for PHR05 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 62%
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * (0.25); %for PHR08 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 25%
    hourlyTab.Effort_Sec(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * (0.17); %for PHR09-12 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 17%
    hourlyTab.Effort_Sec(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Saipan');
    hourlyTab.Effort_Bin(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * (0.13); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal1:endVal1) = hourlyTab.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * (0.25); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal2:endVal2) = hourlyTab.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * (0.83); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal3:endVal3) = hourlyTab.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    hourlyTab.Effort_Bin(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * (0.71); %No continous data so I linear boosted by the duty cycle
    hourlyTab.Effort_Sec(startVal4:endVal4) = hourlyTab.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
else
end

%% Two ways to account for effort
%PROPORTION OF MINUTES WITH CLICKS
hourlyTab.Minutes = hourlyTab.Count_Bin * 5; %convert bins to minutes
hourlyTab.MinProp = hourlyTab.Minutes./(hourlyTab.Effort_Sec ./ (60 * 60)); %proportion of minutes per hour w/clicks

%NORMALIZING BIN COUNT ACCORDING TO EFFORT
hourlyTab.NormEffort_Bin = hourlyTab.Effort_Bin./hourlyTab.MaxEffort_Bin; %what proportion of the hour (min) was there effort
hourlyTab.NormEffort_Sec = hourlyTab.Effort_Sec./hourlyTab.MaxEffort_Sec; %what proportion of the hour (sec) was there effort
hourlyTab.NormBin = hourlyTab.Count_Bin ./ hourlyTab.NormEffort_Bin; %what would the normalized bin count be given the amount of effort
hourlyTab.NormClick = hourlyTab.Count_Click ./ hourlyTab.NormEffort_Sec; %what would be the normalized click count given the amount of effort
hourlyTab.MinNorm = (hourlyTab.NormBin * 5); %convert the number of 5-min bins to minutes
%% Add seasons + year + day
%Winter starts on January (closest to the real thing, which is Dec. 21st)
hourlyTab.Season = zeros(p,1);
hourlyTab.month = month(hourlyTab.tbin);
summeridxD = (hourlyTab.month == 7  | hourlyTab.month == 8 | hourlyTab.month == 9);
fallidxD = (hourlyTab.month == 10  | hourlyTab.month == 11 | hourlyTab.month == 12);
winteridxD = (hourlyTab.month == 1  | hourlyTab.month == 2 | hourlyTab.month == 3);
springidxD = (hourlyTab.month == 4  | hourlyTab.month == 5 | hourlyTab.month == 6);

%adds the season according to the month the data was collected
hourlyTab.Season(summeridxD) = 1;
hourlyTab.Season(fallidxD) = 2;
hourlyTab.Season(winteridxD) = 3;
hourlyTab.Season(springidxD) = 4;

%add year and day to data
hourlyTab.Year = year(hourlyTab.tbin); 
hourlyTab.day = day(hourlyTab.tbin,'dayofyear');
%% Replace NANs where there was recording effort but no detections (this usually happens in the beginning or end of the recordings)
NANidx = ismissing(hourlyTab(:,{'NormBin'}));
hourlyTab{:,{'Count_Click'}}(NANidx) = 0; %if there was effort, but no detections change the Count_Click column to zero
hourlyTab{:,{'Count_Bin'}}(NANidx) = 0; %if there was effort, but no detections change the Count_Bin column to zero
hourlyTab{:,{'NormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
hourlyTab{:,{'NormClick'}}(NANidx) = 0; %if there was effort, but no detections change the NormClick column to zero
hourlyTab{:,{'MinNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero
hourlyTab{:,{'MinProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero

%save day table to csv for modeling in R
writetable(timetable2table(hourlyTab),[saveDir,'\',siteabrev,'_binData_forGAMGEE.csv']); %save table to .csv to continue stats in R
%% Accounting for the duty cycle and effort (day data)
[p,~]=size(dayTable);
dayTable.MaxEffort_Bin = ones(p,1)*(288); %total number of bins possible in one day
dayTable.MaxEffort_Sec = ones(p,1)*(86400); %seconds in one day

if DutyCy == 1
    DutyCycleIdxStart =  dayTable.tbin < startTime;
    startVal = find(DutyCycleIdxStart == 0,1);
    DutyCycleIdxEnd =  dayTable.tbin > endTime;
    endVal = find(DutyCycleIdxEnd == 1,1);
elseif DutyCy == 2
    DutyCycleIdxStart1 =  dayTable.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  dayTable.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  dayTable.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  dayTable.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
elseif DutyCy == 3
    DutyCycleIdxStart1 =  dayTable.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  dayTable.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  dayTable.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  dayTable.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  dayTable.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  dayTable.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
elseif DutyCy == 4
    DutyCycleIdxStart1 =  dayTable.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  dayTable.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  dayTable.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  dayTable.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  dayTable.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  dayTable.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
    DutyCycleIdxStart4 =  dayTable.tbin < startTime4;
    startVal4 = find(DutyCycleIdxStart4 == 0,1);
    DutyCycleIdxEnd4 =  dayTable.tbin > endTime4;
    endVal4 = find(DutyCycleIdxEnd4 == 1,1);
elseif DutyCy == 5
    DutyCycleIdxStart1 =  dayTable.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  dayTable.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  dayTable.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  dayTable.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  dayTable.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  dayTable.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
    DutyCycleIdxStart4 =  dayTable.tbin < startTime4;
    startVal4 = find(DutyCycleIdxStart4 == 0,1);
    DutyCycleIdxEnd4 =  dayTable.tbin > endTime4;
    endVal4 = find(DutyCycleIdxEnd4 == 1,1);
    DutyCycleIdxStart5 =  dayTable.tbin < startTime5;
    startVal5 = find(DutyCycleIdxStart5 == 0,1);
    DutyCycleIdxEnd5 =  dayTable.tbin > endTime5;
    endVal5 = find(DutyCycleIdxEnd5 == 1,1);
elseif DutyCy == 6
    DutyCycleIdxStart1 =  dayTable.tbin < startTime1;
    startVal1 = find(DutyCycleIdxStart1 == 0,1);
    DutyCycleIdxEnd1 =  dayTable.tbin > endTime1;
    endVal1 = find(DutyCycleIdxEnd1 == 1,1);
    DutyCycleIdxStart2 =  dayTable.tbin < startTime2;
    startVal2 = find(DutyCycleIdxStart2 == 0,1);
    DutyCycleIdxEnd2 =  dayTable.tbin > endTime2;
    endVal2 = find(DutyCycleIdxEnd2 == 1,1);
    DutyCycleIdxStart3 =  dayTable.tbin < startTime3;
    startVal3 = find(DutyCycleIdxStart3 == 0,1);
    DutyCycleIdxEnd3 =  dayTable.tbin > endTime3;
    endVal3 = find(DutyCycleIdxEnd3 == 1,1);
    DutyCycleIdxStart4 =  dayTable.tbin < startTime4;
    startVal4 = find(DutyCycleIdxStart4 == 0,1);
    DutyCycleIdxEnd4 =  dayTable.tbin > endTime4;
    endVal4 = find(DutyCycleIdxEnd4 == 1,1);
    DutyCycleIdxStart5 =  dayTable.tbin < startTime5;
    startVal5 = find(DutyCycleIdxStart5 == 0,1);
    DutyCycleIdxEnd5 =  dayTable.tbin > endTime5;
    endVal5 = find(DutyCycleIdxEnd5 == 1,1);
    DutyCycleIdxStart6 =  dayTable.tbin < startTime6;
    startVal6 = find(DutyCycleIdxStart6 == 0,1);
    DutyCycleIdxEnd6 =  dayTable.tbin > endTime6;
    endVal6 = find(DutyCycleIdxEnd6 == 1,1);
end

% Dealing with duty cycled data for each site
if strcmp(siteabrev,'CSM');
    dayTable.Effort_Bin(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * (0.2); %for Cross_01 and 02 only, 5 on 20 off (25 minute cycle)-- meaning you're recording 20% (0.2) of the time
    dayTable.Effort_Sec(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'LSM'); %all duty cycled
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .5); %for LSM only, 5 on 5 off (10 minute cycle)-- meaning you're recording 50% (0.5) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Pagan'); %all duty cycled
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .33); %for Pagan_01 only, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'BD');
    dayTable.Effort_Bin(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * (0.49); %for BD02,  I evaluated the duty cycle using continous deployments and I should adjust by 49%
    dayTable.Effort_Sec(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CB');
    dayTable.Effort_Bin(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * (0.83); %for CB02,  I evaluated the duty cycle using continous deployments and I should adjust by %
    dayTable.Effort_Sec(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'HOKE'); %all duty cycled
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * (5/35)); %for Pagan_01 only, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CORC'); %all duty cycled
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .5); %for CORC_01 and 02 only, 15 on 15 off (30 minute cycle)-- meaning you're recording 50% (0.5) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'CA'); %all duty cycled
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * .333); %for GofCA10 and 11, 5 on 10 off (15 minute cycle)-- meaning you're recording 33% (0.33) of the time
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'QC');
    dayTable.Effort_Bin(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * (0.33); %for QC06 I evaluated the duty cycle using continous deployments and I should adjust by 33%
    dayTable.Effort_Sec(startVal:endVal) = dayTable.Effort_Bin(startVal:endVal) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Kauai');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.26); %for Kauai01 I evaluated the duty cycle using Kauai02 (continous deployment) and I should adjust by 26%
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.71); %for Kaua05 I evaluated the duty cycle using Kauai02 (continous deployment) and I should adjust by 71%
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PS1');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.33); %for PS1_01 and 02 I evaluated the duty cycle using PS 12 (continous deployment) and I should adjust by 33%
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.25); %for PS1_03 and 04 I evaluated the duty cycle using PS 12 (continous deployment) and I should adjust by 25%
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * (0.5); %for PS1_13 I evaluated the duty cycle using PS 12 (continous deployment) and I should adjust by 50%
    dayTable.Effort_Sec(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Palmyra'); %all duty cycled
    dayTable.Effort_Bin = floor(dayTable.Effort_Bin * 0.25); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec = dayTable.Effort_Bin * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Tinian');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.25); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.83); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * (0.71); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Wake');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.5); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.83); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * (0.17); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Kona');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.33); %for Hawaii05 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 33%
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.42); %for Hawaii08 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 42%
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * (0.20); %for Hawaii09 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 20%
    dayTable.Effort_Sec(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * (0.63); %for Hawaii10-11 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 63%
    dayTable.Effort_Sec(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal5:endVal5) = dayTable.Effort_Bin(startVal5:endVal5) * (0.26); %for Hawaii16 I evaluated the duty cycle using Hawaii17-30 (continous deployment) and I should adjust by 26%
    dayTable.Effort_Sec(startVal5:endVal5) = dayTable.Effort_Bin(startVal5:endVal5) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PS2');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.33); %for PS2_05 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 33%
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.5); %for PS2_06 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 25%
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * (0.33); %for PS2_07 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 50%
    dayTable.Effort_Sec(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * (0.5); %for PS2_08 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 50%
    dayTable.Effort_Sec(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal5:endVal5) = dayTable.Effort_Bin(startVal5:endVal5) * (0.2); %for PS2_09 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 50%
    dayTable.Effort_Sec(startVal5:endVal5) = dayTable.Effort_Bin(startVal5:endVal5) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal6:endVal6) = dayTable.Effort_Bin(startVal6:endVal6) * (0.5); %for PS2_11 I evaluated the duty cycle using PS 14-15 (continous deployment) and I should adjust by 50%
    dayTable.Effort_Sec(startVal6:endVal6) = dayTable.Effort_Bin(startVal6:endVal6) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'PHR');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.25); %for PHR01 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 25%
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.62); %%for PHR05 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 62%
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * (0.25); %for PHR08 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 25%
    dayTable.Effort_Sec(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * (0.17); %for PHR09-12 I evaluated the duty cycle using PHR02 and 04 (continous deployment) and I should adjust by 17%
    dayTable.Effort_Sec(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
elseif strcmp(siteabrev,'Saipan');
    dayTable.Effort_Bin(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * (0.13); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal1:endVal1) = dayTable.Effort_Bin(startVal1:endVal1) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * (0.25); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal2:endVal2) = dayTable.Effort_Bin(startVal2:endVal2) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * (0.83); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal3:endVal3) = dayTable.Effort_Bin(startVal3:endVal3) * 5 * 60; %convert from bins into efforts in seconds per day
    dayTable.Effort_Bin(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * (0.71); %No continous data so I linear boosted by the duty cycle
    dayTable.Effort_Sec(startVal4:endVal4) = dayTable.Effort_Bin(startVal4:endVal4) * 5 * 60; %convert from bins into efforts in seconds per day
else    
end
%% Two ways to account for effort
%PROPORTION OF MINUTES WITH CLICKS
dayTable.Minutes = dayTable.Count_Bin * 5; %convert bins to minutes
dayTable.Hours = (dayTable.Count_Bin * 5) ./ 60; %convert the number of bins sperm whales were detected in to hours per day
dayTable.HoursProp = dayTable.Hours./(dayTable.Effort_Sec ./ (60 * 60)); %proportion of hours per day w/clicks

%NORMALIZING BIN COUNT BASED ON EFFORT
dayTable.NormEffort_Bin = dayTable.Effort_Bin./dayTable.MaxEffort_Bin; %what proportion of the day was there effort
dayTable.NormEffort_Sec = dayTable.Effort_Sec./dayTable.MaxEffort_Sec; %what proportion of the day was there effort
dayTable.NormBin = dayTable.Count_Bin ./ dayTable.NormEffort_Bin; %what would the normalized bin count be given the amount of effort
dayTable.NormClick = dayTable.Count_Click ./ dayTable.NormEffort_Sec; %what would be the normalized click count given the amount of effort
dayTable.HoursNorm = (dayTable.NormBin * 5) ./ 60; %convert the number of 5-min bins per day to hours
%% Add seasons + year + day
%Winter starts on January (closest to the real thing, which is Dec. 21st)
dayTable.Season = zeros(p,1);
dayTable.month = month(dayTable.tbin);
summeridxD = (dayTable.month == 7  | dayTable.month == 8 | dayTable.month == 9);
fallidxD = (dayTable.month == 10  | dayTable.month == 11 | dayTable.month == 12);
winteridxD = (dayTable.month == 1  | dayTable.month == 2 | dayTable.month == 3);
springidxD = (dayTable.month == 4  | dayTable.month == 5 | dayTable.month == 6);

%adds the season according to the month the data was collected
dayTable.Season(summeridxD) = 1;
dayTable.Season(fallidxD) = 2;
dayTable.Season(winteridxD) = 3;
dayTable.Season(springidxD) = 4;

%add year and day to data
dayTable.Year = year(dayTable.tbin); 
dayTable.day = day(dayTable.tbin,'dayofyear');
%% Replace NANs where there was recording effort but no detections (this usually happens in the beginning or end of the recordings)
NANidx = ismissing(dayTable(:,{'NormBin'}));
dayTable{:,{'Count_Click'}}(NANidx) = 0; %if there was effort, but no detections change the Count_Click column to zero
dayTable{:,{'Count_Bin'}}(NANidx) = 0; %if there was effort, but no detections change the Count_Bin column to zero
dayTable{:,{'NormBin'}}(NANidx) = 0; %if there was effort, but no detections change the NormBin column to zero
dayTable{:,{'NormClick'}}(NANidx) = 0; %if there was effort, but no detections change the NormClick column to zero
dayTable{:,{'HoursNorm'}}(NANidx) = 0; %if there was effort, but no detections change the HoursNorm column to zero
dayTable{:,{'HoursProp'}}(NANidx) = 0; %if there was effort, but no detections change the HoursProp column to zero

%save day table to csv for R
writetable(timetable2table(dayTable),[saveDir,'\',siteabrev,'_dayData_forGLMR125.csv']); %save table to .csv to continue stats in R
%% Day table with days grouped together (summed and averaged)
[MD,~] = findgroups(dayTable.day);

if length(MD) < 365
    meantab365 = table(dayTable.day(:), dayTable.HoursProp(:));
    meantab365.Properties.VariableNames = {'day' 'HoursProp'};
    sumtab365 = meantab365;
    meantab365.HoursProp(isnan(meantab365.HoursProp)) = 0;
    dayTable.HoursProp(isnan(dayTable.HoursProp)) = 0; 
else
dayTable.day = categorical(dayTable.day);
meantab365 = grpstats(timetable2table(dayTable),'day',{'mean','sem','std','var','range'},'DataVars',{'HoursProp','HoursNorm','NormBin','NormClick'}); %takes the mean of each day of the year
end

writetable(meantab365, [saveDir,'\',siteabrev,'_days365GroupedMean_forGLMR125.csv']); %table with the mean for each day of the year
%% Save workspace variable
save([saveDir,'\',siteabrev,'_workspaceStep2.mat']);