clearvars
close all

%This code was modified from the CANARC duty cycle code and was written on 
%07/21/2021 to evaluate the effect of duty cycle on
%the GofAK_CB02 and ALEUT03BD data sets.

%This code was modified again to work with CCE data on 7/11/2022.

%NP

%% Parameters defined by user
filePrefix = 'PS'; % File name to match 
siteabrev = 'PS2'; %abbreviation of site.
sp = 'Pm'; % your species code
srate = 200; % sample rate
region = 'CCE';
GDrive = 'I';
NumSamples = 100; %Number of random duty cycles you want it to test
DutyCont = 0; %If the Duty cycle is at the beginning or end of a deployment
%and it's continuous (1) if only 1 deployment isn't duty cycled and that's 
%the one you want to test on (0)
MinRec = 5; %How many minutes did the duty cycle record for
MinPer = 25; %what was the cycle interval
tpwsPath = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\TPWS_125\',siteabrev]; %directory of TPWS files
dir = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %seasonality analysis directory
effortXls = [GDrive,':\My Drive\',region,'_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\Pm_Effort.xlsx']; % specify excel file with effort times
saveDir = [GDrive,':\My Drive\GofAK_TPWS_metadataReduced\Plots\',siteabrev]; %specify directory to save files
load([dir,'\',siteabrev,'_workspace125.mat']); %load workspace from sumPPICIbin_seasonality code
%% Which deployments are duty cycled
clearvars -except vTT tbin TTall PPall effort er p binEffort DutyCont MinRec MinPer NumSamples
%If the duty cycled data is one after the other
if DutyCont == 1 %duty cycle is continous and you want to remove that single deployment
    startTime = datetime(2006,10,03);
    endTime = datetime(2012,09,13); 
else
    %only one deployment IS NOT duty cycled and that's the one you want to
    %test on
    startTime = datetime(2018,11,14);
    endTime = datetime(2020,1,25); 
end
SecRec = MinRec *60;
SecPer = MinPer * 60;
%% Evaluating the duty cycle by shifting the 15 minute listening period by 1 minute - THIS IS THE ONE I ENDED UP USING
%within the entire 20 minute cycle. This will result in 20 samples.
%group data by 1 second bins
tbin = datetime(vTT);
tbin = dateshift(tbin, 'start','second');
data = timetable(tbin,TTall,PPall);
clear TTall PPall
SecData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
SecData.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
SecData.max_PPall = [];

%Calculate 1-second bin effort
if er > 1
    MinbinEffort = intervalTo1SecBinTimetable(effort.Start,effort.End,p); % convert intervals in bins when there is multiple lines of effort
else
    %binEffort = intervalToBinTimetable_Only1RowEffort(effort.Start,effort.End,p); % convert intervals in bins when there is only one line of effort
end

clickTable = synchronize(SecData,MinbinEffort);
clear MinbinEffort

%Remove duty cycled portion
if DutyCont == 1
    SecData_1 = clickTable;
    SecData_2 = clickTable;
    SecData_1(SecData_1.tbin > startTime,:) = [];
    SecData_2(SecData_2.tbin < endTime,:)=[];
    SecData_cont = [SecData_1;SecData_2];
    [SDC,~] = size(SecData_cont);
    clear SecData_1
    clear SecData_2
else
    SecData_1 = clickTable;
    SecData_1(SecData_1.tbin < startTime | SecData_1.tbin > endTime,:) = [];
    SecData_cont = SecData_1;
    [SDC,~] = size(SecData_cont);
    clear SecData_1
end

%make it divisble by the duty cycle
divis = floor(SDC/SecPer);
ROWS = divis * SecPer;
SecData_contRound = SecData_cont(1:ROWS,:);
[SDCR,~] = size(SecData_contRound);
clear SecData_cont

%choose 100 random samples
R = randperm(SecPer,100); %choose 100 random numbers

%clear extra things to make room
clear SecData clickTable data tbin vTT

All_Clicks = [];
for j = 1:NumSamples
    %make an array of zeros that's the length of one duty cycle 
    Z_array = zeros(SecPer,1); %array of zeros the length of 1 period
    Z_array(R(j),1) = 1; %replace a random value with a 1
    cycle_skeleton = repmat(Z_array, divis,1);%repeat the length of the dataset
    zidx = find(cycle_skeleton);
    PositionsToFill = SecRec;
    for na = 1:numel(zidx)
        anchor = zidx(na);
        cycle_skeleton(anchor:anchor+PositionsToFill-1) = 1;
    end
    cycle_skeleton = cycle_skeleton(1:ROWS,:); %make the duty cycled table equal to the original data table
    
    cycle_skeleton(cycle_skeleton==1)=SecData_contRound.Count(cycle_skeleton == 1);
    cycle_skeleton = array2table(cycle_skeleton);
    cycle_skeleton.Properties.VariableNames{'cycle_skeleton'} = ['Count_Sub',num2str(j)];
    if j > 1
        All_Clicks = [All_Clicks,cycle_skeleton];
    else
        All_Clicks = [SecData_contRound,cycle_skeleton];
    end
end

%Group the duty cycled data back into 5 min bins
All_Clicks_Bin = retime(All_Clicks, 'minutely','sum');
vTT = datevec(All_Clicks_Bin.tbin);
All_Clicks_Bin.tbin = datetime([vTT(:,1:4), floor(vTT(:,5)/p.binDur)*p.binDur, ...
    zeros(length(vTT),1)]);
binEffort_1 = binEffort;
binEffort_2 = binEffort;
binEffort_1(binEffort_1.tbin > startTime,:) = [];
binEffort_2(binEffort_2.tbin < endTime,:)=[];
binEffort_cont = [binEffort_1;binEffort_2];
All_Clicks_Bin_Effort = synchronize(All_Clicks_Bin,binEffort_cont);

%remove rows with no effort
All_Clicks_Bin_Effort(All_Clicks_Bin_Effort.effortSec == 0, :) = [];

%Average # of bins with sperm whales
TableLength = NumSamples+2;
All_Clicks_Bin_Effort.DutyAvg = mean(All_Clicks_Bin_Effort{:,3:TableLength},2);
All_Clicks_Bin_Effort.Diff = All_Clicks_Bin_Effort.Count - All_Clicks_Bin_Effort.DutyAvg; %average number of missed clicks in each bin
%what to multiply the duty cycled data by to supplement so it can look like the continuous data
All_Clicks_Bin_Effort.Supp = All_Clicks_Bin_Effort.Count./All_Clicks_Bin_Effort.DutyAvg; 
%what does the duty cycle percent look like, compared to the actual duty cycle which was recording 43% of the time
All_Clicks_Bin_Effort.DutyPercent = All_Clicks_Bin_Effort.DutyAvg./All_Clicks_Bin_Effort.Count; 
%Multiply the duty cycled average number of clicks in each bin by the
%'supplement' so you can see what number of clicks you'd have if you
%recorded continuously
All_Clicks_Bin_Effort.Adj = All_Clicks_Bin_Effort.DutyAvg .* All_Clicks_Bin_Effort.Supp;
All_Clicks_Bin_Effort.DiffAdj = All_Clicks_Bin_Effort.Count - All_Clicks_Bin_Effort.Adj; %sanity check, what's the difference between the 
%actual recorded number of clicks and the 'adjusted' number of clicks based
%on the supplemented duty cycled data

%Missed Clicks in 5-Minute bins
figure
idx = All_Clicks_Bin_Effort.Diff > 0;
hist(All_Clicks_Bin_Effort.Diff(idx))
title('Histogram of Missed Clicks in 5-Minute Bins')
xlabel('# of Missed Clicks in Each 5-Min Bin')
ylabel('Count')

%Average Duty cycle
MeanBin = nanmean(All_Clicks_Bin_Effort.DutyPercent);
Avg_DutyCycle = ['The average duty cycle for this site was ',num2str(MeanBin)];
disp(Avg_DutyCycle)

%Average # of days with sperm whales 
%retime bin table for daily
columnIndices2Delete = [2 TableLength+1:TableLength+8];
All_BinsINT = All_Clicks_Bin_Effort;
All_BinsINT(:,columnIndices2Delete) = [];
All_Clicks_Bin_Effort_Days = retime(All_BinsINT,'daily','sum');

%recalculate all columns
%All_2016_Days{:,2:end}(All_2016_Days{:,2:end} == 0) = NaN;
All_Clicks_Bin_Effort_Days.DutyAvg = mean(All_Clicks_Bin_Effort_Days{:,2:end},2); %average number of clicks in each bin
All_Clicks_Bin_Effort_Days.Diff = All_Clicks_Bin_Effort_Days.Count - All_Clicks_Bin_Effort_Days.DutyAvg; %average number of missed clicks in each bin
%what to multiply the duty cycled data by to supplement so it can look like the continuous data
All_Clicks_Bin_Effort_Days.Supp = All_Clicks_Bin_Effort_Days.Count./All_Clicks_Bin_Effort_Days.DutyAvg; 
%what does the duty cycle percent look like, compared to the actual duty cycle which was recording 43% of the time
All_Clicks_Bin_Effort_Days.DutyPercent = All_Clicks_Bin_Effort_Days.DutyAvg./All_Clicks_Bin_Effort_Days.Count; 
%Multiply the duty cycled average number of clicks in each bin by the
%'supplement' so you can see what number of clicks you'd have if you
%recorded continuously
All_Clicks_Bin_Effort_Days.Adj = All_Clicks_Bin_Effort_Days.DutyAvg .* All_Clicks_Bin_Effort_Days.Supp;
All_Clicks_Bin_Effort_Days.DiffAdj = All_Clicks_Bin_Effort_Days.Count - All_Clicks_Bin_Effort_Days.Adj; %sanity check, what's the difference between the 
%actual recorded number of clicks and the 'adjusted' number of clicks based
%on the supplemented duty cycled data

%Missed 5-Minute bins each day
figure
idx = All_Clicks_Bin_Effort_Days.Diff > 0;
hist(All_Clicks_Bin_Effort_Days.Diff(idx))
title('Histogram of Missed 5-Minute Bins each day')
xlabel('# of Missed 5-Min Bins Each Day')
ylabel('Count')

%Average Duty cycle
Day_2016 = nanmean(All_Clicks_Bin_Effort_Days.DutyPercent);
Avg_DutyCycle = ['The average duty cycle was ',num2str(Day_2016)];
disp(Avg_DutyCycle)

%QC (5/35) - 0.33
%PS1 (5/15) - 0.33
%PS1 (5/20) - 0.25
%PS1 (5/10) - 0.499 or 0.50
%PS2 (5/15) - 0.33
%PS2 (5/10) - 0.50
%PS2 (5/25) - 0.20

