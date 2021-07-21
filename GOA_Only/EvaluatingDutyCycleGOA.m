clearvars
close all

%This code was modified from the CANARC duty cycle code and was written on 
%07/21/2021 to evaluate the effect of duty cycle on
%the GofAK_CB02 and ALEUT03BD data sets.

%NP

%% Parameters defined by user
filePrefix = 'GofAK_CB'; % File name to match 
siteabrev = 'CB'; %abbreviation of site.
sp = 'Pm'; % your species code
itnum = '2'; % which iteration you are looking for
srate = 200; % sample rate
tpwsPath = 'I:\My Drive\GofAK_TPWS_metadataReduced\TPWS_125\CB'; %directory of TPWS files
dir = 'I:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\CB'; %seasonality analysis directory
effortXls = 'I:\My Drive\GofAK_TPWS_metadataReduced\ICIgrams\Pm_Effort_CB.xlsx'; % specify excel file with effort times
saveDir = 'I:\My Drive\GofAK_TPWS_metadataReduced\Plots\CB'; %specify directory to save files
load([dir,'\',siteabrev,'_workspace125.mat']); %load workspace from sumPPICIbin_seasonality code
%% group data by 5min bins, days, weeks, and seasons
%group data by 5 minute bins
binDataIDX = (binData.Count < 5); %remove anything with less than 5 clicks in a bin
binData.Count(binDataIDX) = 0;
binTable = synchronize(binData,binEffort);
binTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
binTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
binTable.maxPP = [];
%binidx1 = (binTable.Count < 5);
%binTable.Count(binidx1) = 0;
[y,~]=size(binTable);
binTable.PreAbs = zeros(y,1);
binidx2 = (binTable.Count >= 5);
binTable.PreAbs(binidx2) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 5min bin
%no effort bins are excluded 
binTable.Year = year(binTable.tbin); %add year
binTable.Minutes = minute(binTable.tbin); %add minutes
binTable.Day = day(binTable.tbin, 'dayofyear'); %add julian day

Click = retime(binData(:,1),'daily','sum'); % #click per day
binDataIDX_zeros = binData.Count == 0;
binData(~binData.Count,:) = [];
Bin = retime(binData(:,1),'daily','count'); % #bin per day

%group data by day
dayData = synchronize(Click,Bin);
dayEffort = retime(binEffort,'daily','sum');
dayTab = synchronize(dayData,dayEffort);
dayTable = synchronize(dayData,dayEffort);
dayTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
dayTable.Properties.VariableNames{'sec'} = 'Effort_Sec';
dayTableZeros = dayTable;
[y,~]=size(dayTable);
dayTable.PreAbs = zeros(y,1);
dayTableidx1 = (dayTable.Count_Click >=5);
dayTable.PreAbs(dayTableidx1) = 1; %table with 0 for no presence in 5min bin and 1 with presence in 5min bin
dayTable(~dayTable.Effort_Bin,:)=[]; %removes days with no effort, NOT days with no presence
dayTable.Year = year(dayTable.tbin); % add year
%% Evaluating the duty cycle by shifting the 15 minute listening period by 1 minute - THIS IS THE ONE I ENDED UP USING
%within the entire 20 minute cycle. This will result in 20 samples.
%group data by 1 second bins
tbin = datetime([vTT(:,1:4), floor(vTT(:,5)), ...
    zeros(length(vTT),1)]); %round to the nearest minute
data = timetable(tbin,TTall,PPall);
MinData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
MinData.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
MinData.Properties.VariableNames{'max_PPall'} = 'maxPP';

%Calculate 1-minute bin effort
if er > 1
    MinbinEffort = intervalTo1MinBinTimetable(effort.Start,effort.End,p); % convert intervals in bins when there is multiple lines of effort
    %binEffort.sec = binEffort.bin*(p.binDur*60);
else
    %binEffort = intervalToBinTimetable_Only1RowEffort(effort.Start,effort.End,p); % convert intervals in bins when there is only one line of effort
    %binEffort.sec = binEffort.bin*(p.binDur*60);
end

clickTable = synchronize(MinData,MinbinEffort); %table with clicks per 1-min bin
clickTable.Properties.VariableNames{'effortBin'} = 'Effort_Bin';
clickTable.Properties.VariableNames{'effortSec'} = 'Effort_Sec';
clickTable.maxPP = [];

clickTable = timetable2table(clickTable);
clickTable.Year = year(clickTable.tbin);

%2016
clickTable2016 = clickTable(find(clickTable.Year == 2016,1,'first'):find(clickTable.Year == 2016,1,'last'),:);
clickTable2016.Effort_Bin = [];
clickTable2016.Effort_Sec = [];
clickTable2016.Year = [];
[cT16,~] = size(clickTable2016);

clickTable2016 = table2timetable(clickTable2016);

%group data in 1 second bins
tbin = datetime(vTT);
tbin = dateshift(tbin, 'start','second');
data = timetable(tbin,TTall,PPall);
SecData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
SecData.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
SecData.Properties.VariableNames{'max_PPall'} = 'maxPP';

%Calculate 1-second bin effort
if er > 1
    MinbinEffort = intervalTo1SecBinTimetable(effort.Start,effort.End,p); % convert intervals in bins when there is multiple lines of effort
    %binEffort.sec = binEffort.bin*(p.binDur*60);
else
    %binEffort = intervalToBinTimetable_Only1RowEffort(effort.Start,effort.End,p); % convert intervals in bins when there is only one line of effort
    %binEffort.sec = binEffort.bin*(p.binDur*60);
end

clickTable = synchronize(SecData,MinbinEffort); %table with clicks per 1-min bin
clickTable.Properties.VariableNames{'effortBin'} = 'Effort_Bin';
clickTable.Properties.VariableNames{'effortSec'} = 'Effort_Sec';
clickTable.maxPP = [];

clickTable = timetable2table(clickTable);
clickTable.Year = year(clickTable.tbin);
clickTable.Effort_Bin = [];
clickTable.Effort_Sec = [];

%2016
clickTable2016 = clickTable(find(clickTable.Year == 2016,1,'first'):find(clickTable.Year == 2016,1,'last'),:);
clickTable2016.Effort_Bin = [];
clickTable2016.Effort_Sec = [];
clickTable2016.Year = [];
[cT16,~] = size(clickTable2016);

clickTable2016 = table2timetable(clickTable2016);

All_2016_Clicks = [];
for j = 1:20
    Sub_2016 = [];
    for i = j:20:cT16-2
        cycleRange = i:i+19;
        columnsToDelete = cycleRange > 312481;
        cycleRange(columnsToDelete) = [];
        dataRange = clickTable2016(cycleRange,:);
        [xx,~] = size(dataRange);
        if xx < 15
            SubS = dataRange;
        else
            SubS = dataRange(1:15,:);
        end
        Sub_2016 = [Sub_2016;SubS];
    end
    Sub_2016.Properties.VariableNames{'Count'} = ['Count_Sub',num2str(j)];
    if j > 1
        All_2016_Clicks = synchronize(All_2016_Clicks,Sub_2016);
    else
        All_2016_Clicks = synchronize(clickTable2016,Sub_2016);
    end
end

binEffort.Year = year(binEffort.tbin); 
binEffort16 = binEffort(find(binEffort.Year == 2016,1,'first'):find(binEffort.Year == 2016,1,'last'),:);

All_2016_Bins = synchronize(binEffort16,All_2016_Clicks,'regular','sum','TimeStep',minutes(5));

%remove bins with less than 5 clicks for continuous data
Count2016IDX = (All_2016_Bins.Count < 5);
All_2016_Bins.Count(Count2016IDX) = 0;

for i = 1:20
    variableName = ['Count_Sub',num2str(i)];
    All2016IDX = (All_2016_Bins.(variableName) < 5); %remove anything with less than 5 clicks in a bin
    All_2016_Bins.(variableName)(All2016IDX) = 0;
end

%Average # of bins with sperm whales
%All_2016_Bins{:,2:end}(All_2016_Bins{:,2:end} == 0) = NaN;
%All_2016_Bins.DutyAvg = mean(All_2016_Bins{:,3:end},2,'omitnan'); %average number of clicks in each bin
All_2016_Bins.DutyAvg = mean(All_2016_Bins{:,3:end},2);
All_2016_Bins.Diff = All_2016_Bins.Count - All_2016_Bins.DutyAvg; %average number of missed clicks in each bin
%what to multiply the duty cycled data by to supplement so it can look like the continuous data
All_2016_Bins.Supp = All_2016_Bins.Count./All_2016_Bins.DutyAvg; 
%what does the duty cycle percent look like, compared to the actual duty cycle which was recording 43% of the time
All_2016_Bins.DutyPercent = All_2016_Bins.DutyAvg./All_2016_Bins.Count; 
%Multiply the duty cycled average number of clicks in each bin by the
%'supplement' so you can see what number of clicks you'd have if you
%recorded continuously
All_2016_Bins.Adj = All_2016_Bins.DutyAvg .* All_2016_Bins.Supp;
All_2016_Bins.DiffAdj = All_2016_Bins.Count - All_2016_Bins.Adj; %sanity check, what's the difference between the 
%actual recorded number of clicks and the 'adjusted' number of clicks based
%on the supplemented duty cycled data

%Missed Clicks in 5-Minute bins
figure
idx = All_2016_Bins.Diff > 0;
hist(All_2016_Bins.Diff(idx))
title('Histogram of Missed Clicks in 5-Minute Bins in 2016')
xlabel('# of Missed Clicks in Each 5-Min Bin')
ylabel('Count')

%Average Duty cycle
Mean_2016 = nanmean(All_2016_Bins.DutyPercent);
Avg2016_DutyCycle = ['The average duty cycle for 2016 was ',num2str(Mean_2016)];
disp(Avg2016_DutyCycle)

%Average # of days with sperm whales 
%retime bin table for daily
columnIndices2Delete = [1 38 39 40 41 42 43];
All_2016_BinsINT = All_2016_Bins;
All_2016_BinsINT(:,columnIndices2Delete) = [];
All_2016_Days = retime(All_2016_BinsINT,'daily','sum');

%recalculate all columns
%All_2016_Days{:,2:end}(All_2016_Days{:,2:end} == 0) = NaN;
All_2016_Days.DutyAvg = mean(All_2016_Days{:,3:end},2); %average number of clicks in each bin
All_2016_Days.Diff = All_2016_Days.Count - All_2016_Days.DutyAvg; %average number of missed clicks in each bin
%what to multiply the duty cycled data by to supplement so it can look like the continuous data
All_2016_Days.Supp = All_2016_Days.Count./All_2016_Days.DutyAvg; 
%what does the duty cycle percent look like, compared to the actual duty cycle which was recording 43% of the time
All_2016_Days.DutyPercent = All_2016_Days.DutyAvg./All_2016_Days.Count; 
%Multiply the duty cycled average number of clicks in each bin by the
%'supplement' so you can see what number of clicks you'd have if you
%recorded continuously
All_2016_Days.Adj = All_2016_Days.DutyAvg .* All_2016_Days.Supp;
All_2016_Days.DiffAdj = All_2016_Days.Count - All_2016_Days.Adj; %sanity check, what's the difference between the 
%actual recorded number of clicks and the 'adjusted' number of clicks based
%on the supplemented duty cycled data

%Missed 5-Minute bins each day
figure
idx = All_2016_Days.Diff > 0;
hist(All_2016_Days.Diff(idx))
title('Histogram of Missed 5-Minute Bins each day in 2016')
xlabel('# of Missed 5-Min Bins Each Day')
ylabel('Count')

%Average Duty cycle
Mean_2016 = nanmean(All_2016_Days.DutyPercent);
Avg2016_DutyCycle = ['The average duty cycle for 2016 was ',num2str(Mean_2016)];
disp(Avg2016_DutyCycle)

%2018/2019 deployment
clickTable20189 = clickTable(find(clickTable.Year == 2018,1,'first'):find(clickTable.Year == 2019,1,'last'),:);
clickTable20189.Effort_Bin = [];
clickTable20189.Effort_Sec = [];
clickTable20189.Year = [];
[cT189,~] = size(clickTable20189);

clickTable20189 = table2timetable(clickTable20189);

All_20189_Clicks = [];
for j = 1:20
    Sub_20189 = [];
    for i = j:20:cT189-2
        cycleRange = i:i+19;
        columnsToDelete = cycleRange > 549204;
        cycleRange(columnsToDelete) = [];
        dataRange = clickTable20189(cycleRange,:);
        [xx,~] = size(dataRange);
        if xx < 15
            SubS = dataRange;
        else
            SubS = dataRange(1:15,:);
        end
        Sub_20189 = [Sub_20189;SubS];
    end
    Sub_20189.Properties.VariableNames{'Count'} = ['Count_Sub',num2str(j)];
    if j > 1
        All_20189_Clicks = synchronize(All_20189_Clicks,Sub_20189);
    else
        All_20189_Clicks = synchronize(clickTable20189,Sub_20189);
    end
end

binEffort.Year = year(binEffort.tbin); 
binEffort189 = binEffort(find(binEffort.Year == 2018,1,'first'):find(binEffort.Year == 2019,1,'last'),:);

All_20189_Bins = synchronize(binEffort189,All_20189_Clicks,'regular','sum','TimeStep',minutes(5));

%remove bins with less than 5 clicks for continuous data
Count20189IDX = (All_20189_Bins.Count < 5);
All_20189_Bins.Count(Count20189IDX) = 0;

for i = 1:20
    variableName = ['Count_Sub',num2str(i)];
    All20189IDX = (All_20189_Bins.(variableName) < 5); %remove anything with less than 5 clicks in a bin
    All_20189_Bins.(variableName)(All20189IDX) = 0;
end

%Average # of bins with sperm whales
All_20189_Bins.DutyAvg = mean(All_20189_Bins{:,3:end},2);
All_20189_Bins.Diff = All_20189_Bins.Count - All_20189_Bins.DutyAvg; %average number of missed clicks in each bin
%what to multiply the duty cycled data by to supplement so it can look like the continuous data
All_20189_Bins.Supp = All_20189_Bins.Count./All_20189_Bins.DutyAvg; 
%what does the duty cycle percent look like, compared to the actual duty cycle which was recording 43% of the time
All_20189_Bins.DutyPercent = All_20189_Bins.DutyAvg./All_20189_Bins.Count; 
%Multiply the duty cycled average number of clicks in each bin by the
%'supplement' so you can see what number of clicks you'd have if you
%recorded continuously
All_20189_Bins.Adj = All_20189_Bins.DutyAvg .* All_20189_Bins.Supp;
All_20189_Bins.DiffAdj = All_20189_Bins.Count - All_20189_Bins.Adj; %sanity check, what's the difference between the 
%actual recorded number of clicks and the 'adjusted' number of clicks based
%on the supplemented duty cycled data

%Missed Clicks in 5-Minute bins
figure
idx = All_20189_Bins.Diff > 0;
hist(All_20189_Bins.Diff(idx))
title('Histogram of Missed Clicks in 5-Minute Bins in 2018/2019')
xlabel('# of Missed Clicks in Each 5-Min Bin')
ylabel('Count')

%Average Duty cycle
Mean_20189 = nanmean(All_20189_Bins.DutyPercent);
Avg20189_DutyCycle = ['The average duty cycle for 2018/2019 was ',num2str(Mean_20189)];
disp(Avg20189_DutyCycle)

%Average # of days with sperm whales 
%retime bin table for daily
columnIndices2Delete = [1 38 39 40 41 42 43];
All_20189_BinsINT = All_20189_Bins;
All_20189_BinsINT(:,columnIndices2Delete) = [];
All_20189_Days = retime(All_20189_BinsINT,'daily','sum');

%recalculate all columns
%All_2016_Days{:,2:end}(All_2016_Days{:,2:end} == 0) = NaN;
All_20189_Days.DutyAvg = mean(All_20189_Days{:,3:end},2); %average number of clicks in each bin
All_20189_Days.Diff = All_20189_Days.Count - All_20189_Days.DutyAvg; %average number of missed clicks in each bin
%what to multiply the duty cycled data by to supplement so it can look like the continuous data
All_20189_Days.Supp = All_20189_Days.Count./All_20189_Days.DutyAvg; 
%what does the duty cycle percent look like, compared to the actual duty cycle which was recording 43% of the time
All_20189_Days.DutyPercent = All_20189_Days.DutyAvg./All_20189_Days.Count; 
%Multiply the duty cycled average number of clicks in each bin by the
%'supplement' so you can see what number of clicks you'd have if you
%recorded continuously
All_20189_Days.Adj = All_20189_Days.DutyAvg .* All_20189_Days.Supp;
All_20189_Days.DiffAdj = All_20189_Days.Count - All_20189_Days.Adj; %sanity check, what's the difference between the 
%actual recorded number of clicks and the 'adjusted' number of clicks based
%on the supplemented duty cycled data

%Missed 5-Minute bins each day
figure
idx = All_20189_Days.Diff > 0;
hist(All_20189_Days.Diff(idx))
title('Histogram of Missed 5-Minute Bins each day in 2018/2019')
xlabel('# of Missed 5-Min Bins Each Day')
ylabel('Count')

%Average Duty cycle
Mean_20189 = nanmean(All_20189_Days.DutyPercent);
Avg20189_DutyCycle = ['The average duty cycle for 2018/2019 was ',num2str(Mean_20189)];
disp(Avg20189_DutyCycle)

