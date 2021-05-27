clearvars
close all
%% Parameters defined by user
filePrefix = 'Baja_GI'; % File name to match. 
siteabrev = 'GI'; %abbreviation of site.
titleNAME = 'Baja California - Guadalupe Island';
sp = 'Pm'; % your species code
tpwsPath = ['G:\',siteabrev,'\TPWS_125\TPWS2\TPWS3\']; %directory of TPWS files
%% load data from step 3
filename = ['G:\Baja\Seasonality\',siteabrev,'_workspaceStep3'];
load(filename);
saveDir = 'G:\Baja\Seasonality'; %specify directory to save files
%% combine tables into one
meanTAB = array2table([meantabFE365.Day meantabFE365.HoursPropFE meantabJU365.HoursPropJU meantabMA365.HoursPropMA]);
meanTAB.Properties.VariableNames = {'Day' 'HoursPropF' 'HoursPropJ' 'HoursPropM'};
meanTAB.Fem = meanTAB.HoursPropF > 0;
meanTAB.Juv = meanTAB.HoursPropJ > 0;
meanTAB.Mal = meanTAB.HoursPropM > 0;
%% Calculate A, AB, B, BC, C, CA, ABC
F = sum(meanTAB.Fem); %A
J = sum(meanTAB.Juv); %B
M = sum(meanTAB.Mal); %C

FJ = F + J; %AB
FM = F + M; %AC
JM = J + M; %BC
FJM = F + J + M; %ABC

data = [F FJ J JM M FM FJM];
vennX(data,0.01)


