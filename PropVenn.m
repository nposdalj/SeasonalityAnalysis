clearvars
close all
%% Parameters defined by user
filePrefix = 'ALEUT01KS'; % File name to match. 
siteabrev = 'KS'; %abbreviation of site.
titleNAME = 'Aleutians - Kiska Island';
sp = 'Pm'; % your species code
tpwsPath = ['D:\My Drive\GofAK_TPWS_metadataReduced\TPWS_125',siteabrev,'\TPWS_125\TPWS2\TPWS3\']; %directory of TPWS files
%% load data from step 3
filename = ['D:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev,'\',siteabrev,'_workspaceStep3'];
load(filename);
saveDir = ['D:\My Drive\GofAK_TPWS_metadataReduced\SeasonalityAnalysis\',siteabrev]; %specify directory to save files
%% combine tables into one
if strcmp(siteabrev,'KOA') == 1 || strcmp(siteabrev,'KS')
    meanTAB = array2table([meantab365.Day meantab365.HoursPropJU meantab365.HoursPropMA]);
meanTAB.Properties.VariableNames = {'Day' 'HoursPropJ' 'HoursPropM'};
meanTAB.JuvMal = meanTAB.HoursPropJ >0 & meanTAB.HoursPropM > 0;
meanTAB.Juv = meanTAB.HoursPropJ > 0 & meanTAB.JuvMal == 0;
meanTAB.Mal = meanTAB.HoursPropM > 0 & meanTAB.JuvMal == 0;
elseif strcmp(siteabrev,'AB') ==1
meanTAB = array2table([meantab365.Day meantab365.HoursPropFE meantab365.HoursPropJU meantab365.HoursPropMA]);
meanTAB.Properties.VariableNames = {'Day' 'HoursPropF' 'HoursPropJ' 'HoursPropM'};
meanTAB.FemJuv = meanTAB.HoursPropF > 0 & meanTAB.HoursPropJ > 0; 
meanTAB.FemMal = meanTAB.HoursPropF > 0 & meanTAB.HoursPropM > 0;
meanTAB.JuvMal = meanTAB.HoursPropJ >0 & meanTAB.HoursPropM > 0;
meanTAB.FemJuvMal = meanTAB.HoursPropF > 0 & meanTAB.HoursPropJ > 0 & meanTAB.HoursPropM > 0;
meanTAB.Fem = meanTAB.HoursPropF > 0 & meanTAB.FemJuv == 0 & meanTAB.FemMal == 0 & meanTAB.FemJuvMal == 0;
meanTAB.Juv = meanTAB.HoursPropJ > 0 & meanTAB.FemJuv == 0 & meanTAB.JuvMal == 0 & meanTAB.FemJuvMal == 0;
meanTAB.Mal = meanTAB.HoursPropM > 0 & meanTAB.FemMal == 0 & meanTAB.JuvMal == 0 & meanTAB.FemJuvMal == 0;
else
    meanTAB = array2table([meantabFE365.Day meantabFE365.HoursPropFE meantabJU365.HoursPropJU meantabMA365.HoursPropMA]);
meanTAB.Properties.VariableNames = {'Day' 'HoursPropF' 'HoursPropJ' 'HoursPropM'};
meanTAB.FemJuv = meanTAB.HoursPropF > 0 & meanTAB.HoursPropJ > 0; 
meanTAB.FemMal = meanTAB.HoursPropF > 0 & meanTAB.HoursPropM > 0;
meanTAB.JuvMal = meanTAB.HoursPropJ >0 & meanTAB.HoursPropM > 0;
meanTAB.FemJuvMal = meanTAB.HoursPropF > 0 & meanTAB.HoursPropJ > 0 & meanTAB.HoursPropM > 0;
meanTAB.Fem = meanTAB.HoursPropF > 0 & meanTAB.FemJuv == 0 & meanTAB.FemMal == 0 & meanTAB.FemJuvMal == 0;
meanTAB.Juv = meanTAB.HoursPropJ > 0 & meanTAB.FemJuv == 0 & meanTAB.JuvMal == 0 & meanTAB.FemJuvMal == 0;
meanTAB.Mal = meanTAB.HoursPropM > 0 & meanTAB.FemMal == 0 & meanTAB.JuvMal == 0 & meanTAB.FemJuvMal == 0;
end
%% Calculate A, AB, B, BC, C, CA, ABC
if strcmp(siteabrev,'KOA') == 1 || strcmp(siteabrev,'KS')
J = sum(meanTAB.Juv); %B
M = sum(meanTAB.Mal); %C
JM = sum(meanTAB.JuvMal); %BC
data = [J JM M];
else
F = sum(meanTAB.Fem); %A
J = sum(meanTAB.Juv); %B
M = sum(meanTAB.Mal); %C
FJ = sum(meanTAB.FemJuv); %AB
FM = sum(meanTAB.FemMal); %AC
JM = sum(meanTAB.JuvMal); %BC
FJM = sum(meanTAB.FemJuvMal); %ABC
data = [F FJ J JM M FM FJM];
end

vennX(data,0.01)
title(['Proportion of Classes in the ',titleNAME])

%save figure
export_fig([saveDir,'\',siteabrev,'PropVenn.png'],'-png','-transparent');



