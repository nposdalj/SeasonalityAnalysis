%% 
close all;clear all;clc;
% This script is made to compare two transfer functions (i.e. wind and Rev B)
% at two peak frequencies and make plots.
% NP 08/08/2022
%% User Definied Variables
GDrive = 'I';
Freq = {9000,10000,11000}; %frequency of comparing interest in kHz
% REV B
MBARC_TF = [GDrive,':\Shared drives\MBARC_TF'];

% WIND TF
Wind_TF = [GDrive,':\Shared drives\Wind_deltaTF\pub_TF\TF_Wind'];

saveDIR = [GDrive,':\My Drive\TestTFs']; %directory where to save outputs
HARPsum = [saveDIR,'\HARPdataSummary.xlsx']; %HARP data summary sheet
TFdata = [saveDIR,'\SiteTransferFunctions.xlsx']; %All the TFs I want to test
%% Load Transfer Functions to test
TFtable = readtable(TFdata);
dtable = readtable(HARPsum);
stxt = size(dtable); 
%% Prepare tables
sz = [107 4];
varNames = {'Phone','WindDiff','Adjustment','RevBDiff'};
T = array2table(zeros(sz));
T.Properties.VariableNames = varNames;
T.Adjustment = num2cell(T.Adjustment);
T.Phone = TFtable.TF;
%% Loop through HARP data summary sheet and find matching sites and TFs
for allTF = 1:height(TFtable)
Phone = num2str(TFtable.TF(allTF));
disp(['Beginning Analysis for Hydrophone ' Phone])
ifoundx = 0;
for itab = 1 : stxt(1)
    ifound = strfind(dtable.PreAmp(itab),Phone);
    if cell2mat(ifound) >0
        ifoundx = ifoundx + 1;
        Site(ifoundx) = dtable.Data_ID(itab);
        WindSite(ifoundx) = dtable.Data_ID_Wind(itab);
    end
end
%% Find and load TF files
% Rev B
[Vals,RevBTF] = getRevB(Site(1),MBARC_TF,str2double(Phone));

%Wind
[~,qq] = size(Site);
windTF = cell(1,qq);
Valss = cell(1,qq);
AllWindSite = cell(1,qq);
[Valss,windTF,AllWindSite] = getWindTFALLL(Wind_TF,Phone,Valss,AllWindSite,windTF);
%Delete empty cells
if isa(AllWindSite,'double')
else
    windTF = windTF(~cellfun('isempty',windTF));
    Valss = Valss(~cellfun('isempty',Valss));
    AllWindSite = AllWindSite(~cellfun('isempty',AllWindSite));
end
%% Test if the wind TFs are different from site to site
[~,qqq] = size(AllWindSite);
diff = [];
if qqq > 1
    for itrD = 2:qqq
        diffInt = windTF{1,1} - windTF{1,qqq};
        diff{itrD} = max(diffInt);
    end
    diff = diff(~cellfun('isempty',diff));
    diffMAT = cell2mat(diff);
    GreaterOned = diffMAT > 1;
    LessOned = diffMAT < -1;
    if sum(sum(GreaterOned)) >= 1 || sum(sum(LessOned)) >= 1
        disp(['The wind TF for Hydrophone ',Phone,' have differences greater than 1 dB.'])
        T.WindDiff(allTF) = 1;
    end
end
%% Look for differences between each wind TF and Rev B
adjustTF = [];
if ~isa(AllWindSite,'double') && isa(RevBTF,'double')
    adjustTF = getAdjustment(windTF,Valss,RevBTF,Vals,Freq,Phone,AllWindSite);
    T.Adjustment(allTF) = {adjustTF};
    adjustTFMAT = cell2mat(adjustTF);
    GreaterOne = adjustTFMAT > 1;
    LessOne = adjustTFMAT < -1;
    if sum(sum(GreaterOne)) >= 1 || sum(sum(LessOne)) >= 1
        T.RevBDiff(allTF) = 1;
    end
end
%% Plots
if isa(AllWindSite,'double')
    disp(['No Wind TFs to Plot for ',Phone])
else
figure
semilogx(Vals,RevBTF,'r-','LineWidth',2)
hold on
semilogx(Valss{1},windTF{1},'LineWidth',2)
if qqq>1
    for itrP = 2:qqq
        semilogx(Valss{itrP},windTF{itrP},'LineWidth',2)
    end
end
hold off
grid on
xlabel('Frequency [Hz]')
ylabel('Inverse Sensitivity [dB re \muPa//Count]')
axis([1 3e5 10 100])
legend(['RevB',AllWindSite],'Location','best')
title(['All Sites - ',Phone])

% save plots
plotName = [saveDIR,'\',Phone,'_AllSites'];
ylim([30 100])
xlim([1 200000])
saveas(gcf,[plotName,'.fig'])
saveas(gcf,[plotName,'.png'])

set(gca,'XScale','linear')
xlim([1000 10000])
ylim([40 90])
saveas(gcf,[plotName,'_windRange.png'])

xlim([5000 100000])
ylim([25 90])
saveas(gcf,[plotName,'_5-95kHz.png'])

xlim([8500 11500])
ylim([25 90])
saveas(gcf,[plotName,'_9-11kHz.png'])
close all;
end
end
%% Save table as .mat file
save([saveDIR,'\SummaryTable.mat'],'T');