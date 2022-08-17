%% 
close all;clear all;clc;
% This script is made to compare two transfer functions (i.e. wind and Rev B)
% at two peak frequencies and make plots.
% NP 08/08/2022
%% User Definied Variables
GDrive = 'I';

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
%% Loop through HARP data summary sheet and find matching sites and TFs
for allTF = 1:height(TFtable)
Phone = num2str(TFtable.TF(allTF));
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
end
%% Plots
if isa(AllWindSite,'double')
    disp('No Wind TFs to Plot')
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
end
end