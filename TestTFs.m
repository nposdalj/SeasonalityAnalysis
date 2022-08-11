%% 
close all;clear all;clc;
% This script is made to compare two transfer functions (i.e. wind and Rev B)
% at two peak frequencies and make plots.
% NP 08/08/2022
%% User Definied Variables
GDrive = 'I';
Freq = {9000,10000,11000}; %frequency of comparing interest in kHz
site = 'PAL';
region = '';
if ~isempty(region)
fullSite = [region,'_',site];
fullSitePlots = [region,'\_',site];
else
fullSite = site;
fullSitePlots = site;
end
dpns = {'02','03','05','06','07','08'};

% REV B
MBARC_TF = [GDrive,':\Shared drives\MBARC_TF'];
%TF = 560; %Which TF are you looking for
%Date = '090612'; %date for TF

% WIND TF
Wind_TF = [GDrive,':\Shared drives\Wind_deltaTF\pub_TF\TF_Wind'];
%Site = 'SOCALN35'; %site and deployment for wind TF

saveDIR = [GDrive,':\My Drive\TestTFs']; %directory where to save outputs
HARPsum = [saveDIR,'\HARPdataSummary.xlsx']; %HARP data summary sheet
%% Loop through HARP data summary sheet and find matching sites and TFs
dtable = readtable(HARPsum);
stxt = size(dtable); 
tfnum = [];
ifoundx = 0;
for itab = 1 : stxt(1)
    ifound = strfind(dtable.Data_ID(itab),fullSite);
    if cell2mat(ifound) >0
        ifoundx = ifoundx + 1;
        tfnum(ifoundx) = str2double(dtable.PreAmp(itab));
    end
end
%% Loop through TFs to look for TF adjustment
[~,q] = size(Freq);
[~,qq] = size(tfnum);
adjustTF = zeros(q,qq); %preallocate vecctor

for itf = 1:qq
    for ifr = 1:q
        Freqq = Freq{ifr};
        adjustTF(ifr,itf) = getTF_BvsWind(site,MBARC_TF,Wind_TF,tfnum(itf),Freqq);
    end
end

%max adjustment for each TF (regardless of peak frequency)
maxAdjust = max(adjustTF);

%save adjustments as text files
save([saveDIR,'\',fullSite,'_Adjustments.mat'],'tfnum','Freq','adjustTF','maxAdjust');
%% Find and load TF files
% Rev B
RevBTF = [];
Vals = [];
for itf = 1:qq
    [Vals{itf},RevBTF{itf}] = getRevB(site,MBARC_TF,tfnum(itf));
end

%Wind
windTF = [];
Valss = [];
for itf = 1:qq
    [Valss{itf},windTF{itf}] = getWindTF(site,Wind_TF,tfnum(itf));
end
%% Plots
for itrP = 1:qq
tf_now = num2str(tfnum(itrP));
figure
semilogx(Vals{itrP},RevBTF{itrP},'r-','LineWidth',2)
hold on
semilogx(Valss{itrP},windTF{itrP},'k--','LineWidth',2) %comment out if there is no wind TF
grid on
xlabel('Frequency [Hz]')
ylabel('Inverse Sensitivity [dB re \muPa//Count]')
axis([1 3e5 10 100])
legend('RevB','WindTF','Location','best')
title([fullSitePlots,'-',tf_now])

% save plots
plotName = [saveDIR,'\',fullSite,'_',tf_now];
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