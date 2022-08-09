%% 
close all;clear all;clc;
% This script is made to compare two transfer functions (i.e. wind and Rev B)
% at two peak frequencies and make plots.
% NP 08/08/2022
%% User Definied Variables
GDrive = 'I';
%Freq = {9000,10000,11000}; %frequency of comparing interest in kHz
Freq = 10000;
site = 'NC';
region = 'WAT';
fullSite = [region,'_',site];
dpns = {'01','02','03','04'};

% REV B
MBARC_TF = [GDrive,':\Shared drives\MBARC_TF'];
TF = 560; %Which TF are you looking for
Date = '090612'; %date for TF

% WIND TF
Wind_TF = [GDrive,':\Shared drives\Wind_deltaTF\pub_TF\TF_Wind'];
Site = 'SOCALN35'; %site and deployment for wind TF

saveDIR = [GDrive,':\My Drive\TestTFs']; %directory where to save outputs
HARPsum = [saveDIR,'\HARPdataSummary.xlsx']; %HARP data summary sheet
%% Determine Series based on TF
if TF >= 100 && TF <= 399
    Series = '100-399';
    elseif TF >= 400 && TF <= 499
        Series = '400-499';
        elseif TF >= 500 && TF <= 599
            Series = '500-599';
            elseif TF >= 600 && TF <= 699
               Series = '600-699';
               elseif TF >= 700 && TF <= 799
                   Series = '700-799';
                   elseif TF >= 800 && TF <= 899
                       Series = '800-899';
                       elseif TF >= 900 && TF <= 999
                           Series = '900-999';
end
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
%% Loop through TFs to compare wind vs RevB
adjustTF = [];
for itf = 1:length(tfnum)
    adjustTF(itf) = getTF_BvsWind(site,MBARC_TF,Wind_TF,tfnum(itf),Freq);
end
%% Find and load TF files
% Rev B
RevB_DIR = [MBARC_TF,Series,'\',num2str(TF),'\',num2str(TF),'_',Date,'_B_HARP.tf'];
RevB = fopen(RevB_DIR);
[calB,~] = fscanf(RevB,'%f %f',[2,inf]);
fclose(RevB);

%Wind
WindTF_DIR = [Wind_TF,num2str(TF),'_',Site,'_TFnew.tf'];
WindTF = fopen(WindTF_DIR);
[calW,~] = fscanf(WindTF,'%f %f',[2,inf]);
fclose(WindTF);
%% Compare at frequencies of interest

%% Plots
figure
semilogx(calB(1,:),calB(2,:),'r-','LineWidth',2)
hold on
semilogx(calW(1,:),calW(2,:),'k--','LineWidth',2) %comment out if there is no wind TF
grid on
xlabel('Frequency [Hz]')
ylabel('Inverse Sensitivity [dB re \muPa//Count]')
axis([1 3e5 10 100])
legend('RevB','WindTF','Location','best')
title(Site)

% save plots
plotName = [saveDIR,'\',num2str(TF),'_',Site];
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