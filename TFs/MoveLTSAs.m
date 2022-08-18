%% This code was created to move all LTSAs in a directory to another specificed folder
%NP 08182022
close all;clear all;clc;
%% Set Directories
StartDIR = '\\snowman.ucsd.edu\Atlantic_Region_Decimated_3\WAT_HZ_01'; %Directory where the LTSAs currently live
EndDIR = '\\frosty.ucsd.edu\LTSA\WAT\HZ'; %Directory where you want to move the LTSAs
%% Find all the LTSAs in a given directory and move the files into the specified directory
dirinfo = dir(StartDIR);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
dirinfo(1:2) = [];
subdirinfo = cell(length(dirinfo),1);
for K = 1 : length(dirinfo)
  thisdir = [dirinfo(K).folder,'\',dirinfo(K).name];
  subdirinfo{K} = dir(fullfile(thisdir, '*.LTSA'));
  [~,szsu] = size(subdirinfo{K,1});
  if szsu > 1
  for KK = 1:height(subdirinfo{K,1})
    disp(['Moving ',subdirinfo{K,1}(KK).name])
    copyfile([subdirinfo{K,1}(KK).folder,'\',subdirinfo{K,1}(KK).name],EndDIR)
  end
  else
    disp(['Moving ',subdirinfo{K,1}.name])
    copyfile([subdirinfo{K,1}.folder,'\',subdirinfo{K,1}.name],EndDIR)
  end
end

disp(['Done Moving LTSAs'])