function [Vals,RevBTF] = getRevB(site,tfbdir,tfnum)
% find TF number from XLS and then read in TF"B" and TF "Wind" to find dif
%Search TFs folder for the appropriate preamp
if ~isnan(tfnum)
if exist(tfbdir)
    tfdstruct = dir(tfbdir);
    tfdcell = struct2cell(tfdstruct);

    tfd = floor(tfnum(1)/100)*100;
    tfMatch = [];
    iTF = 3;
     if tfd <= 300
        suggestedTFPath = [tfbdir,'\100-399'];
    else
    while isempty(tfMatch) && iTF<=length(tfdstruct)
        tfMatch = strfind(tfdcell{1,iTF},num2str(tfd));
        iTF = iTF+1;
    end
    if ~isempty(tfMatch)
        tfMatchIdx = iTF-1;
        suggestedTFPath = cell2mat(fullfile(tfdcell(2,tfMatchIdx),tfdcell(1,tfMatchIdx)));
    end
     end
end

stfnum = num2str(tfnum);
tffile = fullfile(suggestedTFPath,stfnum,[stfnum,'*B_HARP.tf']);
tffilename = dir(tffile);
if ~isempty(tffilename)
    [ax,~]=size(tffilename);
    if ax > 1
     tffile = fullfile(suggestedTFPath,stfnum,tffilename(1).name);
    [Vals,RevBTF] = loadTF(tffile); % open and read B Transfer Function
        else
    tffile = fullfile(suggestedTFPath,stfnum,tffilename.name);
    [Vals,RevBTF] = loadTF(tffile); % open and read B Transfer Function
    end
    else
    tffile = fullfile(suggestedTFPath,stfnum,[stfnum,'*A_HARP.tf']);
    tffilename = dir(tffile);
    if ~isempty(tffilename)
    tffile = fullfile(suggestedTFPath,stfnum,tffilename.name);
    [Vals,RevBTF] = loadTF(tffile); % open and read A Transfer Function
    else
        disp([num2str(tfnum),' TF does not exist'])
        Vals = NaN;
        RevBTF = NaN;
    end
end
else
    disp('Transfer Function does not exist; no plots made')
    Vals = NaN;
    RevBTF = NaN;
end
end

