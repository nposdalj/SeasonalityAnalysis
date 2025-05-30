function adjustTF = getTF_BvsWind(site,tfbdir,tfwdir,tfnum,Freq)
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
    [ax,~] = size(tffilename);
    if ax > 1
    tffile = fullfile(suggestedTFPath,stfnum,tffilename(1).name);
    [freqB,uppcB] = loadTF(tffile); % open and read B Transfer Function
    else
    tffile = fullfile(suggestedTFPath,stfnum,tffilename.name);
    [freqB,uppcB] = loadTF(tffile); % open and read B Transfer Function
    end
    else
    tffile = fullfile(suggestedTFPath,stfnum,[stfnum,'*A_HARP.tf']);
    tffilename = dir(tffile);
    if ~isempty(tffilename)
    tffile = fullfile(suggestedTFPath,stfnum,tffilename.name);
    [freqB,uppcB] = loadTF(tffile); % open and read A Transfer Function
    else
    disp([num2str(stfnum),' TF does not exist'])
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(tfwdir)
    tfdstructw = dir(tfwdir);
    tfdcellw = struct2cell(tfdstructw);

    tfM1 = []; tfM2 = []; tfM3 = [];
    tfMatch = [];
    for iTF = 3 : length(tfdstructw)
        tfM1 = strfind(tfdcellw{1,iTF},num2str(tfnum));
        tfM2 = strfind(tfdcellw{1,iTF},site);
        tfM3 = strfind(tfdcellw{1,iTF},'.tf');
        if ~isempty(tfM1) && ~isempty(tfM2) && ~isempty(tfM3)
            tfMatch =  iTF;
        else
            tfM1 = []; tfM2 = []; tfM3 = [];
        end
    end
end   

    if ~isempty(tfMatch)
        tffilew = cell2mat(fullfile(tfdcellw(2,tfMatch),tfdcellw(1,tfMatch)));
        [freqW,uppcW] = loadTF(tffilew); % open and
        iB = find(freqB == Freq);
        iW = find(freqW == Freq);
        if ~isempty(iB)
            adjustTF = uppcW(iW) - uppcB(iB) ;
            disp([num2str(tfnum),'  Wind TF : ',num2str(uppcW(iW)),'  B TF : ',num2str(uppcB(iB)),...
            ' AdjustTF = ',num2str(adjustTF)])
        else
            disp([num2str(Freq),' dBs does not exist for wind TF'])
            adjustTF = NaN;
        end
    else
        disp(['No wind TF exists for ', num2str(tfnum)])
        adjustTF = NaN;
        
    end
else
    disp(['Transfer Function does not exist; double check data summary sheet'])
    adjustTF = NaN;
end
end

