function [Valss,windTF] = getWindTF(site,tfwdir,tfnum)
% find TF number from XLS and then read in TF"B" and TF "Wind" to find dif
%Search TFs folder for the appropriate preamp
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
        [Valss,windTF] = loadTF(tffilew); % open and
    else
        Valss = NaN;
        windTF = NaN;
    end
end

