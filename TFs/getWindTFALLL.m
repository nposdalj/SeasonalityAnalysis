function [Valss,windTF,AllWindSite] = getWindTFALL(tfwdir,tfnum,Valss,windTF,AllWindSite)
if exist(tfwdir)
    tfdstructw = dir(tfwdir);
    tfdcellw = struct2cell(tfdstructw);
    tfM1 = [];
    tfM2 = [];
    tfMatch = [];
    
    for iTF = 3 : length(tfdstructw)
        tfM1 = strfind(tfdcellw{1,iTF},num2str(tfnum));
        tfM2 = strfind(tfdcellw{1,iTF},'.tf');
        if ~isempty(tfM1) && ~isempty(tfM2)
            tfMatch =  iTF;
            tffilew = cell2mat(fullfile(tfdcellw(2,tfMatch),tfdcellw(1,tfMatch)));
            newStr = extractAfter(tffilew,[tfnum,'_']);
            finalStr = char(extractBefore(newStr,'_'));
            AllWindSite{find(cellfun(@isempty,AllWindSite),1)} = finalStr;
            [Valss{find(cellfun(@isempty,Valss),1)},windTF{find(cellfun(@isempty,windTF),1)}] = loadTF(tffilew);
        else
            tfM1 = []; tfM2 = [];
        end 
    end

    if isempty(tfMatch)
        Valss = 0;
        windTF = 0;
        AllWindSite = 0;
    end
end

