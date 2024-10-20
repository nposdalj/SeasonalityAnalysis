function [Valss,windTF,AllWindSite] = getWindTFALL(site,tfwdir,tfnum,sitenum,AllWindSite,Valss,windTF)
% find TF number from XLS and then read in TF"B" and TF "Wind" to find dif
%Search TFs folder for the appropriate preamp
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
            finalStr = extractBefore(newStr,'_');
            AllWindSite{1,sitenum} = finalStr;
            [Valss{1,sitenum},windTF{1,sitenum}] = loadTF(tffilew);
        else
            tfM1 = []; tfM2 = [];
        end 
    end

    %if exist('tfMatch','var')
     %   if ~isempty(tfMatch)
%
 %       else
  %      Valss = NaN;
   %     windTF = NaN;
    %    end
    %else
     %   Valss = NaN;
      %  windTF = NaN;
       % AllWindSite = NaN;
    %end
    %end
end

