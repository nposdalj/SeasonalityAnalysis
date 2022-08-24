function adjustTF = getAdjustment(windTF,Vals,RevBTF,Valss,Freq,Phone,AllWindSite)
adjustTF = [];
for ww = 1:length(windTF)
for ff = 1:length(Freq)
        iB = find(Valss == Freq{ff});
        iW = find(Vals{ww} == Freq{ff});
        if ~isempty(iB)
            adjustTFint = windTF{ww}(iW) - RevBTF(iB);
            disp(['Site:',AllWindSite{ww},' Hydrophone:',Phone,'  Wind TF: ',num2str(windTF{ww}(iW)),'  B TF: ',num2str(RevBTF(iB)),...
            ' AdjustTF = ',num2str(adjustTFint)])
        else
            disp([num2str(Freq{ff}),' dBs does not exist for wind TF'])
            adjustTFint = NaN;
        end
        adjustTF{ww,ff} = adjustTFint;
end
end    
end

