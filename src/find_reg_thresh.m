function [regThresh,regSizes] = find_reg_thresh(bIm,labelIm,allowedGapLength)

%     regSizes = cellfun(@(x) size(x,1),bIm);
    regSizes = arrayfun(@(x) sum(labelIm(:)==x),1:length(bIm));

    uniqueSizeIm = unique(sort([0 regSizes]));
    morethanGap = find(diff(uniqueSizeIm)>allowedGapLength,1,'first');
    regThresh = uniqueSizeIm(morethanGap);


    
end

