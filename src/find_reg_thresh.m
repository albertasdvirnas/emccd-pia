function [regThresh,regSizes] = find_reg_thresh(bIm,allowedGapLength)

    regSizes = cellfun(@(x) size(x,1),bIm);
    uniqueSizeIm = unique(sort(regSizes));
    morethanGap = find(diff(uniqueSizeIm)>allowedGapLength,1,'first');
    regThresh = uniqueSizeIm(morethanGap);

    
end

