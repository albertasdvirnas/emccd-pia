function [chipParsSpec] = chippars_specific_gain(chipPars,gainIdx)

    chipParsSpec.gain = chipPars.gain(gainIdx);
    chipParsSpec.adFactor = chipPars.adFactor;
    chipParsSpec.countOffset = chipPars.countOffset(gainIdx);
    chipParsSpec.roNoise = chipPars.roNoise(gainIdx);
end

