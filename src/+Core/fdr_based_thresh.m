

function [nOutliers,hasOutliers,stats,fdrEst] =  fdr_based_thresh(intensities, lambdaBg, gain, adFactor, countOffset, roNoise, qStar, U, stopIntensity)
    % calculate FDR based threshold

    % CDF - 1-pvalue (alternatively 3rd output of p_values_emccd_sorted)
    import Core.p_values_emccd_sorted;
    cdfSorted =  @(x)  1-p_values_emccd_sorted(x+0.5,lambdaBg,gain,adFactor,countOffset,roNoise);

    % Now calculate the estimate of FDR and FOR:
    [stats,allVals ] = estimate_stats(intensities, ...
        cdfSorted, U, stopIntensity);

    elt = find(stats.fdr<qStar,1,'first');
    if ~isempty(elt)
        fdrEst = allVals(elt);%find(cdfFun>0.5,1,'first');% find(cdfFun>1,0.9,'first');
    else
        fdrEst = allVals(1);
    end

    outliers = intensities > fdrEst;
    if ~isempty(outliers)
        nOutliers = sum(outliers);
        hasOutliers = 1;
    else
        hasOutliers = 0;
        nOutliers = 0;
    end
end