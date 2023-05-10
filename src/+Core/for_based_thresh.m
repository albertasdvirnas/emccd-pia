
function [nOutliers,hasOutliers,stats,forEst] =  for_based_thresh(intensities, lambdaBg,gain,adFactor,countOffset,roNoise,alphaStar, U, stopIntensity)
        import Core.p_values_emccd_sorted;
        cdfSorted =  @(x)  1-p_values_emccd_sorted(x+0.5,lambdaBg,gain,adFactor,countOffset,roNoise);
        % Now calculate the estimate of FDR and FOR:
        [stats,allVals ] = estimate_stats(intensities, ...
             cdfSorted,U, stopIntensity);

        elt = find(stats.FOR < alphaStar,1,'last');
        stats.allVals = allVals;
        if ~isempty(elt)
            forEst = allVals(elt);%find(cdfFun>0.5,1,'first');% find(cdfFun>1,0.9,'first');
            stats.finalFor = stats.FOR(elt);
        else
            forEst = allVals(1);
            stats.finalFor = stats.FOR(1);
        end

        outliers = intensities > forEst;
        if ~isempty(outliers)
            nOutliers = sum(outliers);
            hasOutliers = 1;
        else
            hasOutliers = 0;
            nOutliers = 0;
        end
end
