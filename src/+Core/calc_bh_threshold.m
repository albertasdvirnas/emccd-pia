
function [threshold, nOutliers, outliers,hasOutliers,estThresh ] = calc_bh_threshold(sortI, lambdaBg, gain, adFactor, countOffset, roNoise,qStar)
        % Calculate p-values. Anything above mean+6STD from the dist automatically assumed to be   outliers      
        import Core.p_values_emccd_sorted;
        [pVals, pdfUnique, cdfsUnique, intUnique] = p_values_emccd_sorted(sortI,lambdaBg, gain, adFactor, countOffset, roNoise);   

%         [pVals, pdfUnique, cdfsUnique, intUnique] = p_values_emccd_continuous(sortI,lambdaBg, gain, adFactor, countOffset, roNoise);        
        pValsFlipped = fliplr(pVals);

        m = length(sortI);

        % Find outliers https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-6/The-positive-false-discovery-rate--a-Bayesian-interpretation-and/10.1214/aos/1074290335.full
        threshold = ((1:m)./(m)).*qStar; % m-number of remaining pixels % BH procedure

        outliers = find(pValsFlipped < threshold,1,'last');

        estThresh = sortI(end-outliers);
        if ~isempty(outliers)
            nOutliers = outliers(end);
            hasOutliers = 1;
        else
            hasOutliers = 0;
            nOutliers = 0;
        end
        
        %         figure,histogram(sortTruncI,'normalization','pdf')
        %         hold on
        %         plot(intUnique,pdfUnique)
        %         xlim([0 300])
        %         
        %         histogram(sortTruncI,'normalization','cdf','DisplayStyle','stairs')
        %         hold on
        %         plot(intUnique,cdfsUnique)
end