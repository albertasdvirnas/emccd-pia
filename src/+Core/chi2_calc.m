function [chi2Score] = chi2_calc(histAll, pdf,cdf,Nthresh,LU,nQuantiles)

    % truncated cdf
    pdf = pdf/(cdf(min(end,Nthresh))); % should sum to 1
    %     nansum(pdf(1:min(end,Nthresh)))

%     cdf = cdf/cdf(min(end,Nthresh));

    % two alternatives: even num of quantiles, or min bin size
    if nargin < 6
        minBinSize = 200;
        binEdges = find(histAll(1:end-1) >= minBinSize);
    else
        % find quantiles
        binEdges = zeros(1,nQuantiles+1);
        for i=1:nQuantiles-1
            % find the value where pvalue=1-cdf > pValThresh
            binEdges(i+1) = find(cdf > i/nQuantiles,1,'first');
        end
        %
        binEdges(1) = LU(1);
        binEdges(end) = min(LU(2),Nthresh);
    end
   
    
    numCounts = sum(histAll(1:Nthresh)); % number of counts in histogram
    % bin counts
    binCountsFit = zeros(1,length(binEdges)-1);
    histAllGroup = zeros(1,length(binEdges)-1);
    for i=1:length(binEdges)-1
          binCountsFit(i) = numCounts*sum(pdf(binEdges(i):binEdges(i+1)-1));
          histAllGroup(i) = sum(histAll(binEdges(i):binEdges(i+1)-1));
    end

    
%         figure,plot((histAllGroup - binCountsFit).^2./binCountsFit)

    chi2Score = sum ( (histAllGroup - binCountsFit).^2./binCountsFit);
end

