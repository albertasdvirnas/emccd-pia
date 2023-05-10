function [chi2Score,pdf,cdf,intensitiesU] = chi2_fun(sortTruncI,Nthresh,L,U,lambdaBg,gain,adFactor,offset,roNoise,nQuantiles)
    % intensities we integrate over (thresh can't be outside of these
    intensitiesU = [max(1,ceil(L)):min(Nthresh,floor(U))]; % take edge point of each box, so +0.5
    
%     Nthresh = 50;
    % cdf, in this case do not need to be truncated
    cdf = nan(1,max(intensitiesU));
    pdf = nan(1,max(intensitiesU));
% 
%     gain = 20;
%     lambdaBg = 38;
%     adFactor = 36;
%     roNoise = 1.44;
%     offset = 27;

    % calc cdf/pdf for all of these
    [pdf(intensitiesU), cdf(intensitiesU)] = pdf_cdf_from_characteristic_fun(intensitiesU,lambdaBg,gain,adFactor,offset,roNoise);
     cdf(cdf<0) = 0;

    % calc cdf for thresh value // DO NOT! WILL GIVE FALSE CDF!
%     [pdfEnd, cdfEnd] = pdf_cdf_from_characteristic_fun(Nthresh,lambdaBg,gain,adFactor,offset,roNoise);
%     [pdfStart, cdfStart] = pdf_cdf_from_characteristic_fun(ceil(L),lambdaBg,gain,adFactor,offset,roNoise);

    % truncated cdf
    pdf = pdf/(cdf(min(end,Nthresh)));
%     nansum(pdf(1:min(end,Nthresh)))

    cdf = cdf/cdf(min(end,Nthresh));
%         [pdf, cdf] = pdf_cdf_from_characteristic_fun([18:60],lambdaBg,gain,adFactor,offset,roNoise);

%     [histAll,edges] = histcounts(sortTruncI,100,'Normalization','cdf');

    % q = (2:9)/nQuantiles;
    
    binEdges = zeros(1,nQuantiles+1);
    for i=1:nQuantiles-1
        % find the value where pvalue=1-cdf > pValThresh
        binEdges(i+1) = find(cdf > i/nQuantiles,1,'first');
    end
    %
    binEdges(1) = intensitiesU(1);
    binEdges(end) = min(intensitiesU(end),Nthresh);
    
%     [pdf, cdf] = pdf_cdf_from_characteristic_fun( binEdges(2)+1,lambdaBg,gain,adFactor,offset,roNoise);

    % Use CDF of histcounts instead
    [histAll,edges] = histcounts(sortTruncI,binEdges-0.5);

    
% [~, cdfvals] = pdf_cdf_from_characteristic_fun(binEdges,lambdaBg,gain,adFactor,offset,roNoise);

    %
    binCountsFit = zeros(1,length(binEdges)-1);
    for i=1:length(binEdges)-1
          binCountsFit(i) = length(sortTruncI)*sum(pdf(binEdges(i):binEdges(i+1)-1));
%         [pdfvals, cdfvals] = pdf_cdf_from_characteristic_fun([binEdges(i):binEdges(i+1)],lambdaBg,gain,adFactor,offset,roNoise);
%        binCountsFit(i) = sum(length(sortTruncI)*pdfvals);
    end
%     cdfCenters = binEdges(1:end-1)-0.5 + diff(binEdges-0.5)/2;

%     cdf
    %
    % cdfEmccd = cdf(binEdges);
%     [pdfEmccd, cdfEmccd] = pdf_cdf_from_characteristic_fun(cdfCenters,lambdaBg,gain,adFactor,offset,roNoise);
%     binCountsFit = length(sortTruncI)*cdfEmccd;
    
%         figure,plot((histAll - binCountsFit).^2./binCountsFit)

    chi2Score = sum ( (histAll - binCountsFit).^2./histAll);
end

