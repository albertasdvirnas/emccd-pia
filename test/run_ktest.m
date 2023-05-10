function [threshVal,binEdgesAll,posMax,a,threshPos] = run_ktest(sortTruncI,lambda,gain, adFactor, offset, roNoise)


    import Core.calc_bounds;
    [L, U, EX, STD] = calc_bounds(lambda,gain,adFactor,offset,roNoise);

    
% Get bin edges
binEdgesAll = [ceil(L):1:floor(U)]-0.5; % bin edges shifted by half  

binEdgesThresh = quantile(sortTruncI,0.05)-0.5; % bin edges shifted by half  
stPos = find(binEdgesThresh==binEdgesAll)
%
posMax = nan(1,length(binEdgesAll));
for ii=stPos:length(binEdgesAll)
    binEdges = binEdgesAll(1:ii);
 
    histAll = histcounts(sortTruncI(sortTruncI<binEdges(end)),binEdges)';
    histAllcdf = histcounts(sortTruncI(sortTruncI<binEdges(end)),binEdges,'Normalization','cdf')';

    % center position for each bin
    binPos = binEdges(1:end-1) + diff(binEdges)/2;
   
    import Core.pdf_cdf_emccd;
    [pdfEmccd, cdfEmccd] = pdf_cdf_emccd(binPos,lambda,gain, adFactor, offset, roNoise,L,U);

    [pdfEmccd,cdfEmccdEnd] = pdf_cdf_emccd(binPos(end),lambda,gain, adFactor, offset, roNoise,L,U);
    
    cdfEmccd = cdfEmccd/cdfEmccdEnd;

%     figure;plot(histAllcdf);hold on;
%     plot(cdfEmccd)
%     [~ ,cdfEmccdEnd] = pdf_cdf_emccd(min(binEdges(end),max(sortTruncI)+0.5),lambda,gain, adFactor, offset, roNoise,L,U);
    numBins = sum(histAll);
    numBins = length(binPos);
    posMax(ii) = max(abs(cdfEmccd-histAllcdf)*sqrt(numBins));
end
% 
[a,b] = min(posMax);
threshVal = binEdgesAll(b);
threshPos = stPos:length(binEdgesAll);

end

