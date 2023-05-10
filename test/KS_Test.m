[f,x] = ecdf(sortTruncI);
hold on
plot(binPos,cdfEmccd)


histAllcdf = histcounts(sortTruncI(sortTruncI<binEdges(end)),binEdges,'Normalization','cdf')';

figure;
plot(binPos,cdfEmccd)
hold on
plot(binPos,histAllcdf)

figure,plot(abs(cdfEmccd-histAllcdf)/sqrt(length(binPos)))

                
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

%     [pdfEmccd,cdfEmccd] = pdf_cdf_emccd(binPos,lambda,chipPars,N);
    
    [~ ,cdfEmccdEnd] = pdf_cdf_emccd(min(binEdges(end),max(sortTruncI)+0.5),lambda,gain, adFactor, offset, roNoise,L,U);
    posMax(ii) = max(abs(cdfEmccd-histAllcdf)/sqrt(length(binPos)));
end
% 
[a,b] = min(posMax);

figure,plot(binEdgesAll,posMax)
