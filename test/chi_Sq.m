%%

% 1) Simulate data
% 2) Estimate lambda based on Nthresh
% 3) Calculate pdf/cdf for intensities in range
% 4) Create bins for 
dims = [64 64];
numRuns = 100;
    %% Chi square dist generation
%     truncationIdx = 1;      % truncation is done at truncationIdx*100 % of the data  
%     nRandNumbersBg = 1E3;     % number of random numbers generated
%     nRuns = 100;                % number of simulation runs
    nQuantiles = 5;          % number of quantiles used for binning
    Nthresh = 60;
    import Core.chi2_calc;

chi2Score = zeros(1,numRuns);
for runIdx = 1:numRuns;
    tic
    SNRVals = [0];  
    [filenames,chipPars,data] = simulate_random_beads_full(100, 20, 100, 0,SNRVals, 1, [],dims); % zoom 100x
    
    gain = chipPars.gain;
    adFactor = chipPars.adFactor;
    roNoise = chipPars.roNoise;
    offset = chipPars.countOffset;
    r = gain/adFactor;
    
    
    %% Estimate lambdaBg
    intensities = data{1}{1}.image;
%     sortI = sort(data{1}{1}.image,'ascend');
    %%
    opt = statset('MaxIter',200,'MaxFunEvals',400,'FunValCheck','on');


    lowestIntThresh = ceil(quantile(intensities,0.25)); 
    structRes.lowestIntThresh = lowestIntThresh;
    %
    sortI = sort(intensities);
    S = numel(sortI);
%     Nthresh = sortI(round(S/2));

%     m = find(sortI>Nthresh,1,'First')-1; 
%     sortTruncI = sortI(1:m);
    % Fit lambda based on data
    lamGuess = abs((sortI(round(end/2)) - offset)/(gain/adFactor));

    % Prep data for truncated fit
    import Core.calc_bounds;
    [L, U, EX, STD] = calc_bounds(lamGuess, gain, adFactor, offset, roNoise, 6);
    structRes.LU = [max(1,ceil(L)) floor(U)];


    % Get bin edges
    binCenters = [max(1,ceil(L)):1:floor(U)+1 inf]; % bin edges shifted by half  

    % histogram for intensities. Calculated only once!
    histAll = zeros(1,binCenters(end-1));
    histAll(binCenters(1:end-1)) = histcounts(sortI,binCenters-0.5)';
    
    %%
    [lambdaBg, pci, pdf, cdf] = est_lambda(histAll,lamGuess, Nthresh, gain, adFactor, offset, roNoise,structRes.LU,opt);

%     [lambdaBg, pci, L, U, binEdges, binPos,sortTruncI] = est_lambda(sortI, Nthresh, gain, adFactor, offset, roNoise, r);
%     lambdaBg
        chi2Score(runIdx) = chi2_calc(histAll, pdf, cdf, Nthresh,structRes.LU, nQuantiles);

%         [ chi2Score(runIdx)] = chi2_fun(sortTruncI,Nthresh,L,U,lambdaBg,gain,adFactor,offset,roNoise,nQuantiles);

    
%     % intensities we integrate over (thresh can't be outside of these
%     intensitiesU = [max(1,ceil(L)):min(Nthresh,floor(U))]; % take edge point of each box, so +0.5
%     
% %     Nthresh = 50;
%     % cdf, in this case do not need to be truncated
%     cdf = nan(1,max(intensitiesU));
%     pdf = nan(1,max(intensitiesU));
% % 
% %     gain = 20;
% %     lambdaBg = 38;
% %     adFactor = 36;
% %     roNoise = 1.44;
% %     offset = 27;
% 
%     % calc cdf/pdf for all of these
%     [pdf(intensitiesU), cdf(intensitiesU)] = pdf_cdf_from_characteristic_fun(intensitiesU,lambdaBg,gain,adFactor,offset,roNoise);
%      cdf(cdf<0) = 0;
% 
%     % calc cdf for thresh value // DO NOT! WILL GIVE FALSE CDF!
% %     [pdfEnd, cdfEnd] = pdf_cdf_from_characteristic_fun(Nthresh,lambdaBg,gain,adFactor,offset,roNoise);
% %     [pdfStart, cdfStart] = pdf_cdf_from_characteristic_fun(ceil(L),lambdaBg,gain,adFactor,offset,roNoise);
% 
%     % truncated cdf
%     pdf = pdf/(cdf(min(end,Nthresh)));
%     nansum(pdf(1:min(end,Nthresh)))
% 
%     cdf = cdf/cdf(min(end,Nthresh));
% %         [pdf, cdf] = pdf_cdf_from_characteristic_fun([18:60],lambdaBg,gain,adFactor,offset,roNoise);
% 
% %     [histAll,edges] = histcounts(sortTruncI,100,'Normalization','cdf');
% 
%     % q = (2:9)/nQuantiles;
%     
%     binEdges = zeros(1,nQuantiles+1);
%     for i=1:nQuantiles-1
%         % find the value where pvalue=1-cdf > pValThresh
%         binEdges(i+1) = find(cdf > i/nQuantiles,1,'first');
%     end
%     %
%     binEdges(1) = intensitiesU(1);
%     binEdges(end) = min(intensitiesU(end),Nthresh);
%     
% %     [pdf, cdf] = pdf_cdf_from_characteristic_fun( binEdges(2)+1,lambdaBg,gain,adFactor,offset,roNoise);
% 
%     % Use CDF of histcounts instead
%     [histAll,edges] = histcounts(sortTruncI,binEdges-0.5);
% 
%     
% % [~, cdfvals] = pdf_cdf_from_characteristic_fun(binEdges,lambdaBg,gain,adFactor,offset,roNoise);
% 
%     %
%     binCountsFit = zeros(1,length(binEdges)-1);
%     for i=1:length(binEdges)-1
%           binCountsFit(i) = length(sortTruncI)*sum(pdf(binEdges(i):binEdges(i+1)-1));
% %         [pdfvals, cdfvals] = pdf_cdf_from_characteristic_fun([binEdges(i):binEdges(i+1)],lambdaBg,gain,adFactor,offset,roNoise);
% %        binCountsFit(i) = sum(length(sortTruncI)*pdfvals);
%     end
% %     cdfCenters = binEdges(1:end-1)-0.5 + diff(binEdges-0.5)/2;
% 
% %     cdf
%     %
%     % cdfEmccd = cdf(binEdges);
% %     [pdfEmccd, cdfEmccd] = pdf_cdf_from_characteristic_fun(cdfCenters,lambdaBg,gain,adFactor,offset,roNoise);
% %     binCountsFit = length(sortTruncI)*cdfEmccd;
%     
% %         figure,plot((histAll - binCountsFit).^2./binCountsFit)
% 
%     chi2Score(runIdx) = sum ( (histAll - binCountsFit).^2./binCountsFit);
    toc
end

%%
% % plot
 figure
[countsChi2,binEdges] = histcounts(chi2Score);
 binPos = binEdges(1:end-1) + diff(binEdges)/2;
bar(binPos,countsChi2); hold on;
xlabel('chi2-score')
ylabel('counts')
title(['Sample size =',num2str(dims(1)*dims(2)), ' Nthresh = ', num2str(Nthresh)])

% % chi2 - distribution 
x = 0:0.03:30;
y = chi2pdf(x,nQuantiles-2 );  
% yTr = chi2cdf(Nthresh,nQuantiles - 2)
               % Should be NQuantiles -1 degree of freedom if lambdaBg is
               % known.
               % In case NQuantiles is ML estimated using the same histogram
               % as the one that is used for the chi2 test then
               % one should use nQuantiles - 2 degree of freedom.

plot(x,numRuns*y*(binPos(2)-binPos(1)),'k-','linewidth',3)
legend({'Sample fit','$\chi^2$ dist'},'Interpreter','latex')


