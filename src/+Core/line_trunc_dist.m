function [logL] = line_trunc_dist(sortTruncI,lambda,...
                             gain, adFactor, offset, roNoise, stopIntensity)
    % Calculates the  log-likelihood for the truncated EMCCD-PIA distribution
    % 
    % Args: 
    % 
    %   sortTruncI = sorted and truncated intensity values 
    %   lambda = Poisson parameter
    %   gain, adFactor, countOffset, roNoise - chipPars
    %
    % Returns:
    % 
    % logL = log likelihood
    %
    % Comment: 
    % The truncated PDF is
    %      PDF_trunc = pdfEmccd(I)/cdfEmccd(I_trunc) for I <= I_trunc
    %      PDF_trunc = 0 elsewhere
    % Here, I_trunc is the truncation intensity.
    %  
    import Core.calc_bounds;
    [L, U, EX, STD] = calc_bounds(lambda,gain,adFactor,offset,roNoise);

    %     U = min(max(sortTruncI),U); % limit to U for truncated case ?
                 
    % Get bin edges
    binEdges = [ceil(L):1:floor(U)]-0.5; % bin edges shifted by half  
    % histogram for intensities
    histAll = histcounts(sortTruncI,binEdges)';
    % center position for each bin
%     binPos = binEdges(1:end-1)+ diff(binEdges)/2;
    binPos = binEdges(2:end);
   
    import Core.pdf_cdf_emccd;
    [pdfEmccd, cdfEmccd] = pdf_cdf_emccd(binPos,lambda,gain, adFactor, offset, roNoise,L,U);

%     [pdfEmccd,cdfEmccd] = pdf_cdf_emccd(binPos,lambda,chipPars,N);
    
    [~ ,cdfEmccdEnd] = pdf_cdf_emccd(min(binEdges(end),max(sortTruncI)+0.5),lambda,gain, adFactor, offset, roNoise,L,U);
%     [~ ,cdfEmccdStart] = pdf_cdf_emccd(30,lambda,gain, adFactor, countOffset, roNoise,L,U);

    % log-likelihood
%     logL = sum(histAll.*log(pdfEmccd)) - sum(histAll)*log(cdfEmccdEnd);

   

    [~ ,cdfThresh] = pdf_cdf_emccd(stopIntensity,lambda,gain, adFactor, offset, roNoise,L,U);
%     cdfEmccd(cdfEmccd<0.25) = nan;
    cdfEmccd(cdfEmccd>cdfThresh) = nan;

%     nbg =   arrayfun(@(x) sum(sortTruncI<=binPos(x))/cdfEmccd(x), 1:length(binPos));

    cdfVals = arrayfun(@(x) cdfEmccd(x), 1:length(binPos));
    nVals = arrayfun(@(x) sum(sortTruncI<=binPos(x)), 1:length(binPos));
%     figure,plot(cdfVals,nVals)

    nVals = nVals(~isnan(cdfVals));
    cdfVals =cdfVals(~isnan(cdfVals));

%     cdfVals(cdfVals(find(~isnan(nbg))))
%     nbg= nbg(find(~isnan(nbg)));
%     lineVals = polyfit(cdfVals,nVals,1);

    lineVals = cdfVals' \ nVals';

%     lineVals = sum(nVals)/sum(cdfVals);

%     lineVals = polyfit(1:length(nbg),nbg,1);
% figure,plot(cdfVals,cdfVals*lineVals(1)+lineVals(2))
% hold on
% plot(cdfVals,nVals)
% 

    logL = 1/length(cdfVals)*abs(nansum((cdfVals*lineVals(1)-nVals).^2));

%     logL = -1/length(cdfVals)*abs(nansum((cdfVals*lineVals(1)+lineVals(2)-nVals).^2));
%     logL = -abs(lineVals(1));

%               % Now calculate the estimate of FDR and FOR:
%         [forEst, stats ] = estimate_stats(intensities, ...
%              sortI(m),cdfSorted,qStar,U);

  
end