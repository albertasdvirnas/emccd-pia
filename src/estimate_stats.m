function  [stats,allVals ] = estimate_stats(intensities, cdfFun, U, stopIntensity)
    % estimate_stats
    %   
    % Estimation of p(black|bg) and p(black|signal), as well as other stats without split histogram
    %
    % Input:
    %
    % intensities = matrix or vector with pixel intensity values
    % intThreshBlackWhite = choice of intensity threshold (integer) 
    % intThreshBg = intensity threshold (integer) below which there 
    %               should be only bg pixels.
    % lambdaBg = lambda-parameter for background Poisson distribution
    % chipPars = struct containing the chip parameters
    %     cdfFun - cdf function
    %   
    % Output:
    %
    % pBlackBg = probability that a "black" pixel (at threshold set by inThreshBg) 
    %          is a background pixel [pWhiteBg = 1-pBlackBg]
    % pBlackSignal = probability that a "black" pixel is a signal pixel
    %              [pWhiteSignal = 1-pBlackSignal] 
    % misClassRate = estimate of the miclassification rate at the chosen
    %                intensity threshold
    %
    %  stats - FPR, FNR, FDR, FRR
    %
    % Dependencies: pdf_cdf_from_characteristic_fun
    %               emccd_distribution/pdf_cdf_emccd.m
    %
    % J. D. M. Rennie, L. Shih, J. Teevan, and D. R. Karger, “Tackling
    % the poor assumptions of naive bayes text classifiers,
    
    % Hard-coded variables
    alpha = 1;    % pseudocount used to avoid zeros 
                  %  when calculating probabilities
                  
    if nargin < 5
        U = max(intensities);
    end
%     cdfFun = @(x)  pdf_cdf_from_characteristic_fun(x+0.5,lambdaBg,gain,adFactor,offset,roNoise);
    % Pre-processing  
    intensityVec = intensities(:);
    uniqueVals = unique(intensityVec);
    nPixels = numel(intensityVec);
    allVals = stopIntensity:min(U,uniqueVals(end));
    
%     intThresh = max(intensityVec)
    % Estimate number of bg pixels nBg   (eq 12) manuscript
%     [cdfval] = cdfFun([uniqueVals(1):intThresh]); % cdf at these intensities
%     cdfPercentage = 0.25; % 1st inter-quartile
%     val = find(cdfval>cdfPercentage,1,'first');%find(cdfFun>0.5,1,'first');% find(cdfFun>1,0.9,'first');

    % first estimate nBG (up to cdfThresh which gives 0.25 of intensities)
%     [cdfval] = cdfFun([uniqueVals(1):stopIntensity]); % cdf at these intensities
    %% nbg - fitting line, not accurate
        cdfIntensities = uniqueVals(1):stopIntensity-0.5;

    stats.yval =    arrayfun(@(x) sum(intensities<=x), cdfIntensities-1);
    stats.xval =    arrayfun(@(x) cdfFun(x), cdfIntensities-1);
    nBg =  stats.xval'\stats.yval';
    stats.nBg = nBg;
 

    cdfThresh = cdfFun(stopIntensity-0.5);

%     nBg = 100000;
    stats.xval =    arrayfun(@(x) cdfFun(x), cdfIntensities-1);


%     [~,cdfEmccdNew] = pdf_cdf_emccd(binEdges,lambdaBg,gain,adFactor,countOffset,roNoise,L,U);
difVals = diff(stats.xval);
binCountsFit = nBg*diff(stats.xval);
binPos = cdfIntensities(1:end-1);


% figure

% histAll = histcounts(intensities(intensities<max(cdfIntensities-0.5)),cdfIntensities-0.5);    
% h1 = bar(binPos,histAll,1); 
% hold on
% plot(binPos,binCountsFit,'--','Color','black','LineWidth',2)
% figure,plot(binCountsFit-histAll)

% check max for last value
% shifts = binCountsFit-histAll;
% [maxShift,pos] = max(abs(shifts));

% nBg = (binCountsFit(end)-(sign(shifts(pos)))*maxShift)/difVals(end);

% nbgFun = @(x) max(binCountsFit-x*diff(stats.xval));
% figure,plot(10^4:10^4:10^6,arrayfun(@(x) nbgFun(x),10^4:10^4:10^6))
% 
% % nbgFun = @(x) max(abs(histAll-x*diff(stats.xval)));
% nBg = fminunc(nbgFun,nBg);

% figure,
% plot(binPos,nBg*diff(1-stats.xval),'--','Color','black','LineWidth',2)


    % (eq 13) Total number of signal pixels
%     nSignal = max(0,nPixels - nBg);
    
    %  (eq 14) Number of signal pixels which are "black" (below the intensity
    %  threshold)
%     nBlack = length(find(intensityVec <= intThresh));

    % Now we vary intThresh to get the desired FOR/FDR
    fdr = ones(1,length(allVals));
    FOR = zeros(1,length(allVals));

    cdfs =  cdfFun(allVals-0.5);

     for i=1:length(allVals)
        fdr(i) = (nBg*(1-cdfs(i)))/sum(intensityVec>allVals(i));
        FOR(i) = (max(0,sum(intensityVec<=allVals(i))-nBg*cdfs(i)))/(sum(intensityVec<=allVals(i)));
     end
%         if cdfs(i) > cdfThresh % if enough data points

%         else
%             FOR(i) = 0;
%             fdr(i) = 1;
%         end

        % FOR = FN/(FN+TN) 
        % FN+TN = sum(intensityVec<=intThreshBg)
        % FN = sum(intensityVec<=intThreshBg)-nBg*cdfAtBG


        % eq 15
%         n_black_bg = nBg*cdfs(i);
        
        % eq 16 n(black|s)
%         n_black_sig = max(0,nBlack - n_black_bg);
%         n_black_sig = min(n_black_sig,nSignal);
%         
%         % p(black|bg) and p(black|signal)
%         % eq 17
%         p_black_bg = (n_black_bg + alpha)/(nBg + 2*alpha);
%         p_black_sig = (n_black_sig + alpha)/(nSignal + 2*alpha);
%         
%         % Misclassification rate
%         % eq 18?
%         fBg = (nBg + alpha)/(nPixels + 2*alpha);  % ratio of black pixels
%         stats.misClassRate(i) = (1-p_black_bg).*fBg + p_black_sig*(1-fBg); % ra
%         
%         %     [~,cdfAtBG] = pdf_cdf_from_characteristic_fun(intThreshBg,lambdaBg,gain,adFactor,offset,roNoise);
%         
% 
%      end
     stats.fdr = fdr;
     stats.FOR = FOR;


end
