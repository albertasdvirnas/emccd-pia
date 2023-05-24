function [lambdaBg,intThreshBg,structRes ] = ...
                    estimate_lambda(intensities, gain, adFactor, countOffset, roNoise)
    % estimate_lambda 
    %
    % determines lambdaBg (Poisson parameter for background)
    % and the largest intensity cut-off intThreshBg so thate goodness-of-fit chi^2 test
    % is successful

    % Input:
    %
    % intensities = matrix or vector with intensity values
    % chipPars = struct containing the chip parameters
    % qStar = parameter (q^*) which controls the false detection rate (FDR)
    %         [q^* should be "rather" close to 1 to make sure that there are
    %          essentially no pixels below intThreshBg. 
    %          However, if "too close" to 1 the estmated number of
    %          background pixels may be unrealiable (see estimate_n_bg.m)
    %          and, as a consequence, the amplitude in the fit may not be
    %          accurate].
    %
    % Output:
    %
    % lambdaBg = lambda-parameter for Poisson distribution
    % intThreshBg = intensity threshold (integer) below which 
    %               most pixels are background
    % structRes - results structure
    %

    % other parameters in the algorithm:
    %
    %   lowestIntThresh - lower limit on the intThreshBg (25% of data)
    %   sortI - sorted data intensities
    %   Nthresh - initial estimate for intThreshBg for quick estimation of
    %   lambda
    %   sortTruncI - truncated data
    %   lamGuess - lambda guess
    %   histAll - histogram counts for data

    opt = statset('MaxIter',200,'MaxFunEvals',400,'FunValCheck','on');
    % Type: statset('mlecustom')   for additional fields

    if nargin < 7
        method = 'FOR';
    end

    lowestIntThresh = ceil(quantile(intensities,0.25)); 
    structRes.lowestIntThresh = lowestIntThresh;
    %
    sortI = sort(intensities);
    S = numel(sortI);
    Nthresh = sortI(round(S/2));

%     m = find(sortI>Nthresh,1,'First')-1; 
%     sortTruncI = sortI(1:m);
    % Fit lambda based on data
    lamGuess = abs((sortI(round(end/2)) - countOffset)/(gain/adFactor));

    % Prep data for truncated fit
    import Core.calc_bounds;
    [L, U, EX, STD] = calc_bounds(lamGuess, gain, adFactor, countOffset, roNoise, 6);
    structRes.LU = [max(1,ceil(L)) floor(U)];


    % Get bin edges
    binCenters = [max(1,ceil(L)):1:floor(U)+1 inf]; % bin edges shifted by half  

    % histogram for intensities. Calculated only once!
    histAll = zeros(1,binCenters(end-1));
    histAll(binCenters(1:end-1)) = histcounts(sortI,binCenters-0.5)';
    
    structRes.histAll = histAll;

    %% CHI^2 test
    nQuantiles = 5; % number quantiles
    intVals = lowestIntThresh:structRes.LU(2); 
    % calculate chi^2 test to check which bin sizes are ok
    lambdaBgMLE = nan(1,max(intVals));
    chi2Score = nan(1,max(intVals));
    distVals = cell(1,max(intVals));

    import Core.chi2_calc;
    for Nthresh = intVals;
        [lambdaBgMLE(Nthresh), pci, pdf, cdf] = est_lambda(histAll,lamGuess, Nthresh, gain, adFactor, countOffset, roNoise,structRes.LU,opt);
%         [chi2Score(Nthresh)] = chi2_calc(histAll, pdf, cdf, Nthresh);
        [chi2Score(Nthresh)] = chi2_calc(histAll, pdf, cdf, Nthresh,structRes.LU, nQuantiles);

        distVals{Nthresh}.pdf = pdf;
        distVals{Nthresh}.cdf = cdf;

    end
% figure,plot(chi2Score)
% [a,b] = min(chi2Score)

pval = chi2pdf(chi2Score,nQuantiles-2 );%/length(intVals); 

% figure,plot(pval)
% 
% figure,plot(pval>nanmean(pval))



    pvalThresh = pval > 0.01;
    intThreshBg = find(pvalThresh==1,1,'last');
    if isempty(intThreshBg)
        [maxPval,intThreshBg] = min(chi2Score);
        structRes.passthresh = 0;
        warning('pvalues too low');
    else
        structRes.passthresh = 1;
    end
    % 
%     [maxPval,intThreshBg] = min(chi2Score);
    lambdaBg = lambdaBgMLE(intThreshBg);

    structRes.pval = pval;
    structRes.chi2Score = chi2Score;
    structRes.lambdaBgMLE = lambdaBgMLE;
    structRes.pdf = distVals{intThreshBg}.pdf;
    structRes.cdf =  distVals{intThreshBg}.cdf;

% 
%     binPos = binCenters(1:end-1) + diff(binCenters)/2;
% 
% 
% 
%    
%     % Specify likelihood function
%     import Core.log_likelihood_trunc_dist; 
% %     logL = @(data,lambda) log_likelihood_trunc_dist(data , lambda , ...
% %                        gain, adFactor, countOffset, roNoise);
% 
% 
%     import Core.line_trunc_dist; 
%     lineL = @(data,lambda) line_trunc_dist(data , lambda , ...
%                        gain, adFactor, countOffset, roNoise, lowestIntThresh);
%     % specify function for calc for sorted values
%     import Core.p_values_emccd_sorted;
%     import Core.calc_bh_threshold;
% 
% %             vv = 10:0.1:60;
% %             dat=  arrayfun(@(x) logL(sortTruncI,x),vv);
%    
%     r = gain/adFactor;
%     % Sort data
%      m = S;
%     nOutliers = round(m/2); % initial amount of outliers
%     %         figure,plot(vv,dat)
%     %         
% % 
% %     vv = 10:0.1:60;
% %     dat=  arrayfun(@(x) lineL(sortI(sortI<50),x),vv);
% %     figure,plot(vv,dat)
% %     set(gca, 'YScale', 'log')
% % 
%      import Core.for_based_thresh;
%      import Core.fdr_based_thresh;
%      % Recursively reduce the data set until there are no "outliers" left.
%     hasOutliers = 1;
% %     nOutliers = 0;
%     runs = 0;
%     diffLambdas = Inf;
%     lambdaPrev = 0;
%     
%     structRes = [];
%     %     while hasOutliers && runs < 10 && diffLambdas > 0.00001
%     runs = runs + 1;
%     % Remove outliers
%     %         m = S - nOutliers;
% 
% 
%     %            tic
%     %     [threshVal,posMax] = run_ktest(sortI,lamGuess,gain, adFactor, countOffset, roNoise);
%     %         toc
%     threshVal = Nthresh;
% 
%     m = find(sortI>threshVal,1,'First')-1; % find first element greater than Nthresh
%     sortTruncI = sortI(1:m);
%  
% 
%  
% 
%     logL = @(lambda,data,cens,freq,trunc) log_likelihood_trunc_dist(lambda,data,cens,freq,trunc, ...
%                        gain, adFactor, countOffset, roNoise);
%     
% 
% %     [lambdaBg(Nthresh), ~, L, U, ~, ~,sortTruncI] = est_lambda(sortI, Nthresh, gain, adFactor, countOffset, roNoise, r);
% 
%      
% 
%     % use mle with negative loglikelihood (nloglf)
%     tic
%     [lambdaBgMLE, pci] = mle(binPos,'nloglf',logL,'start',lamGuess,'lowerBound',0,'Frequency',histAll,'TruncationBounds',structRes.LU,'Options',opt);
%     toc
% 
%     structRes.intensitiesU = [max(1,ceil(L)):min(Nthresh,floor(U))]; % take edge point of each box, so +0.5
%     structRes.cdf = nan(1,max(structRes.intensitiesU));
%     structRes.pdf = nan(1,max(structRes.intensitiesU));
% 
%     % calc cdf/pdf for all of these
%     [structRes.pdf(structRes.intensitiesU), structRes.cdf(structRes.intensitiesU)] = pdf_cdf_from_characteristic_fun(structRes.intensitiesU,lambdaBgMLE,gain,adFactor,countOffset,roNoise);
% 
% %                 tic
% %         [lambdaBg, pci] = mle(histAll,'logpdf',logL,'start',lamGuess,'lowerBound',0,'Options',opt);
% %         toc
% %         tic % alternative
% %         lambdaBg2 = fminsearch(@(x) lineL(sortTruncI,x),lamGuess)
% %         toc
% %         tic
% %        [lambdaBg, pci] = mle(sortTruncI,'logpdf',lineL,'start',lamGuess,'lowerBound',0,'Options',opt);
% %         toc
%         
%         % calculate nBG (lower bound)
% 
%         diffLambdas = abs((lambdaBgMLE - lambdaPrev)/lambdaPrev);
%         lambdaPrev = lambdaBgMLE;
%          [~,U,~,~] = calc_bounds(lambdaBgMLE,gain,adFactor,countOffset,roNoise);
%         lambdaBg = 38;
%         switch method
%             case 'FOR'
%                 [nOutliers,hasOutliers,stats, Nthresh] =  for_based_thresh(intensities, lambdaBgMLE,gain,adFactor,countOffset,roNoise,qStar,U,lowestIntThresh);
%             case 'FDRest'
%                 [nOutliers,hasOutliers,stats, Nthresh] =  fdr_based_thresh(intensities, lambdaBgMLE,gain,adFactor,countOffset,roNoise,qStar,U,lowestIntThresh);
% 
%             case 'FDR'
%                 [threshold, nOutliers, outliers, hasOutliers,Nthresh ] = calc_bh_threshold(sortI, lambdaBgMLE, gain, adFactor, countOffset, roNoise,qStar);
%                 stats = [];
% %                 forEst = 
%             otherwise
%         end
%         % save outputs
%         structRes.lamGuess(runs) = lamGuess;
%         structRes.nOutliers(runs) =  nOutliers;
%         structRes.lambdaBg(runs) =  lambdaBgMLE;
%         structRes.threshold(runs) =  sortI(S - nOutliers);
%         structRes.Nthresh(runs) = threshVal;
%         structRes.stats{runs} = stats;
% %     end
%     idxBgEstimation = S - nOutliers;
%     intThreshBg = threshVal;
end


%     intThreshBg = sortI(idxBgEstimation);
%     intThreshBg = round(intThreshBg);  % why -1?
%     
%    f= figure,
%     tiledlayout(1,2)
%     nexttile 
%     hold on
%     plot(structRes.nOutliers/S)
%     title('a) Fraction of outliers','Interpreter','latex')
%     xlabel('Number of iterations')
%     nexttile 
%     plot(structRes.lambdaBg)
%     title('$\hat \lambda_{bg}$','Interpreter','latex')
%        print('FigS7.eps','-depsc','-r300')

% %    
%     [~, ~, cdfEnd, ~] = p_values_emccd_sorted(max(sortTruncI)+1,lambdaBg,gain, adFactor, countOffset, roNoise);        
%     figure,histogram(sortTruncI,'normalization','pdf')
%     hold on
%     plot(intUnique,pdfUnique/cdfEnd)

            % subtract 1 to make sure the threshold does not include
            % extra signal intensities
%    f=histogram(pVals,'normalization','count')
%    xlabel('P-value')
%    
%    print('pvaldist.png','-depsc','-r300')
% 


%         figure,histogram(sortTruncI,'normalization','pdf')
%         hold on
%         plot(intUnique,pdfUnique)
%         xlim([0 300])
%         
%         histogram(sortTruncI,'normalization','cdf','DisplayStyle','stairs')
%         hold on
%         plot(intUnique,cdfsUnique)

%% MLE alternative
% % 
% firstIntensity = quantile(intensities,0.1)+0.5; % stop intensity for calculating nbg
% lastIntensity = quantile(intensities,0.9)+0.5; % stop intensity for calculating nbg
% 
% lastIntensity = 85+0.5;
% 
% lval = nan(1,lastIntensity-0.5);
% logvals = nan(1,lastIntensity-0.5);
% 
% [L, U, EX, STD] = calc_bounds(lamGuess,gain,adFactor,countOffset,roNoise);
% 
% for ii=firstIntensity:min(floor(U),lastIntensity)
%     m = find(sortI>ii,1,'First')-1; % find first element greater than Nthresh
%     sortTruncI = sortI(1:m);
%     % Fit lambda
%     lamGuess = abs((sortTruncI(round(m/2)) - countOffset)/r);
% 
%     % Get bin edges
%     binEdges = [ceil(L):1:ii+1]-0.5; % bin edges shifted by half  
%     % histogram for intensities we can have this outside of the function
%     histAll = histcounts(sortTruncI,binEdges)';
%     
%     binPos = binEdges(1:end-1) + diff(binEdges)/2;
% 
% %      [threshVal,posMax] = run_ktest(sortI,lamGuess,gain, adFactor, countOffset, roNoise);
% 
%     tic
%     [lval(ii-0.5), pci] = mle(binPos,'nloglf',logL,'start',lamGuess,'lowerBound',0,'Frequency',histAll,'TruncationBounds',[L U],'Options',opt);
%     toc
% 
% %     [lval(ii-0.5), pci] = mle(sortTruncI,'logpdf',logL,'start',lamGuess,'lowerBound',0,'Options',opt);
%     toc
%     logvals(ii-0.5) = logL(lval(ii-0.5),binPos',[],histAll,[L U]);
% end
% 
% 
% 
% figure,plot(logvals)
% figure,plot(lval)
% 
% import Core.log_likelihood_trunc_dist; 
%  [L, U, EX, STD] = calc_bounds(lamGuess,gain,adFactor,countOffset,roNoise);
%     logL = @(data,lambda) log_likelihood_trunc_dist(data , lambda , ...
%                        gain, adFactor, countOffset, roNoise,L, U);
% 
% 
% vv = 36:0.1:40;
% 
% dat=  arrayfun(@(x) logL(x,binPos',[],histAll,[L U]),vv);
% figure,plot(vv,dat)
% % [a,b] = max(dat)
% 
% vv = 32:0.1:44;
% 
% 
%     import Core.calc_bounds;
%     [L, U, EX, STD] = arrayfun(@(x) calc_bounds(x,gain,adFactor,countOffset,roNoise),vv,'un',true);
% 
% figure,plot(vv,floor(L))
% figure,plot(vv,floor(U))