function [lambdaBg, pci,pdf,cdf] = est_lambda(histAll,lamGuess,Nthresh,gain,adFactor,offset,roNoise, LU, opt)
    %   est_lambda - estimating lambda_bg function
    % 
    %   Returns:
    %       lambdaBg - lambda using MLE
    %       pci - predicted confidence intervals


    binPos = LU(1):Nthresh;
    import Core.log_likelihood_trunc_dist; 
    logL = @(lambda,data,cens,freq,trunc) log_likelihood_trunc_dist(lambda,data,cens,freq,trunc, ...
                       gain, adFactor, offset, roNoise);
            
    % use mle with negative loglikelihood (nloglf)
    [lambdaBg, pci] = mle(binPos,'nloglf',logL,'start',lamGuess,'lowerBound',0,'Frequency',histAll(binPos),'TruncationBounds',LU,'Options',opt);

% %     %   If plot sample values:
%         logL2 = @(x) logL(x,binPos',[],histAll(binPos)',LU);
%         vec = 30:0.1:45;
%         figure,plot(vec,arrayfun(@(x) logL2(x),vec))
%
    cdf = nan(1,length(histAll)-1);
    pdf = nan(1,length(histAll)-1);
    import Core.pdf_cdf_emccd;

    [pdf(LU(1):length(histAll)-1), cdf(LU(1):length(histAll)-1)] = pdf_cdf_emccd(LU(1):length(histAll)-1, lambdaBg, gain, adFactor, offset, roNoise, LU(1), LU(2));


end

