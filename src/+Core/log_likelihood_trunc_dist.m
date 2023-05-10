function [logL] = log_likelihood_trunc_dist(lambda,data,cens,freq,trunc,...
                             gain, adFactor, offset, roNoise)
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

    
%     binEdges
     % center position for each bin
%     binPos = binEdges(1:end-1) + diff(binEdges)/2;
   
    import Core.pdf_cdf_emccd;
    [pdfEmccd, cdfEmccd] = pdf_cdf_emccd(data',lambda,gain, adFactor, offset, roNoise, trunc(1), trunc(2));
        logL = -sum(freq.*log(pdfEmccd)) + sum(freq)*log(cdfEmccd(end));

    % 
%     [~ ,cdfEmccdEnd] = pdf_cdf_emccd(data(end)+0.5,lambda,gain, adFactor, offset, roNoise,trunc(1), trunc(2));

    % https://stats.stackexchange.com/questions/48897/maximum-likelihood-estimators-for-a-truncated-distributionlog-likelihood

%     logL = -sum(freq.*log(pdfEmccd)) + sum(freq)*log(cdfEmccdEnd);
    
 
  
end