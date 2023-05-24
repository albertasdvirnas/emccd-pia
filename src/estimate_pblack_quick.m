function  [pBlackBg , pBlackSignal , misClassRate, FDR, FOR ] = estimate_pblack_quick(intensities, ...
             intThreshBlackWhite, intThreshBg , lambdaBg , gain,adFactor,offset,roNoise   )

    %   
    % Estimation of p(black|bg) and p(black|signal) without split histogram
    %
    % Input:
    %
    % intensities = matrix or vector with pixel intensity values
    % intThreshBlackWhite = choice of intensity threshold (integer) 
    % intThreshBg = intensity threshold (integer) below which there 
    %               should be only bg pixels.
    % lambdaBg = lambda-parameter for background Poisson distribution
    % chipPars = struct containing the chip parameters
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
                  
                  
    % Pre-processing  
    intensityVec = intensities(:);
    nPixels = numel(intensityVec);
    
    % Number of background pixels. For synthetic images we know true
    % numbers
    nBackTruePart = length(find(intensityVec(:) <= intThreshBg));
    [~,cdfAtThreshTrueBg] = pdf_cdf_from_characteristic_fun(intThreshBg+1,lambdaBg,gain,adFactor,offset,roNoise);

    %     (eq 12) manuscript
    nBg = min(round(nBackTruePart/cdfAtThreshTrueBg), nPixels); 


    
    % from test_fdr_pixel
%     val = round(min(intensityVec))+5:round(max(intensityVec));
%     nbgv = zeros(1,length(val));
%     for i=1:length(val)
%         [~,cdfAtThreshTrueBg] = pdf_cdf_from_characteristic_fun(val(i)+1,lambdaBg,gain,adFactor,offset,roNoise);
%         nbgv(i) = sum(intensityVec<=val(i))/cdfAtThreshTrueBg; % can we estimate std for this
%     end
% %     
%     figure,plot(val,nbgv);hold on
%     plot([intThreshBg intThreshBg],[min(nbgv) max(nbgv)],'red-')
%     ylim([min(nbgv) max(nbgv)])
%     xlabel('Intensity')
%     ylabel('Estimate bg pixels')
%     legend({'num of bg estimate','N thresh'})
%     
    % (eq 13) Total number of signal pixels
    nSignal = max(0,nPixels - nBg);
    
    %  (eq 14) Number of signal pixels which are "black" (below the intensity
    %  threshold)
    nBlack = length(find(intensityVec <= intThreshBlackWhite));
%     nBlack = numel(idx);   % number of black pixels
 
                         
     % Number of background pixels which are "black" (<= the intensity
     % threshold)
    [~,cdfAtThreshBlackWhite] = pdf_cdf_from_characteristic_fun(intThreshBlackWhite,lambdaBg,gain,adFactor,offset,roNoise);
    % eq 15
    nBlackBg = nBg*cdfAtThreshBlackWhite;
     
    % eq 16 n(black|s)
    nBlackSignal = max(0,nBlack - nBlackBg);
    nBlackSignal = min(nBlackSignal,nSignal);
   
    % p(black|bg) and p(black|signal)
    % eq 17
    pBlackBg = (nBlackBg + alpha)/(nBg + 2*alpha);
    pBlackSignal = (nBlackSignal + alpha)/(nSignal + 2*alpha);
   
    % Misclassification rate
    % eq 18?
    fBg = (nBg + alpha)/(nPixels + 2*alpha);  % ratio of black pixels
    misClassRate = (1-pBlackBg).*fBg + pBlackSignal*(1-fBg); % ra
    
    [~,cdfAtBG] = pdf_cdf_from_characteristic_fun(intThreshBg,lambdaBg,gain,adFactor,offset,roNoise);

     FDR = (nBg*(1-cdfAtBG))/sum(intensityVec>intThreshBg);

     % FOR = FN/(FN+TN) 
     % FN+TN = sum(intensityVec<=intThreshBg)
     % FN = sum(intensityVec<=intThreshBg)-nBg*cdfAtBG
    FOR = (max(0,sum(intensityVec<=intThreshBg)-nBg*cdfAtBG))/(sum(intensityVec<=intThreshBg));

end
