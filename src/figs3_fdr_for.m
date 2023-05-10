function [] = figs3_fdr_for(imageNames, SNRVals, chipPars, outFig)

    % controls p(white|bg) (p_thresh)
    p_white_bg = 1E-2;

    % input parameter values
     
         
    % output filename
    outputFilename = 'fdr_for_results_temp';
    
    imageName = filenames{1}{1}(1);
    chipPars = chipParsCur;
    % we use single image
    imageName = strrep(imageName,'.tif','.mat');
    data = importdata(imageName{1}); 
    im = data.image;      
    images.imAverage = reshape(double(im),size( data.groundTruthImage));
    images.imageName = imageName;
    snr = data.snr;  % signal-to-noise ratio
    
    % Extract chip, optics parameters, etc from input file
    lambdaBgGroundTruth = data.lambdabg;
    %         lambdaSigGroundTruth = data.lambdasig;
    
    %         groundTruthPositions = data.placements;      
    groundTruthImage = data.groundTruthImage > 0;

    qStar = 0.5;
             % Maximum-likelihood estimation for background strength. Todo
         % here: use estimated values, otherwise this is not
         % realistic
         disp('Estimating lambda_bg.');
         % todo: select for specific gain
%          chipParsCurrent = 
%%

qStar = 0.001:0.001:0.01;

tic
FDR = zeros(1,length(qStar));
FOR = zeros(1,length(qStar));
FDR2 = zeros(1,length(qStar));
FOR2 = zeros(1,length(qStar));
Npar = zeros(1,length(qStar));
gain = chipPars.gain;
adFactor = chipPars.adFactor;
roNoise = chipPars.roNoise;
offset = chipPars.countOffset;
    chipPars = chipParsCur;
% % 
% chipPars.adFactor = 36; % ADU
% chipPars.countOffset = 27; % offset (bias)
% chipPars.roNoise = 1.44; % noise std
% chipPars.gain = 20;
alphaStar = 0.1;
[lambdaBg,intThreshBg,stats] = emccdpia_estimation(chipPars,'out.jpg',images,alphaStar,1);


%     [lambdaBg,intThreshBg] = fig2_calibration(chipPars,[],images,0.8,1);

structRes_1 = [];
structRes_2 = [];

for i = 1:length(qStar)
%     [lambdaBg_1,intThreshBg_1,structRes_1{i}] = fig2_calibration(chipPars,[],images,0.9,1);
%      Npar(i) = intThreshBg_1;
%      % Binarize image
%     [binarizedImage , intThreshBlackWhite ] = binarize_image_pval_thresh2(...
%         images.imAverage, p_white_bg ,lambdaBg_1 , gain, adFactor, offset, roNoise,intThreshBg_1);
%      disp(['intensity threshold value = ',num2str(intThreshBlackWhite)]);
%     
%     [fpr,fnr,tmr,FDR(i),FOR(i)] = compare_regions_to_ground_truth_beads_pixelbased(binarizedImage,groundTruthImage);      


    [lambdaBg_2,intThreshBg_2,structRes_2{i}] = fig2_calibration_ROC(chipPars,[],images,qStar(i),1);

%     [lambdaBg,intThreshBg] = fig2_calibration_finetuned(chipPars,[],images,qStar(i),0);
%     [lambdaBg_2,intThreshBg_2,structRes_2{i}] = fig2_calibration_FOR_control(chipPars,[],images,qStar(i),0);

    Npar2(i) = intThreshBg_2;


     % Binarize image
    [binarizedImage2 , intThreshBlackWhite2 ] = binarize_image_pval_thresh2(...
        images.imAverage, p_white_bg ,lambdaBg_2 , gain, adFactor, offset, roNoise,intThreshBg_2);
     disp(['intensity threshold value = ',num2str(intThreshBlackWhite2)]);
    
    [fpr,fnr,tmr,FDR2(i),FOR2(i)] = compare_regions_to_ground_truth_beads_pixelbased(binarizedImage2,groundTruthImage);      

end

%     % Determine p(black|bg), p(black|signal) and misclassification rate
%     % at the chosen threshold
%     [pBlackBgOptimal , pBlackSignalOptimal , misClassRateOptimal,FDR(i),FOR(i)] = estimate_pblack_quick(images.imAverage, ...
%         intThreshBlackWhite, intThreshBg , lambdaBg , gain,adFactor,offset,roNoise   );
%     disp(['Estimated p(white|bg)  = ',num2str(1-pBlackBgOptimal)]);
%     disp(['Estimated p(black|signal) = ',num2str(pBlackSignalOptimal)]); 
%     disp(['Estimated total misclassification rate = ',num2str(misClassRateOptimal)]);   
%     disp(' ')  
figure,plot(qStar,FOR)
hold on
plot(qStar,FDR)
legend({'FOR','FDR'})
xlabel('q^*')
title('FDR control method 1 (BH procedure)')


figure,plot(qStar,FOR2)
hold on
plot(qStar,FDR2)
% plot(qStar,qStar)
legend({'FOR','FDR'})
xlabel('q^*')
title('FDR control method new')
%%
        % Plot binarized image
% %         figure
% %         imshow(reshape(binarizedImage,size(groundTruthImage)),'InitialMagnification','fit')
% %         titleStr = ['Binarized image, SNR = ',num2str(snr(j))]; 
% %         title(titleStr)
%               
%         % Plot histogram along with fit
%    
% %         % Histogram
% %         figure     
% %         binPos = binEdges(1:end-1) + diff(binEdges)/2;
% %         bar(binPos,histAll); 
% %         hold on
%             
%         % Fit
% %         N=2^8; % number of integration points for CDF evaluation
% %         nIntVals =100;  % number of intensity values at which we plot the fit
%         % Extract relevant information
%         imageInput = images.imAverage;
%         nPixels = numel(imageInput); 
%         % Estimate the number of bg pixels
%         idxThreshTrueBg = length(find(imageInput(:) <= intThreshBg));      
%         [~,cdfAtThreshTrueBg] = pdf_cdf_from_characteristic_fun(intThreshBg+0.5,lambdaBg,gain,adFactor,offset,roNoise);  
%         nBg = min(idxThreshTrueBg/cdfAtThreshTrueBg, nPixels); 
%         % Estimate the expected number of bin counts 
%         % based on the fit parameters
%         [~,cdfEmccdNew] = pdf_cdf_emccd(binEdges,lambdaBg,chipPars,N);
%         binCountsFit = nBg*diff(cdfEmccdNew);
%         binPos = binEdges(1:end-1) + diff(binEdges)/2;
%         plot(binPos,binCountsFit,'-','Color','red','LineWidth',1)
%         
%         % Set some labels
%         xlabel('$n_{ic}$','Interpreter','latex','FontSize',12)
%         ylabel('Counts','Interpreter','latex','FontSize',12)
%         titleStr = ['Image count histogram, SNR = ',num2str(snr(j))]; 
%         title(titleStr)


        % Calculate performance results 
        [fpr,fnr,tmr] = compare_regions_to_ground_truth_beads_pixelbased(binarizedImage,groundTruthImage);      
        ppFPRPixelbased(j) = fpr;
        ppFNRPixelbased(j) = fnr;
        ppTmrPixelbased(j) = tmr;
        
        % Estimates
        fprEstimate(j) = 1 - pBlackBgOptimal;
        fnrEstimate(j) = pBlackSignalOptimal;
        tmrEstimate(j) = misClassRateOptimal;


end

    [lambdaBg_2,intThreshBg_2,structRes_2{i}] = fig2_calibration_ROC(chipPars,[],images,0.1,1);
