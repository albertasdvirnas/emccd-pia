function [] = fig3_probabilitstic_segmentation(imageNames, SNRVals, chipPars, outFig)
    
%     if nargin < 1
%         inputFolder = 'synthSingleFrame\100x\';
% %         calc = 0;
%     end
    
    % input parameter values
     qStar = 0.5; % parameter which controls the false detection rate (FDR)
     pValThreshBinarization = 1E-2;
     
     
  % p-value threshold using the background distribution
     %outputFilename = 'eval_binarization_results';
     outputFilename = 'eval_binarization_results_temp';
%      inputFolder = './';     
                       % folder where we store the images to be analyzed

                       
                       

%     probabilistic_binarization()
%     if calc
% if need to generate simulated data at different SNR ratio's:
% simulate_random_beads_images(chipPars)


    % Input parameters
    
%      SNRVals = ["100","150","200","250","300","350","400","450","500","550", ...
%              "600","650","700","750","800","850","900","950","1000"];
%      SNRVals = ["300","500"];
  
     
    % Assign memory for FDR and FRR, etc.
    nImages = length(imageNames);
    lambdaBgSave = zeros(1,nImages);
    intThreshBgSave = zeros(1,nImages);
    ppFdr = zeros(1,nImages);
    ppFrr = zeros(1,nImages);
    ppFdrPixelbased = zeros(1,nImages);
    ppFrrPixelbased = zeros(1,nImages);
    fdrEstimate = zeros(1,nImages);
    frrEstimate = zeros(1,nImages);
    tmrEstimate = zeros(1,nImages);
    otsuFdr = zeros(1,nImages);
    otsuFrr = zeros(1,nImages);
    otsuFdrPixelbased = zeros(1,nImages);
    otsuFrrPixelbased = zeros(1,nImages);
  


    for j = 1:nImages 
        
        
        disp('-------------')
        disp(['image with ','SNR = ', num2str(SNRVals(j))])
        disp('------------- ')
        
        
        % Pre-processing   
%         imageNames = filenames{1}{1};
        
        imageName = strrep(imageNames{j},'.tif','.mat');
%         path = strcat(inputFolder,imageFile);
        data = importdata(imageName); 
        im = data.image;      
        images.imAverage = reshape(double(im),size( data.groundTruthImage));
        images.imageName = imageName;
        snr(j) = data.snr;  % signal-to-noise ratio

        
        % Extract chip, optics parameters, etc from input file
        lambdaBgGroundTruth = data.lambdabg;
%         lambdaSigGroundTruth = data.lambdasig;

        
%         groundTruthPositions = data.placements;      
        groundTruthImage = data.groundTruthImage > 0;
    
%         % ----- Original image ------
%         figure
%         imshow(mat2gray(im),'InitialMagnification','fit')
%         titleStr = ['Original image, SNR = ',num2str(snr(j))]; 
%         title(titleStr)
%       
%         
%         % ----- Ground truth image ------
%         figure
%         imshow(groundTruthImage,'InitialMagnification','fit')
%         titleStr = ['Ground truth image, SNR = ',num2str(snr(j))]; 
%         title(titleStr)
        
        %
        % ---- Photophysical image binarization  -----
        %
           
     
         % Maximum-likelihood estimation for background strength. Todo
         % here: use estimated values, otherwise this is not
         % realistic
         disp('Estimating lambda_bg.');
         % todo: select for specific gain
%          chipParsCurrent = 
        gain = chipPars.gain;
        adFactor = chipPars.adFactor;
        roNoise = chipPars.roNoise;
        offset = chipPars.countOffset;

%          chipPars = struct();
%          chipPars.gain = data.gain;
%          chipPars.adFactor = data.adFactor;
%          chipPars.roNoise = data.roNoise;
%          chipPars.countOffset = data.offset; 
%         chipPars.gain = 19;
%         chipPars.countOffset = 27;
%         chipPars.roNoise = 1.44;
%         chipPars.adFactor = 36;
        [lambdaBg,intThreshBg] = fig2_calibration(chipPars,[],images,qStar,1);
%          [lambdaBg , intThreshBg] = ...
%                         estimate_bg_params(images.imAverage(:),chipPars,qStar);
         lambdaBgSave(j) = lambdaBg;
         intThreshSave(j) = intThreshBg;
         disp(' ')

         % Binarize image
        disp('Binarizing image');
        [binarizedImage , intThreshBlackWhite ] = binarize_image_pval_thresh2(...
            images.imAverage, pValThreshBinarization ,lambdaBg , gain, adFactor, offset, roNoise);
         disp(['intensity threshold value = ',num2str(intThreshBlackWhite)]);


        % Determine p(black|bg), p(black|signal) and misclassification rate
        % at the chosen threshold
        [pBlackBgOptimal , pBlackSignalOptimal , misClassRateOptimal ] = estimate_pblack_quick(images.imAverage, ...
                intThreshBlackWhite, intThreshBg , lambdaBg , gain,adFactor,offset,roNoise   );
        disp(['Estimated p(white|bg)  = ',num2str(1-pBlackBgOptimal)]);
        disp(['Estimated p(black|signal) = ',num2str(pBlackSignalOptimal)]); 
        disp(['Estimated total misclassification rate = ',num2str(misClassRateOptimal)]);   
        disp(' ')    

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
        [fdr,frr,tmr] = compare_regions_to_ground_truth_beads_pixelbased(binarizedImage,groundTruthImage);      
        ppFdrPixelbased(j) = fdr;
        ppFrrPixelbased(j) = frr;
        ppTmrPixelbased(j) = tmr;
        
        % Estimates
        fdrEstimate(j) = 1 - pBlackBgOptimal;
        frrEstimate(j) = pBlackSignalOptimal;
        tmrEstimate(j) = misClassRateOptimal;


        
        fprintf('False detection rate : %f.\n',ppFdr(j))
        fprintf('False rejection rate : %f.\n',ppFrr(j))     
        fprintf('False detection rate (pixelbased): %f.\n',ppFdrPixelbased(j))
        fprintf('False rejection rate (pixelbased): %f.\n',ppFrrPixelbased(j)) 
        fprintf('Estimated false detection rate: %f.\n',fdrEstimate(j))
        fprintf('Estimated false rejection rate: %f.\n',frrEstimate(j))       
        fprintf('Estimated total misclassification rate: %f.\n',tmrEstimate(j))       


         % -----------------------------------------------
       
        
        %
        % ---- Otsu based binarization and segmentation ----
        % ---- (segmentation result not used) ----- 
        %
        segOutputOtsu = segmentation_otsu(images);
        signalLabelImOtsu = segOutputOtsu.labelIm;
        binarizedImageOtsu = segOutputOtsu.binarizedImage;
        
%         % Plot 
%         figure
%         imshow(binarizedImageOtsu,'InitialMagnification','fit')
%         titleStr = ['Otsu binarization image, SNR = ',num2str(snr(j))]; 
%         title(titleStr)
%       
      
         % Calculate performance results 
        [fdrOtsuTemp,frrOtsuTemp,tmrOtsuTemp] = compare_regions_to_ground_truth_beads_pixelbased(binarizedImageOtsu,groundTruthImage);
      
        otsuFdrPixelbased(j) = fdrOtsuTemp;
        otsuFrrPixelbased(j) = frrOtsuTemp;
        otsuTmrPixelbased(j) = tmrOtsuTemp;
                  
        fprintf('Otsu false detection rate : %f.\n',otsuFdr(j))
        fprintf('Otsu false rejection rate : %f.\n',otsuFrr(j))     
        fprintf('Otsu false detection rate (pixelbased): %f.\n',otsuFdrPixelbased(j))
        fprintf('Otsu false rejection rate (pixelbased): %f.\n',otsuFrrPixelbased(j)) 
   
         % -----------------------------------
       

    end

    % Save results 
    results = struct();
    results.snr = snr;
    results.lambdaBgGroundTruth = lambdaBgGroundTruth;    
    results.lambdaBgEst = lambdaBgSave;
    results.intThreshEst = intThreshBgSave;      
    results.ppFdrPixelbased = ppFdrPixelbased;
    results.ppFrrPixelbased = ppFrrPixelbased;
    results.ppTmrPixelbased = ppTmrPixelbased;
    results.fdrEst = fdrEstimate;
    results.frrEst = frrEstimate;
    results.tmrEst = tmrEstimate;
    results.otsuFdrPixelbased = otsuFdrPixelbased;
    results.otsuFrrPixelbased = otsuFrrPixelbased;
    results.otsuTmrPixelbased = otsuTmrPixelbased;
    save(outputFilename,'results')
    
    plot_binarization_evaluation_results_beads(outFig);

% else
% end

end



function simulate_random_beads_images(chipPars)

    %
    % Script for generating an image with randomly positioned
    % particles (point emitters).
    % 
    % Dependencies: generate_image_randomly_deposited_particles.m
    %
   
   
    % Input parameters
    SNRVals = [1.00 , 1.50 , 2.00 , 2.50 , 3.00 , 3.50 , 4.00 , 4.50 , ...
               5.00 , 5.50 , 6.00 , 6.50 , 7.00 , 7.50 , 8.00 , 8.50 , ...
               9.00 , 9.50 , 10.00 ];             
                              % signal-to-noise ratio  values 
    lambdaBg = 38.0;          % Poisson parameter for background
    circRadius = 20;          % circle radius (in pixels)  
    nRows = 2048;             % number of rows in image
    nCols = 2048;              % number of columns in image
    particleDensity = 1/8000;   % signal particle density
    % estimates for parameters from Mean Variance (MV) callibration
    gain = mean(chipPars.gain{3});
    adFactor = mean(chipPars.adFactor);
    countOffset = mean(chipPars.countOffset{3});
    roNoise = mean(chipPars.roNoise{3});

    for idxSNR = 1:length(SNRVals)
        
        SNR = SNRVals(idxSNR);
        fprintf('SNR: %.3f\n',SNR);
       
        % Translate SNR to lambda for signal regions.
        %  Comment: 
        % The signal-to-noise ratio is defined:
        %       SNR = lambdaSig /sqrt(lambdaSig + lambdaBg)
        % 
        % Inverting these expressions we get: 
        %      lambdaSig = SNR^2/2 + sqrt(SNR^2*lambdaBg + SNR^4/4)
        %
        
        lambdaSignal = SNR*(SNR+sqrt(SNR^2+4*lambdaBg))/2;
        
%         lambdaSignal = SNR^2/2 + sqrt(SNR^2*lambdaBg + SNR^4/4);

        % Generate the image
        [noisyImage,groundTruthImage, placements] =...
            simulate_fluorophores_emccd(lambdaSignal,lambdaBg,circRadius,nRows,nCols,particleDensity,...,
            gain, adFactor, countOffset, roNoise);

%          [noisyImage,groundTruthImage, placements] = ...
%              generate_image_randomly_deposited_beads(lambdaBg,circRadius,nRows,nCols,...
%              particleDensity, gain, adFactor, countOffset, roNoise);

%         % Plot 
%         figure
%         imshow(groundTruthImage)
%         title('Ground truth image');
% 
%         figure()
%         imshow(mat2gray(noisyImage));
%         title('Noisy image')


        % Store image and associated information
        data.image = noisyImage;
        data.groundTruthImage = groundTruthImage;
        data.lambdabg = lambdaBg;
        data.lambdasig = lambdaSignal;
        data.gain = gain;
        data.roNoise = roNoise;
        data.offset = countOffset;
%         data.waveLength = waveLength;
%         data.NA = NA;
%         data.pixelSize = pixelSize;
        data.adFactor = adFactor;
        data.placements = placements;
        data.snr = SNR;
        titleSNR = round(100*SNR);
        saveFileName = sprintf('testImageEMCCDBeadsSNR%i',titleSNR);
        save(saveFileName,'data')
        data.imageName = saveFileName;
        
    end
    

end



function  [pBlackBg , pBlackSignal , misClassRate ] = estimate_pblack_quick(intensities, ...
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
    % Dependencies: binarization_segmentation/estimate_n_bg.m
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
    
     
end
% 
% 
% function [binarizedImage , intThresh ] = binarize_image_pval_thresh(...
%               imageInput, pValThresh ,lambdaBg ,  gain,adFactor,offset,roNoise)
% 
%     % 
%     % Binarizes an image using a p-value threshold
%     %
%     % Input:
%     %
%     % imageInput = image
%     % pValThresh = p-value threshold. Pixels with a p-value 
%     %             below this threshold are turned white ( = 1)
%     %             Pixels with a p-value above this threshold 
%     %             are black ( = 0).
%     % lambdaBg = Poissoin parameter for the background
%     % chipPars = struct containing chip-parameters.
%     %
%     % Output:
%     %
%     % binarizedImage = binarized image 
%     % intThresh = intensity threshold at the specified p-value
%     %
%     % Dependencies: emccd_distribution/inverse_cdf_emccd.m
%     % 
%     
%     % Hard-coded variables
% %     N= 2^8;     % number of integration points for evaluating the EMCCD CDF
% %     tol=1E-5;   % accuracy for the inverse CDF calculation.
%     
% % inverse cdf to get intThresh
%     % first these are quickly calculated bounds
%     [L,U,EX,STD] = calc_bounds(lambdaBg,gain,adFactor,offset,roNoise);
%     
%     % that give intensities to calculate over
%     intensities = ceil(L):floor(U);
%     
%     % cdf, in this case do not need to be truncated
%     [pdf,cdf] = pdf_cdf_from_characteristic_fun(intensities,lambdaBg,gain,adFactor,offset,roNoise);
%     
%     % find the value where pvalue=1-cdf > pValThresh
%     intThresh = find(1-cdf < pValThresh,1,'first');
% %     intensities(find(1-cdf >pValThresh,1,'last'));
% 
%     intThresh = intensities(intThresh);
% 
% %     intThresh = inverse_cdf_emccd( 1-pValThresh , lambdaBg , chipPars , N , tol);
%     intThresh = floor(intThresh)-0.5;  % since intensities are integers 
%                                        % in experimental images, 
%                                        % we set the threshold to be a half-integer
%    
%     % Binarize using intensity threshold
%     binarizedImage = imbinarize(imageInput,intThresh);
% 
% 
% end

function plot_binarization_evaluation_results_beads(outFig3)

    %
    % Script for analyzing the output of the function 
    % evaluate_seg_performance_beads
    %
   
    % Load results from file 
    results = importdata('eval_binarization_results_temp.mat');
    snr = [results.snr];

    ppFdrPixelbased = [results.ppFdrPixelbased];
    ppFrrPixelbased = [results.ppFrrPixelbased];
    ppTmrPixelbased = [results.ppTmrPixelbased];

    fdrEst = [results.fdrEst];
    frrEst = [results.frrEst];
    tmrEst = [results.tmrEst];

    otsuFdrPixelbased = [results.otsuFdrPixelbased];
    otsuFrrPixelbased = [results.otsuFrrPixelbased];
    otsuTmrPixelbased = [results.otsuTmrPixelbased];

    % Sort according to SNR values
    [snr,idx] = sort(snr);

    ppFdrPixelbased = ppFdrPixelbased(idx);
    ppFrrPixelbased = ppFrrPixelbased(idx);
    ppTmrPixelbased = ppTmrPixelbased(idx);

    fdrEst = fdrEst(idx);
    frrEst = frrEst(idx);
    tmrEst = tmrEst(idx);

    otsuFdrPixelbased = otsuFdrPixelbased(idx);
    otsuFrrPixelbased = otsuFrrPixelbased(idx);
    otsuTmrPixelbased = otsuTmrPixelbased(idx);

    % Plot lambda_bg
   % Plot FDR, pixelbased 
    fS = 12;
    figure
    tiledlayout(2,2,'TileSpacing','tight');
 
    nexttile
    plot(snr,results.lambdaBgEst,'red-d','linewidth',2);
    hold on
    y = results.lambdaBgGroundTruth*ones(1,length(snr));
    plot(snr,y,'black--','linewidth',2)
    ylim([37 41])
    text(4,37.5,'Ground-truth')
    text(4,39.5,'Estimated')

    hold off
%     try
%     axis([snr(1) , snr(end) , results.lambdaBgGroundTruth - 1 , ...
%              results.lambdaBgGroundTruth + 1])
%     catch
%     end
%      legend('estimated value','ground truth','Location','northeast','Fontsize',fS)
%     xlabel('Signal to noise ratio','Interpreter','latex','Fontsize',fS)
    ylabel('\lambda_{bg}','Fontsize',12)
    title('a) Poisson parameter','Fontsize',fS)
%     set(gca,'Fontsize',15)
  
    
    nexttile
    plot(snr,ppFdrPixelbased,'-d','linewidth',2);
    hold on
    plot(snr,otsuFdrPixelbased,'--x','linewidth',2);
    plot(snr,fdrEst,'.-','linewidth',2)
    hold off
%     legend(["non-Otsu","Otsu's method","estimate"],'Location','east','Fontsize',15)
%     xlabel('SNR','Interpreter','latex','Fontsize',fS)
    xlim([1 10])
    ylabel('FDR','Interpreter','latex','Fontsize',fS)
    fig = gcf;
    fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 2.5];
    title('b) FDR','Fontsize',fS)
%     set(gca,'Fontsize',15)

    % Plot FRR, pixelbased
    nexttile
    plot(snr,ppFrrPixelbased,'-d','linewidth',2);
    hold on
    plot(snr,otsuFrrPixelbased,'--x','linewidth',2);
    plot(snr,frrEst,'.-','linewidth',2)
    hold off
%     legend("non-Otsu","Otsu's method","estimate",'Fontsize',15)
%     xlabel('SNR','Interpreter','latex','Fontsize',fS)
    xlim([1 10])
    ylabel('FRR','Interpreter','latex','Fontsize',fS)
    fig = gcf;
    fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 2.5];
    title('c) FRR','Fontsize',fS)
%     set(gca,'Fontsize',15)
%    print('C:\Users\Lenovo\postdoc\PAPERS\emccd-paper\draft\Figs\Fig6.eps','-depsc','-r300')
    xlabel('SNR','Interpreter','latex','Fontsize',fS)


    % Plot TMR, pixelbased
%     figure()
    nexttile
    plot(snr,ppTmrPixelbased,'-d','linewidth',2);
    hold on
    plot(snr,otsuTmrPixelbased,'--x','linewidth',2);
    plot(snr,tmrEst,'.-','linewidth',2)
    hold off
    lgnd = legend("EMCCD-PIA","Otsu's method","a priori estimate",'Fontsize',fS)
    lgnd.Layout.Tile = 'north';
    fig = gcf;
    fig.PaperUnits = 'inches';
    xlabel('SNR','Interpreter','latex','Fontsize',fS)
    xlim([1 10])
    ylabel('TMR','Interpreter','latex','Fontsize',fS)

%     fig.PaperPosition = [0 0 4 2.5];
    title('d) TMR','Fontsize',fS)
%     set(gca,'Fontsize',fS)
  
    print(outFig3,'-depsc','-r300')

  
            
            
end
% 
% function probabilistic_binarization()
% 
% end
 