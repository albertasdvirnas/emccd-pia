% fig2_calibration(chipPars)

function [lambdaBg,intThreshBg,structRes] = emccdpia_estimation(chipPars, outFig, images, pValThresh, showfig)
%   emccdpia_estimation - estimation of lambda_bg parameter value and intThreshBg
%   threshold
%   Args:
%       chipPars - chip parameters
%       outFig - output figure to save to
%       images - input image
%       alphaStar - level to control FDR/FOR at
%       showfig - whether to show fig
%
%   Returns:
%       lambdaBg - lambda background,intThreshBg - intensity
%       threshold,structRes - extra settings
%
%
% Plot of the estimated lambda_bg using a beads image
%

% close all

% pixel size in exp
if isfield(chipPars,'pixelsize')
    pixelsize = chipPars.pixelsize;
else
    pixelsize = 160;
end

if nargin < 5
    showfig = 0;
end

% nicer to work with passing individual parameters
gain = chipPars.gain;
adFactor = chipPars.adFactor;
countOffset = chipPars.countOffset;
roNoise = chipPars.roNoise;
% method = chipPars.method;

% Maximum-likelihood estimation for lambda_bg estimation. Should we also
% re-fit other params?
import Core.estimate_lambda;
import Core.pdf_cdf_emccd;

disp('Estimating lambda_bg.');
[lambdaBg ,intThreshBg,structRes] = ...
    estimate_lambda(images.imAverage(:), gain, adFactor, countOffset, roNoise,pValThresh);
disp(' ')

% Estimate statistics at p-value
        % Binarize image
import Core.estimate_statistics;
[stats,binarizedImage] = estimate_statistics(images.imAverage,structRes,pValThresh);
structRes.stats = stats;
structRes.binarizedImage = binarizedImage;
% Show image with contrast set to show noise
if showfig
    figure
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact')
    nexttile

    sampIm = mat2gray(images.imAverage);
    minInt = min(sampIm(:));
    medInt = median(sampIm(:));
    maxInt = max(sampIm(:));
    J = imadjust(sampIm,[minInt min(1,4*medInt)]);
    imshow(J,'InitialMagnification','fit');
    %imshow(images.imAverage/max(images.imAverage(:)))
    hold on    
    % scale bar (ten microns in number of pixels)
    nPixels = 1e4/pixelsize;
    x = [5, 5 + nPixels ];
    y = [0.9*size(sampIm,1) , 0.9*size(sampIm,1)];
    plot(x,y,'Linewidth',8,'Color',[1 1 1])
    text(0,0.05,'10 microns','Fontsize',10,'Color',[1 1 1],'Units','normalized')
    title('(a)','Interpreter','latex')
%     axis equal
%     pbaspect([1 0.8 0.8])

end

if showfig
% Estimate mean and variances for background and for signal
disp('Estimating mean/variances for bg/signal');
[meanBg , varBg, ~, ~ ] = mean_variance_emccd_analytic( lambdaBg, gain, adFactor, countOffset, roNoise );
L = structRes.LU(1);
U = structRes.LU(2);

disp(['Estimated mean of background  = ',num2str(meanBg)]);
disp(['Estimated std dev for background = ',num2str(sqrt(varBg))]); 
% [meanSignal , varSignal ] = estimate_mean_variance_for_signal( images.imAverage(:) ,...
% intThreshBg,  lambdaBg, chipPars );
% disp(['Estimated mean of signal = ',num2str(meanSignal)]);
% disp(['Estimated std dev for signal = ',num2str(sqrt(varSignal))]); 
% disp(' ');
    
% Generate a histogram over all intensities and determine bin edges
disp('Generating histogram over pixel intensities')
% binEdges  =  generate_binedges(images.imAverage, ...
%                   binWidthPerSigmaBg, lambdaBg,gain, adFactor, countOffset, roNoise );

% binEdges = [max(1,ceil(L)):floor(U)]-0.5; 
% histAll = histcounts(images.imAverage(:),binEdges);   
histAll = structRes.histAll(1:end-1); % last idx is all remaining bins
disp(' ')   
% hold on
% Plot histogram
nexttile     
binPos = 1:structRes.LU(2) + 0.5;
[minVal , idx] = min(abs(binPos - intThreshBg));
h1 = bar(binPos(1:idx),histAll(1:idx),1); 
set(h1,'FaceColor',[0.4 0.6 0.9])
set(h1,'EdgeColor','black')
hold on
h2 = bar(binPos(idx+1:end),histAll(idx+1:end),1); 
set(h2,'FaceColor',[1 0.5 0.3])
set(h2,'EdgeColor','black')
hold on


% Overlay histogram with the fitted PDF 

% Number of background pixels. For synthetic images we know true
% import Core.p_values_emccd_sorted;
% cdfFun =  @(x)  1-p_values_emccd_sorted(x,lambdaBg,gain,adFactor,countOffset,roNoise);
% uniqueVals = unique((images.imAverage(:)));
% stopIntensity = quantile(images.imAverage(:),0.25);
% cdfIntensities = uniqueVals(1):stopIntensity;
% stats.yval =    arrayfun(@(x) sum(images.imAverage(:)<=x), cdfIntensities);
% stats.xval =   structRes.cdf(cdfIntensities);
% nBg =  stats.xval'\stats.yval'; % nBg estimated as the slope of this
% nBg

% nPixels = numel(images.imAverage(:));
% Estimate the number of bg pixels
% idxThreshTrueBg = length(find(images.imAverage(:) <= min(U,intThreshBg)));   
% %
% [~,cdfAtThreshTrueBg] = pdf_cdf_emccd(min(U,intThreshBg)+0.5,lambdaBg,gain,adFactor,countOffset,roNoise,L,U);  
% %
% nBg = min(idxThreshTrueBg/cdfAtThreshTrueBg, nPixels); 
% Estimate the expected number of bin counts 
% based on the fit parameters
% [pdfEmccdNew, cdfEmccdNew] = pdf_cdf_emccd(binPos,lambdaBg,gain,adFactor,countOffset,roNoise,L,U);
% pdfEmccdNew = pdfEmccdNew;%/cdfEmccdNew(binPos==structRes.Nthresh);
binCountsFit = stats.nBg.*structRes.pdf;
% binPos = binEdges(1:end-1) + diff(binEdges)/2;
plot(binCountsFit,'--','Color','black','LineWidth',2)

% Set figure labels, etc
xlabel('Image counts','Interpreter','latex')
ylabel('Histogram counts','Interpreter','latex')
% set(gca,'Fontsize',15)
% axis([30 80 0 46000])
title('(b)','Interpreter','latex')
% axis equal
pbaspect([1 0.8 0.8])
legendEntry = strcat(['Fit, $\lambda_{bg} =  ' num2str(lambdaBg,2) ', N_{icr}^{bg}=' num2str(intThreshBg) '$']);
lgnd = legend('Image counts, true background','Image counts, not true background',legendEntry,'Interpreter','latex');
lgnd.Layout.Tile = 'south';
% print('C:\Users\Lenovo\postdoc\PAPERS\emccd-paper\draft\Figs\Fig4.eps','-depsc','-r300')
print(outFig,'-depsc','-r300');


%  % Add error bars
% X = 3; % how many stds away from mean we accept
% pBin = histAll/sum(histAll)
% varBin = sum(histAll)*pBin.*(1-pBin)
% errHigh = X*sqrt(varBin)
% errLow =  X*sqrt(varBin)
% errorbar(binPos,histAll,errLow,errHigh)
end
% false_detection_rate(images,gain, adFactor, countOffset, roNoise);
% print('output\FigS1.eps','-depsc','-r300')

end





% function [L,U] = dist_bounds(sortedInt, lambda, gain, adFactor, countOffset, roNoise)
% 
%     % bounds
% 
%     r = gain/adFactor;
% 
%     % Analytic expressions for the mean and variance
%     EX = lambda*gain/adFactor+countOffset; 
%     STD = sqrt(roNoise^2 + 2*lambda*r^2 + 1/12);  
%     % limits where pdf is nonzero
%     numstds = 6;
%     L = EX-numstds*STD; % mean - 6 std
%     U = EX+numstds*STD;
% end

function [pVals,pdfUnique,cdfsUnique,intUnique] = p_values_emccd_continuous(sortedInt, lambda, gain, adFactor, countOffset, roNoise)

    % Possible modification of p_values_emccd_sorted:
    % Calculates p-values using the EMCCD distribution
    % The input intensity values are unrounded to continuous by adding
    % uniform random number
    %
    % Input:
    %
    % intensities = vector with sorted intensity values
    % chipPars = struct containing the chip parameters
    %
    % Output:
    %
    % pVals = p-values for each input intensity
    %
    % Dependencies: emccd_distribution/pdf_cdf_emccd.m 
    % 
    
    % Hard-coded variable
%     N = 2^8;     % number of integration points when calculating 
                 % PDF through  numerical inverse Fourier-transform
                 % of the characteristic function   

     % can pass it with a top function
    r = gain/adFactor;
    % Analytic expressions for the mean and variance
    EX = lambda*gain/adFactor+countOffset; 
    STD = sqrt(roNoise^2 + 2*lambda*r^2 + 1/12);  
    % limits where pdf is nonzero
    numstds = 6;
    L = EX-numstds*STD; % mean - 6 std
    U = EX+numstds*STD;

   % Turn into a row vector
    sortedInt = sortedInt(:);  
    nAll = length(sortedInt);

   [intUnique,idx]  = unique(sortedInt);

    sortedInt = sortedInt + rand(length(sortedInt),1)-0.5;
    sortedInt = sort(sortedInt);
    idx = 1:length(sortedInt);
    % Find all unique intensity values
%     intUnique  unique(sortedInt)
   
   % take only those smaller than integration range
    idx = idx(sortedInt <= floor(U));
    sortedInt = sortedInt(sortedInt <= floor(U));
    
    nVals = length(idx);

    % Evaluate the CDF only at the unique intensities 

    [pdfUnique,cdfsUnique] = pdf_cdf_emccd(sortedInt',lambda,gain, adFactor, countOffset, roNoise,L,U);
%     [~, cdfsEnd] = pdf_cdf_emccd(min(U,max(intUnique)+1),lambda,gain, adFactor, countOffset, roNoise,L,U);
%     cdfsUnique = cdfsUnique./cdfsEnd;
%     pdfUnique = pdfUnique./cdfsEnd;

    % remove ones out of range
    cdfsUnique = min(cdfsUnique,1);
    cdfsUnique = max(cdfsUnique,0);
% 
%     cdfsUnique(intUnique>floor(U)) = nan; % should not be the case
%     pdfUnique(intUnique>floor(U)) = nan;

%     [~ , cdfsUnique] = pdf_cdf_emccd(intUnique,lambda,chipPars,N);

    % Calculate p-values for all input intensities
    pVals = zeros(1,nAll);
    for k=1:nVals-1   
        pVals(idx(k):idx(k+1)-1) = 1 - cdfsUnique(k);
    end
    pVals(idx(nVals):end) = 1 - cdfsUnique(end);


end


function  [meanAnalytic , varAnalytic,L,U ] = mean_variance_emccd_analytic( lambda,gain, adFactor, countOffset, roNoise )

    %
    % Gives the analytic mean and variance of the EMCCD
    % probability density function   
   r = gain/adFactor;

    % Analytic expressions for the mean and variance
    meanAnalytic = lambda*gain/adFactor+countOffset; 
    varAnalytic = roNoise^2 + 2*lambda*r^2 + 1/12;  
    
   
    %
   % limits where pdf is nonzero
    numstds = 6;
    L = meanAnalytic-numstds*sqrt(varAnalytic); % mean - 6 std
    U = meanAnalytic+numstds*sqrt(varAnalytic);
    

    
end
 

function  [meanSignal , varSignal ] = estimate_mean_variance_for_signal( imageInput ,...
       intThreshBg, lambdaBg, gain, adFactor, countOffset, roNoise)

    %
    % Estimates the mean and variance of the signal regions. 
    % A "subtraction" type method is applied which uses the 
    % the sample moments for all pixels and the analytic estimates
    % of the moments of the background  distributions.
    %
    % Input:
    %
    % imageInput = image (matrix with intensity values)
    % intThreshBg = intensity threshold below which there are only
    %              "true" background intensities
    % lambdaBg = lambda parameter for background Poisson distribution  
    % chipPars = struct containing the chip parameters
    %  
    % Output:
    %
    % meanSignal = mean intensity for signal regions
    % varSignal = variance in intensity for signal regions
    %
    % Dependencies: emccd_distribution/mean_variance_emccd_analytic.m 
    %               binarization_segmentation/estimate_n_bg.m
    %               
    %

    % Hard-coded variable
    N = 2^8;   % number of integration points used when
                % calculating the EMCCD CDF
                
    % All intensities           
    intensityVec = imageInput(:);
     
    % Estimate analytic mean and variance for background
    [meanBgAnalytic , varBgAnalytic ] = mean_variance_emccd_analytic( lambdaBg, gain, adFactor, countOffset, roNoise  );

    % Estimate the fraction, f_bg, of background pixels
    nPixels = numel(intensityVec);
    nBg = estimate_n_bg( intensityVec, intThreshBg,  nPixels, lambdaBg, chipPars );
    fBg = nBg/(nPixels+1);
                           
    % Estimation the mean intensity for the signal pixels
    meanAll = sum(intensityVec)/nPixels;
    meanSignal = ( meanAll - meanBgAnalytic*fBg)/(1-fBg);
    
    % Estimation the variance in intensity for the signal pixels
    sqIntBg = varBgAnalytic + meanBgAnalytic^2;    % <I^2> for background.
    sqIntAll = sum(intensityVec.^2)/nPixels;       % <I^2> for all pixels
    sqIntSignal = (sqIntAll - sqIntBg*fBg)/(1-fBg);
    varSignal = sqIntSignal - meanSignal^2;
   
     
end

% 
% 
% function  nBg = estimate_n_bg( intensities, intThreshBg,  ...
%                                nPixels, lambdaBg, gain, adFactor, countOffset, roNoise  )
% 
%     %
%     % Estimates the number of background pixels in the image
%     %
%     % Input:
%     %
%     % intensities = matrix or vector with intensity values 
%     % intThreshBg = intensity  below which there are only
%     %                   "true" background intensities
%     % nPixels = number of pixels in the image
%     % lambdaBg = lambda parameter for background Poisson distribution  
%     % chipPars = struct containing the chip parameters
%     %  
%     % Output:
%     %
%     % meanSignal = mean intensity for signal regions
%     % varSignal = variance in intensity for signal regions
%     %
%     % Dependencies: emccd_distribution/pdf_cdf_emccd.m 
%     %
% 
%     % Hard-coded variable
% %      N = 2^8;   % number of integration points used when
% %                 % calculating the EMCCD CDF
%                 
%    
%      % Total number of background pixels
%     nTrueBg = length(find(intensities(:) <= intThreshBg));
%     [~,cdfAtThreshTrueBg] = pdf_cdf_emccd(intThreshBg+0.5,lambdaBg,chipPars,N);
%                                % bin edges are located at "half integers"
%     nBg = min(nTrueBg/cdfAtThreshTrueBg, nPixels);    
%     
% end
%  
%  

function  binEdges  = generate_binedges(intensities, ...
                     binWidthPerSigmaBg, lambdaBg, gain, adFactor, countOffset, roNoise )

    %
    % A histogram for all intensities is generated using the 
    % expected variance of the background to determine the bin width. 
    %
    % Input:
    %
    % intensities = matrix or vector with intensity values 
    % binWidthPerSigmaBg = ratio of bin width and the standard deviation
    %                      for the background intensities
    % lambdaBg = lambda-parameter for background Poisson distribution  
    % chipPars = struct containing the chip parameters
    %
    % Output:
    %
    % binEdges = edge positions of the bins 
    %           [there is one more bin edge than there are bins]
    %
    % Comment: to generate a histgram use: 
    %          h = histcounts(intensities(:),binEdges)
    %
    % Dependencies: emccd_distribution/mean_variance_emccd_analytic.m
    
    %

    % Calculate bin width
    [ ~ , varBg ] = mean_variance_emccd_analytic( lambdaBg, gain, adFactor, countOffset, roNoise );
    binWidth = binWidthPerSigmaBg*sqrt(varBg);

    %  Make a histogram, H_all, over all pixel intensities in the image
    intensityVec = intensities(:);
    binEdges = min(intensityVec)-binWidth:binWidth:max(intensityVec) + binWidth;
   
  
end

% maybe move to separate function as fig1sup
function []=false_detection_rate(images,gain, adFactor, countOffset, roNoise)
% Plot of the estimated lambda_bg using a beads image
%
% close all

% Input parameters
% filename = 'beads_low_conc_100x_gain100_lamp100_013.tif';
qStarStart = 0.01; % start-value, parameter which controls the false detection rate (FDR)
qStarStep = 0.01; % step-value
qStarStop = 0.95; % stop-value
    

% Loop over different values for qStar
qStarVec = qStarStart:qStarStep:qStarStop;
lambdaBgVec = zeros(1,length(qStarVec));
intThreshBgVec = zeros(1,length(qStarVec));
disp('Estimating lambda_bg.');
counter = 1;
%
for qStar = qStarStart:qStarStep:qStarStop
    disp(['qStar = ',num2str(qStar)])
    % Maximum-likelihood estimation for background strength 
    [lambdaBg , intThreshBg] = ...
        estimate_bg_params(images.imAverage(:),gain, adFactor, countOffset, roNoise,qStar);
    lambdaBgVec(counter) = lambdaBg;
    intThreshBgVec(counter) = intThreshBg;
    counter = counter + 1;
end
disp(' ')


% Plot
figure

yyaxis left
colorLeft = [0 0.4470 0.7410];
plot(qStarVec,lambdaBgVec,'-','linewidth',2,'Color',colorLeft)
xlabel('q^*','FontSize',12)
ylabel('\lambda_{bg}','FontSize',12)
set(gca,'YColor',colorLeft,'Fontsize',12)
axis([0 1 22 24])

yyaxis right
colorRight = [0.8500 0.3250 0.0980];
plot(qStarVec,intThreshBgVec,'-','linewidth',2,'Color',colorRight)
ylabel('Image count threshold','Color',colorRight)
set(gca,'YColor',colorRight,'Fontsize',12)

legend('estimated \lambda_{bg}','estimated image count threshold','location','southoutside')

% save (for replotting)


end

% % Score should be improving with less data to fit?
% %
% % vals = 0:10000: numel(sortI);
% % lam = [];
% % score = [];
% % for i=1:length(vals)
% %     nOutliers = vals(i);
% % % test calib
% %        % Remove outliers
% %            m = numel(sortI);
% % 
% % m = m - nOutliers;
% % sortTruncI = sortI(1:m);
% % % Fit lambda
% % lamGuess = abs((sortTruncI(round(m/2)) - countOffset)/r);
% % [lambdaBg,pci] = mle(sortTruncI,'logpdf',logL,'start',lamGuess,'lowerBound',0,'Options',opt);
% %       
% % 
% % % 
% % % Calculate p-values. Anything above mean+6STD from the dist automatically assumed to be   outliers      
% % [pVals, pdfUnique, cdfsUnique, intUnique] = p_values_emccd_sorted(sortTruncI,lambdaBg,gain, adFactor, countOffset, roNoise);        
% % pValsFlipped = fliplr(pVals);
% % 
% % % lim = 40;
% % [counts,edges] = histcounts(sortTruncI,[intUnique-0.5;max(intUnique)+0.5],'Normalization','pdf');
% % pos = (edges(2:end)+edges(1:end-1))/2;
% % 
% % % figure,plot(pos,counts);
% % % hold on
% % % plot(intUnique,pdfUnique)
% % counts = counts(~isnan(pdfUnique));
% % cdfsUnique = cdfsUnique(~isnan(pdfUnique));
% % pdfUnique = pdfUnique(~isnan(pdfUnique));
% % score(i) = 1/length(counts)*sqrt(sum((counts-pdfUnique'/cdfsUnique(end)).^2));
% % % 
% % % figure,histogram(sortTruncI,'normalization','pdf')
% % % hold on
% % % plot(intUnique,pdfUnique/cdfsUnique(end))
% % xlim([0 100])
% % lam(i) = lambdaBg; 
% % end
% %     
% % %    
% %  