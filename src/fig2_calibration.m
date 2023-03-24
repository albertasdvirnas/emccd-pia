% fig2_calibration(chipPars)

function [lambdaBg,intThreshBg] = fig2_calibration(chipPars,outFig,images,qStar,showfig)

%
% Plot of the estimated lambda_bg using a beads image
%

% close all

% pixel size in exp
pixelsize = 160;

if nargin < 5
    showfig = 0;
end

if nargin < 3 % in case we run calibration on specific test image for fig2

    % Input parameters
    filename = chipPars.inputImage; %'data\100x_gain100_lamp50_018.tif';
    qStar = 0.5; % parameter which controls the false detection rate (FDR)
    
    % Images and associated information                             
    im = imread(filename); 
    images.imAverage = double(im);
    images.registeredIm{1} = double(im);
    images.imageName = 'BeadsOnSurface';
    images.imageNumber = 1;
    images.runNo = 1;
    showfig = 1;

    % estimates for parameters from Mean Variance (MV) callibration
    gain = mean(chipPars.gainQ(3,2)); % based on gain in the image.. 50 100 300..
    adFactor = mean(chipPars.aduQ(2));
    countOffset = mean(chipPars.deltaQ(3,2));
    roNoise = mean(chipPars.sigmaQ(3,2));

%     gain = mean(chipPars.gain{3});
%     adFactor = mean(chipPars.adFactor);
%     countOffset = mean(chipPars.countOffset{3});
%     roNoise = mean(chipPars.roNoise{3});
else
    gain = chipPars.gain;
    adFactor = chipPars.adFactor;
    countOffset = chipPars.countOffset;
    roNoise = chipPars.roNoise;
%     showfig = 0;

end

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

% Maximum-likelihood estimation for lambda_bg estimation. Should we also
% re-fit other params?
disp('Estimating lambda_bg.');
[lambdaBg , intThreshBg] = ...
    estimate_bg_params(images.imAverage(:), gain, adFactor, countOffset, roNoise, qStar);
disp(' ')

if showfig
% Estimate mean and variances for background and for signal
disp('Estimating mean/variances for bg/signal');
[meanBg , varBg, L, U ] = mean_variance_emccd_analytic( lambdaBg, gain, adFactor, countOffset, roNoise );

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

binEdges = ceil(L):floor(U); 
histAll = histcounts(images.imAverage(:),binEdges);    
disp(' ')   
% hold on
% Plot histogram
nexttile     
binPos = binEdges(1:end-1) + diff(binEdges)/2;
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
% N=2^10; % number of integration points for CDF evaluation
% nIntVals =100;  % number of intensity values at which we plot the fit
nPixels = numel(images.imAverage(:));
% Estimate the number of bg pixels
idxThreshTrueBg = length(find(images.imAverage(:) <= intThreshBg));   
%
[~,cdfAtThreshTrueBg] = pdf_cdf_emccd(intThreshBg+1.5,lambdaBg,gain,adFactor,countOffset,roNoise,L,U);  
%
nBg = min(idxThreshTrueBg/cdfAtThreshTrueBg, nPixels); 
% Estimate the expected number of bin counts 
% based on the fit parameters
[~,cdfEmccdNew] = pdf_cdf_emccd(binEdges,lambdaBg,gain,adFactor,countOffset,roNoise,L,U);
binCountsFit = nBg*diff(cdfEmccdNew);
binPos = binEdges(1:end-1) + diff(binEdges)/2;
plot(binPos,binCountsFit,'--','Color','black','LineWidth',2)

% Set figure labels, etc
xlabel('Image counts','Interpreter','latex')
ylabel('Histogram counts','Interpreter','latex')
% set(gca,'Fontsize',15)
% axis([30 80 0 46000])
title('(b)','Interpreter','latex')
% axis equal
pbaspect([1 0.8 0.8])
legendEntry = strcat(['Fit, $\lambda_{bg} =  ' num2str(lambdaBg,2) ', N_{icr}^{bg}=' num2str(intThreshBg) '$']);
lgnd = legend('Image counts, true background','Image counts, not true background',legendEntry,'Interpreter','latex')
lgnd.Layout.Tile = 'south';
% print('C:\Users\Lenovo\postdoc\PAPERS\emccd-paper\draft\Figs\Fig4.eps','-depsc','-r300')
print(outFig,'-depsc','-r300')

end
% false_detection_rate(images,gain, adFactor, countOffset, roNoise);
% print('output\FigS1.eps','-depsc','-r300')

end


function [lambdaBg,intThreshBg ] = ...
                    estimate_bg_params(intensities,gain, adFactor, countOffset, roNoise,qStar)

    %
    % Determines lambda_bg (Poisson parameter for background),
    % and the intensity cut-off below which essentially all pixels
    % are background, using a recursive method.
    %
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
    %
    % Dependencies: emccd_distribution/log_likelihood_trunc_dist.m
    %               emccd_distribution/p_values_emccd.m
    %
    %

    % Hard-coded variables
    % number of points. no need, we integrate up to pi because of
    % discreteness
%     N = 2^8;     % number of integration points when calculating 
                 % PDF through  numerical inverse Fourier-transform
                 % of the characteristic function 
    opt = statset('MaxIter',200,'MaxFunEvals',400,'FunValCheck','on');
                % here one can add input to the 
                % maximum-likelihood estimation (MLE) routine.
                % Type: statset('mlecustom') at the matlab prompt
                % to see what fields are accessible. 
                
    % Extract chip parameters

    
    r = gain/adFactor;
    
    %% Pre-process. If data has two or more wide peaks, algo would fail unless they are well separated.
    
%     intensityVec = intensities(:);        
% 
%     [a,b] = histcounts(intensities);
%     [pk,pos] = findpeaks(imgaussfilt(a,3));
%     if length(pk) > 1
%         intensityVec = intensities(intensities<b(pos(2)));
%     else
%         intensityVec = intensities(:);        
%     end

    % Sort data
    intensityVec = intensities(:);
    sortI = sort(intensities);
    S = numel(sortI);
    
    m = S;
    nOutliers = round(m/2);
    

    

    % Specify likelihood function
    logL = @(data,lambda) log_likelihood_trunc_dist(data , lambda , ...
                       gain, adFactor, countOffset, roNoise); 
                         
     % Recursively reduce the data set until there are no "outliers" left.
    hasOutliers = 1;
%     nOutliers = 0;
    runs = 0;
    diffLambdas = Inf;
    lambdaPrev = 0;
    
    structRes = [];
    while hasOutliers && runs < 10 && diffLambdas > 0.00001
        runs = runs + 1;
        % Remove outliers
        m = S - nOutliers;
        sortTruncI = sortI(1:m);
        % Fit lambda
        lamGuess = abs((sortTruncI(round(m/2)) - countOffset)/r);
               
%    direct check (without mle
%         vv = 10:0.1:60;
%         dat=  arrayfun(@(x) logL(sortTruncI,x),vv);
%         figure,plot(vv,dat)
%         
        try 
            [lambdaBg,pci] = mle(sortTruncI,'logpdf',logL,'start',lamGuess,'lowerBound',0,'Options',opt);
                % in the first iteration (when there are outliers) 
                % the MLE routine may not converge too well and
                % hence will pass on a warning. Here we suppress this warning.
        catch
        end

        diffLambdas = abs((lambdaBg - lambdaPrev)/lambdaPrev);
        lambdaPrev = lambdaBg;
        % 
        % Calculate p-values. Anything above mean+6STD from the dist automatically assumed to be   outliers      
        [pVals, pdfUnique, cdfsUnique, intUnique] = p_values_emccd_sorted(sortI,lambdaBg, gain, adFactor, countOffset, roNoise);        
        pValsFlipped = fliplr(pVals);

% %         
%         [pVals, pdfUnique, cdfsUnique, intUnique] = p_values_emccd_sorted(sortTruncI,20,gain, adFactor, countOffset, roNoise);        
% % 
%         figure,histogram(sortTruncI,'normalization','pdf')
%         hold on
%         plot(intUnique,pdfUnique)
%         xlim([0 300])
%         
%         histogram(sortTruncI,'normalization','cdf','DisplayStyle','stairs')
%         hold on
%         plot(intUnique,cdfsUnique)

        % Find outliers https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-6/The-positive-false-discovery-rate--a-Bayesian-interpretation-and/10.1214/aos/1074290335.full
        threshold = ((1:S)./(S)).*qStar; % m-number of remaining pixels 
%                 threshold = ((1:m)./(S)).*qStar; % m-number of remaining pixels 

        outliers = find(pValsFlipped < threshold,1,'last');
        if ~isempty(outliers)
            nOutliers = outliers(end);
            hasOutliers = 1;
        else
            hasOutliers = 0;
        end
%         outliers
        structRes.nOutliers(runs+1) =  nOutliers;
        structRes.lambdaBg(runs+1) =  lambdaBg;
%                 structRes.lambdaBg(runs+1) =  lambdaBg;

    end
    idxBgEstimation = m;
    intThreshBg = sortI(idxBgEstimation);
    intThreshBg = round(intThreshBg) - 1;  % why -1?
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

end

function [logL] = log_likelihood_trunc_dist(sortTruncI,lambda,...
                             gain, adFactor, countOffset, roNoise)

    % Calculates the  log-likelihood for the truncated EMCCD distribution . 
    % 
    % Input: 
    % 
    % sortTruncI = sorted and truncated intensity values 
    % lambda = Poisson parameter
    % chipPars = struct containing the chip parameters
    % N = number of integration points when calculating PDF 
    %     through numerical inverse Fourier-transform
    %     of the characteristic function 
    %
    % Output:
    % 
    % logL = log likelihood
    %
    % Comment: 
    % The truncated PDF is
    %      PDF_trunc = pdfEmccd(I)/cdfEmccd(I_trunc) for I <= I_trunc
    %      PDF_trunc = 0 elsewhere
    % Here, I_trunc is the truncation intensity.
    %  
    % Dependencies: emccd_distribution/pdf_cdf_emccd.m
    %
    r = gain/adFactor;

        % Analytic expressions for the mean and variance
    EX = lambda*gain/adFactor+countOffset; 
    STD = sqrt(roNoise^2 + 2*lambda*r^2 + 1/12);  
    
    % Analytic expression for the characteristic function 
    % for the EMCCD distribution
%     cfAnaly = @(t) exp(-t.^2*roNoise^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*offset)*2*sin(t/2)/t;

    %
   % limits where pdf is nonzero
    numstds = 6;
    L = EX-numstds*STD; % mean - 6 std
    U = EX+numstds*STD;
    
%     U = min(max(sortTruncI),U); % limit to U for truncated case

    
    % Hard-coded variable
    binWidth = 1;     % if this parameter is = 1 and the intensities are integers, 
                      % then the log-likelihood calculation is exact.
                     
    % Get bin edges
    binEdges = ceil(L):binWidth:floor(U);
%     binEdges = min(sortTruncI)-binWidth/2:binWidth:max(sortTruncI) - binWidth/2;
%     binEdges = [binEdges , max(sortTruncI) + binWidth/2];   
    histAll = histcounts(sortTruncI,binEdges)';
    binPos = binEdges(1:end-1) + diff(binEdges)/2;
   
    [pdfEmccd,cdfEmccd] = pdf_cdf_emccd(binPos,lambda,gain, adFactor, countOffset, roNoise,L,U);

%     [pdfEmccd,cdfEmccd] = pdf_cdf_emccd(binPos,lambda,chipPars,N);
    
    [~ ,cdfEmccdEnd] = pdf_cdf_emccd(min(binEdges(end),max(sortTruncI)+1),lambda,gain, adFactor, countOffset, roNoise,L,U);
%     [~ ,cdfEmccdStart] = pdf_cdf_emccd(30,lambda,gain, adFactor, countOffset, roNoise,L,U);

    % log-likelihood
    logL = sum(histAll.*log(pdfEmccd)) - sum(histAll)*log(cdfEmccdEnd);
    
 
  
end



function [pdfEmccd,cdfEmccd] = pdf_cdf_emccd(intensities,lambda,gain,adFactor,offset,roNoise,L,U)

    % Generates EMCCD probability density function (PDF) and 
    % cumulative distribution functin (CDF) by numerical inversion 
    % of the characteristic function.
    %
    % Input:
    % 
    % intensities = vector (or matrix) with intensity values
    % lambda = Poisson parameter
    % chipPars = struct containing the chip parameters
    % N = number of integration points.
    %
    % Output:
    % 
    % pdfEmccd = EMCCD probability density function 
    % cdfEmccd = cumulative distribution function 
    %
    % Refs: V. Witkovský, "Numerical inversion of a characteristic function: 
    % An alternative tool to form the probability distribution of 
    % output quantity in linear measurement models.", 
    % Acta IMEKO 5.3 (2016): 32-44, see Eqs. (8) and (9).
    %
    %
    
    % Hard-coded variables:
%     pdfMin = 1E-14;   % smallest allow value for PDF 
                      % (need to be > 0 to avoid errors in 
                      % log-likelihood calculations)
%     cdfDelta = 1E-14; % smallest allowed CDF value is cdfDelta,
                      % and largest allowed CDF value is 1-cdfDelta.
      
     % Extract chip parameters
%     gain = chipPars.gain;
%     adFactor = chipPars.adFactor;
%     offset = chipPars.countOffset;
%     roNoise = chipPars.roNoise;
    r = gain/adFactor;
      
%     % Analytic expressions for the mean and variance
    EX = lambda*gain/adFactor+offset; 
%     STD = sqrt(roNoise^2 + 2*lambda*r^2 + 1/12);  
%     
%     % Analytic expression for the characteristic function 
%     % for the EMCCD distribution
% %     cfAnaly = @(t) exp(-t.^2*roNoise^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*offset)*2*sin(t/2)/t;
% 
%     %
%    % limits where pdf is nonzero
%     numstds = 6;
%     L = EX-numstds*STD; % mean - 6 std
%     U = EX+numstds*STD;
%     U = min(max(intensities),U); % limit to U for truncated case


	% optimal value for step parameter
    dt = 2*pi/(U-L);

    % For discrete, integral is -pi..pi, because the output variable is
    % discretized
    N = pi/dt;
    
    
    % Estimate step size, dt, for numerical integration
    t = (1:1:N)' * dt;  

    cf = char_fun(t , roNoise,lambda,r,offset);

        % y is the grid for our pdf (from L to U)
    y = intensities;
    
    
    % calculate main integral
    pdfEmccd = trapezoidal_pdf(y,dt,t,cf);
    cdfEmccd = trapezoidal_cdf(y,dt,t,cf,EX);

   
    
end


function cfCombined = char_fun(t , roNoise,lambda,r,offset)
%

    cfAnaly = exp(-t.^2*roNoise^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*offset);
    cfROUND = 2*sin(t/2)./t;
    cfROUND(t==0) = 1;


    
    cfCombined = cfAnaly.*cfROUND;

end

%
function pdf = trapezoidal_pdf(y,dt,t,cf)
    w = ones(length(t),1);
    w(end) = 1/2; % last coef is 1/2
       
    pdf = dt/pi*(1/2 +cos(t*y)'*(real(cf).*w)+sin(t*y)'*(imag(cf).*w));
end
%
function cdf = trapezoidal_cdf(y,dt,t,cf,ex)
    w = ones(length(t),1);
    w(end)=1/2; % last coef is 1/2
    cdf = 1/2 - dt/pi*(1/2*(ex-y') +cos(t*y)'*(imag(cf./t).*w)-sin(t*y)'*(real(cf./t).*w));
end


function [pVals,pdfUnique,cdfsUnique,intUnique] = p_values_emccd_sorted(sortedInt, lambda, gain, adFactor, countOffset, roNoise)

    %
    % Calculates p-values using the EMCCD distribution
    % The input intensity values are assumed to be sorted.
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
    
    % Find all unique intensity values
   [intUnique,idx]  = unique(sortedInt);
   
   % take only those smaller than integration range
    idx = idx(intUnique <= floor(U));
    intUnique = intUnique(intUnique <= floor(U));
    
    nVals = length(idx);

    % Evaluate the CDF only at the unique intensities 

    [pdfUnique,cdfsUnique] = pdf_cdf_emccd(intUnique',lambda,gain, adFactor, countOffset, roNoise,L,U);
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
    pVals = zeros(1,length(sortedInt));
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