function [] = fig4_segmentation(imageFilenames, chipParsCur,outFig,pValThresh,pValThreshBinarization)

if nargin < 1
    imageFilenames = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif', 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV1_DAPI\20x_gain300_lamp100_001.tif'};
    
end

    ff = figure;
    t=tiledlayout(1,2);
    axtile =[];
    axtile{1}= nexttile(t);
    axtile{2} = nexttile(t);
% Input PARAMETERS 
if nargin < 4
pValThresh = 0.01; % for chi2
pValThreshBinarization = 1E-2; % for binarization
end
T = 64; %  tile size

allowedGapLength = 1;        % allowed gap length. In sorted list of 
                             % region sizes this, this gives us how large gaps, 
                             % starting at region size = 1, we allow for the 
                             % regions to be said to belong to 
                             % the "noise regions" cluster. 
images.runNo = 1;

%
% Script for generating a segmentated image (white regions = signal, 
% and black regions = background) from  an input fluroescence image. 

    % Hard-coded variables
%     fontSizePlot = 12;         % fontsizes in plots
%     lineWidthPlot = 2;         % linewidth in plots
lineWidthBoundaries = 1;   % in images, this gives the 
                               % size of the overlaid boundaries
%     nMedians = 4;              % upper threshold (in terms of number of medians)
                               % when setting the constrast
     

% Actions to be performed on the segmentation output
%     NOTE: No error handling of actions associated 
%     with missing required fields in segOutput.
actions.plotBinarizedImage = 1;      % Plot binarized image
actions.showOriginalRegions = 0;     % Show plots of boundaries of original regions 
                                     % obtained from binarized image
                                     % and of associated label
                                     % matrices
actions.adjustContrast = 1;          % Adjust contrast when displaying the original image 
actions.plotUsingBwboundaries = 0;   % When plotting region boundaries,
                                     % detect the boundaries using the matlab-function
                                     % nbwboundaries function. 
                                     % Else a custom-written function is used. 
actions.showMergedRegions = 1;		 % Show plots of boundaries of final regions 
actions.showGroundTruthMarkers = 0;  % set = 1 if we want to show markers 
                                     % at "ground" truth locations  
actions.plotHistAll = 1;             % plot a histogram over all pixel intensites
actions.plotOverlaidPdf = 1;         % overlay the fitted PDF on top of the 
                                     % histogram for all intensities
actions.plotHistBlackRegions = 0;    % plot a histogram of intensities over black  
                                     % regions in the final segmentation 
                                     
                                     
for i=1:length(imageFilenames)
    imageFilename = imageFilenames{i};

    % Images and associated information                             
    im = imread(imageFilename); % Specify target image
    images.imAverage = double(im);
    images.registeredIm{1} = double(im);
    images.imageNumber = i;

    titles = {'a)','b)'};

    % Perform image segmentration
    t = cputime;  

    chipPars = chipParsCur{i};

    matrixBlocks = reshape(permute(reshape( double(images.imAverage),T,size( images.imAverage,1)/T,T,[]),[1,3,2,4]),T,T,[]);
    
    N = size(matrixBlocks,3);
    scores = nan(1,size(matrixBlocks,3));
    lambdaPars = nan(1,size(matrixBlocks,3));
    intthreshPars = nan(1,size(matrixBlocks,3));
    
    statsAll = cell(1,size(matrixBlocks,3));
    for idx = 1:size(matrixBlocks,3);
        image2.imAverage = matrixBlocks(:,:,idx);
    
        tic
        [lambdaBg, intThreshBg, stats] = emccdpia_estimation(chipPars,[],image2,pValThresh,0);
        scores(idx) = min(stats.chi2Score);
        lambdaPars(idx) = lambdaBg;
        intthreshPars(idx) = intThreshBg;
        statsAll{idx} = stats;
    
        toc
    end
    gain = chipPars.gain;
    adFactor = chipPars.adFactor;
    roNoise = chipPars.roNoise;
    offset = chipPars.countOffset;

    binarizedImageT = zeros(size(matrixBlocks));
    disp('Binarizing image');
    for j=1:size(binarizedImageT,3)
        [binarizedImageT(:,:,j) , intThreshBlackWhite(j) ] = binarize_image_pval_thresh2(...
            matrixBlocks(:,:,j), pValThreshBinarization , lambdaPars(j) , gain, adFactor, offset, roNoise);
    end
    %
%     binarizedImage = reshape(permute(reshape( double(binarizedImageT),T,size( binarizedImageT,1)/T,T,[]),[1,3,2,4]),T,T,[]);
     binarizedImage = imtile(pagetranspose(binarizedImageT),'thumbnailsize',[T T])';
%     [lambdaBg,intThreshBg] = fig2_calibration(chipPars,outFig{i,2},images,qStar,1);
%

        
%      disp('Binarizing image');
%         [binarizedImage , intThreshBlackWhite ] = binarize_image_pval_thresh2(...
%             images.imAverage, pValThreshBinarization ,lambdaBg , gain, adFactor, offset, roNoise);
%%
  
    % Find regions in the binarized image
    disp('Finding regions in the binarized image')

    % w_white
    [bIm, labelIm] = bwboundaries(binarizedImage, 'noholes','CONN',4);
    labelFilled = imfill(labelIm,'holes');

    [regSizeThreshWhite, regSizesImg] = find_reg_thresh(bIm,labelIm,allowedGapLength);
    
    % w_black
    flipIm = 1-binarizedImage;
%     flipIm(:, [1 end]) = 0;
%     flipIm([1 end],:) = 0;

    [bBg, labelBg] = bwboundaries(flipIm, 'noholes','CONN',4);
%     bBg{1} = []; % ignore first component (likely only background)
    [regSizeThreshBg, regSizesBg] = find_reg_thresh(bBg,labelBg,allowedGapLength);

    bgFilled  = imfill(labelBg,'holes');

    regSizesBg(1) = 0;
    disp(['Region size threshold for black regions  = ',num2str(regSizeThreshBg)]);
    disp(['Region size threshold for white regions  = ',num2str(regSizeThreshWhite)]);  
    segOutput.regSizeThreshBlack = regSizeThreshBg;
    segOutput.regSizeThreshWhite = regSizeThreshWhite;
    disp(' ')
    
    % filter out all regions which has size smaller than
    bIm = bIm(regSizesImg>regSizeThreshWhite);
    labeLocs = find(regSizesImg>regSizeThreshWhite);


    % maybe discard first one as artefact
    bBg = bBg(regSizesBg>regSizeThreshBg);
    regSizesLoc = find(regSizesBg>regSizeThreshBg);
    
    % flip large bg region in a signal image 
    
    % Calculate p-values for final regions using summed intensities
    disp('Calculating p-values for the final segmented regions') 
    
    nReg = length(bIm);
    summedInt = zeros(1,nReg);
    cdfSummedInt =  zeros(1,nReg);

    img = zeros(size(images.imAverage));
    for regIdx=1:nReg
        img(labelFilled==labeLocs(regIdx)) = labeLocs(regIdx);
    end

    for regIdx=1:length(regSizesLoc)
        img(bgFilled==regSizesLoc(regIdx))= 0;
    end
% end

    for regIdx=1:nReg
        pixelvals = images.imAverage(img==labeLocs(regIdx));
        summedInt(regIdx) = sum(pixelvals);
%         summedInt(regIdx) = sum(arrayfun(@(x) images.imAverage(bIm{regIdx}(x,1),bIm{regIdx}(x,2)),1:length(bIm{regIdx})));
%         M = regSizesImg(regIdx);
        M = length(pixelvals);
        
       [~,cdfTemp,~,~] = pdf_cdf_from_characteristic_fun_sum( summedInt(regIdx), lambdaBg,gain,adFactor, offset,roNoise,M); 
%                 [pdfEmccd,cdfTemp,L,U] = pdf_cdf_from_characteristic_fun( summedInt(regIdx), M*lambdaBg,gain,adFactor, M*offset,sqrt(M)*roNoise); 

        cdfSummedInt(regIdx) = min(1,cdfTemp);
    end
    pValsMerged = 1 - cdfSummedInt;
    segOutput.pValsMerged = pValsMerged;
    
    colorValU = 1 - segOutput.pValsMerged;  
    imageInput = images.imAverage;


    %% Plot original image with adjusted contrast
%     figure
% 
%     tiledlayout(1,2,'TileSpacing','tight')
%     ax1=nexttile;
 %     imshow(J,'InitialMagnification','fit');
    axes(axtile{i})
   plot_res(images,img,colorValU,lineWidthBoundaries)

%     figure
%    plot_out(images,bIm,nReg,colorValU,lineWidthBoundaries)
   title(titles{i},'Interpreter','latex')
        
end
    
    print(ff,outFig{1,1},'-depsc','-r300')

end



function plot_res(images,img,colorValU,lineWidthBoundaries)
   sampIm = mat2gray(images.imAverage);
    minInt = min(sampIm(:));
    medInt = median(sampIm(:));
    maxInt = max(sampIm(:));
    J = imadjust(sampIm,[minInt min(1,4*medInt)]);
%     nexttile
%         clf
% %         ax1 = axes;
% figure
        imshowpair(J,img>0,'ColorChannels','green-magenta');
%         parMap = parula;
% %         ax1.Visible = 'off';
% %         ax1.XTick = [];
% %         ax1.YTick = [];
%         hold on    
%         imshow(img==1)
% 
%         % Plot boundaries contours
%         parMap = parula;
%          for regIdx = 1:nReg  % loop over regions
%              boundariesReg = boundariesCellArray{regIdx};
%              colorVal = 1 + floor(colorValU(regIdx)*255); % in the range [1,256]
%              lineColor = parMap(colorVal,:);    
%                  plot(boundariesReg(:,2),boundariesReg(:,1),'LineWidth',lineWidthBoundaries,'Color',lineColor)
%         end   
    


end





function plot_out(images,boundariesCellArray,nReg,colorValU,lineWidthBoundaries)
   sampIm = mat2gray(images.imAverage);
    minInt = min(sampIm(:));
    medInt = median(sampIm(:));
    maxInt = max(sampIm(:));
    J = imadjust(sampIm,[minInt min(1,4*medInt)]);
%     nexttile
%         clf
%         ax1 = axes;
        imshow(J,'InitialMagnification','fit');
        parMap = parula;
%         ax1.Visible = 'off';
%         ax1.XTick = [];
%         ax1.YTick = [];
        hold on    

        % Plot boundaries contours
        parMap = parula;
         for regIdx = 1:nReg  % loop over regions
             boundariesReg = boundariesCellArray{regIdx};
             colorVal = 1 + floor(colorValU(regIdx)*255); % in the range [1,256]
             lineColor = parMap(colorVal,:);    
                 plot(boundariesReg(:,2),boundariesReg(:,1),'LineWidth',lineWidthBoundaries,'Color',lineColor)
        end   
    


end


