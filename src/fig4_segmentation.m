function [] = fig4_segmentation(imageFilenames, chipParsCur,outFig)

if nargin < 1
    imageFilenames = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif', 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV1_DAPI\20x_gain300_lamp100_001.tif'};
    
end

    ff = figure;
    t=tiledlayout(1,2);
    axtile =[];
    axtile{1}= nexttile(t);
    axtile{2} = nexttile(t);
% Input PARAMETERS 
binWidthPerSigmaBg = 0.3;    % bin width / sigma for background 
qStar = 0.5;                 % parameter which controls the 
                             % false detection rate (FDR)
pValThreshBinarization = 1E-2;
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
%     images.imageName = 'BeadsHighConc';
    images.imageNumber = i;



%% WHYY?
% Crop out a movie in the center of the image
% cropLength = 150
% imageTemp = images.imAverage ;
% [nRows , nCols ] = size(imageTemp);
% startRow = round(nRows/2 - cropLength/2) + 0;
% endRow = startRow + cropLength - 1;
% startCol = round(nCols/2 - cropLength/2) + 0 ;
% endCol = startCol + cropLength - 1;
% images.imAverage = imageTemp(startRow:endRow,startCol:endCol);



% Add path to relevant folders

% addpath('../../photophys_image_analysis/util/binarization_segmentation')
% addpath('../../photophys_image_analysis/util/emccd_distribution')
% addpath('../../photophys_image_analysis/util/read_write')



% Perform image segmentration
t = cputime;  

    chipPars = chipParsCur{i};

    [lambdaBg,intThreshBg] = fig2_calibration(chipPars,outFig{i,2},images,qStar,1);
%

    gain = chipPars.gain;
    adFactor = chipPars.adFactor;
    roNoise = chipPars.roNoise;
    offset = chipPars.countOffset;

        
        
     disp('Binarizing image');
        [binarizedImage , intThreshBlackWhite ] = binarize_image_pval_thresh2(...
            images.imAverage, pValThreshBinarization ,lambdaBg , gain, adFactor, offset, roNoise);
%%
  
    % Find regions in the binarized image
    disp('Finding regions in the binarized image')
    [labelIm, isRegBlack ] = ...
             find_regions_in_binarized_image( binarizedImage );
    segOutput.labelIm = labelIm;
    segOutput.isRegBlack = isRegBlack;
    disp(' ');       
          
   
    % Determine region size thresholds
     disp('Determining region size thresholds')
    regSizes = calc_region_sizes(labelIm);    
    [regSizeThreshBlack , regSizeThreshWhite ] =  ...
        calc_region_size_threshold_using_clustering( regSizes , isRegBlack , allowedGapLength);
     disp(['Region size threshold for black regions  = ',num2str(regSizeThreshBlack)]);
    disp(['Region size threshold for white regions  = ',num2str(regSizeThreshWhite)]);  
    segOutput.regSizeThreshBlack = regSizeThreshBlack;
    segOutput.regSizeThreshWhite = regSizeThreshWhite;
    disp(' ')
    
         
    % Generate the final segmentation by removing small black and white
    % regions
    disp('Generating the final segmented image')
    [mergedLabelIm , ~ ] = ...
    generate_final_segmented_regions(labelIm, isRegBlack, regSizeThreshBlack , regSizeThreshWhite );
    fprintf('Found %i merged regions .\n',max(mergedLabelIm(:)))
    segOutput.mergedLabelIm = mergedLabelIm;
    disp(' ')
    
%     [L,U,EX,STD] = calc_bounds(lambdaBg,gain,adFactor,offset,roNoise);
%                 [~,cdfTemp] = pdf_cdf_from_characteristic_fun(L:U,lambdaBg,gain,adFactor, regSizes(regIdx)*offset,sqrt(regSizes(regIdx))*roNoise);

%     % that give intensities to calculate over
%     intensities = ceil(L):floor(U);
%     [pdf,cdf] = pdf_cdf_from_characteristic_fun(intensities,lambdaBg,gain,adFactor,offset,roNoise);
% 
%     % lets say we have 10 images
%     M = 10;
%         [L,U,EX,STD] = calc_bounds(lambdaBg*M,gain,adFactor,M*offset,sqrt(M)*roNoise);
% 
%     intensities = ceil(L):floor(U); % to be more accurate, 
%     
%     [pdf,cdf] = pdf_cdf_from_characteristic_fun(intensities,M*lambdaBg,gain,adFactor,M*offset,sqrt(M)*roNoise);

    % cdf, in this case do not need to be truncated
    
    
    % Calculate p-values for final regions using summed intensities
    disp('Calculating p-values for the final segmented regions')     
    nReg = max(mergedLabelIm(:));  % number of regions
    regSizes = calc_region_sizes( mergedLabelIm  ); % number pixels in each region
    summedInt = calc_region_summed_intensities( images.imAverage , mergedLabelIm  );
                        % summed intensities for each region.
    cdfSummedInt = zeros(1,nReg); 
    for regIdx = 1:nReg
        % redo using new pdf cdf emcc
        M = regSizes(regIdx);
        [L,U,EX,STD] = calc_bounds(lambdaBg*M,gain,adFactor, M*offset,sqrt(M)*roNoise);
%         intensities = ceil(L):floor(U); % to be more accurate, 

        [~,cdfTemp] = pdf_cdf_from_characteristic_fun( summedInt(regIdx), M*lambdaBg,gain,adFactor, M*offset,sqrt(M)*roNoise);
%             [~,cdfTemp] = pdf_cdf_from_characteristic_fun( summedInt(regIdx), regSizes(regIdx)*lambdaBg,gain,adFactor, regSizes(regIdx)*offset,sqrt(regSizes(regIdx))*roNoise);

%          [ ~ , cdfTemp ] = pdf_cdf_emccd_summed_intensities( ...
%            summedInt(regIdx), regSizes(regIdx), lambdaBg , chipPars , N );   
        cdfSummedInt(regIdx) = min(1,cdfTemp);

    end
    pValsMerged = 1 - cdfSummedInt;
    segOutput.pValsMerged= pValsMerged;
           
        
        
%         
% segOutput = segmentation_photophys(images,experiment, ...
%     binWidthPerSigmaBg , qStar , pValThreshBinarization , allowedGapLength);
% e = cputime - t;


 imageInput = images.imAverage;
 % Choose colors [in the range 0 to 1] for the boundaries
 % based on region probbilities of p-values
 nReg = max(mergedLabelIm(:));  % number of regions 
 if isfield(segOutput,'mergedProbSignal')
     colorValU = segOutput.mergedProbSignal;
 elseif isfield(segOutput,'pValsMerged')
     colorValU = 1 - segOutput.pValsMerged;       
 else
     colorValU = 0.9*ones(1,nReg);
 end

% Trace out the boundaries of all regions
% if actions.plotUsingBwboundaries 
%    boundariesCellArray = trace_boundaries_bwboundaries( mergedLabelIm  );
% else
   boundariesCellArray = trace_boundaries( mergedLabelIm );
% end

    %% Plot original image with adjusted contrast
%     figure
% 
%     tiledlayout(1,2,'TileSpacing','tight')
%     ax1=nexttile;
 %     imshow(J,'InitialMagnification','fit');
    axes(axtile{i})

   plot_out(images,boundariesCellArray,nReg,colorValU,lineWidthBoundaries)
%         
%         clabeltext = text(475,320,'$1- p_-value$','Color','black','Interpreter','latex');
%         set(clabeltext,'Rotation',-90) 
       
        
end
    
    print(ff,outFig{i,1},'-depsc','-r300')


% Plot output
% plot_segmentation_photophys_result( segOutput, images , experiment , actions)



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
             for boundaryIdx = 1:length(boundariesReg)
                 contour = boundariesReg{boundaryIdx};
                 plot(contour(:,2),contour(:,1),'LineWidth',lineWidthBoundaries,'Color',lineColor)
             end 

        end   
    
%         % Plot fround truth positions (if available)
%         if (~isempty(groundTruthPositions) & actions.showGroundTruthMarkers)
%             scatter(groundTruthPositions(:,2),groundTruthPositions(:,1),100,'x','MarkerEdgeColor','white')
%         end 
%         hold off
    
        % Set title
        %if isfield(images,'imageName')
        %    strTitle = strcat('Photophysical segmentation, image = ',images.imageName);  
        %    title(strTitle)
        %end
        
%         % Some axis adjustments etc
%         ax2 = axes;
%         linkaxes([ax1,ax2])
%         ax2.Visible = 'off';
%         ax2.XTick = [];
%         ax2.YTick = [];
%         cb = colorbar('Ticks',[0,.2,0.4,0.6,0.8,1],...
%              'TickLabels',{0,0.2,0.4,0.6,0.8,1},'Position',[.85 .11 .055 .815]); %.815
%         set(gca,'ColorScale','log')
%         cb.Label.String ='1- p';

end
