
calibFold = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\';
imagefiles = {'Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp50_018.tif',...
    '2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif'};


i=1;
    imageFilename = fullfile(calibFold,imagefiles{i});

    % Images and associated information                             
    im = imread(imageFilename); % Specify target image
    im = double(im);

% Otsu
imageMax = max(im(:));
normImage = im/imageMax;
normThreshold = graythresh(normImage);
intThresh= normThreshold*imageMax;
threshIm = im>intThresh;

% adaptive thresh
T = adaptthresh(normImage);
BW = imbinarize(normImage,T);


%% Otsu
lineWidthBoundaries = 1;   % in images, this gives the 

sampIm = mat2gray(im);
minInt = min(sampIm(:));
medInt = median(sampIm(:));
maxInt = max(sampIm(:));
J = imadjust(sampIm,[minInt min(1,4*medInt)]);

[boundariesCellArray, labelIm] = bwboundaries(threshIm);

nReg = length(boundariesCellArray);

ff=figure;tiledlayout(2,2)
nexttile([2 1])

imshow(J,'InitialMagnification','fit');
hold on    

% Plot boundaries contours
parMap = parula;
 for regIdx = 1:nReg  % loop over regions
     boundariesReg = boundariesCellArray{regIdx}; 
     plot(boundariesReg(:,2),boundariesReg(:,1),'LineWidth',lineWidthBoundaries,'Color','blue')
 end   

pixelsize = 160;
nPixels = 1e4/pixelsize;
x = [5, 5 + nPixels ];
y = [0.9*size(im,1) , 0.9*size( im,1)];
plot(x,y,'Linewidth',4,'Color','white')
text(0,0.05,'10 microns','Fontsize',10,'Color','white','Units','normalized')
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'color','white');
    



% imshowpair(J,threshIm,'ColorChannels','green-magenta');

title('a) Otsu segmentation','Interpreter','latex')

nexttile([2 1])
[boundariesCellArray, labelIm] = bwboundaries(BW);


boundariesCellArray(cellfun(@(x) size(x,1)<4,boundariesCellArray)) = [];

nReg = length(boundariesCellArray);

imshow(J,'InitialMagnification','fit');
hold on    

% Plot boundaries contours
parMap = parula;
 for regIdx = 1:nReg  % loop over regions
     boundariesReg = boundariesCellArray{regIdx}; 
     plot(boundariesReg(:,2),boundariesReg(:,1),'LineWidth',lineWidthBoundaries,'Color','blue')
 end   

pixelsize = 160;
nPixels = 1e4/pixelsize;
x = [5, 5 + nPixels ];
y = [0.9*size(im,1) , 0.9*size( im,1)];
plot(x,y,'Linewidth',4,'Color','white')
text(0,0.05,'10 microns','Fontsize',10,'Color','white','Units','normalized')
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'color','white');

   title('b) Adaptive thresholding','Interpreter','latex')


    print(ff,'output/FigS9.eps','-depsc','-r300')

