nn=10;
lambdaBg = zeros(1,nn);
for i=1:nn
SNRVals = [10.00 ];  
filenames = simulate_random_beads_full(100, 20, 100, 1/1000,SNRVals, 1, 'synthSingleFrame22'); % zoom 100x

% controls p(white|bg) (p_thresh)
p_white_bg = 1E-2;

% input parameter values
 
     
% output filename
outputFilename = 'fdr_for_results_temp';

imageName = filenames{1}{1}(end);
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

qStar = 0.9;



chipPars.adFactor = 36; % ADU
chipPars.countOffset = 27; % offset (bias)
chipPars.roNoise = 1.44; % noise std
chipPars.gain = 20;
[lambdaBg(i),intThreshBg] = fig2_calibration(chipPars,[],images,qStar,1);
end