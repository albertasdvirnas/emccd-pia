%% Generate/load image data
% load generated image
SNRVals = [2.00];  
filenames = simulate_random_beads_full(100, 20, 100, 1/1000,SNRVals, 1, 'synthSingleFrameNew'); % zoom 100x
imageName = filenames{1}{1}{1};
chipParsCur.gain  = 20;
chipParsCur.countOffset  = 27;
chipParsCur.roNoise  = 1.44;
chipParsCur.adFactor  = 36;
p_white_bg = 1E-2;

chipPars = chipParsCur;
% we use single image
imageName = strrep(imageName,'.tif','.mat');
data = importdata(imageName); 
im = data.image;      
images.imAverage = reshape(double(im),size( data.groundTruthImage));
images.imageName = imageName;
groundTruthImage = data.groundTruthImage > 0;
sortI = sort(images.imAverage(:));

% find Nthresh for each intensity
intensRange = unique(images.imAverage(:));
%
FDRGT = zeros(1,length(intensRange));
FORGT =  zeros(1,length(intensRange));

% GT
for i=1:length(intensRange)
    binarizedImage2 = imbinarize(images.imAverage,intensRange(i));
    [fpr,fnr,tmr,FDRGT(i),FORGT(i)] = compare_regions_to_ground_truth_beads_pixelbased(binarizedImage2,groundTruthImage);      
end

% GT threshold
qStar=0.1:0.1:0.9;

nGT = zeros(1,length(qStar));
fdrGT = zeros(1,length(qStar));

for i=1:length(qStar)
    idxQstar = find(FDRGT < qStar(i),1,'first');
    if isempty(idxQstar)
        nGT(i) = intensRange(end);
        fdrGT(i) = FDRGT(end);
    else
        nGT(i) = intensRange(idxQstar);
        fdrGT(i) = FDRGT(idxQstar);

    end
end


%% Calculate FDR/FOR for different alphaStar/qStar
% FDR  method comparison
Nthresh1 = zeros(1,length(qStar));
Nthresh2 = zeros(1,length(qStar));
% Nthresh3 = zeros(1,length(qStar));

import Core.fdr_based_thresh;
import Core.for_based_thresh;
import Core.calc_bounds;
import Core.calc_bh_threshold;

lambdaBg = 39; % lambda estimate
[~,U,~,~] = calc_bounds(lambdaBg,chipParsCur.gain,chipParsCur.adFactor,chipParsCur.countOffset,chipParsCur.roNoise); % max intensity for integration
stopIntensity = quantile(sortI,0.25)+0.5; % intensity for nbg estimation


for i=1:length(qStar)
    [nOutliers,hasOutliers,stats, Nthresh1(i)] =  fdr_based_thresh(sortI, ...
        lambdaBg,chipParsCur.gain,chipParsCur.adFactor,chipParsCur.countOffset,chipParsCur.roNoise,qStar(i),U,stopIntensity);
    [threshold, nOutliers, outliers, hasOutliers,Nthresh2(i) ] = calc_bh_threshold(sortI, lambdaBg, chipParsCur.gain,chipParsCur.adFactor,chipParsCur.countOffset,chipParsCur.roNoise,qStar(i));
%     [nOutliers,hasOutliers,stats, Nthresh3(i)] =  for_based_thresh(intensities,sortI, m, lambdaBg,gain,adFactor,countOffset,roNoise,qStar(i),U,stopIntensity);

end

figure,plot(qStar,Nthresh1-0.5)
hold on
plot(qStar,Nthresh2)
plot(qStar,nGT)

xlabel('qStar')
ylabel('Nthresh')
legend({'FDRest','BH','Ground Truth'})

%%





% FDR  method comparison
alphaStar=0.001:0.01:0.2;

nGTalpha = zeros(1,length(alphaStar));
forGTalpha = zeros(1,length(alphaStar));

for i=1:length(alphaStar)
    idxQstar = find(FORGT< alphaStar(i),1,'last');
    if isempty(idxQstar)
        nGTalpha(i) = intensRange(end);
        forGTalpha(i) = FORGT(end);
    else
        nGTalpha(i) = intensRange(idxQstar);
        forGTalpha(i) = FORGT(idxQstar);

    end
end

% load('tempm.mat')

Nthresh3 = zeros(1,length(alphaStar));
FOR2 =  zeros(1,length(alphaStar));


for i=1:length(alphaStar)
    [nOutliers,hasOutliers,stats, Nthresh3(i)] =  for_based_thresh(sortI, lambdaBg,chipParsCur.gain,chipParsCur.adFactor,chipParsCur.countOffset,chipParsCur.roNoise,alphaStar(i),U,stopIntensity);
    FOR2(i) = stats.finalFor;
end


figure; plot(alphaStar,Nthresh3-0.5)
hold on
plot(alphaStar,nGTalpha)

xlabel('alphastar')
ylabel('Nthresh')
legend({'FORest','FOR GT'})


figure; plot(alphaStar,FOR2)
hold on
plot(alphaStar,forGTalpha)
xlabel('alphastar')
ylabel('FOR')
legend({'FOR estimate','FOR GT'})

% %%
% % v = VideoWriter('myFile','Archival');
% % open(v)
% % writeVideo(v,images.imAverage )
% % close(v)
% 
% % T = T / max(T(:));
% VidObj = VideoWriter('movie.avi', 'Uncompressed AVI'); %set your file name and video compression
% VidObj.FrameRate = 5; %set your frame rate
% open(VidObj);
% for f = unique(sortI)'
%     writeVideo(VidObj, double(images.imAverage>f));
% end
% close(VidObj);