%% output folder for results
figFold = 'output';

%% Simulate synthetic / for 0 gain we zoom in and take 100x image
zoom = 100;
simulate_random_beads_full(zoom,1,0,1/300); % zoom 100x
%
zoom = 20;
gain = [1 12 20 46];
gainNames = [0 50 100 300];
particleDensity = 1/300;
simulate_random_beads_full(zoom, gain, gainNames, particleDensity); % zoom 20x, particles are 5 times smaller


% Calculate moments for synthetic data
calculate_moments('data\100x\*.tif'); 
calculate_moments('data\20x\*.tif');

% For real data
foldRealData = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 beads high conc';
calculate_moments(fullfile(foldRealData,'100x\*.tif'));
foldRealData = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics';
calculate_moments(fullfile(foldRealData,'20x\*.tif'));

%% FigS1.eps - synthetic data.
no_gain_mat = 'synth100x.mat';
gain_mat = 'synth20x.mat';

outFig = fullfile(figFold,'FigS0.eps');
chipParsS = fig1_calibration(no_gain_mat, gain_mat,outFig);

%% FigS2.eps / on synthetic image


%% Fig1.eps
no_gain_mat = '100x.mat';
gain_mat = '20x.mat';

outFig = fullfile(figFold,'Fig1.eps');
chipPars = fig1_calibration(no_gain_mat, gain_mat,outFig);

%% Fig2.eps
outFig2 = fullfile(figFold,'Fig2.eps');
chipPars.inputImage = 'data\100x_gain100_lamp50_018.tif';
[lambdaBg,intThreshBg] = fig2_calibration(chipPars,outFig2);

%% Fig3.eps
    SNRVals = [3.00 , 3.50 , 4.00 , 4.50 , ...
               5.00 , 5.50 , 6.00 , 6.50 , 7.00 , 7.50 , 8.00 , 8.50 , ...
               9.00 , 9.50 , 10.00 ];  
filenames = simulate_random_beads_full(100, 20, 100, 1/1000,SNRVals, 1, 'synthSingleFrame'); % zoom 100x

outFig3 = 'output\Fig3.eps';
chipParsCur = struct();
chipParsCur.gain = mean(chipParsS.gainQ(3,2)); % based on gain in the image.. 50 100 300..
chipParsCur.adFactor = mean(chipParsS.aduQ(2));
chipParsCur.countOffset = mean(chipParsS.deltaQ(3,2));
chipParsCur.roNoise = mean(chipParsS.sigmaQ(3,2));
% chipParsCur.gain  = mean(chipParsS.gain{3});
% chipParsCur.countOffset  = mean(chipParsS.countOffset{3});
% chipParsCur.roNoise  = mean(chipParsS.roNoise{3});
% chipParsCur.adFactor  = mean(chipParsS.adFactor);
% 
fig3_probabilitstic_segmentation(filenames{1}{1}(1:end),SNRVals,chipParsCur,outFig3)

% fig3_probabilitstic_segmentation(filenames{1}{1}([1 2]),SNRVals,chipParsCur,outFig3)
% [lambdaBg,intThreshBg] = fig2_calibration(chipPars,outFig2);
%% Fig4.eps

% for the beads
imagefiles = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif', 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV1_DAPI\20x_gain300_lamp100_001.tif'};
imagefiles = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp50_018.tif',...
    'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif'};
% imagefiles = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif',...
%     'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif'};


chipParsCur1.gain = mean(chipPars.gainQ(3,2)); % based on gain in the image.. 50 100 300..
chipParsCur1.adFactor = mean(chipPars.aduQ(2));
chipParsCur1.countOffset = mean(chipPars.deltaQ(3,2));
chipParsCur1.roNoise = mean(chipPars.sigmaQ(3,2));
    
% chipParsCur1.gain = mean(chipPars.gainQ(4,2)); % based on gain in the image.. 50 100 300..
% chipParsCur1.adFactor = mean(chipPars.aduQ(2));
% chipParsCur1.countOffset = mean(chipPars.deltaQ(4,2));
% chipParsCur1.roNoise = mean(chipPars.sigmaQ(4,2));



chipParsFig4 = {chipParsCur1,chipParsCur1}; % the two experiments will have slightly varying parameters because of different Gain setting

outFig4 = { fullfile(figFold,'Fig4a.eps'),fullfile(figFold,'FigS4a.eps');fullfile(figFold,'Fig4b.eps'),fullfile(figFold,'FigS4b.eps')};
% outFig3 = 'C:\Users\Lenovo\postdoc\PAPERS\emccd-paper\draft\Figs\Fig4a.eps';


fig4_segmentation(imagefiles,chipParsFig4,outFig4);
% fig4_segmentation(chipPars,'beads_low_conc_100x_gain100_lamp100_013.tif')