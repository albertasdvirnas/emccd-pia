%% output folder for results
figFold = 'output';

%% Simulate synthetic / for 0 gain we zoom in and take 100x image

gain = [1 12 20 46];
gainNames = [0 50 100 300];
particleDensity = 1/300;

zoom = 100;
simulate_random_beads_full(zoom,gain,gainNames,particleDensity); % zoom 100x
%
zoom = 20;
simulate_random_beads_full(zoom, gain, gainNames, particleDensity); % zoom 20x, particles are 5 times smaller


% Calculate moments for synthetic data
calculate_moments('data\100x\*.tif'); 
calculate_moments('data\20x\*.tif');

% For real data
limits = [100 200 260 360];
foldRealData = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 beads high conc';
calculate_moments(fullfile(foldRealData,'100x\*.tif'));
foldRealData = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics';
calculate_moments(fullfile(foldRealData,'20x\*.tif'));

%% FigS0.eps - synthetic data.
no_gain_mat = 'synth100x.mat';
gain_mat = 'synth100x.mat';

rng('default')

outFig = fullfile(figFold,'FigS2.eps');
[chipParsS,chipParsSAll] = fig1_calibration(no_gain_mat, gain_mat,outFig);

%% FigS1.eps - real data
rng('default')
no_gain_mat = '100x.mat';
gain_mat = '20x.mat';

outFig = fullfile(figFold,'FigS1.eps');
[chipPars,chipParsAll] = fig1_calibration(no_gain_mat, gain_mat,outFig);

table_1

mdl = fitlm(gainNames(2:end),chipParsAll.gain(2:end))
p = signrank(gainNames(2:end),chipParsAll.gain(2:end))
%% FigS2.eps / on synthetic image

%% Fig1.eps - real experiment
figures_paper_fig2 %- Fig1 in the paper
% individuak block
outFig2 = fullfile(figFold,'Fig2.eps');
filename = 'data\100x_gain100_lamp50_018.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp10_022.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif';

chipPars.inputImage = filename;
[image] = Core.load_image(filename);

T = 128;
matrixBlocks = reshape(permute(reshape( double(image.imAverage),T,size( image.imAverage,1)/T,T,[]),[1,3,2,4]),T,T,[]);

idx = 12;
image2.imAverage = matrixBlocks(:,:,idx);


% image.imAverage=image.imAverage(100:200,260:360);
% image.imAverage=image.imAverage(100:200,350:420);
% image.imAverage=image.imAverage(200:300,200:300);

% image.imAverage=image.imAverage(150:250,300:400);

% image.imAverage=image.imAverage(400:500,400:500);

chipParsSpec = Core.chippars_specific_gain(chipParsAll,3);
alphaStar = 0.01; % tnr control
chipParsSpec.pixelsize = 160;
tic
[lambdaBg, intThreshBg, stats] = emccdpia_estimation(chipParsSpec,outFig2,image2,alphaStar,1);
min(stats.chi2Score)
toc
%% Fig2.eps
rng('default')

SNRVals =3:0.5:10;  
filenames = simulate_random_beads_full(100, 20, 100, 1/1000,SNRVals, 1, 'synthSingleFrameNew'); % zoom 100x

outFig3 = 'output\Fig2.eps';
% chipParsCur = struct();
% chipParsCur.gain = mean(chipParsS.gainQ(3,2)); % based on gain in the image.. 50 100 300..
% chipParsCur.adFactor = mean(chipParsS.aduQ(2));
% chipParsCur.countOffset = mean(chipParsS.deltaQ(3,2));
% chipParsCur.roNoise = mean(chipParsS.sigmaQ(3,2));
chipParsCur = Core.chippars_specific_gain(chipParsSAll,3);


chipParsCur.pval = 0.01;
% emccdpia_thresholding(filenames{1}{1}(6),SNRVals,chipParsCur,outFig3,chipParsCur.pval,1);

[results,stats] = emccdpia_thresholding(filenames{1}{1}(1:end),SNRVals,chipParsCur,outFig3,chipParsCur.pval,0);


figure_sup; % supplementary figure S4

% chipParsCur.gain  = mean(chipParsS.gain{3});
% chipParsCur.countOffset  = mean(chipParsS.countOffset{3});
% chipParsCur.roNoise  = mean(chipParsS.roNoise{3});
% chipParsCur.adFactor  = mean(chipParsS.adFactor);
%     
% chipParsCur.method = 'FOR';
% chipParsCur.alphaStar = 0.05;
%  [results,stats] = emccdpia_thresholding(filenames{1}{1}(5),SNRVals,chipParsCur,outFig3,chipParsCur.pval,1);

% chipParsCur.method = 'FDR';
% chipParsCur.alphaStar = 0.9;
% emccdpia_thresholding(filenames{1}{1}(1:end),SNRVals,chipParsCur,outFig3,chipParsCur.alphaStar,1);
% 
% 
% chipParsCur.gain  = 20;
% chipParsCur.countOffset  = 27;
% chipParsCur.roNoise  = 1.44;
% chipParsCur.adFactor  = 36;
%  [results,stats] = emccdpia_thresholding(filenames{1}{1}(1:end),SNRVals,chipParsCur,outFig3,0.1,0);
%  [results,stats] = emccdpia_thresholding(filenames{1}{1}(1),SNRVals,chipParsCur,outFig3,chipParsCur.pval,1);

% fig3_probabilitstic_segmentation(filenames{1}{1}(1),SNRVals,chipParsCur,outFig3)

% fig3_probabilitstic_segmentation(filenames{1}{1}([1 2]),SNRVals,chipParsCur,outFig3)
% [lambdaBg,intThreshBg] = fig2_calibration(chipPars,outFig2);


%% Fig3.eps

% for the beads
% imagefiles = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif', 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV1_DAPI\20x_gain300_lamp100_001.tif'};
imagefiles = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp50_018.tif',...
    'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif'};
% imagefiles = {'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif',...
%     'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif'};


% chipParsCur1.gain = mean(chipPars.gainQ(3,2)); % based on gain in the image.. 50 100 300..
% chipParsCur1.adFactor = mean(chipPars.aduQ(2));
% chipParsCur1.countOffset = mean(chipPars.deltaQ(3,2));
% chipParsCur1.roNoise = mean(chipPars.sigmaQ(3,2));
chipParsCur1     = Core.chippars_specific_gain(chipParsAll,3);
% chipParsCur1.gain = mean(chipPars.gainQ(4,2)); % based on gain in the image.. 50 100 300..
% chipParsCur1.adFactor = mean(chipPars.aduQ(2));
% chipParsCur1.countOffset = mean(chipPars.deltaQ(4,2));
% chipParsCur1.roNoise = mean(chipPars.sigmaQ(4,2));



chipParsFig4 = {chipParsCur1,chipParsCur1}; % the two experiments will have slightly varying parameters because of different Gain setting

outFig4 = { fullfile(figFold,'Fig4a.eps'),fullfile(figFold,'FigS4a.eps');fullfile(figFold,'Fig4b.eps'),fullfile(figFold,'FigS4b.eps')};
% outFig3 = 'C:\Users\Lenovo\postdoc\PAPERS\emccd-paper\draft\Figs\Fig4a.eps';


outres = fig4_segmentation(imagefiles, chipParsFig4, outFig4, 0.01,0.01);
% fig4_segmentation(chipPars,'beads_low_conc_100x_gain100_lamp100_013.tif')


%% Otsu vs Adaptive FigS9

fig_adaptive

%% Segmentation another example FigS10
figS10_segmentation