outFig2 = fullfile(figFold,'Fig2.eps');
filename = 'data\100x_gain100_lamp50_018.tif';
filename = filenames{1}{1}{5};
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp10_022.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif';

chipPars.inputImage = filename;
[image] = Core.load_image(filename);

chipParsSpec = Core.chippars_specific_gain(chipParsSAll,3);
alphaStar = 0.01; % tnr control
chipParsSpec.pixelsize = 160;

T = 64;
matrixBlocks = reshape(permute(reshape( double(image.imAverage),T,size( image.imAverage,1)/T,T,[]),[1,3,2,4]),T,T,[]);

N = size(matrixBlocks,3);
scores = nan(1,size(matrixBlocks,3));
lambdaPars = nan(1,size(matrixBlocks,3));
intthreshPars = nan(1,size(matrixBlocks,3));

statsAll = cell(1,size(matrixBlocks,3));
for idx = 1:size(matrixBlocks,3);
    image2.imAverage = matrixBlocks(:,:,idx);

    tic
    [lambdaBg, intThreshBg, stats] = emccdpia_estimation(chipParsSpec,outFig2,image2,alphaStar,0);
    scores(idx) = min(stats.chi2Score);
    lambdaPars(idx) = lambdaBg;
    intthreshPars(idx) = intThreshBg;
    statsAll{idx} = stats;

    toc
end

figure,nexttile;imagesc(reshape(scores,[sqrt(N) sqrt(N)])');colorbar;colormap gray
title(['$\chi^2$ scores for blocks of size = ', num2str(T)],'Interpreter','latex')
nexttile
imagesc(reshape(lambdaPars,[sqrt(N) sqrt(N)])');colorbar;colormap gray
title(['$\lambda_{bg}$ scores for blocks of size = ', num2str(T)],'Interpreter','latex')

nexttile
imagesc(reshape(intthreshPars,[sqrt(N) sqrt(N)])');colorbar;colormap gray
title(['$N_{thresh}$ scores'],'Interpreter','latex')


nexttile
imagesc(reshape(cellfun(@(x) x.passthresh,statsAll),[sqrt(N) sqrt(N)])');colorbar;colormap gray
title(['Scores passing p-val thresh'],'Interpreter','latex')
