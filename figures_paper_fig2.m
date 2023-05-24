%% NOTES
% todo
% make sure there are orange bars in the histogram i
% Small thing in Fig. 1, panel b) change" lambda_bg score for blocks of size 64" to “Estimates of lambda_bg”. In panel c) change to “Fit to truncated EMCCD distribution for tile (i,j)” [and replace I and j with relevant numbers)
% % In the figure caption we say that each tile is labeled by its row index (i) and column index (j), in the form (i,j). 

% Move panel B to figure 2

outFig2 = fullfile(figFold,'Fig1.eps');
% filename = 'data\100x_gain100_lamp50_018.tif';
filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\20x_gain100_lamp100_004.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp10_022.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif';
% filename = 'C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\100x_gain100_lamp100_013.tif';

chipPars.inputImage = filename;
[image] = Core.load_image(filename);

chipParsSpec = Core.chippars_specific_gain(chipParsAll,3);
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
title(['Estimate of $\lambda_{bg}$],'Interpreter','latex')

nexttile
imagesc(reshape(intthreshPars,[sqrt(N) sqrt(N)])');colorbar;colormap gray
title(['$N_{thresh}$ scores'],'Interpreter','latex')


nexttile
imagesc(reshape(cellfun(@(x) x.passthresh,statsAll),[sqrt(N) sqrt(N)])');colorbar;colormap gray
title(['Scores passing p-val thresh'],'Interpreter','latex')

%% REDO FIG2 with tiles

%
% 
idx = 15;

    figure
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
    nexttile([2 1])
    pixelsize = 160;

    sampIm = mat2gray(image.imAverage);
    minInt = min(sampIm(:));
    medInt = median(sampIm(:));
    maxInt = max(sampIm(:));
    J = imadjust(sampIm,[minInt min(1,4*medInt)]);
    matrixBlocksJ = reshape(permute(reshape( double(J),T,size( J,1)/T,T,[]),[1,3,2,4]),T,T,[]);
%     out = imtile(matrixBlocksJ,'thumbnailsize',[64 64]);
%     figure,imshow(out)

    IM3 = padarray(matrixBlocksJ,[2 2],1,'both');
    out = imtile(pagetranspose(IM3),'thumbnailsize',[64 64]);
%     figure,imshow(out')

    imshow(out','InitialMagnification','fit');
    %imshow(images.imAverage/max(images.imAverage(:)))
    hold on    
    % scale bar (ten microns in number of pixels)
    nPixels = 1e4/pixelsize;
    x = [5, 5 + nPixels ];
    y = [0.9*size(sampIm,1) , 0.9*size(sampIm,1)];
    plot(x,y,'Linewidth',8,'Color',[1 1 1])
    text(0,0.05,'10 microns','Fontsize',10,'Color',[1 1 1],'Units','normalized')
    title('(a)','Interpreter','latex')

    nexttile
    outL = imtile(lambdaPars);
%     figure,imshow(outL)
    imshow(reshape(lambdaPars,[sqrt(N) sqrt(N)]),[min(lambdaPars) max(lambdaPars)]);colorbar;colormap gray
title('(b) Estimates of $\lambda_{bg}$' ,'Interpreter','latex')


% Plot single tile
structRes = statsAll{idx};
intThreshBg =   intthreshPars(idx);
histAll = statsAll{idx}.histAll;
stats = statsAll{idx}.stats;
    nexttile    
binPos = 1:structRes.LU(2) + 0.5;
[minVal , idx] = min(abs(binPos - intThreshBg));
h1 = bar(binPos(1:idx),histAll(1:idx),1); 
set(h1,'FaceColor',[0.4 0.6 0.9])
set(h1,'EdgeColor','none')
hold on
h2 = bar(binPos(idx+1:end),histAll(idx+1:end-1),1); 
set(h2,'FaceColor',[1 0.5 0.3])
set(h2,'EdgeColor','none')
hold on
binCountsFit = stats.nBg.*structRes.pdf;
% binPos = binEdges(1:end-1) + diff(binEdges)/2;
plot(binCountsFit,'--','Color','black','LineWidth',2)

% Set figure labels, etc
xlabel('Image counts','Interpreter','latex')
ylabel('Histogram counts','Interpreter','latex')
% set(gca,'Fontsize',15)
% axis([30 80 0 46000])
[a,b] = ind2sub([8 8],15)
title(['(c) Fit for tile \{',num2str(a),',',num2str(b) , '\}'],'Interpreter','latex')
% axis equal
pbaspect([1 0.8 0.8])
legendEntry = strcat(['Fit, $\lambda_{bg} =  ' num2str(lambdaBg,2) ', N_{icr}^{bg}=' num2str(intThreshBg) '$']);
lgnd = legend('Image counts, true background','Image counts, not true background',legendEntry,'Interpreter','latex');
lgnd.Layout.Tile = 'south';
% print('C:\Users\Lenovo\postdoc\PAPERS\emccd-paper\draft\Figs\Fig4.eps','-depsc','-r300')
print(outFig2,'-depsc','-r300');

