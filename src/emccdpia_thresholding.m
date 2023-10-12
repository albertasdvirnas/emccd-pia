function [results,stats] = emccdpia_thresholding(imageNames, SNRVals, chipPars, outFig,pvalthresh,ifplot)

    % p-value threshold using the background distribution
     outputFilename = 'eval_binarization_results_temp';
 
     nImages = length(imageNames);    
    lambdaBgSave = zeros(1,nImages);
    intThreshBgSave = zeros(1,nImages);
    
    results = struct();

    results.fprEst = zeros(1,nImages);
    results.fnrEst = zeros(1,nImages);
    results.tmrEst = zeros(1,nImages);
    results.otsuFdr = zeros(1,nImages);
    results.otsuFrr = zeros(1,nImages);
  
    stats = cell(1,nImages);

    % Estimate lambda_bg for each image
    for j = 1:nImages        
        disp('-------------')
        disp(['image with ','SNR = ', num2str(SNRVals(j))])
        disp('------------- ')           
        % Pre-processing      
        imageName = strrep(imageNames{j},'.tif','.mat');
        data = importdata(imageName); 
        im = data.image;      
        images.imAverage = reshape(double(im),size( data.groundTruthImage));
        images.imageName = imageName;
%         snr(j) = data.snr;  % signal-to-noise ratio

        % Extract chip, optics parameters, etc from input file
        lambdaBgGroundTruth = data.lambdabg;
        groundTruthImage = data.groundTruthImage > 0;

        disp('Estimating lambda_bg.'); 
        [lambdaBg, intThreshBg, stats{j}] =  emccdpia_estimation(chipPars,[],images,pvalthresh,ifplot);
        stats{j}.intThreshBg = intThreshBg;
        lambdaBgSave(j) = lambdaBg;
        % Calculate rates when compared to GT
        [results.fprGT(j),results.fnrGT(j),results.tmrGT(j),results.fdrGT(j),results.forGT(j)] = compare_regions_to_ground_truth_beads_pixelbased(stats{j}.binarizedImage,groundTruthImage);      


        % ---- Otsu based binarization and segmentation ----
        % ---- (segmentation result not used) ----- 
        %
        segOutputOtsu = segmentation_otsu(images);
        signalLabelImOtsu = segOutputOtsu.labelIm;
        binarizedImageOtsu = segOutputOtsu.binarizedImage;


         % Calculate performance results 
        [results.fprOtsu(j),results.fnrOtsu(j),results.tmrOtsu(j),results.fdrOtsu(j),results.forOtsu(j)] = compare_regions_to_ground_truth_beads_pixelbased(binarizedImageOtsu,groundTruthImage);
      
        results.intThreshEst(j) = intThreshBg;
        results.lambdaBgEst(j) = lambdaBg;
        results.lambdaBgGroundTruth(j) = lambdaBgGroundTruth;
        results.snr(j) =  data.snr; 

    end
    
    results.fprEst = cellfun(@(x) x.stats.fpr, stats);
    results.fnrEst = cellfun(@(x) x.stats.fnr, stats);
    results.tmrEst = cellfun(@(x) x.stats.misClassRate, stats);
    results.forEst = cellfun(@(x) x.stats.for, stats);
    results.fdrEst = cellfun(@(x) x.stats.fdr, stats);

    % Save results 
%     results = struct();
%     results.lambdaBgGroundTruth = lambdaBgGroundTruth;    
%     results.lambdaBgEst = lambdaBgSave;
%     results.intThreshEst = intThreshBgSave;     
    results.stats = stats;

    save(outputFilename,'results')
    
    plot_binarization_evaluation_results_beads(outFig,results);
end




function plot_binarization_evaluation_results_beads(outFig3,results)

    %
    % Script for analyzing the output of the function 
    % evaluate_seg_performance_beads
    %
   
    % Load results from file 
%     results = importdata('eval_binarization_results_temp.mat');
    snr = [results.snr];

     % Sort according to SNR values
    [snr,idx] = sort(snr);

%%     fS = 12;
    figure
    tiledlayout(3,2,'TileSpacing','tight');
 
    nexttile
    plot(snr,results.lambdaBgEst,'m-o','linewidth',2);
    hold on
    y = results.lambdaBgGroundTruth.*ones(1,length(snr));
    plot(snr,y,'black--','linewidth',2)
    ylim([min([results.lambdaBgEst y])-0.1 max([y results.lambdaBgEst])+0.1])
%     text(4,37.5,'Ground-truth')
%     text(4,39.5,'Estimated')
    lgnd2 = legend({'Estimated','Ground-truth'})
    lgnd2.Layout.Tile = 'north';

    hold off
    ylabel('\lambda_{bg}')
    title('a) Estimated $\lambda_{bg}$ vs. GT','Interpreter','latex')
  
    nexttile
    plot(snr,results.fprGT,'-d','linewidth',2);
    hold on
    plot(snr,results.fprOtsu,'--x','linewidth',2);
    plot(snr,results.fprEst,'.-','linewidth',2)
    hold off
%     legend(["non-Otsu","Otsu's method","estimate"],'Location','east','Fontsize',15)
%     xlabel('SNR','Interpreter','latex','Fontsize',fS)
    xlim([1 10])
    ylabel('FPR','Interpreter','latex')
    fig = gcf;
    fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 2.5];
    title('b) FPR','Interpreter','latex')
%     set(gca,'Fontsize',15)

    % Plot FNR, pixelbased
    nexttile
    plot(snr,results.fnrGT,'-d','linewidth',2);
    hold on
    plot(snr,results.fnrOtsu,'--x','linewidth',2);
    plot(snr,results.fnrEst,'.-','linewidth',2)
    hold off
%     legend("non-Otsu","Otsu's method","estimate",'Fontsize',15)
%     xlabel('SNR','Interpreter','latex','Fontsize',fS)
    xlim([1 10])
    ylabel('FNR','Interpreter','latex')
    fig = gcf;
    fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 2.5];
    title('c) FNR','Interpreter','latex')
    xlabel('SNR','Interpreter','latex')


    % Plot TMR, pixelbased
%     figure()
    nexttile
    plot(snr,results.tmrGT,'-d','linewidth',2);
    hold on
    plot(snr,results.tmrOtsu,'--x','linewidth',2);
    plot(snr,results.tmrEst,'.-','linewidth',2)
    hold off
    lgnd = legend("EMCCD-PIA binarization","Otsu's method","EMCCD-PIA a priori estimate")
    lgnd.Layout.Tile = 'north';
    fig = gcf;
    fig.PaperUnits = 'inches';
    xlabel('SNR','Interpreter','latex')
    xlim([1 10])
    ylabel('1-ACC','Interpreter','latex')

%     fig.PaperPosition = [0 0 4 2.5];
    title('d) 1-ACC','Interpreter','latex')

      nexttile
    plot(snr,results.fdrGT,'-d','linewidth',2);
    hold on
    plot(snr,results.fdrOtsu,'--x','linewidth',2);
    plot(snr,results.fdrEst,'.-','linewidth',2)
    hold off
    xlim([1 10])
    ylabel('FDR','Interpreter','latex')
    fig = gcf;
    fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 2.5];
    title('e) FDR','Interpreter','latex')
    xlabel('SNR','Interpreter','latex')

    nexttile
    plot(snr,results.forGT,'-d','linewidth',2);
    hold on
    plot(snr,results.forOtsu,'--x','linewidth',2);
    plot(snr,results.forEst,'.-','linewidth',2)
    hold off
    xlim([1 10])
    ylabel('FOR','Interpreter','latex')
    fig = gcf;
    fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 2.5];
    title('f) FOR','Interpreter','latex')
    xlabel('SNR','Interpreter','latex')

    print(outFig3,'-depsc','-r300')

  
            
            
end
% 
% function probabilistic_binarization()
% 
% end
 