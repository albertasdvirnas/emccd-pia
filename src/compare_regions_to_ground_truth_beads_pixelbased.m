function [fpr,fnr,tmr,FDR,FOR] = compare_regions_to_ground_truth_beads_pixelbased(signalLabelIm,groundTruthImage)

    %
    % Compares regions obtained using an image segmentation to 
    % ground truth results. 
    %
    % Args:
    %
    % signalLabelIm = matrix of the same size as original image, where the
    %           value of each matrix element gives a white region label
    %           (a '0' represents a black region)
    % groundTruthImage = ground truth image, where a '1' corresponds to a
    %                    a pixel where there are fluorescent molecules (signal),
    %                    and a '0' corresponds to a pixel where there are
    %                    no fluroescent molecules (background).
    %
    % Returns:
    %
    % fpr = false positive rate  FP/(FP+TN)
    % fnr = false negative rate FN/(FN+TP)
    % tmr = "error rate" 1-accuracy (1-ACC)
    % FDR
    % FOR
    %
    % Dependencies: None.
    %
    
    % Make a binary image based on the image segmentation result
    [nRows,nCols] = size(signalLabelIm); 
    binarySegImage = zeros(nRows,nCols);
    idx = find(signalLabelIm > 0);
    binarySegImage(idx) = 1; 
    noOfPixels = numel(signalLabelIm);
   
    
    % Compare the binary image to the ground truth
    
    % False negative rate (FNR) FN/(FN+TP)
    n_white_sig = length (find(groundTruthImage(:)== 1 & binarySegImage(:) == 1));
    n_sig = length( find(groundTruthImage == 1 ));
    fnr = 1 - n_white_sig/n_sig;
    
    % False Positive rate p(white|bg)= FPR = FP/(FP+TN)
    n_black_bg = length(find(groundTruthImage == 0 & binarySegImage == 0));
    n_bg = length( find(groundTruthImage == 0 ));
    fpr = 1 - n_black_bg/n_bg;
   
    % 1-ACCURACY (ACC)
    noOfCorrect = n_white_sig + n_black_bg;
    tmr = 1 - noOfCorrect/noOfPixels;

    % FDR = FP/(FP+TN)
    n_white_bg = length (find(groundTruthImage(:)== 0 & binarySegImage(:) == 1));
    n_white = length( find(binarySegImage(:) == 1));

    FDR = n_white_bg/n_white;

    % FOR = FN/(FN+TN)
    n_black_sig = length (find(groundTruthImage(:)== 1 & binarySegImage(:) == 0));
    n_black = length( find(binarySegImage == 0 ));
    FOR = n_black_sig/n_black;
    
end