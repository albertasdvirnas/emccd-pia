function [fdr,frr,tmr] = compare_regions_to_ground_truth_beads_pixelbased(signalLabelIm,groundTruthImage)

    %
    % Compares regions obtained using an image segmentation to 
    % ground truth results. 
    %
    % Input:
    %
    % signalLabelIm = matrix of the same size as original image, where the
    %           value of each matrix element gives a signal region label
    %           (a '0' represents a bg region)
    % groundTruthImage = ground truth image, where a '1' corresponds to a
    %                    a pixel where there are fluorescent molecules (signal),
    %                    and a '0' corresponds to a pixel where there are
    %                    no fluroescent molecules (background).
    %
    % Output:
    %
    % fdr = false detection rate
    % frr = false rejection rate
    % tmr = total misclassification rate
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
    
    % False rejection rate
    noOfCorrectSignal = length (find(groundTruthImage(:)== 1 & binarySegImage(:) == 1));
    noOfSignal = length( find(groundTruthImage == 1 ));
    frr = 1 - noOfCorrectSignal/noOfSignal;
    
    % False detection rate
    noOfCorrectBg = length(find(groundTruthImage == 0 & binarySegImage == 0));
    noOfBg = length( find(groundTruthImage == 0 ));
    fdr = 1 - noOfCorrectBg/noOfBg;
   
    % Total misclassification rate
    noOfCorrect = noOfCorrectSignal + noOfCorrectBg;
    tmr = 1 - noOfCorrect/noOfPixels;
    
    
end