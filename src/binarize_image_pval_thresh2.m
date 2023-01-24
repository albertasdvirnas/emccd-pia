

function [binarizedImage , intThresh ] = binarize_image_pval_thresh(...
              imageInput, pValThresh ,lambdaBg ,  gain,adFactor,offset,roNoise)

    % 
    % Binarizes an image using a p-value threshold
    %
    % Input:
    %
    % imageInput = image
    % pValThresh = p-value threshold. Pixels with a p-value 
    %             below this threshold are turned white ( = 1)
    %             Pixels with a p-value above this threshold 
    %             are black ( = 0).
    % lambdaBg = Poissoin parameter for the background
    % chipPars = struct containing chip-parameters.
    %
    % Output:
    %
    % binarizedImage = binarized image 
    % intThresh = intensity threshold at the specified p-value
    %
    % Dependencies: emccd_distribution/inverse_cdf_emccd.m
    % 
    
    % Hard-coded variables
%     N= 2^8;     % number of integration points for evaluating the EMCCD CDF
%     tol=1E-5;   % accuracy for the inverse CDF calculation.
    
% inverse cdf to get intThresh
    % first these are quickly calculated bounds
    [L,U,EX,STD] = calc_bounds(lambdaBg,gain,adFactor,offset,roNoise);
    
    % that give intensities to calculate over
    intensities = ceil(L):floor(U);
    
    % cdf, in this case do not need to be truncated
    [pdf,cdf] = pdf_cdf_from_characteristic_fun(intensities,lambdaBg,gain,adFactor,offset,roNoise);
    
    % find the value where pvalue=1-cdf > pValThresh
    intThresh = find(1-cdf < pValThresh,1,'first');
%     intensities(find(1-cdf >pValThresh,1,'last'));

    intThresh = intensities(intThresh);

%     intThresh = inverse_cdf_emccd( 1-pValThresh , lambdaBg , chipPars , N , tol);
    intThresh = floor(intThresh)-0.5;  % since intensities are integers 
                                       % in experimental images, 
                                       % we set the threshold to be a half-integer
   
    % Binarize using intensity threshold
    binarizedImage = imbinarize(imageInput,intThresh);


end