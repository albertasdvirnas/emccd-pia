% Simulate random beads images - full example. We simulate for both 0 gain
% and for selected intensities
function [filenames,chipPars,data] = simulate_random_beads_full(zooms,gainF,gainName,particle_density,SNR,numFrames,outFold,dims,gtImage, placements)

if nargin < 1
    zooms = [20 100];
    gainF = [ 1 12 20 46];
    gainName = [0 50 100 300];
    particle_density = [1/400 1/800*5];
end
% chipParameters
chipPars.adFactor = 36; % ADU
chipPars.countOffset = 27; % offset (bias)
chipPars.roNoise = 1.44; % noise std
chipPars.numFrames = 50;
if nargin < 8
    chipPars.nRows = 512; % num rows
    chipPars.nCols = 512; % num columns
else
    chipPars.nRows = dims(1); % num rows
    chipPars.nCols = dims(2); % num columns
end
chipPars.maxRadius = 10;
chipPars.lambdaBg = 38;
chipPars.lambdaSigBase = 1;

if nargin < 5
    chipPars.lampPower = 10*[0:10:100];
else
%     chipPars.lampPower = lampPower;
    chipPars.lambdaSigBase = 1;
    chipPars.lampPower = SNR.*(SNR+sqrt(SNR.^2+4*chipPars.lambdaBg))/2;
    chipPars.numFrames = numFrames; 

end

if nargin < 7
    outFold = 'data';
end
    
% filenamesAll = cell(1,length(gainF));
filenames = cell(1,length(zooms));
data = cell(1,length(zooms));

for jj=1:length(zooms)
    chipPars.zoom = zooms(jj); % zoom
    chipPars.particleDensity = particle_density(jj); % how many particles





    % fluorophore positions are fixed and the same for all experiments
    circRadius = chipPars.maxRadius*chipPars.zoom/100;          % circle radius (in pixels)  

    nRows = chipPars.nRows;             % number of rows in image
    nCols = chipPars.nCols;              % number of columns in image
    particleDensity = chipPars.particleDensity;   % signal particle density

%     if gainF == 1
%         gtImage = ones(nRows,nCols);
%     else
    if nargin < 9
        [gtImage, placements] = fixed_beads(particleDensity ,nRows, nCols, circRadius);
    end
%     end
    for j=1:length(gainF)
        chipPars.gain = gainF(j); % different gain settings
        chipPars.gainName = gainName(j);
        [filenames{jj}{j},data{jj}{j} ]= simulate_random_beads_images(chipPars, gtImage, placements, outFold);
    end
end

end
 
function [filenames,data] = simulate_random_beads_images(chipPars, gtImage, placements, outFold)

    %
    % Script for generating an image with randomly positioned
    % particles (point emitters).
    % 
    % Dependencies: generate_image_randomly_deposited_particles.m
    %
   
    % lamp intensity: how much photons / instead of SNR vals
	lampPower = chipPars.lampPower;
    
    lambdaBg = chipPars.lambdaBg;          % Poisson parameter for background
  
    % estimates for parameters from Mean Variance (MV) callibration
    gain = chipPars.gain;
    adFactor = chipPars.adFactor;
    countOffset = chipPars.countOffset;
    roNoise = chipPars.roNoise;
    lambdaSigBase = chipPars.lambdaSigBase;
    
    nRows = chipPars.nRows;             % number of rows in image
    nCols = chipPars.nCols;              % number of columns in image
    
    numFrames = chipPars.numFrames;
    
    if ~isempty(outFold)
        subFold = strcat(num2str(chipPars.zoom),'x');
        [~,~] = mkdir(outFold,subFold);
    end
    
    filenames = cell(1,length(lampPower));
    for idxSNR = 1:length(lampPower)
        
        lambdaSig = lambdaSigBase*lampPower(idxSNR);
        fprintf('Intensity: %.3f\n',lambdaSig);
        
        % photons image
        photonsPlaced = gtImage.*lambdaSig + lambdaBg;
        
        if ~isempty(outFold)
            filename = fullfile(outFold,subFold,strcat(['synth' num2str(chipPars.zoom) 'x_gain' num2str(chipPars.gainName) '_lamp' num2str(lambdaSig) '.tif'] ));
            name = strcat(['synth' num2str(chipPars.zoom) 'x_gain' num2str(chipPars.gainName) '_lamp' num2str(lambdaSig)]);
            delete(filename);
            filenames{idxSNR} = filename;
        end

        for jj=1:numFrames
            noisyImage = noise_model_emccd(photonsPlaced(:)', gain, adFactor, countOffset, roNoise);
        
            finalImage = reshape(noisyImage,nRows,nCols);
            if ~isempty(outFold)
                imwrite(uint16(finalImage),filename,'WriteMode','append');
            end
        end
        % Translate SNR to lambda for signal regions.
        %  Comment: 
        % The signal-to-noise ratio is defined:
        %       SNR = lambdaSig /sqrt(lambdaSig + lambdaBg)
        % 
        % Inverting these expressions we get: 
        %      lambdaSig = SNR^2/2 + sqrt(SNR^2*lambdaBg + SNR^4/4)
        %
        lambdaSig = lampPower(idxSNR);
        SNR = lambdaSig /sqrt(lambdaSig + lambdaBg);
%         lambdaSignal = SNR*(SNR+sqrt(SNR^2+4*lambdaBg))/2;
        


    % save image as tif file
        % Store image and associated information
        data.image = noisyImage;
        data.groundTruthImage = gtImage;
        data.lambdabg = lambdaBg;
        data.lambdasig = lambdaSig;
        data.gain = gain;
        data.roNoise = roNoise;
        data.offset = countOffset;
        data.adFactor = adFactor;
        data.placements = placements;
        data.snr = SNR;

        if ~isempty(outFold)
            data.imageName = strcat(name,'.mat');
            %         titleSNR = round(100*SNR);
            %         saveFileName = sprintf('testImageEMCCDBeadsSNR%i',titleSNR);
            save(fullfile(outFold,subFold,data.imageName),'data');
        end
    end
    

end


% noise image emccd
function noisyImage = noise_model_emccd(photonsPlaced, gain, adFactor, countOffset, roNoise)
    %
    % from Super-resolution fight club: assessment  of 2D and 3D single-molecule 
    %localization microscopy software
    %     gain = chipPars.gain;
    %     roNoise = chipPars.roNoise;
    %     adFactor = chipPars.adFactor;
    %     offset = chipPars.offset;
    %     c = chipPars.c; 
    %     QE = chipPars.QE;
 
    % Poisson value
    poissVals = poissrnd(photonsPlaced);
    
    % Generate amplified electric signal
    if gain == 1
        elecImage = poissVals;
    else
        elecImage = gamrnd(poissVals,gain);
    end
    % Perform readout on each pixel     
    noisyImage = (roNoise) * randn(1,length(photonsPlaced)) + countOffset + elecImage/adFactor;
    
    % round (uniformly distributed noise)
    noisyImage = round(noisyImage);

end


function [gtImage, placements] = fixed_beads(particleDensity ,nRows, nCols, circRadius)
    
    %   Args:
    %       particleDensity - particular density in the image
    %       nRows - number of rows in image
    %       nCols - number of columns in the image
    %       circRadius - radius of the image
    %
    %   Returns:
    %       gtImage - ground truth image without lambdaSig and lambdaBg
    %       placements - locations of centers of beads

    % Generate black and white signal image (randomly deposited circles)
    nCircles = round(particleDensity*nRows*nCols);
%     signalBwImage = zeros(nRows,nCols);
    iMin = ceil(circRadius) + 1;
    iMax = nRows - ceil(circRadius);
    circleCenterRows = randi([ iMin, iMax],1,nCircles);
    circleCenterCols = randi([ iMin ,iMax],1,nCircles);
    gtImage = zeros(nRows,nCols);
    placements = [circleCenterRows' , circleCenterCols'];

    % placement of circles on the image
    for circleIdx = 1:nCircles
        for row = 1:nRows
            for col = 1:nCols
                if ( (row - circleCenterRows(circleIdx))^2 + ...
                        (col - circleCenterCols(circleIdx))^2) < circRadius^2
%                     signalBwImage(row,col) = lambdaSig;
                    gtImage(row,col) = gtImage(row,col)+1; % still might be two circles at the same place if image is dense
                end
            end
        end
    end
%    signalBwImage = signalBwImage + lambdaBg;
%     
end