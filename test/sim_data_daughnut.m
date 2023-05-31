
N = 256;
circleCenterRows = N/2;
circleCenterCols = N/2;
circRadius = 50;
gt = zeros(N,N);
for row = 1:N
    for col = 1:N
        if ( (row - circleCenterRows(1))^2 + ...
                (col - circleCenterCols(1))^2) < circRadius^2
%                     signalBwImage(row,col) = lambdaSig;
            gt(row,col) = gt(row,col)+1; % still might be two circles at the same place if image is dense
        end
    end
end

circRadius = 25;
% gt = zeros(N,N);
for row = 1:N
    for col = 1:N
        if ( (row - circleCenterRows(1))^2 + ...
                (col - circleCenterCols(1))^2) < circRadius^2
%                     signalBwImage(row,col) = lambdaSig;
            gt(row,col) = 0; % still might be two circles at the same place if image is dense
        end
    end
end


SNRVals =[3 4];  
filenames = simulate_random_beads_full(100, 20, 100, 1/1000,SNRVals, 1, 'synthSingleFrameNew',[N N],gt,[]); % zoom 100x


chipParsCur = struct();
chipParsCur.gain = mean(chipParsS.gainQ(3,2)); % based on gain in the image.. 50 100 300..
chipParsCur.adFactor = mean(chipParsS.aduQ(2));
chipParsCur.countOffset = mean(chipParsS.deltaQ(3,2));
chipParsCur.roNoise = mean(chipParsS.sigmaQ(3,2));

chipParsCur.pval = 0.01;



chipParsFig4 = {chipParsCur,chipParsCur}; % the two experiments will have slightly varying parameters because of different Gain setting

outFig4 = { fullfile(figFold,'FigDonut.eps')};
% outFig3 = 'C:\Users\Lenovo\postdoc\PAPERS\emccd-paper\draft\Figs\Fig4a.eps';

pValThresh = 0.01; % for chi2
pValThreshBinarization = 0.1; % for binarization

fig4_segmentation(filenames{1}{1},chipParsFig4,outFig4,pValThresh,pValThreshBinarization);

