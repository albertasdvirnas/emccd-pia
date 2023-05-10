   
minVal = quantile(sortI,0.25);
maxVal = quantile(sortI,0.9);
nvvals = minVal:maxVal;

% nvvals = 35:60;
lambdaBg = [];
for Nthresh=nvvals
    m = find(sortI>Nthresh,1,'First')-1; % find first element greater than Nthresh
    sortTruncI = sortI(1:m);
    % Fit lambda
    lamGuess = abs((sortTruncI(round(m/2)) - countOffset)/r);
    
    %         tic
    [lbg, pci] = mle(sortTruncI,'logpdf',logL,'start',lamGuess,'lowerBound',0,'Options',opt);
    lambdaBg =  [lambdaBg lbg];
end

figure,plot(nvvals,lambdaBg)