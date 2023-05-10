function [lambdaBg, pci,L, U,binEdges,binPos,sortTruncI] = est_lambda(sortI,Nthresh,gain,adFactor,offset,roNoise,r)
    opt = statset('MaxIter',200,'MaxFunEvals',400,'FunValCheck','on');

    m = find(sortI>Nthresh,1,'First')-1; % find first element greater than Nthresh
    if isempty(m)
        m = length(sortI);
    end
    sortTruncI = sortI(1:m);
    % Fit lambda
    lamGuess = abs((sortTruncI(round(m/2)) - offset)/r);
    
    % Prep data for truncated fit
    import Core.calc_bounds;
    [L, U, EX, STD] = calc_bounds(lamGuess,gain,adFactor,offset,roNoise);
    % Get bin edges
    binEdges = [max(1,ceil(L)):1:floor(U)]-0.5; % bin edges shifted by half  
    % histogram for intensities we can have this outside of the function
    histAll = histcounts(sortTruncI,binEdges)';
            
    binPos = binEdges(1:end-1) + diff(binEdges)/2;
    import Core.log_likelihood_trunc_dist; 
    logL = @(lambda,data,cens,freq,trunc) log_likelihood_trunc_dist(lambda,data,cens,freq,trunc, ...
                       gain, adFactor, offset, roNoise);
            
   
%     logLValues = @(x) logL(x,data,cens,freq,trunc);
    
    % use mle with negative loglikelihood (nloglf)
%     tic
    [lambdaBg, pci] = mle(binPos,'nloglf',logL,'start',lamGuess,'lowerBound',0,'Frequency',histAll,'TruncationBounds',[L U],'Options',opt);
%     toc
    

end

