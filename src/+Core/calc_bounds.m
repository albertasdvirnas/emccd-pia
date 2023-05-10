function [L,U,EX,STD] = calc_bounds(lambda,gain,adFactor,offset,roNoise)
    % reference
    r = gain/adFactor;
      
%     % Analytic expressions for the mean and variance
    EX = lambda*gain/adFactor+offset; 
% 
% if nargin < 7
    STD = sqrt(roNoise^2 + 2*lambda*r^2 + 1/12);  
    numstds = 6;
    L = EX-numstds*STD; % mean - 6 std
    U = EX+numstds*STD;
end
% end

