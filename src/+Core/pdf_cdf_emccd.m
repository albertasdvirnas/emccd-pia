function [pdfEmccd,cdfEmccd] = pdf_cdf_emccd(intensities,lambda,gain,adFactor,offset,roNoise,L,U)

    % Generates EMCCD probability density function (PDF) and 
    % cumulative distribution functin (CDF) by numerical inversion 
    % of the characteristic function.
    %
    % Input:
    % 
    % intensities = vector (or matrix) with intensity values
    % lambda = Poisson parameter
    % chipPars = struct containing the chip parameters
    % N = number of integration points.
    %
    % Output:
    % 
    % pdfEmccd = EMCCD probability density function 
    % cdfEmccd = cumulative distribution function 
    %
    % Refs: V. Witkovský, "Numerical inversion of a characteristic function: 
    % An alternative tool to form the probability distribution of 
    % output quantity in linear measurement models.", 
    % Acta IMEKO 5.3 (2016): 32-44, see Eqs. (8) and (9).
    %
    %
    

    import Core.calc_bounds;

    if nargin < 7
        [L, U, EX, STD] = calc_bounds(lambda, gain, adFactor, offset, roNoise);
    end
      
     % Extract chip parameters
%     gain = chipPars.gain;
%     adFactor = chipPars.adFactor;
%     offset = chipPars.countOffset;
%     roNoise = chipPars.roNoise;
    r = gain/adFactor;
      
%     % Analytic expressions for the mean and variance
%     EX = lambda*gain/adFactor+offset; 
%     STD = sqrt(roNoise^2 + 2*lambda*r^2 + 1/12);  
%     
%     % Analytic expression for the characteristic function 
%     % for the EMCCD distribution
% %     cfAnaly = @(t) exp(-t.^2*roNoise^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*offset)*2*sin(t/2)/t;
% 
%     %
%    % limits where pdf is nonzero
%     numstds = 6;
%     L = EX-numstds*STD; % mean - 6 std
%     U = EX+numstds*STD;
%     U = min(max(intensities),U); % limit to U for truncated case


	% optimal value for step parameter
    dt = 2*pi/(U-L);

    % For discrete, integral is -pi..pi, because the output variable is
    % discretized
    N = pi/dt;
    
    
    % Estimate step size, dt, for numerical integration
    t = (1:1:N)' * dt;  

    cf = char_fun(t , roNoise,lambda,r,offset);

        % y is the grid for our pdf (from L to U)
    y = intensities;
    
    
    % calculate main integral
    pdfEmccd = trapezoidal_pdf(y,dt,t,cf);
    cdfEmccd = cumsum(pdfEmccd);

%     cdfEmccd = trapezoidal_cdf(y,dt,t,cf,EX);

   
    
end


function cfCombined = char_fun(t , roNoise,lambda,r,offset)
%

    cfAnaly = exp(-t.^2*roNoise^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*offset);
    cfROUND = 2*sin(t/2)./t;
    cfROUND(t==0) = 1;


    
    cfCombined = cfAnaly.*cfROUND;

end

%
function pdf = trapezoidal_pdf(y,dt,t,cf)
    w = ones(length(t),1);
    w(end) = 1/2; % last coef is 1/2
       
    pdf = dt/pi*(1/2 +cos(t*y)'*(real(cf).*w)+sin(t*y)'*(imag(cf).*w));
end

% Not needed to calculate since we have PMF
% function cdf = trapezoidal_cdf(y,dt,t,cf,ex)
%     w = ones(length(t),1);
%     w(end)=1/2; % last coef is 1/2
%     cdf = 1/2 - dt/pi*(1/2*(ex-y') +cos(t*y)'*(imag(cf./t).*w)-sin(t*y)'*(real(cf./t).*w));
% end

