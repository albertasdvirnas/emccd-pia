function char_fun_test(lamGuess, gainGuess, adFactorGuess, offsetGuess, roNoiseGuess )
if nargin < 1
    offsetGuess = 100;
    roNoiseGuess = 1.44;
    adFactorGuess = 20;
    gainGuess = 100;
    lamGuess = 10;
end
% Tests how CDF/PDF are calculated for given parameters

% useful papers:
% https://www.tandfonline.com/doi/pdf/10.1080/00949655.2021.1982942

% tvals = -10:0.001:10;
% [phi] = char_fun_to_dist_test(tvals);
% figure,plot(tvals,imag(phi))
%
%  Nr = char_fun_uniform(tvals);
% calc_gaussian()

% calc_poisson()
% %% EXAMPLE UNIFORM
% calc_uniform()
% %
% calc_logistic()
%
% calc_combined()

% parameter chouces


run_main(lamGuess, gainGuess, adFactorGuess, offsetGuess, roNoiseGuess)

end

function []=run_main(lamGuess, gainGuess, adFactorGuess, offsetGuess, roNoiseGuess)

    
    % Prep data for truncated fit
    import Core.calc_bounds;
    [L, U, EX, STD] = calc_bounds(lamGuess, gainGuess, adFactorGuess, offsetGuess, roNoiseGuess);

    r = gainGuess/adFactorGuess;
  
	% optimal value for step parameter
    dt = 2*pi/(U-L);

    % For discrete, integral is -pi..pi, because the output variable is
    % discretized
    N = pi/dt;
    
    
    % Estimate step size, dt, for numerical integration
    t = (1:1:N)' * dt;  

    cf = char_fun_full(t , roNoiseGuess,lamGuess,r,offsetGuess);

        % y is the grid for our pdf (from L to U)
    y = ceil(L):floor(U);
    y = 304:floor(U);

    
    % calculate main integral
    pdfEmccd = trapezoidal_pdf(y,dt,t,cf);
    cdfEmccd = cumsum(pdfEmccd);


        figure,
    tiledlayout(2,1)
    nexttile
    plot(y,pdfEmccd)
    hold on
%     cdf = trapezoidal_cdf(y,dt,t,cf,EX);
    nexttile
    hold on
    plot(y,cdfEmccd)
    title('CDF')
  



end

function cfCombined = char_fun_full(t , roNoise,lambda,r,offset)
%

    cfAnaly = exp(-t.^2*roNoise^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*offset);
    cfROUND = 2*sin(t/2)./t;
    cfROUND(t==0) = 1;


    
    cfCombined = cfAnaly.*cfROUND;

end



 function []=calc_uniform()

    EX=0;
    numstds = 6;
    STD = sqrt(1/12);
    
    L = EX-numstds*STD; % mean - 6 std
    U = EX+numstds*STD;

%     sF = 0.01;
%     N = (U-L)/sF;
% %     %      T = 1000;
% %     %      dt = T/N; % alt
%     dt = pi/(U-L)*sF;
%     
%      N = 2^12;
% %      T = 1000;
% %      dt = T/N; % alt
%      dt = 2*pi/(U-L);
%      t = (1:N)' * dt;  
%      cf = char_fun_uniform(t);
     
         N = 2^22;
        dt = 2*pi/(U-L);

    t = (1:1:N)' * dt;  
    cf = char_fun_uniform(t);
%     eps = 10^(-8);
%     val = abs(cf./t);
    % https://arxiv.org/pdf/1701.08299.pdf

%     potStop = find(val<eps,1,'first');
%     t = (1:1:potStop)' * dt;  
%     cf = char_fun_uniform(t);

    
    % figure,plot(t,real(cf))
%     h = 0.000001;
%     mEst = 1/(12*1i*h)*(char_fun_uniform(-2*h)-8*char_fun_uniform(-h)+8*char_fun_uniform(-h)-8*char_fun_uniform(2*h));
    figure,plot(t,imag(exp(-1i*t*0.2).*cf))

    y = -0.4:0.01:0.4;
    % y  = L:0.01:U;
    pdf=trapezoidal_pdf(y,dt,t,cf);
    figure,
    tiledlayout(3,1)
    nexttile
    plot(y,pdf)
    hold on
    trueC = zeros(length(y),1);
    trueC(logical((y>=-0.5)+(y<=0.5)-1))=1;
    plot(y,trueC)
    title('PDF')
    legend({'Fitted','True'},'location','eastoutside')
    nexttile
    plot(y,sqrt((pdf-trueC).^2))
    xlabel('RMSE')

    cdf = trapezoidal_cdf(y,dt,t,cf,EX);
    nexttile
    hold on
    plot(y,cdf)
    title('CDF')
 end
% % 

% 
 %%
 
 function []=calc_gaussian()
 %% EXAMPLE UNIFORM
EX=0;
numstds = 6;
STD = sqrt(1);
L = EX-numstds*STD; % mean - 6 std
U = EX+numstds*STD;
 
 N = 2;
%  T = 1000;
%  dt = T/N; % alt
 dt = 2*pi/(U-L);
 t = (1:N)' * dt;  
 cf = char_fun_gaussian(t);
% figure,plot(t,real(cf))

y = -5:0.1:5;
% y  = L:0.01:U;
pdf=trapezoidal_pdf(y,dt,t,cf);
figure,
tiledlayout(3,1)
nexttile
plot(y,pdf)
hold on
trueC  = 1/sqrt(2*pi)*exp(-1/2*y.^2);
trueC = trueC';
plot(y,trueC)
title('PDF')
legend({'Fitted','True'},'location','eastoutside')
nexttile
plot(y,sqrt((pdf-trueC).^2))
xlabel('RMSE')

cdf = trapezoidal_cdf(y,dt,t,cf,EX);
nexttile
hold on
plot(y,cdf)
title('CDF')

 
 end
 
 
 
function calc_logistic()
    s = 50;
    EX = 0;
    numstds = 6;
    STD = s*pi/sqrt(3);
    L = EX-numstds*STD; % mean - 6 std
    U = EX+numstds*STD;
    
    dt = 2*pi/(U-L);

    % what should N be?
    N = 2^16;
    
    t = (1:1:N)' * dt;  
%     cy = char_fun_logistic(t,s)
    cf = char_fun_logistic(t,s);
    eps = 10^(-12);
    val = abs(cf./t);
    % https://arxiv.org/pdf/1701.08299.pdf

    potStop = find(val<eps,1,'first');
    t = (1:1:potStop)' * dt;  
%     cf = char_fun_poisson(t,g,lambdabg);

% val = abs(imag(exp(-1i*t*EX).*cf./t))

%     sF = 1;
% 	N = (U-L)/sF;
%      T = 1000;
%      dt = T/N; % alt
%     dt = pi/(U-L)*sF;
%     t = (1:1:N)' * dt;  
    cf = char_fun_logistic(t,s);
%     sum(imag(exp(-1i*t*3).*cf))
%     figure,plot(t,real(exp(-1i*t*30).*cf))
    % 
%     movsum(imag(exp(-1i*t*3).*cf),100)
%     val = imag(exp(-1i*t*3).*cf);
%     locs=find(val<min(val)+0.01);
    
    y = round(L:1:U);
    % y  = L:0.01:U;
    pdf=trapezoidal_pdf(y,dt,t,cf);
    figure,
    tiledlayout(3,1)
    nexttile
    plot(y,pdf)
    hold on
    trueC = exp(-y/s)./(s*(1+exp(-y/s)).^2);
    trueC = trueC';
    % trueC = zeros(length(y),1);
    % trueC(logical((y>=-0.5)+(y<=0.5)-1))=1;
    plot(y,trueC)
    title('PDF')
    legend({'Fitted','True'},'location','eastoutside')
    nexttile
    plot(y,sqrt((pdf-trueC).^2))
    xlabel('RMSE')

    cdf = trapezoidal_cdf(y,dt,t,cf,EX);
    nexttile
    hold on
    plot(y,cdf)
    title('CDF')
  
end
 
 
%  end
 
function calc_poisson()
    g = 1;
    lambdabg = 20;

    EX = lambdabg;
    numstds = 6;
    STD = sqrt(lambdabg);
    L = round(max(1,EX-numstds*STD)); % mean - 6 std
    U = round(EX+numstds*STD);
    
%     Nfixed = 256;
    
    dt = 1*2*pi/(U-L);

    % For discrete, N*dt = pi, so N = 2/(U-L)
    N = pi/dt;
    
    t = (1:1:N)' * dt;  
    cf = char_fun_poisson(t,g,lambdabg);
%     eps = 10^(-15);
%     val = abs(cf./t);
    % https://arxiv.org/pdf/1701.08299.pdf

%     potStop = find(val<eps,1,'first');
%     t = (1:1:potStop)' * dt;  
%     cf = char_fun_poisson(t,g,lambdabg);

% val = abs(imag(exp(-1i*t*EX).*cf./t))

%     sF = 1;
% 	N = (U-L)/sF;
%      T = 1000;
%      dt = T/N; % alt
%     dt = pi/(U-L)*sF;
%     t = (1:1:N)' * dt;  
%     cf = char_fun_poisson(t,g,lambdabg);
%     sum(imag(exp(-1i*t*3).*cf))
%     figure,plot(t,real(exp(-1i*t*30).*cf))
    % 
%     movsum(imag(exp(-1i*t*3).*cf),100)
%     val = imag(exp(-1i*t*3).*cf);
%     locs=find(val<min(val)+0.01);
    
    y = round(L:1:U);
    % y  = L:0.01:U;
    pdf=trapezoidal_pdf(y,dt,t,cf);
    figure,
    tiledlayout(3,1)
    nexttile
    plot(y,pdf)
    hold on
    trueC =g*lambdabg.^y*exp(-lambdabg)./factorial(y);
    trueC = trueC';
    % trueC = zeros(length(y),1);
    % trueC(logical((y>=-0.5)+(y<=0.5)-1))=1;
    plot(y,trueC)
    title('PDF')
    legend({'Fitted','True'},'location','eastoutside')
    nexttile
    plot(y,sqrt((pdf-trueC).^2))
    xlabel('RMSE')

    cdf = trapezoidal_cdf(y,dt,t,cf,EX);
    nexttile
    hold on
    plot(y,cdf)
    title('CDF')
  
end
 
function calc_combined()
    % parameters
    s = 20; % scale
    g = 100; % gain
    lambdabg = 4; % lambda bg
    o = 364; % offset.
%     eps = 10^(-12); % epsilon, small number to stop calculation at.

    % find analytic mean and standard deviation
    EX = g*lambdabg+o;
    STD = sqrt(g^2*lambdabg+s^2*pi^2/3+1/12);

    % limits where pdf is nonzero
    numstds = 6;
    L = EX-numstds*STD; % mean - 6 std
    U = EX+numstds*STD;
    
    % optimal value for step parameter
    dt = 1*2*pi/(U-L);

    % For discrete, N*dt = pi, so N = 2/(U-L)
    N = pi/dt;
    
    
    % https://arxiv.org/pdf/1701.08299.pdf
    % estimate N, based on where  abs(cf./t) starts being smaller than eps
    t = (1:1:N)' * dt;  
%     cf = char_fun(t, g, o, lambdabg, s);
%     val = abs(cf./t);
% 
%     potStop = find(val < eps,1,'first'); % assume val to be decreasing/periodic
%     
%     % define steps
%     t = (1:1:potStop)' * dt;  
    cf = char_fun(t, g, o, lambdabg, s);

   
    % y is the grid for our pdf (from L to U)
    y = round(L:1:U);
    
    
    % calculate main integral
    pdf=trapezoidal_pdf(y,dt,t,cf);
    %
    
    %
    figure,
    tiledlayout(2,1)
    nexttile
    plot(y,pdf)
    hold on
    cdf = trapezoidal_cdf(y,dt,t,cf,EX);
    nexttile
    hold on
    plot(y,cdf)
    title('CDF')
  
end

function [phi] = char_fun_to_dist_test(t)

    g = 20; % gain
    o = 100; % ofset
    lambda_bg = 10; % lambda bg
    s = 4;      % scale
    
%     t=0:0.1:10;
    phi=char_fun(t , g, o, lambda_bg, s);
    
end


function cfCombined = char_fun(t , g, o, lambda_bg, s)
%
cfGX =  exp(lambda_bg*(exp(1i*g*t)-1));

% t = 0:0.01:100;
cfNREAD = pi*s*t./(sinh(pi*s*t));

 cfROUND = 2*sin(t/2)./t;
 cfROUND(t==0) = 1;
% if t~=0
%    
% else
%     Nr = 1;
% end
cfO = exp(1i*t*o);

cfCombined = cfGX.*cfNREAD.*cfROUND.*cfO;

% how to calculate very small values?
%phi =  exp(lambda_bg*(exp(1i*g*t-1)))*pi*s*t./(sinh(pi*s*t))* 2*sin(t/2)./t*exp(1i*t*o);

end


function Nr = char_fun_uniform(t)
 Nr = 2*sin(t/2)./t;
 Nr(t==0) = 1;
end
%
function Nr = char_fun_gaussian(t)
 Nr = exp(-1/2*t.^2);
end
%
function gX = char_fun_poisson(t,g,lambda_bg)
 gX = exp(lambda_bg*(exp(1i*g*t)-1));
end
%
function cy = char_fun_logistic(t,s)
 cy = pi*s*t./(sinh(pi*s*t));

end
%
%
function pdf = trapezoidal_pdf(y,dt,t,cf)
    w = ones(length(t),1);
    w(end) = 1/2; % last coef is 1/2
       
    pdf = dt/pi*(1/2 +cos(t*y)'*(real(cf).*w)+sin(t*y)'*(imag(cf).*w));
end
%
function cdf = trapezoidal_cdf(y,dt,t,cf,ex)
    w = ones(length(t),1);
    w(end)=1/2; % last coef is 1/2
    cdf = 1/2 - dt/pi*(1/2*(ex-y') +cos(t*y)'*(imag(cf./t).*w)-sin(t*y)'*(real(cf./t).*w));
end

