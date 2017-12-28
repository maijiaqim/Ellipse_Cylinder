%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear hil;

    %
    char = 'decay';
    flag = -1;
    % Set region that apply Hilbert Transformation to extract envelope
    % amplitude.
    hilbert_region = 250000:length(time);
    t_hil = time(hilbert_region);
    [hil,temp] = envelope(dragy(hilbert_region),10000,'peak');  
    % (Alternative) Hilbert Transformation: = abs(hilbert(decay(hilbert_region,i)));
    
% Check the envelope amplitude    
%   for i = 1:col
%       figure
%       plot(t,hil);hold on
%       plot(time,dragy);
%   end

    % Set region that apply LM algorithm for nonlinear fitting.
    fit_region = 30000:(length(t_hil)-35000);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       linear fitting to determine
%       initial guess:
%           1. sigma  
%           2. saturated value CL_square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_fit = t_hil(fit_region);
    hil_fit = hil(fit_region);
    
    S = (hil((fit_region+1)) - hil((fit_region-1)))./(t_hil(fit_region+1)-t_hil(fit_region-1))./hil(fit_region);
    R = hil(fit_region).^2;
        % figure; plot(R,S);

    P = polyfit(R(100000:end),S(100000:end),1);
    sigma = P(2);
    k = -P(1);
    CL_square = -P(2)/P(1);
    t0 = log(flag*(sigma/hil_fit(1)^2-k))/(2*sigma)+t_fit(1);
    f0 = sqrt(sigma./(k+flag*exp(-2*sigma*(t_fit-t0))));

    x0 = [sigma;k;t0];

    [J,r] = jacobian_SL(x0,t_fit,hil_fit,char);
    res0 = norm(r);
        % figure
        % plot(t_fit,f);hold on
        % plot(t_fit,hil_fit)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Levenberg–Marquardt Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [x,r] = LM(x0,t_fit,hil_fit,300);
    [x,r] = LM(x,t_fit,hil_fit,300);
%% 
    [x,r] = LM(x,t_fit,hil_fit,300);
    sgm = x(1);
    k = x(2);
    t0 = x(3);
    f = sqrt(sgm./(k+flag*exp(-2*sgm*(t_fit-t0))));
    reltv_chg = (x - x0)./x0
    res = norm(r);
    res_reltv_chg = (res-res0)/res0
%% 
    
% Check plot
    figure
    plot(t_fit,f0);hold on
    plot(t_fit,hil_fit);%hold on
        % plot(time(hilbert_region),dragy(hilbert_region))
