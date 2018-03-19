function [param,yfit,resnorm,xx0] = SL_equation_fit(time,ydata,char)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Nonlinear fitting to Stuart-Landau Equation
%
%
%   INPUT:
%       time    : observed points
%       ydata   : response data
%       char    : switches
%           'grow'  => dA/dt > 0
%           'decay' => dA/dt < 0
%   OUTPUT:
%       param   : fitted parameters (6) -> sr si lr li t0 phi
%       yfit    : fitted ydata
%       resnorm : norm of residual
%       xx0     : initial guess (5) -> sr_0 si_0 lr_0 t0_0 phi_0
%
%   DEPENDENCY:
%       findfrequency.m
%       SL_cos.m
%       SL_analytical.m
%       
% Written by Mai, Jiaqi, on 03/19/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch char
        case 'grow'
            flag = 1;
        case 'decay'
            flag = -1;
        otherwise
            flag = -1;
    end
    
    SL_pre = @(x,t)SL_cos(x,t,char);
    SL_corr = @(x,t)SL_analytical(x,t,char);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Generate initial Guesses
    %       1. \omega, \phi
    %       2. \sigma, k, t0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1. Initial guess of \omega. 
    %   \phi is set arbitrarily.
    [amp,freq] = findfrequency(time,ydata);
    [M,ind] = max(amp);
    omega = 2*pi*freq(ind);
    phi = -1;
    
    % 2. Initial guess of \sigma, k and t0
    A = envelope(ydata,round(length(ydata)/100),'peak');
    S = (A(3:end)-A(1:end-2))./(time(3:end)-time(1:end-2))./A(2:end-1);
    R = A(2:end-1).^2;
    P = polyfit(R,S,1);
    sigma = P(2);
    k = -P(1);
    t0 = log(flag*(sigma/A(1)^2-k))/(2*sigma)+time(1);

    x0 = [sigma;k;t0;omega;phi];
    
    % Adjust initial guess of \omega and \phi by fitting data to SL_cos.
    options = optimoptions('lsqcurvefit','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',1e5,'MaxIter',50,'DiffMinChange',1e-10,'TolPCG',1e-5,'FinDiffType','central','Display','iter');
    lb = x0;
    ub = x0;
    lb(4:5) = -Inf;
    ub(4:5) = Inf;
    x = lsqcurvefit(SL_pre,x0,time,ydata,lb,ub,options);
    
    % 3. L is set (reasonably) arbitrarily.
    L = 0.04;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Nonlinear Fitting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xx0 = [x;L];
    lb = [];
    ub = [];
    [xx,rr] = lsqcurvefit(SL_corr,xx0,time,ydata,lb,ub,options);
    
    yfit = SL_corr(xx,time);
   
    param = xx;
    resnorm = rr;
end