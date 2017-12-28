function [J,r] = jacobian_SL(x,t,y,char)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Jacobian of r_i = y_i - A(t_i),
%   where A(t) is the solution of Stuart-Landau Equation
%       1  dA                
%      --- --- = sigma - k A^2 
%       A  dt           
%                 /            sigma              \
%    A(t) = Sqrt | ------------------------------- |
%                 \   k +(-) exp[-2 sigma (t-t0)] /
%    grow : dA/dt > 0, '+'
%    decay: dA/dt < 0, '-'
%
%   Input:
%       x   : parameters
%               x(1) = sigma
%               x(2) = k
%               x(3) = t0
%       t   : observed points
%       y   : response data
%       char: switches
%           'grow'  => dA/dt > 0
%           'decay' => dA/dt < 0
%   Output:
%       J   : Jacobian of residual function
%       r   : residuals
%
% Written by Mai, Jiaqi, on 12/28/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = length(t);
    m = length(x);
    sigma = x(1);
    k = x(2);
    t0 = x(3);
    J = zeros(n,m);
    r = zeros(n,1);
    flag = -1;
    switch char
        case 'grow'
            flag = 1;
        case 'decay'
            flag = -1;
        otherwise
            flag = -1;
    end
    f = sqrt(sigma./(k+flag*exp(-2*sigma*(t-t0))));
    i = 1:n;
        J(i,1) = 0.5.*(k+(1+2*(t-t0))*flag.*exp(-2*sigma*(t-t0)))./(k+flag*exp(-2*sigma*(t-t0))).^2./f;
        J(i,2) = -0.5*sigma./(k+flag*exp(-2*sigma*(t-t0))).^2./f;
        J(i,3) = -0.5*sigma.*(2*sigma*flag*exp(-2*sigma*(t-t0)))./(k+flag*exp(-2*sigma*(t-t0))).^2./f;
     r = y-f;
     J = -J;
end