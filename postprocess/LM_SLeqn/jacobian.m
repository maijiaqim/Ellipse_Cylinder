function [J,r] = jacobian(x,t,y,model_func)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Jacobian of r_i = y_i - A(t_i), required for LM algorithm
%
%   Input:
%       x   : parameters
%       t   : observed points
%       y   : response data
%       model_func: function handle of model function f(x,t)
%   Output:
%       J   : Jacobian of residual function
%       r   : residuals
%
%   Written by Mai, Jiaqi, on 12/29/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = length(t);
    m = length(x);
    J = zeros(n,m);
    r = zeros(n,1);
    
    f = model_func(x,t);        
    
    %
    % Compute Jacobian (first derivative) -- Complex number trick
    %             1
    %   f'(x) ~ ----- Im[ f(x + i eps) ]
    %            eps
    %
    h = sqrt(eps);
    ii = sqrt(-1);
    for j = 1:m
        xh = x;
        xh(j) = x(j)*(1+ii*h);
        fh = model_func(xh,t);
        J(:,j) = imag(fh)/x(j)/h;
    end
    r = y-f;
    J = -J;
end
