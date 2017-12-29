function [x,res] = LM(x0,t,y,model_func,iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Levenberg-Marquardt Algorithm
%   Input:
%       x0  : initial guess - parameters, in column-form
%       t   : t_i - observed points
%       y   : y_i - response data
%       model_func: function handle of model function f(x,t)
%       iter: number of iterations
%   Output:
%       x   : fitted parameters
%       res : residual
%   
%   Required: 
%       jacobian.m
%
% Written by Mai, Jiaqi, on 12/29/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = length(x0);
    x = x0;
    I = eye(m);
    for k = 1:iter
        [J,r] = jacobian(x,t,y,model_func);
        smu = k^(-0.5);%%norm(r) %;exp(-k);
        
        curr_res = norm(r)
        
        L = [J ; smu*I];
        r0 = [r ; zeros(m,1)];
        s = - L \ r0;
        
%         L = [J;smu*I];
%         r0 = [r;zeros(m,1)];
%         [Q,R] = qr(L,0);
%         s = linsolve(R,-Q'*r0);

        x  = x+s;
        
    end
    res = r;
end
