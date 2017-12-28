function [x,res] = LM(x0,t,y,iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Levenberg-Marquardt Algorithm
%   Input:
%       x0  : initial guess - parameters, in column-form
%       t   : t_i - observed points
%       y   : y_i - response data
%       iter: number of iterations
%   Required:
%       Jacobian function: 
%       [J,r] = jacobian(x,t,y)
%   Output:
%       x   : fitted parameters
%       res : residual
%   
% Written by Mai, Jiaqi, on 12/28/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = length(x0);
    x = x0;
    I = eye(m);
%     figure
    for k = 1:iter
        [J,r] = jacobian_SL(x,t,y,'decay');
        smu = k^(-1);%norm(r) %;exp(-k);
        
        curr_res = norm(r)
        
%         L = [J;smu*I];
%         r0 = [r;zeros(m,1)];
%         [Q,R] = qr(L,0);
%         s = linsolve(R,-Q'*r0)';
        L = [J ; smu*I];
        r0 = [r ; zeros(m,1)];
        s = - L \ r0;
        x  = x+s;
        
%         flag = -1;
%         sgm = x(1);
%         k = x(2);
%         t0 = x(3);
%         f = sqrt(sgm./(k+flag*exp(-2*sgm*(t-t0))));
%         plot(t,f); pause(0.01)
    end
    res = r;
end
