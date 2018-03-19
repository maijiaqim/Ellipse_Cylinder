function output = SL_analytical(x,t,char)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The solution of Stuart-Landau Equation
%       1  dA                
%      --- --- = sigma - l |A|^2 
%       A  dt           
%
%   sigma = sr + si
%   l     = lr + li
%
%                      /            sr               \     /          li                                        \    
%    Real[A(t)] = Sqrt | --------------------------- |  cos| si*t - ------ ln(lr*exp[2 sr (t-t0)] +(-) 1) + phi |
%                      \   lr +(-) exp[-2 sr (t-t0)] /     \         2 lr                                       /
%    grow : dA/dt > 0, '+'
%    decay: dA/dt < 0, '-'
%
%   INPUT:
%       x   : parameters
%               x(1) = sr
%               x(2) = lr
%               x(3) = t0
%               x(4) = si
%               x(5) = phi
%               x(6) = li
%       t   : observed points
%       y   : response data
%       char: switches
%           'grow'  => dA/dt > 0
%           'decay' => dA/dt < 0
%   OUTPUT:
%       A(t)
%
%   DEPENDENCY:
%       none
%
% Written by Mai, Jiaqi, on 03/19/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sr  = x(1);
    lr  = x(2);
    t0  = x(3);
    si  = x(4);
    phi = x(5);
    li  = x(6);
    switch char
        case 'grow'
            flag = 1;
        case 'decay'
            flag = -1;
        otherwise
            flag = -1;
    end
    X = exp(2*sr*(t-t0));
    I = li/(2*lr)*log(lr*X+flag);
    output = sqrt(sr./(lr+flag./X)).*cos(si*t+phi-I);
end