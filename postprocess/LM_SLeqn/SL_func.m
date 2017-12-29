function output = myfunc(x,t,char)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The solution of Stuart-Landau Equation
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
%       A(t)
%
% Written by Mai, Jiaqi, on 12/28/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigma = x(1);
    k = x(2);
    t0 = x(3);
    switch char
        case 'grow'
            flag = 1;
        case 'decay'
            flag = -1;
        otherwise
            flag = -1;
    end
    output = sqrt(sigma./(k+flag*exp(-2*sigma*(t-t0))));
end