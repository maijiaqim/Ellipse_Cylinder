close all;
clear hil;

hilbert_region = 250000:length(time);
t = time(hilbert_region);
fit_region = 30000:(length(t)-35000);
col = size(decay,2);

for i = 1:col
    [hil(:,i),temp] = envelope(decay(hilbert_region,i),10000,'peak');  %abs(hilbert(decay(hilbert_region,i)));
end
    
% Check the envelope amplitude    
%   for i = 1:col
%       figure
%       plot(t,hil(:,i));hold on
%       plot(time,decay(:,i));
%   end


%% 

% Re = 46.1:0.1:47.0;

dt = time(2)-time(1);
for i = 1:col
    S(:,i) = (hil((fit_region+1),i) - hil((fit_region-1),i))./(2.0*dt)./hil(fit_region,i);
    R(:,i) = hil(fit_region,i).^2;
%     figure
%     plot(t(fit_region),S(:,i))
%     figure
%     plot(t(fit_region),R(:,i))
end
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       linear fitting to determine 
%       sigma and 
%       saturated value CL_square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:col
    P = polyfit(R(100000:end,i),S(100000:end,i),1);
    sigma(i) = P(2);
    CL_square(i) = -P(2)/P(1);
%     figure
%     plot(R(100000:end,i),S(100000:end,i))
end

P1 = polyfit(Re,sigma,1);
Re_crit_1 = -P1(2)/P1(1)
figure
plot(Re,sigma,'o'); hold on
plot(Re,P1(1)*Re+P1(2),'-'); hold on
plot(Re_crit_1,0,'x')

P2 = polyfit(Re,CL_square,1);
Re_crit_2 = -P2(2)/P2(1)
figure
plot(Re,CL_square,'o'); hold on
plot(Re,P2(1)*Re+P2(2),'-'); hold on
plot(Re_crit_2,0,'x')

