close all;
clear hil col dragy_temp hilbert_region t_hil lift fit_region t_fit hil_fit lift_fit x0 x r f sigma;

char = 'grow';
    
col = size(decay,2);
for i = 1:col
    dragy_temp = decay(:,i);
    hilbert_region = round(0.4*length(time)):length(time);
    t_hil = time(hilbert_region);
    lift  = dragy_temp(hilbert_region);
    [hil,temp] = envelope(dragy_temp(hilbert_region),10000,'peak');
    
    fit_region = round(0.1*length(t_hil)):round(0.9*length(t_hil));
    t_fit = t_hil(fit_region);
    hil_fit = hil(fit_region);
    lift_fit = lift(fit_region);
    
    [x(:,i),f,r,x0(:,i)] = SL_equation_fit(t_fit,lift_fit,char);
    sigma(i) = x(1,i);
    
    figure
    plot(t_fit,f);hold on
    plot(t_fit,lift_fit);
end
 
P1 = polyfit(Re,sigma,1);
Re_crit = -P1(2)/P1(1)
figure
pp1 = plot(Re,sigma,'ro'); hold on
plot(Re,P1(1)*Re+P1(2),'r-');hold on
plot([Re_crit,Re],[0,P1(1)*Re+P1(2)],'r--'); hold on
pp2 = plot(Re_crit,0,'rx');hold on
    
%     P1 = polyfit(Re,sigma2,1);
%     Re_crit = -P1(2)/P1(1)
%     pp3 = plot(Re,sigma2,'bo'); hold on
%     plot(Re,P1(1)*Re+P1(2),'b-');hold on
%     plot([Re_crit,Re],[0,P1(1)*Re+P1(2)],'b--'); hold on
%     pp4 = plot(Re_crit,0,'bx');hold off
%     legend([pp1,pp2,pp3,pp4]);
    
 