close all;
clear hil dragy_temp hilbert_region t_hil lift fit_region t_fit hil_fit lift_fit x0 x r f sigma;

char = 'grow';
    
dragy_temp = dragy;
hilbert_region = round(0.5*length(time)):length(time);
t_hil = time(hilbert_region);
lift  = dragy_temp(hilbert_region);
[hil,temp] = envelope(dragy_temp(hilbert_region),10000,'peak');
    
fit_region = round(0.1*length(t_hil)):round(0.9*length(t_hil));
t_fit = t_hil(fit_region);
hil_fit = hil(fit_region);
lift_fit = lift(fit_region);
   
[x,f,r,x0] = SL_equation_fit(t_fit,lift_fit,char);
sigma = x(1);
    %% 
    
figure
plot(t_fit,lift_fit);hold on
plot(t_fit,f)
    


    
 