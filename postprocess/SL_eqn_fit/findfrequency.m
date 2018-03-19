function [amp,freq] = findfrequency(t,x)
    L = length(t);
    dt = t(2) - t(1);
    y = fft(x);
    P2 = abs(y/L);
    P1 = P2(1:round(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    frequency = (0:round(L/2))/dt/L;
    [amp,ind] = findpeaks(P1);
    freq = frequency(ind);
end