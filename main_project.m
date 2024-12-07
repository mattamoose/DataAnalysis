clc; close all; clear  

data = readtable('146400-09-25-11-25.csv');

chA = table2array(data(:,2)); chB = table2array(data(:,3));
s = '';
%%
%Data Span 9/25 00:01 - 11/24 23:59  
dt =  2; % mins
n = numel(chA);
df = 1/dt; % cylces/min

t = 0:dt:2*n-1;
figure(1)
plot(t,chA)
hold on
plot(t,chB)
title('Raw PM2.5 vs time')
xlabel('time (min)')
ylabel('PM2.5 Mass Concentration (ug/m^3)')
legend('Channel A','Channel B')

% Average of Both Channels

y = (chA + chB)./2;

figure(2)
plot(t,y)
title('Mean of Channels')
xlabel('time (min)')
ylabel('PM2.5 Mass Concentration (ug/m^3)')
A = [ones(size(t))'  t'];
x = A\y;
fprintf([' \n The linear trend of the data is as follows: \n %10.1s ' ...
    'Intercept = %6.4f \n %10.1s Slope = %6.6f \n\n'],s,x(1),s,x(2))
ym = A*x;
hold on
plot(t,ym)

y = y - ym;
figure(3)
plot(t,y)
title('Detrended Data Set')
dty = mean(y);
fprintf('\n The mean of the detrended data set = %8.8f',dty)
figure(4)

[f,Pyy] = psd(n,df,y);

plot(f,Pyy)
xscale('log')

dom_freq = f(Pyy == max(Pyy));
freq_days = dom_freq*60*24;
fprintf(' \n \n The dominant frequency of PM2.5 is %4.5f cycles/day \n',freq_days)

fprintf('Or a period of %4.2f days \n',1/freq_days)
%% Band-Pass Filter

band = [dom_freq-1e-5*df/n,dom_freq+1e-5*df/n];
 
f_bd = fir1(20,band,'bandpass');

y_l = filter(f_bd,1,y);
figure(6)
plot(t,y_l)
title('Band Pass Filtered Series')

[f,Py_l] = psd(n,df,y_l);
figure(4)
hold on
plot(f,Py_l)
title('PSD of Band-Pass Filter')
xscale('log')
legend('Psd unfiltered','psd filtered')

ex = mean(y_l);
fprintf(' \n The mean contribution is %8.8f',ex)

%% Isolate Commuter traffic band

bandc = 1./[13,9]/60; % cycles/min
f_bc = fir1(20,bandc,'bandpass');

y_c = filter(f_bc,1,y);

Pyy(f < bandc(1)) = 0; Pyy(f > bandc(2));
yc = invpsd(n,Pyy,'symmetric');
I = (yc>30);
yc(I) = [];
tc = t; tc(I) = [];

figure(7)
plot(tc,yc)
hold on
plot(t,y_c)
legend('Crude','FIR')
title('Commuter Band Pass')

cc = mean(yc);

fprintf('\n The mean contribution from commuter traffic is %8.8f',cc)

figure(8)

plot(f,Pyy)
hold on
[f,pyyc] = psd(n,df,y_c);
plot(f,pyyc)
xscale('log')
legend('Psd unfiltered','psd filtered')

figure(9)

spectrogram(y_c)

figure(10)
cwt(y_c)
figure(11)
cwt(y)

figure(12)
[ps,fp] = periodogram(y,[],[],df);
dom = fp(ps == max(ps))*60*24;
fprintf('\nUsing periodogram, the dominant frequency has a period of %4.4f',1/dom)

function [f,Pxx] = psd(n,df,y)
    ffte = fft(y);
    Py = abs(ffte).^2/n;
    Pyy = Py(1:floor(n/2)+1);
    Pxx = [Pyy(1);2*Pyy(2:end)];
    f = 0:df/n:df/2; % cycles/time
end

function [yr] = invpsd(n,Poa,invtype)
    Pxx = [sqrt(Poa(1)*n); sqrt(Poa(2:end)/2*n); sqrt(fliplr(conj(Poa(2:end)/2*n)))];
    yr = ifft(Pxx,n,invtype);

end