clc; close all; clear  

data = readtable('146400 2024-04-01 2024-10-31 0-Minute Average.csv');

chA = table2array(data(:,2)); chB = table2array(data(:,3));

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

y = detrend(chB);

mean = mean(y);
disp(mean)

figure(3)
plot(t,y)
title('Detrended Data Set')

figure(4)

[f,Pyy] = psd(n,df,y);

plot(f,Pyy)
xscale('log')

dom_freq = f(Pyy == max(Pyy));
freq_days = dom_freq*60*24;
fprintf('The dominant frequency of PM2.5 is %4.5f cycles/day \n\n',freq_days)

fprintf('Or a period of %4.2f days \n',1/freq_days)

dom_freq2 = 9.513e-5;
freq_days = dom_freq2*60*24;
fprintf('The third dominant frequency of PM2.5 is %4.5f cycles/day \n\n',freq_days)

fprintf('Or a period of %4.2f days \n',1/freq_days)
%% High-Pass Filter

cutoff = dom_freq+10*df/n;
 
f_l = fir1(50,cutoff,"high");

y_l = filter(f_l,1,y);
figure(6)
plot(t,y_l)
hold on
title('Low Pass Filtered Series')

[f,Py_l] = psd(n,df,y_l);
figure(7)
plot(f,Py_l)
title('PSD of Low-Pass Filter')
xscale('log')


function [f,Pxx] = psd(n,df,y)
    ffte = fft(y);
    Py = abs(ffte).^2/n;
    Pyy = Py(1:floor(n/2)+1);
    Pxx = [Pyy(1);2*Pyy(2:end)];
    f = 0:df/n:df/2; % cycles/time
end