clc; close all; clear  

data = readtable('146400-09-25-11-25.csv');

chA = table2array(data(:,2)); chB = table2array(data(:,3));

%%
%Data Span 9/25 00:01 - 11/24 23:59  
dt =  2; % mins
n = numel(chA);

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

figure(3)
periodogram(y)