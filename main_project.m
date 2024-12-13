clc; close all;

data = readtable('146400-09-25-11-25.csv');
mkdir('Images') % creates directory Images in current working folder

chA = table2array(data(:,2)); chB = table2array(data(:,3)); % converts tables to arrays
time_0 = datetime([2024,9,25,0,1,0]); % intial date time of time series
s = ''; % empty string array for easy white spacing in fprintf
%%
%Data Span 9/25 00:01 - 11/24 23:59  
dt =  2; % mins
n = numel(chA);
df = 1/dt; % cylces/min

t = 0:dt:2*n-1; %time span

%plots raw data channels A & B
figure(1)
plot(t,chA)
hold on
plot(t,chB)
title('Raw PM2.5 vs time')
xlabel('time (min)')
ylabel('PM2.5 Mass Concentration (ug/m^3)')
legend('Channel A','Channel B')
svf('chAB')

% Average Both Channels include plot

y = (chA + chB)./2;

figure(2)
plot(t,y)
title('Mean of Channels')
xlabel('time (min)')
ylabel('PM2.5 Mass Concentration (ug/m^3)')

% Calculate the linear trend of the mean of the channels
A = [ones(size(t))'  t'];
x = A\y;
fprintf([' \n The linear trend of the data is as follows: \n %10.1s ' ...
    'Intercept = %6.4f \n %10.1s Slope = %6.6f \n\n'],s,x(1),s,x(2))
ym = A*x;
hold on
plot(t,ym)
svf('Mean_trend')

% Calculate and plot the detrended data set
y = y - ym;
figure(3)
plot(t,y)
title('Detrended Data Set')
xlabel('time (min)')
ylabel('PM2.5 Mass Concentration (ug/m^3)')
dty = mean(y);
fprintf('\n The mean of the detrended data set = %8.8f',dty)


figure(4)
% Calculate and plot the Power spectral density of the detrended data 
[f,Pyy] = psd(n,df,y);

plot(f,Pyy)
xscale('log')
title('PSD')
xlabel('Frequency (cycles/min)')
svf('PSD')

% Return the dominant frequency of the PSD
dom_freq = f(Pyy == max(Pyy));
fprintf('\n The dominant frequency is %6.6f Hz \n',dom_freq)
freq_days = dom_freq*60*24; % convert to cycles/day
fprintf(' \n \n The dominant frequency of PM2.5 is %4.5f cycles/day \n',freq_days)

fprintf('Or a period of %4.2f days \n',1/freq_days) % period of dom freq


%% Isolate Commuter traffic band

bandc = 1./[14,8]/3600; % cycles/min

% Builds a vector of the dominant frequencies specified in the function
% input, also returns the modified PSD 
[Pyy2,f_iso] = freqbuild(Pyy,30,f,2);

fiso_srt = sort(f_iso); % wsorts the isolated frequencies by magnitude

f_lp = fir1(50,max(f_iso)+0*df/n,"low"); % builds the low pass filter
y_flp = filter(f_lp,1,y); % filters the detrended data
fprintf('\n Mean of Filtered Data set: %6.5f',mean(y_flp))

[M, f_ind] = vand(fiso_srt,t); % builds the Vandermonde matrix
coeff = M\y_flp; % calculates the least square solution for the filtered data regression

% regression model to be used for bootstrapping
b = regress(y_flp,M); 
ym = M*b;
resid = y_flp - ym;

% Bootstrapped standard error of the regresssion ceofficients
SE = std(bootstrp(1000,@(bootr)regress(ym+bootr,M),resid)); 

% Bootstrapped confidence intervals of the regression coefficients
CI = bootci(1000,{@(bootr)regress(ym+bootr,M),resid},'type','normal');

pCHI2 = chi2gof(resid,'Alpha',0.05); % Chi square distribution test of residuals
pZ = ztest(resid,mean(resid),std(resid)); % non-zero mean normal distribution test of residuals

% calculates min and max of CI for regression coefficients
bmin = b - CI(1,:)';
bmax = CI(2,:)' - b;

% plots the SE of the regression coefficients
figure(13)
plot(1:numel(SE),SE'./b)
xlim([0,numel(SE)+1])
xlabel('Coefficent')
ylabel('Percent Standard Error')
svf('SE')

% plots the error bar plot of CI of regression coefficients
figure(5)
errorbar(1:numel(b),b,bmin,bmax)
xlim([0,numel(b) + 1])
ylabel('Confidence Interval')
xlabel('Coefficent')

svf('CI')

% plots the filtered time series and saves before adding the weekend
% patches
figure(6)
plot(t,y_flp)
hold on
plot(t,M*coeff)
xlabel('time (minutes)')
ylabel('Mass Concentration (ug/m^3)')
svf('MdlFilt')
%creates the date time for the first weekend in the series
sat = datetime([2024,9,27,6,0,0]);
sun = datetime([2024,9,29,23,0,0]);

%translates the first weekend to the minutes of the series
t_0 = minutes(sat - time_0);
t_f = minutes(sun - sat);
dwend = 7*24*60; %minutes/week 

% plots the weekend patches on the filtered sereis and model
day = sat;
X = [t_0, t_0, t_0 + t_f, t_0 + t_f]; Y = [-20 50 50 -20];
while day < datetime([2024,11,25,0,0,0])
    patch(X,Y,'g',FaceAlpha=0.2)
    X = X + dwend;
    day = day + minutes(dwend);
end
svf('mdlfilwknd')

% Manual entry of the minutes of spikes of PM
times = [9190,16630,18910,21560,26890,60150,64028,75294,78792,83116,84602,];
times = minutes(times);
hi_resid = time_0 + times;

fprintf('\n\n Dates of High Residuals from the frequency analysis:\n')
disp(hi_resid)

sat = datetime([2024,9,27,6,0,0]);
sun = datetime([2024,9,29,23,0,0]);


% plot of residuals with weekend patches
figure(7)
plot(t,resid,'r')
xlabel('time (minutes)')
ylabel('Residual (ug/m^3)')
hold on
day = sat;
X = [t_0, t_0, t_0 + t_f, t_0 + t_f]; Y = [-20 50 50 -20];
while day < datetime([2024,11,25,0,0,0])
    patch(X,Y,'g',FaceAlpha=0.05)
    X = X + dwend;
    day = day + minutes(dwend);
end
%creates and plots the patches of the elk fire on the residuals
ts11 = minutes(datetime([2024,10,6,0,0,0]) - time_0);   % Definitely Elk Fire
ts12 = minutes(datetime([2024,10,10,9,59,59]) -time_0);
Xs = [ts11, ts11, ts12, ts12];
patch(Xs,Y,[0.7 0.7 0.7],FaceAlpha = 0.3)


% wind event plotted as patch on residuals
tw11 = minutes(datetime([2024,11,18,12,0,0]) - time_0); % Steady winds around 20 mph gusts to 35 mph.
tw12 = minutes(datetime([2024,11,18,21,0,0]) - time_0);
Xw = [tw11 tw11 tw12 tw12]; 
patch(Xw,Y,'b',FaceAlpha = 0.2)

%smoke period patch
tws11 = minutes(datetime([2024,9,30,21,0,0])-time_0);
tws12 = minutes(datetime([2024,10,2,23,59,59]) - time_0); % Elk Fire end of period of southern winds.

% Weather event patch
tws21 = minutes(datetime([2024,11,5,0,0,0]) - time_0);
tws22 = minutes(datetime([2024,11,6,0,0,0]) - time_0);
Xws2 = [tws21 tws21 tws22 tws22];

Xws = [tws11 tws11 tws12 tws12];
patch(Xws,Y,[0.7 0.7 0.7],FaceAlpha = 0.3) %plots smoke event
patch(Xws2,Y,'b',FaceAlpha = 0.3) %plots weather event
svf('physex')

% computes and plots the the continuous wavelet transform of the filtered
% data
figure(11)
[cfs,fw] = cwt(y_flp,df/60);
imagesc(t,fliplr(fw),abs(cfs)) % changes cwt units to not be noramlized
ylabel('frequency (Hz)')
xlabel('time (min)')
svf('cwt')
hold on

% plots weekend patches on the cwt, image not included in report
yyaxis right
day = sat;
X = [t_0, t_0, t_0 + t_f, t_0 + t_f]; Y = [-20 50 50 -20];
while day < datetime([2024,11,25,0,0,0])
    patch(X,Y,'r',FaceAlpha=0.01)
    X = X + dwend;
    day = day + minutes(dwend);
end
svf('cwtwknd')


%% User Defined Functions

function [f,Pxx] = psd(n,df,y)
% calculates the power spectral density of y, with n points and df =fs
% returns the PSD Pxx and frequency vector f
    ffte = fft(y);
    Py = abs(ffte).^2/n;
    Pyy = Py(1:floor(n/2)+1);
    Pxx = [Pyy(1);2*Pyy(2:end)];
    f = 0:df/n:df/2; % cycles/time
end


function [pxx,omegas] = freqbuild(Pxx,n,f,ndf)
% builds the set of dominant frequencies
% n = number of isolated frequencies
% f = frequency vector
% Pxx = PSD 
% ndf = number of df step sizes to zero out around a dominant frequency
%       - helps to avoid double counting
   omegas = NaN(n,1);
   df = f(2) - f(1); % f step size
   for ii = 1:n
       fnan = NaN(size(f));
       index = ((Pxx == max(Pxx)) == 1);
       jj = find(index);
       omegas(ii) = f(index);
       %conditions to avoid runtime errors if dominant frequency is within
       %ndf of ends of vector
       if jj - ndf <= 0
            a = 1;
            b = jj + ndf;
            s1 = ndf - abs(jj -ndf -1);
            s2 = ndf;
       elseif jj + ndf > numel(f)
            a = jj - ndf;
            b = numel(f);
            s2 = ndf - (jj +ndf -numel(f));
            s1 = ndf;
       else 
           a = jj -ndf;
           b = jj + ndf;
           s1 = ndf;
           s2 = ndf;
       end
       fspan = f(index) - s1*df:df: f(index) + s2*df; %span of frequencies to zero
       fnan(1,(a:b)) = fspan; 
       Pxx(f == fnan) = 0; % zeros out the the PSD for identified frequencies
   end
   pxx = Pxx;    % returns the PSD with zeroed out frequencies

end

function [V,identity] = vand(f,t) 
% constructs the Vandermonde matrix
        SIN = @(w,t) sin(w*t);
        COS = @(w,t) cos(w*t);
        
        S = NaN(1,numel(f));
        C = NaN(size(S));
        
        fr = f*2*pi; % converts f from cycles/min -> rad/min

        V = ones(numel(t),(2*numel(f)+2));
        V(:,2) = t';
        identity = NaN(size(V,2),1);
        for ii = 1:numel(f)
           S = SIN(fr(ii),t');
           C = COS(fr(ii),t');
           V(:,(2*ii +1 :2*ii+2)) = [S, C];

           % creates a vector that can be used to easily indentify and
           % extract specific coefficents
           identity((2*ii-1:2*ii)) = [f(ii);f(ii)];
        end
        
        

end

function svf(name)
% saves figures with the name specified
    loc = strcat('Images/',name,'.png'); % provides path and name 
    saveas(gcf,loc) % saves current figure at specified path and name
end