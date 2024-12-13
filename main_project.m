clc; close all;

data = readtable('146400-09-25-11-25.csv');
mkdir('Images')

chA = table2array(data(:,2)); chB = table2array(data(:,3));
time_0 = datetime([2024,9,25,0,1,0]);
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
svf('chAB')

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
svf('Mean_trend')


y = y - ym;
figure(3)
plot(t,y)
title('Detrended Data Set')
xlabel('time (min)')
ylabel('PM2.5 Mass Concentration (ug/m^3)')
dty = mean(y);
fprintf('\n The mean of the detrended data set = %8.8f',dty)


figure(4)

[f,Pyy] = psd(n,df,y);

plot(f,Pyy)
xscale('log')
title('PSD')
xlabel('Frequency (cycles/min)')
svf('PSD')
dom_freq = f(Pyy == max(Pyy));
fprintf('\n The dominant frequency is %6.6f Hz \n',dom_freq)
freq_days = dom_freq*60*24;
fprintf(' \n \n The dominant frequency of PM2.5 is %4.5f cycles/day \n',freq_days)

fprintf('Or a period of %4.2f days \n',1/freq_days)


%% Isolate Commuter traffic band

bandc = 1./[14,8]/3600; % cycles/min

[Pyy2,f_iso] = freqbuild(Pyy,30,f,2);

fiso_srt = sort(f_iso);

f_lp = fir1(50,max(f_iso)+0*df/n,"low");
y_flp = filter(f_lp,1,y);
fprintf('\n Mean of Filtered Data set: %6.5f',mean(y_flp))

%y_flp = y_flp - mean(y_flp);

[M, f_ind] = vand(fiso_srt,t);
coeff = M\y_flp;

b = regress(y_flp,M);
ym = M*b;
resid = y_flp - ym;

% SE = std(bootstrp(1000,@(bootr)regress(ym+bootr,M),resid));
% 
% CI = bootci(1000,{@(bootr)regress(ym+bootr,M),resid},'type','normal');

pCHI2 = chi2gof(resid,'Alpha',0.05);
pZ = ztest(resid,mean(resid),std(resid));

bmin = b - CI(1,:)';
bmax = CI(2,:)' - b;
figure(13)
plot(1:numel(SE),SE'./b)
xlim([0,numel(SE)+1])
xlabel('Coefficent')
ylabel('Percent Standard Error')
svf('SE')

figure(5)
errorbar(1:numel(b),b,bmin,bmax)
xlim([0,numel(b) + 1])
ylabel('Confidence Interval')
xlabel('Coefficent')

svf('CI')


figure(6)
plot(t,y_flp)
hold on
plot(t,M*coeff)
xlabel('time (minutes)')
ylabel('Mass Concentration (ug/m^3)')
svf('MdlFilt')
sat = datetime([2024,9,27,6,0,0]);
sun = datetime([2024,9,29,23,0,0]);

t_0 = minutes(sat - time_0);
t_f = minutes(sun - sat);
dwend = 7*24*60; %minutes/week

day = sat;
X = [t_0, t_0, t_0 + t_f, t_0 + t_f]; Y = [-20 50 50 -20];
while day < datetime([2024,11,25,0,0,0])
    patch(X,Y,'g',FaceAlpha=0.2)
    X = X + dwend;
    day = day + minutes(dwend);
end
svf('mdlfilwknd')
times = [9190,16630,18910,21560,26890,60150,64028,75294,78792,83116,84602,];
times = minutes(times);
hi_resid = time_0 + times;

fprintf('\n\n Dates of High Residuals from the frequency analysis:\n')
disp(hi_resid)

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
ts11 = minutes(datetime([2024,10,6,0,0,0]) - time_0);   % Definitely Elk Fire
ts12 = minutes(datetime([2024,10,10,9,59,59]) -time_0);
Xs = [ts11, ts11, ts12, ts12];
patch(Xs,Y,[0.7 0.7 0.7],FaceAlpha = 0.3)

tw11 = minutes(datetime([2024,11,18,12,0,0]) - time_0); % Steady winds around 20 mph gusts to 35 mph.
tw12 = minutes(datetime([2024,11,18,21,0,0]) - time_0);
Xw = [tw11 tw11 tw12 tw12]; 
patch(Xw,Y,'b',FaceAlpha = 0.2)

tws11 = minutes(datetime([2024,9,30,21,0,0])-time_0);
tws12 = minutes(datetime([2024,10,2,23,59,59]) - time_0); % Elk Fire end of period of southern winds.

tws21 = minutes(datetime([2024,11,5,0,0,0]) - time_0);
tws22 = minutes(datetime([2024,11,6,0,0,0]) - time_0);
Xws2 = [tws21 tws21 tws22 tws22];

Xws = [tws11 tws11 tws12 tws12];
patch(Xws,Y,[0.7 0.7 0.7],FaceAlpha = 0.3)
patch(Xws2,Y,'b',FaceAlpha = 0.3)
svf('physex')

figure(11)

[cfs,fw] = cwt(y_flp,df/60);
imagesc(t,fliplr(fw),abs(cfs))
ylabel('frequency (Hz)')
xlabel('time (min)')
svf('cwt')
hold on
yyaxis right
day = sat;
X = [t_0, t_0, t_0 + t_f, t_0 + t_f]; Y = [-20 50 50 -20];
while day < datetime([2024,11,25,0,0,0])
    patch(X,Y,'r',FaceAlpha=0.01)
    X = X + dwend;
    day = day + minutes(dwend);
end
svf('cwtwknd')

figure(12)
tiledlayout (1,2)
nexttile
histogram(resid)
nexttile
probplot(resid)


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


function [pxx,omegas] = freqbuild(Pxx,n,f,ndf)
   omegas = NaN(n,1);
   df = f(2) - f(1);
   for ii = 1:n
       fnan = NaN(size(f));
       index = ((Pxx == max(Pxx)) == 1);
       jj = find(index);
       omegas(ii) = f(index);
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
       fspan = f(index) - s1*df:df: f(index) + s2*df;
       fnan(1,(a:b)) = fspan;
       Pxx(f == fnan) = 0;
   end
   pxx = Pxx;    

end

function [V,identity] = vand(f,t)
        SIN = @(w,t) sin(w*t);
        COS = @(w,t) cos(w*t);
        
        S = NaN(1,numel(f));
        C = NaN(size(S));
        
        fr = f*2*pi;

        V = ones(numel(t),(2*numel(f)+2));
        V(:,2) = t';
        identity = NaN(size(V,2),1);
        for ii = 1:numel(f)
           S = SIN(fr(ii),t');
           C = COS(fr(ii),t');
           V(:,(2*ii +1 :2*ii+2)) = [S, C];
           identity((2*ii-1:2*ii)) = [f(ii);f(ii)];
        end
        
        

end

function svf(name)
    loc = strcat('Images/',name,'.png');
    saveas(gcf,loc)
end