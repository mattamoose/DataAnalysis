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


%% Isolate Commuter traffic band

bandc = 1./[14,8]/60; % cycles/min

% f_bc = fir1(20,bandc,'bandpass');
% 
% y_c = filter(f_bc,1,y);
% 
% Pyy(f < bandc(1)) = 0; Pyy(f > bandc(2));
% yc = invpsd(n,Pyy,'symmetric');
% I = (yc>30);
% yc(I) = [];
% tc = t; tc(I) = [];
% 
% figure(7)
% plot(tc,yc)
% hold on
% plot(t,y_c)
% legend('Crude','FIR')
% title('Commuter Band Pass')
% 
% cc = mean(yc);
% 
% fprintf('\n The mean contribution from commuter traffic is %8.8f',cc)
% 
% figure(8)
% 
% plot(f,Pyy)
% hold on
% [f,pyyc] = psd(n,df,y_c);
% plot(f,pyyc)
% xscale('log')
% legend('Psd unfiltered','psd filtered')

[Pyy2,f_iso] = freqbuild(Pyy,32,f,2);

fiso_srt = sort(f_iso);
[M, f_ind] = vand(fiso_srt,t);
coeff = M\y;

b = regress(y,M);
ym = M*b;
resid = y - ym;

% SE = std(bootstrp(1000,@(bootr)regress(ym+bootr,M),resid));
% 
% CI = bootci(1000,{@(bootr)regress(ym+bootr,M),resid},'type','normal');

f_lp = fir1(20,max(f_iso)+0*df/n,"low");
y_flp = filter(f_lp,1,y);

p = chi2gof(resid,'Alpha',0.05);


figure(5)
plot(t,y_flp)
hold on
plot(t,M*coeff)
legend('lp filtered','20-Dom freq regress')

figure(4)
hold on
plot(f,Pyy2)
set(gca,'XScale','log')
legend('PSD', 'Reduced PSD')

figure(6)
y_dom = b(1)*sin(f_iso(1)*2*pi*t') + b(2)*cos(f_iso(1)*2*pi*t');
plot(t,y_dom)
hold on
plot(t,y_flp)
% figure(9)
% 
% spectrogram(y_c)

% figure(10)
% cwt(y_c)
figure(11)
cwt(y)

figure(12)
tiledlayout (1,2)
histogram(resid)
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

        V = NaN(numel(t),(2*numel(f)));
        % V(:,2) = t';
        identity = NaN(size(V,2),1);
        for ii = 1:numel(f)
           S = SIN(fr(ii),t');
           C = COS(fr(ii),t');
           V(:,(2*ii - 1 :2*ii)) = [S, C];
           identity((2*ii-1:2*ii)) = [f(ii);f(ii)];
        end
        
        

end