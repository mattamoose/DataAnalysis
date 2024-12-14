clear; clc; close all;

%% Problem 1
% frequency response
f = linspace(0,1/2,1000);

H_derived = @(w) 0.5*(exp(1i*2*pi*w)-exp(-1i*2*pi*w));

[H,fz] = freqz([1/2,0,-1/2],1,length(f),1);

figure(1)
plot(f,abs(H_derived(f)),LineWidth=1)
hold on
plot(fz,abs(H))
legend('Analytical','freqz')
svf('freqresp')

%% Problem 2
fprintf('\n Problem 2 \n')
n = 100;
x_bar = 0;
SE = 1;

binedge = [-inf,-1,0,1,inf];
counts = [15,38,35,12];

figure(2)
histogram('BinEdges',binedge,'BinCounts',counts)

pd = makedist('Normal');
rng default;  % for reproducibility
x = random(pd,100,1);

norm_counts = histcounts(x,binedge);

s = SE*sqrt(100);
chi2 = sum((norm_counts - counts).^2/norm_counts);
fprintf('The Chi Square test statistic = %6.4f',chi2)
dof = 3;

%% Problem 3
fprintf('\n Problem 3 \n')

w = 1/10; %rad\s
f = @(t) cos(2*pi*w*t);

k = 0:100;
dt = [8.2,9.0,9.8];
s1 = f(k*dt(1));
s2 = f(k*dt(2));
s3 = f(k*dt(3));

t = dt'*k;

figure(3)
plot(t(1,:),s1)
svf('dteq8_2')
figure(4)
plot(t(2,:),s2)
svf('dteq9')
figure(5)
plot(t(3,:),s3)
svf('dteq9_8')

%% Problem 4
clc;close all

d4 = load('finalq4.mat');
x4 = d4.x; y4 = d4.y;

figure(6)
plot(x4,y4)

A = [ones(size(x4))' x4'];
% Ac = b

c = regress(y4',A);
fprintf('\n The regression coefficients:\n')
fprintf('%14.6f \n',c)
ym = A*c;
resid = y4' - ym;

CI = bootci(1000,{@(r) regress(ym+r,A),resid},'type', 'normal');
fprintf('\n CI of intercept:\n')
fprintf('%12.4f %12.4f',CI(:,1))
fprintf('\n CI of slope:\n')
fprintf('%12.4f %12.4f',CI(:,2))
cl = c - CI(1,:)';
cu = c + CI(2,:)';

yl = @(x) cl(1) + cl(2)*x;
yu = @(x) cu(2) + cu(2)*x;
y = @(x) c(1) + c(2)*x;
t_test = [100 200 300 1000];

y_test = y(t_test);
yl_test = yl(t_test);
yu_test = yu(t_test);

p4 = chi2gof(resid);

fprintf('\n Lower Bound of predicted values:\n')
fprintf('%12.4f \n',yl_test)
fprintf('\n Upper Bound of predicted values:\n')
fprintf('%12.4f \n',yu_test)
fprintf('\n Predicted Values:\n')
fprintf('%12.4f \n',y_test)
fprintf('\n\n Result of Chi2 test:\n')
fprintf('%12.4f \n',p4)

%% Problem 5

fprintf('\n\n Problem 5')
d5 = load('finalq5.mat');
x5 = d5.x';y5 = d5.y';
dt = mean(diff(x5));

[f,Pxx] = psd(numel(y5),1/dt,y5);

figure(get(gcf,'number') +1)
plot(f,Pxx)
xlabel('frequency (Hz)')
xscale('log')
xlim([0,max(f)+1])
svf('psd51')

freq = [0.00999 0.0999 0.999 3.996];

M = vand(freq,x5);
coeff = M\y5; % coeff pattern sine cosine
ym5 = M*coeff;

fprintf(' \n Coefficients for 0.0099:\n')
fprintf(' %12.4f %12.4f',coeff(1:2))
fprintf(' \n Coefficients for 0.0999:\n')
fprintf(' %12.4f %12.4f',coeff(3:4))
fprintf(' \n Coefficients for 0.999:\n')
fprintf(' %12.4f %12.4f',coeff(5:6))
fprintf(' \n Coefficients for 3.996:\n')
fprintf(' %12.4f %12.4f',coeff(7:8))

yr = y5 - ym5;
[f,Pxx] = psd(numel(y5),1/dt,yr);
figure(get(gcf,'Number')+1)
plot(x5,yr)
svf('p5rmvddomfreq')
figure(get(gcf,'Number')+1)
plot(f,Pxx)
xscale('log')
svf('psd52')

A = [ones(size(x5)) x5 x5.^2 x5.^3];
coeff2 = A\yr;

fprintf('\nRemaining variance estimate:\n')
fprintf('%12.4f \n',coeff2)

figure(get(gcf,'Number')+1)
plot(x5,yr-A*coeff2)
svf('p5cubicrmvd')

%% Problem 6
close all
fprintf('\n\n Problem 6\n')
% Low pass filter:
cutoff = 0.3;
b = fir1(50,cutoff,'low'); 
fprintf('Cutoff Frequency:\n %12.4f Hz\n',cutoff)
yf = filter(b,1,y5);

[f,Pxx] = psd(numel(y5),1/dt,yf);
figure(get(gcf,'Number')+1)
plot(f,Pxx)
xscale('log')
svf('lowpassfir_psd')

band =[0.09 0.11];

bb = fir1(39,band,'bandpass',hann(40));
ybf = filter(bb,1,yf);


[f,Pxx] = psd(numel(y5),1/dt,ybf);
figure(get(gcf,'Number')+1)
plot(f,Pxx)
xscale('log')
svf('bandpassfir_psd')

figure(get(gcf,'Number')+1)
plot(x5,ybf)
svf('fbpy5')

%% Problem 7
d7 = load('finalq7.mat');

time = d7.time;
data = d7.data;
MO = NaN(592,51*41);
ll =1;
for ii = 1:51
    for jj = 1:41
      k = data(ii,jj,:);
      MO(:,ll) = k;
      ll = ll +1;
    end
end

if isempty(find(isnan(MO)))
    disp('good to continue')
end

% Center the data
for jj = 1:2091
    v = MO(:,jj);
    MO(:,jj) = v - mean(v);
end

% SVD
[u,s,v] = svd(MO,'econ');

lambda = s.^2;
vtot = sum(lambda,'all');
lambda = 100*lambda/vtot;

disp(lambda(1,1))
disp(lambda(2,2))
disp(lambda(3,3))
disp(lambda(4,4))
disp(lambda(5,5))
disp(lambda(1,1) + lambda(2,2)+lambda(3,3)+lambda(4,4)+lambda(5,5)+lambda(6,6)+lambda(7,7)+lambda(8,8)+lambda(9,9))
V1 = v*s;

for i = 1:size(v,2)
    Vmean(i) = mean(V1(:,i));
end

for j = 1:size(u,2)
    PCTS(:,j) = u(:,j)*Vmean(j);
end

num_modes = 4;
LVm = zeros(size(v));

for k = 1:size(v,2)
    LVm(:,k) = V1(:,k)./Vmean(k);
end

EOF = PCTS(:,9)*LVm(:,9)';

for i = 1:num_modes
    titl = strcat('PCTSmode',string(i));
    figure(get(gcf,'Number')+1)
    plot(time,PCTS(:,i))
    svf(titl)
    figure(get(gcf,'Number')+1)
    titl = strcat('LVm',string(i));
    plot(time,LVm(i,:))
end

% Re-Center the reconstructed data
for jj = 1:2091
    v = MO(:,jj);
    EOF(:,jj) = EOF(:,jj) + mean(v);
end

% restructure the data
new_dat = NaN(size(data));
ll = 1;
for ii = 1:51
    for jj = 1:41
        new_dat(ii,jj,:) = EOF(:,ll);
        ll = ll+1;
    end
end
[X,Y] = meshgrid(1:41,1:51);
for ii = 1:4
    figure(fig())
    tiledlayout (1,2)
    nexttile
    surf(X,Y,new_dat(:,:,ii*100))
    title('EOF')
    nexttile
    surf(X,Y,data(:,:,ii*100))
    title('Original')
    filename = strcat('EOFteq',strcat(ii*100));
    svf(filename)
end
%% Functions
function svf(name)
% saves figures with the name specified
    loc = strcat('Images/',name,'.png'); % provides path and name 
    saveas(gcf,loc) % saves current figure at specified path and name
end

function [f,Pxx] = psd(n,df,y)
% calculates the power spectral density of y, with n points and df =fs
% returns the PSD Pxx and frequency vector f
    ffte = fft(y);
    Py = abs(ffte).^2/n;
    Pyy = Py(1:floor(n/2)+1);
    Pxx = [Pyy(1);2*Pyy(2:end)];
    f = 0:df/n:df/2; % cycles/time
end

function [V,identity] = vand(f,t) 
% constructs the Vandermonde matrix
        SIN = @(w,t) sin(w*t);
        COS = @(w,t) cos(w*t);
        
        S = NaN(1,numel(f));
        C = NaN(size(S));
        
        fr = f*2*pi; % converts f from cycles/time -> rad/time

        V = ones(numel(t),(2*numel(f)));
        %V(:,2) = t';
        identity = NaN(size(V,2),1);
        for ii = 1:numel(f)
           S = SIN(fr(ii),t);
           C = COS(fr(ii),t);
           V(:,(2*ii-1 :2*ii)) = [S, C];

           % creates a vector that can be used to easily indentify and
           % extract specific coefficents
           identity((2*ii-1:2*ii)) = [f(ii);f(ii)];
        end
end

function n=fig()
    n = get(gcf,'Number') +1;
end