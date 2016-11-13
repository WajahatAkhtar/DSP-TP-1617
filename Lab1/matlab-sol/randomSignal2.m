clc, close all, clear all 

%% Generating Gaussian noise 
N = 1000; 
x_g = randn(1,N);
x_mean = mean(x_g);
x_std = std(x_g);
[h, xh] = hist(x_g);
dh = xh(2)-xh(1);
h = h/(sum(h)*dh);
figure(1), plot(xh,h);
x_gaussian = exp(-0.5*((xh-x_mean)/x_std).^2);
x_gaussian = x_gaussian * (1 / (x_std * sqrt(2*pi)));
hold on; 
plot(xh, x_gaussian, 'r--');
hold off; 

%% Gnerating white noise 
N = 1000; 
x_u = rand(1,N);
x_mean = mean(x_u);
x_std = std(x_u);
[h, xh] = hist(x_u);
dh = xh(2)-xh(1);
h = h/(sum(h)*dh);
figure(2), plot(xh,h);
x_uniform(1:length(xh)) = 1*max(h); 
hold on ; plot(xh, x_uniform, 'r--');
hold off; 

%% cross correlation of a gaussian noise 
gaussian_xcorr = xcorr(x_g,'biased');
uniform_xcorr = xcorr(x_u, 'biased');


figure(4); plot(gaussian_xcorr);
title('white noise');
figure(5); plot(uniform_xcorr)
title('not a white noise');


%% 
s1=round(rand(1,50));
s2=round(rand(1,50));
s3=round(rand(1,50)); 
s(101:150)=s1;
s(191:240)=s2;
s(301:350)=s3;
%sum(s1)
xc1=xcorr(s,s1); 
figure(6), plot(xc1);
xc2=xcorr(s,s2); 
figure(7), plot(xc2);
xc3=xcorr(s,s3); 
figure(8), plot(xc3);