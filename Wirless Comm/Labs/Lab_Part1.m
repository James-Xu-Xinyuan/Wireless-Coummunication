% wireless comm lab 1 - 2020Jan20

% #1
clear all
STDdB = 8;
std = 10^(STDdB/20); %sigma square
SdB = randn(1000,1)*std; % this can be smaller than 0
S = db2pow(SdB); % so there would be imaginary part
% var(SdB,0,1)
nbin = 30;
% histogram(S,nbin)
% histogram(SdB,nbin)

% #2
clear all
h = 1/(sqrt(2))*(randn(1000,1)+j*randn(1000,1));
h_mag = abs(h);
h_pha = angle(h);
h_sqrmag = h.*(conj(h));
nbin = 30;
histogram(h_mag,nbin)
figure
histogram(h_pha,nbin)
figure
histogram(h_sqrmag,nbin)

% #3
clear all
a = 2;
h_bar = exp(j*a);
h_tilt = 1/(sqrt(2))*(randn(1000,1)+j*randn(1000,1));
K = 10; % 0.1,1,10,100
h = sqrt(K/(1+K))*h_bar+sqrt(1/(1+K))*h_tilt;
h_mag = abs(h);
histogram(h_mag,20)

% #4
clear all
n =1000;
SSM = zeros(n,1);
for i=1:n
    h = 1/(sqrt(2*n))*(randn(n,1)+j*randn(n,1));
    ssm = sum(abs(h).^2);% sum of sqaured magnitude
    SSM(i)=ssm;
end
histogram(SSM,40)
% mean(SSM); var(SSM);