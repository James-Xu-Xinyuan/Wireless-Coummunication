% Xinyuan Xu
% Wireless Commm 2020
% Task 2

clear all
N = 10^4;
% drop a large number of users
LT_SINR = zeros(N,1);

rng('default')
for k = 1:N
    % user coordinate
    theta_i = rand(1,1)*2*pi;
    ri = (rand*(250-35)+35)/1000;
    x = ri*cos(theta_i);
    y = ri*sin(theta_i);
    di = ri; %in km
    
    % interfering base station location
    theta_j = [0,1,2,3,4,5]*(2*pi/6);
    x_j = 500/1000*cos(theta_j);
    y_j = 500/1000*sin(theta_j);
    dj = zeros(6,1);
    
    for j = 1:6
        dj(j) = sqrt((x_j(j)-x)^2+(y_j(j)-y)^2);
    end
    
    % pretend A doesnt have the horizontal bar - path loss
    A0j_dB = 128.1+37.6*log10(dj); 
    Sj_dB = randn(6,1)*8; % 6 interfering cells     % std of 8dB
    Aj = db2pow(A0j_dB+Sj_dB);% pass loss + shadowing
    
    A0i_dB = 128.1+37.6*log10(di);
    Si_dB = randn*8; % user cell shadowing     % std of 8dB
    Ai = db2pow(A0i_dB+Si_dB); % pass loss + shadowing
    
    Es = 1*10^(-3)*10^(46/10);      % Es = 46dBm
    Pn = 1*10^(-3)*10^(-174/10);    % noise variance = -174dBm
    % the noise level is super low ?

    LT_SINR(k)=SINR_LT(Ai,Aj,Es, Pn);
end

LT_SINR_dB = pow2db(LT_SINR);

cdfplot(LT_SINR_dB)
title("CDF of Long Term SINR")
xlabel("x - Long Term SINR")
ylabel("F(x) - CDF")
xlim([-50,50])
% normal distribution of long term SINR

mean(LT_SINR_dB)
std(LT_SINR_dB)
% One-sample Kolmogorov-Smirnov test
h = kstest(LT_SINR) % logical 1

