% wirless comm coursework 1
% xx316, Xinyuan Xu, 2020 Jan
% uncoded QPSK transmission over a SISO Rayleigh fading channel
clear all;
%define constellation diagram
c00 = 1/sqrt(2)*(+1+1i); 
c01 = 1/sqrt(2)*(-1+1i); 
c11 = 1/sqrt(2)*(-1-1i); 
c10 = 1/sqrt(2)*(+1-1i);

L = 10^6; 
BitStream = round(rand(L,1)); 
% a random bit stream of 0 and 1 

c = zeros(L/2,1);
for j =1:length(c)
    if BitStream(2*j-1)==0
        if BitStream(2*j)==0
            c(j)=c00;
        else
            c(j)=c01;
        end
    else
        if BitStream(2*j)==0
            c(j)=c10;
        else
            c(j)=c11;
        end
    end
end

%assume normalized noise: variance of n is 1
%average SNR = Es/var(n) = Es numerically
BER = zeros(201,1);
BER_N = zeros(201,1);
SNRdB = 0:0.1:20;
%iterate for each SNR value
for j = 1:201
    Es = db2pow(SNRdB(j)); 
    % generate iid channel and noise
    h = 1/(sqrt(2))*(randn(L/2,1)+1i*randn(L/2,1));
    n = 1/(sqrt(2))*(randn(L/2,1)+1i*randn(L/2,1));

    y = sqrt(Es)*h.*c+n; % simulated received signal
    y_N = sqrt(Es).*c+n;%gaussain channel

    z = y.*conj(h);%assume the channel h is perfectly known at receiver

    ReceivedBits = zeros(L,1);
    % demodulation
    for k = 1:length(y)
        if imag(z(k))>0
            ReceivedBits(k*2-1) = 0;
        else
            ReceivedBits(k*2-1) = 1;
        end
        if real(z(k))>0
            ReceivedBits(2*k) = 0;
        else
            ReceivedBits(2*k) = 1;
        end
    end
    
    ReceivedBits_N = zeros(L,1);%gaussain channel
    for k = 1:length(y_N)
        if imag(y_N(k))>0
            ReceivedBits_N(k*2-1) = 0;
        else
            ReceivedBits_N(k*2-1) = 1;
        end
        if real(y_N(k))>0
            ReceivedBits_N(2*k) = 0;
        else
            ReceivedBits_N(2*k) = 1;
        end
    end
    
    
    %calculate bit error rate
    diff = mod(BitStream+ReceivedBits,2);
    ber = sum(diff)/L;
    BER(j) = ber;
    diff_N = mod(BitStream+ReceivedBits_N,2);
    ber_N = sum(diff_N)/L;
    BER_N(j) = ber_N;
end

plot(SNRdB, log10(BER))
hold on
Pe = 1-sqrt(db2pow(SNRdB)./(1.+db2pow(SNRdB))); % theoretical at high SNR
plot(SNRdB,log10(Pe))
plot(SNRdB, log10(BER_N))
legend('BER in Experiment','BER in Theory','BER of Gaussian Channel')
xlabel('SNR in dB')
ylabel('log10(BER) (Bit Error Rate)')
title('Uncoded QPSK - SISO Rayleigh fading channel')

SISO_BER = BER;
clear BER;
load SIMO.mat
SIMO_BER = BER;
SIMO_SNRout = SNRout;
clear BER;  clear SNRout;
load MISO_MRT.mat
MISO_MRT_BER = BER;
MISO_SNRout = SNRout;
clear BER;  clear SNRout;
load Alamouti.mat
Alamouti_BER = BER;
Alamouti_SNRout = SNRout;
clear BER;  clear SNRout;

figure
plot(SNRdB, log10(SISO_BER),'r')
hold on
plot(SNRdB, log10(SIMO_BER),'b')
plot(SNRdB, log10(MISO_MRT_BER),'g')
plot(SNRdB, log10(Alamouti_BER),'k')
legend('SISO','SIMO-MRC','MISO-MRT','MISO-Alamouti')
title('Bit Error Rates of different channels')
xlabel('SNRdB')
ylabel('log10(BER)')

figure
plot(SNRdB, 10*log10(SIMO_SNRout),'b')
hold on
plot(SNRdB, 10*log10(MISO_SNRout),'g')
plot(SNRdB, 10*log10(Alamouti_SNRout),'r')
legend('SIMO-MRC','MISO-MRT','MISO-Alamouti')
title('Output SNR of different channels')
xlabel('SNRdB')
ylabel('SNRout in dB')
