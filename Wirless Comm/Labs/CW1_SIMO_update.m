% wirless comm coursework 1
% xx316, Xinyuan Xu, 2020 Jan
% uncoded QPSK transmission over a SIMO i.i.d Rayleigh fading channel
% with MRC and two receive antennas
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
BER = zeros(1,201);
SNRdB = 0:0.1:20;
SNRout = zeros(1,201);
for j = 1:201
    Es = db2pow(SNRdB(j)); 
    h = 1/(sqrt(2))*(randn(2,L/2)+1i*randn(2,L/2));
    n = 1/(sqrt(2))*(randn(2,L/2)+1i*randn(2,L/2));
    
%     y = sqrt(Es)*h.*c+n; % simulated received signal
    y = zeros(2,L/2);
    SignalOut = zeros(1,L/2);
    NoiseOut = zeros(1,L/2);
    g = h';  %combiner: assume the channel h is perfectly known at receiver
    z = zeros(L/2,1);
    for m = 1:L/2
        y(:,m)= sqrt(Es)*c(m).*h(:,m)+n(:,m);
        z(m)= g(m,:)*y(:,m);
        SignalOut(m) = sqrt(Es)*g(m,:)*(c(m).*h(:,m));
        NoiseOut(m) = g(m,:)*n(:,m);
    end
    SNRout(j)= var(real(SignalOut))/var(real(NoiseOut)); % the factor of two on both sides cancel out

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

    
    diff = mod(BitStream+ReceivedBits,2);
    ber = sum(diff)/L;
    BER(j) = ber;
end

% plot(SNRdB, log10(BER))
% hold on
% Pe = 1-2*sqrt(db2pow(SNRdB)./(1.+db2pow(SNRdB)))+sqrt(db2pow(SNRdB)./(2.+db2pow(SNRdB))); % theoretical at high SNR
% plot(SNRdB,log10(Pe))
% hold off
% legend('BER in Experiment','BER in Theory')
% xlabel('SNR in dB')
% ylabel('log10(SER) (Symbol Error Rate)')
% title('Uncoded QPSK - SIMO iid Rayleigh fading channel with MRC, Nr=2')
% 
% figure
% gd = -log2(BER)./log2(db2pow(SNRdB));
% plot(SNRdB,gd)
% hold on
% xlabel('SNR in dB')
% ylabel('g_d(SNR)')
% title('Diversity Gain; SIMO MRC')

figure
plot(SNRdB,log10(SNRout))
hold on
SNR = db2pow(SNRdB);
theoretical_SNRout =SNR*(1+pi/4); 
plot(SNRdB,log10(theoretical_SNRout))
xlabel('SNR in dB')
ylabel('SNRout in dB')
title('Array Gain; SIMO MRC')
legend('Experiment','Theory')
xlim([0,20])
ylim([3,300])

SignalPower = 0;
for i = 1:L/2
    SignalPower = SignalPower + Es*(abs(h(1,i))^2+abs(h(2,i))^2);
end

SNRout=mean( (abs(h(1,:))+abs(h(2,:))).^2 )*Es / 2*1;


