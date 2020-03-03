% wirless comm coursework 1
% xx316, Xinyuan Xu, 2020 Jan
% uncoded QPSK transmission over a MISO i.i.d Rayleigh fading channel
% with MRT and two transmit antennas
clear all;
%define constellation diagram
c00 = 1/sqrt(2)*(+1+1i); 
c01 = 1/sqrt(2)*(-1+1i); 
c11 = 1/sqrt(2)*(-1-1i); 
c10 = 1/sqrt(2)*(+1-1i);

L = 10^6; 
BitStream = round(rand(L,1)); 
% a random bit stream of 0 and 1 
%modulation
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
y = zeros(L/2,1);
SNRout = zeros(1,201);
%iterate for each SNR value
for j = 1:201
    Es = db2pow(SNRdB(j)); 
    % generate iid channel and noise
    h = 1/(sqrt(2))*(randn(L/2,2)+1i*randn(L/2,2));
    n = 1/(sqrt(2))*(randn(L/2,1)+1i*randn(L/2,1));
    h_mag = zeros(1,L/2);
    for i=1:L/2
        hi = h(i,:);
         %find magnitude of channel for later SNRout calculation
        h_mag(i)=sqrt(hi*hi');
        % precoder; maximal ratio transmission
        wi = hi'/h_mag(i);
        y(i) = sqrt(Es)*hi*wi*c(i)+n(i);
    end
   
    ReceivedBits = zeros(L,1);
    % demodulation
    for k = 1:length(y)
        if imag(y(k))>0
            ReceivedBits(k*2-1) = 0;
        else
            ReceivedBits(k*2-1) = 1;
        end
        if real(y(k))>0
            ReceivedBits(2*k) = 0;
        else
            ReceivedBits(2*k) = 1;
        end
    end
    %calculate bit error rate
    diff = mod(BitStream+ReceivedBits,2);
    ber = sum(diff)/L;
    BER(j) = ber;
    %calcuate output SNR
    SNRout(j)=Es*mean(h_mag.^2)/1;
end

figure
plot(SNRdB, log10(BER))
xlabel('SNR in dB')
ylabel('log10(BER) (Bit Error Rate)')
title('Uncoded QPSK - MISO i.i.d Rayleigh fading channel, MRT, Nt=2')
save('MISO_MRT', 'BER','SNRout')

figure
gd = -log2(BER)./log2(db2pow(SNRdB));
plot(SNRdB,gd)
hold on
xlabel('SNR in dB')
ylabel('g_d(SNR)')
title('Diversity Gain; MISO MRC')

figure
SNR = db2pow(SNRdB);
plot(SNRdB,10*log10(SNRout))
hold on
Nt=2;
SNRout_theoretical = SNR*Nt;
plot(SNRdB,10*log10(SNRout_theoretical))
legend('Experiment','Theory')
xlabel('SNR in dB')
ylabel('SNRout in dB')
title('SNRout; MISO MRT')

figure
ArrayGain = SNRout./SNR;
ArrayGain_theoretical = ones(201,1)*Nt;
plot(SNRdB,ArrayGain)
hold on
plot(SNRdB,ArrayGain_theoretical)
ylim([1.99,2.01])
xlabel('SNR in dB')
ylabel('Array Gain')
legend('Experiment','Theory')
title('Array Gain; MISO MRT')
