% wirless comm coursework 1
% xx316, Xinyuan Xu, 2020 Jan
% uncoded QPSK transmission over a MISO i.i.d Rayleigh fading channel
% with Space-Time Coding and two transmit antennas
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
%iterate for each SNR value
for j = 1:201
    Es = db2pow(SNRdB(j)); 
    % generate iid channel and noise
    h = 1/(sqrt(2))*(randn(L/4,2)+1i*randn(L/4,2));
    n = 1/(sqrt(2))*(randn(L/2,1)+1i*randn(L/2,1));
    h_mag = zeros(1,L/4);
    y = zeros(L/2,1);
    for i=1:L/4
        hi = h(i,:);
        h_mag(i)=sqrt(hi*hi');
        
        %space time coding
        y(2*i-1)=sqrt(Es)*hi(1)*c(2*i-1)/sqrt(2)+sqrt(Es)*hi(2)*c(2*i)/sqrt(2)+n(2*i-1);
        y(2*i)=-sqrt(Es)*hi(1)*conj(c(2*i))/sqrt(2)+sqrt(Es)*hi(2)*conj(c(2*i-1))/sqrt(2)+n(2*i);
        yi = [y(2*i-1);conj(y(2*i))];
        
        Heff = [hi(1),hi(2);conj(hi(2)),-conj(hi(1))];
        % assume channel is perfectly know at receiver 
        % and using maximum likelihood detection
        zi = Heff'*yi;
        z(2*i-1) = zi(1);
        z(2*i)   = zi(2);
    end
   
    ReceivedBits = zeros(L,1);
    % demodulation
    for k = 1:length(z)
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
    %calculate bit error rate
    diff = mod(BitStream+ReceivedBits,2);
    ber = sum(diff)/L;
    BER(j) = ber;
    %calcuate output SNR
    SNRout(j)=Es/1*mean(((h_mag.^2).^2)./(2*h_mag.^2));
end

figure
plot(SNRdB, log10(BER))
xlabel('SNR in dB')
ylabel('log10(BER) (Bit Error Rate)')
title('Uncoded QPSK - MISO i.i.d Rayleigh fading channel, Alamouti, Nt=2')

figure % diversity is slightly low?
gd = -log2(BER)./log2(db2pow(SNRdB));
plot(SNRdB,gd)
hold on
ref = 2*ones(1,201);
plot(SNRdB,ref)
xlabel('SNR in dB')
ylabel('g_d(SNR)')
legend('Experiment','y=2')
title('Diversity Gain; MISO Alamouti')
ylim([1.5, 3])

figure
SNR = db2pow(SNRdB);
plot(SNRdB,10*log10(SNRout))
hold on
Nt=2;
SNRout_theoretical = SNR;
plot(SNRdB,10*log10(SNRout_theoretical))
legend('Experiment','Theory')
xlabel('SNR in dB')
ylabel('SNRout in dB')
title('SNRout; MISO Alamouti')

figure
ArrayGain = SNRout./SNR;
ArrayGain_theoretical = ones(201,1);
plot(SNRdB,ArrayGain)
hold on
plot(SNRdB,ArrayGain_theoretical)
ylim([0.99,1.01])
xlabel('SNR in dB')
ylabel('Array Gain')
legend('Experiment','Theory')
title('Array Gain; MISO Alamouti')

save('Alamouti', 'BER', 'SNRout')
