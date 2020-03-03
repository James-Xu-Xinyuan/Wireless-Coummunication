% next time L = 10^8 and SNR from 20 to 40+

clear all;

c00 = 1/sqrt(2)*(+1+1i); 
c01 = 1/sqrt(2)*(-1+1i); 
c11 = 1/sqrt(2)*(-1-1i); 
c10 = 1/sqrt(2)*(+1-1i);

L = 10^9; 
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
BER = zeros(1,16);
SNRdB = 30:45;
SNRout = zeros(1,16);
for j = 1:16
    Es = db2pow(SNRdB(j)); 
    h = 1/(sqrt(2))*(randn(L/4,2)+1i*randn(L/4,2));
    n = 1/(sqrt(2))*(randn(L/2,1)+1i*randn(L/2,1));
    h_mag = zeros(1,L/4);
    y = zeros(L/2,1);
    for i=1:L/4
        hi = h(i,:);
        
        y(2*i-1)=sqrt(Es)*hi(1)*c(2*i-1)/sqrt(2)+sqrt(Es)*hi(2)*c(2*i)/sqrt(2)+n(2*i-1);
        y(2*i)=-sqrt(Es)*hi(1)*conj(c(2*i))/sqrt(2)+sqrt(Es)*hi(2)*conj(c(2*i-1))/sqrt(2)+n(2*i);
        
        yi = [y(2*i-1);conj(y(2*i))];
        Heff = [hi(1),hi(2);conj(hi(2)),-conj(hi(1))];
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

    diff = mod(BitStream+ReceivedBits,2);
    ber = sum(diff)/L;
    BER(j) = ber;
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

%'Out of memory.'




