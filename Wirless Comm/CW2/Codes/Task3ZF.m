% wirless comm coursework 2 Task 3 ZF
% xx316, Xinyuan Xu, 2020 Feb
clear all;
c00 = 1/sqrt(2)*(+1+1i); c01 = 1/sqrt(2)*(-1+1i); % define QPSK constellation
c11 = 1/sqrt(2)*(-1-1i); c10 = 1/sqrt(2)*(+1-1i);

L = 3*10^5;
BitStream = round(rand(L,1)); 
nt =2; nr=2;
SymbolStream = zeros(L/2,1);
for j =1:length(SymbolStream)
    if BitStream(2*j-1)==0
        if BitStream(2*j)==0
            SymbolStream(j)=c00;
        else
            SymbolStream(j)=c01;
        end
    else
        if BitStream(2*j)==0
            SymbolStream(j)=c10;
        else
            SymbolStream(j)=c11;
        end
    end
end
c =  reshape(SymbolStream,nr,1,length(SymbolStream)/nr);

%assume normalized noise: variance of n is 1
%average SNR = Es/var(n) = Es numerically
SNRdB = 0:0.2:20;
SNR = db2pow(SNRdB);
BER = zeros(length(SNR),1);

for i=1:length(SNR)
    Es = SNR(i);
    H = 1/(sqrt(2))*(randn(nr,nt,L/4)+1i*randn(nr,nt,L/4));
    n = 1/(sqrt(2))*(randn(nr,1,L/4)+1i*randn(nr,1,L/4));
    c_est = zeros(size(c));
    
    y = zeros(nr,1,L/4);
    for t =1:L/4
        y(:,:,t)= sqrt(Es)*H(:,:,t)*c(:,:,t)+n(:,:,t);
        %c is the symbol stream
    end
    
    z = zeros(size(y));
    for t=1:L/4 % zero forcing
        Ht=H(:,:,t);
        Gt = sqrt(nt/Es)*inv(Ht'*Ht)*Ht'; % zero forcing combiner
        z(:,:,t)=Gt*y(:,:,t); % decoupled channels and symbols
    end
    
    z = reshape(z,length(SymbolStream),1,1);
    
    ReceivedBits = zeros(length(BitStream),1);
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
    BER(i) = ber;
end    

plot(SNRdB, log10(BER))
title('BER: ZF receiver')
xlabel('SNR(dB)')
ylabel('log10(BER)')

gd = -log2(BER)./log2(SNR');
figure
plot(SNRdB, gd)
title('ZF receiver: diversity gain (ignore low SNR part)')
xlabel('SNR(dB)')
ylabel('-log2(BER)/log2(SNR)')

save('ZF', 'SNRdB','BER','gd')