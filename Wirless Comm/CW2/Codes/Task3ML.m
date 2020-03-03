% wirless comm coursework 2 Task 3 ML
% xx316, Xinyuan Xu, 2020 Feb
clear all;
c00 = 1/sqrt(2)*(+1+1i); c01 = 1/sqrt(2)*(-1+1i); % define QPSK constellation
c11 = 1/sqrt(2)*(-1-1i); c10 = 1/sqrt(2)*(+1-1i);

C(:,:,1)=[c00;c00]; C(:,:,2)=[c00;c01];
C(:,:,3)=[c00;c11]; C(:,:,4)=[c00;c10];
C(:,:,5)=[c01;c00]; C(:,:,6)=[c01;c01];
C(:,:,7)=[c01;c11]; C(:,:,8)=[c01;c10];
C(:,:,9)=[c11;c00]; C(:,:,10)=[c11;c01];
C(:,:,11)=[c11;c11]; C(:,:,12)=[c11;c10];
C(:,:,13)=[c10;c00]; C(:,:,14)=[c10;c01];
C(:,:,15)=[c10;c11]; C(:,:,16)=[c10;c10];

L = 10^6;
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
c = reshape(SymbolStream,nr,1,length(SymbolStream)/nr);

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
    
    for t=1:L/4 %Maximum Likelihood detection
        SquredError = zeros(16,1);
        for index = 1:16
            error = y(:,:,t)-sqrt(Es)*H(:,:,t)*C(:,:,index);
            % C is the codebook
            SquredError(index)=error'*error; 
        end
        [M,I] = min(SquredError);
        c_est(:,:,t)=C(:,:,I);
    end
    c_est = reshape(c_est,length(SymbolStream),1,1);
    
    ReceivedBits = zeros(length(BitStream),1);
    for k = 1:length(c_est)
        if imag(c_est(k))>0
            ReceivedBits(k*2-1) = 0;
        else
            ReceivedBits(k*2-1) = 1;
        end
        if real(c_est(k))>0
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
title('BER: ML receiver')
xlabel('SNR(dB)')
ylabel('log10(BER)')

gd = -log2(BER)./log2(SNR');
figure
plot(SNRdB, gd)
title('ML receiver: diversity gain (ignore low SNR part)')
xlabel('SNR(dB)')
ylabel('-log2(BER)/log2(SNR)')

save('ML', 'SNRdB','BER','gd')
