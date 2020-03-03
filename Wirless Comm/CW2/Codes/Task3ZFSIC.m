% wirless comm coursework 2, Task 3 SIC
% xx316, Xinyuan Xu, 2020 Feb
clear all;
c00 = 1/sqrt(2)*(+1+1i); c01 = 1/sqrt(2)*(-1+1i); % define QPSK constellation
c11 = 1/sqrt(2)*(-1-1i); c10 = 1/sqrt(2)*(+1-1i);

C = [c00;c01;c11;c10];

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
    for t=1:L/4 % zero forcing with sussesive interference cancellation
        Ht=H(:,:,t);
        Gt = sqrt(nt/Es)*inv(Ht'*Ht)*Ht'; % zero forcing combiner
        
        % unordered SIC
        g1 = Gt(1,:);
        h1 = H(:,1,t);
        z1 = g1*y(:,:,t);
        for index = 1:4 %ML decode c1 first
            error = z1-sqrt(Es)*g1*h1*C(index);
            % C is the codebook
            SquredError(index)=error'*error; 
        end
        [M,I] = min(SquredError);
        c1 = C(I);
        
        y2 = y(:,:,t)-sqrt(Es)*h1*c1;
        % assume perfect interfernce cancellation
        h2 = H(:,2,t);
        for index = 1:4 %ML detection
            error = y2-sqrt(Es)*h2*C(index);
            SquredError(index)=error'*error; 
        end
        [M,I] = min(SquredError);
        c2 = C(I);
        
        c_est(:,:,t)=[c1;c2];
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

figure
plot(SNRdB, log10(BER))
title('BER: ZF-SIC receiver')
xlabel('SNR(dB)')
ylabel('log10(BER)')

gd = -log2(BER)./log2(SNR');
figure
plot(SNRdB, gd)
title('ZF-SIC receiver: diversity gain (ignore low SNR part)')
xlabel('SNR(dB)')
ylabel('-log2(BER)/log2(SNR)')
% diversity 1.2 explanable

save('ZF_SIC', 'SNRdB','BER','gd')