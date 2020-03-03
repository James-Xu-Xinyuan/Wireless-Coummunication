% wirless comm coursework 2, Task 2
% xx316, Xinyuan Xu, 2020 Feb
clear all;
SNRdB= [-10:0.1:20];
SNR = db2pow(SNRdB);

%perfect CSIT
ErgodicCapacity22CSIT = zeros(1, length(SNR)); % ergodic capacity of 2*2 channel
ErgodicCapacity42CSIT = zeros(1, length(SNR));
ErgodicCapacity24CSIT = zeros(1, length(SNR));
for i=1:length(SNR)
    EC22 = 0; % ergodic capacity
    EC42 = 0;
    EC24 = 0;
    size = 10000;
    for n=1:size
        H22 = 1/sqrt(2)*(randn(2,2)+1i*randn(2,2));
        sigma22 = svd(H22);  % singular values of channel matrix
        % the output vector should be already sorted
        lemda22= sigma22.^2;
        s22= WaterFilling(lemda22,SNR(i)); % power allocation
        C22 = log2(1+SNR(i)*s22(1)*lemda22(1))+ log2(1+SNR(i)*s22(2)*lemda22(2));
        % capacity of this random channel realization
        EC22 = EC22+C22/size; % average
        
        H42 = 1/sqrt(2)*(randn(4,2)+1i*randn(4,2));
        sigma42 = svd(H42);
        lemda42= sigma42.^2;
        s42= WaterFilling(lemda42,SNR(i)); 
        C42 = log2(1+SNR(i)*s42(1)*lemda42(1))+ log2(1+SNR(i)*s42(2)*lemda42(2));
        EC42 = EC42+C42/size;
        
        H24 = 1/sqrt(2)*(randn(2,4)+1i*randn(2,4));
        sigma24 = svd(H24);
        lemda24= sigma24.^2;
        s24= WaterFilling(lemda24,SNR(i)); 
        C24 = log2(1+SNR(i)*s24(1)*lemda24(1))+ log2(1+SNR(i)*s24(2)*lemda24(2));
        EC24 = EC24+C24/size;
    end
    ErgodicCapacity22CSIT(i)= EC22;
    ErgodicCapacity42CSIT(i)= EC42;
    ErgodicCapacity24CSIT(i)= EC24;
end

% partial CSIT
ErgodicCapacity22CDIT = zeros(1, length(SNR));
ErgodicCapacity42CDIT = zeros(1, length(SNR));
ErgodicCapacity24CDIT = zeros(1, length(SNR));
for i=1:length(SNR)
    EC22 = 0; 
    EC42 = 0;
    EC24 = 0;
    size = 10000;
    for n=1:size
        H22 = 1/sqrt(2)*(randn(2,2)+1i*randn(2,2));
        sigma22 = svd(H22);  
        lemda22= sigma22.^2;
        nt = 2; % number of transmit antennas
        C22 = log2(1+SNR(i)/nt*lemda22(1))+ log2(1+SNR(i)/nt*lemda22(2));
        EC22 = EC22+C22/size; 
        
        H24 = 1/sqrt(2)*(randn(2,4)+1i*randn(2,4));
        sigma24 = svd(H24);  
        lemda24= sigma24.^2;
        nt = 4; 
        C24 = log2(1+SNR(i)/nt*lemda24(1))+ log2(1+SNR(i)/nt*lemda24(2));
        EC24 = EC24+C24/size;
        
        H42 = 1/sqrt(2)*(randn(4,2)+1i*randn(4,2));
        sigma42 = svd(H42);  
        lemda42= sigma42.^2;
        nt = 2; 
        C42 = log2(1+SNR(i)/nt*lemda42(1))+ log2(1+SNR(i)/nt*lemda42(2));
        EC42 = EC42+C42/size;
    end
    ErgodicCapacity22CDIT(i)= EC22;
    ErgodicCapacity42CDIT(i)= EC42;
    ErgodicCapacity24CDIT(i)= EC24;

end

plot(SNRdB, ErgodicCapacity22CSIT,'b')
hold on
plot(SNRdB, ErgodicCapacity42CSIT,'r')
plot(SNRdB, ErgodicCapacity24CSIT,'k')
plot(SNRdB, ErgodicCapacity22CDIT,'b--')
plot(SNRdB, ErgodicCapacity42CDIT,'r--')
plot(SNRdB, ErgodicCapacity24CDIT,'k--')
legend('2x2CSIT','4x2CSIT','2x4CSIT','2x2CDIT','4x2CDIT','2x4CDIT')
xlabel('SNR(dB)','FontSize',18)
ylabel('Ergodic Capacity: bits/s/Hz','FontSize',18)
title('Capacity of MIMO channel with perfect CSIT or CDIT','FontSize',22)

