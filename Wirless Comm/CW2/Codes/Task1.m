% wirless comm coursework 2, Task1
% xx316, Xinyuan Xu, 2020 Feb
clear all;

H1= [1, 1;1, 1];
[U1,S1,V1] = svd(H1);

H2 = [sqrt(2),0;0,sqrt(2)];
[U2,S2,V2] = svd(H2);

SNRdB = [-10:0.1:20];
SNR = db2pow(SNRdB);
C1= log2(1+SNR*4);
C2= 2*log2(1+SNR);
plot(SNRdB, C1,'b' )
hold on
plot(SNRdB, C2, 'r')
legend('C(H_1)','C(H_2)')
title('Task1, MIMO capacity, perfect CSIT')
xlabel('SNR (dB)')
ylabel('bits/s/Hz')
