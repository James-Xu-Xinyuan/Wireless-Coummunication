% wirless comm coursework 2, Task 3 graph
% xx316, Xinyuan Xu, 2020 Feb
% run the 3 other task 3 scripts before using this one
clear all
load ML.mat 
BER_ML = BER;
gd_ML = gd;
load ZF.mat 
BER_ZF = BER;
gd_ZF = gd;
load ZF_SIC.mat 
BER_SIC = BER;
gd_SIC = gd;

SNRdB = SNRdB';

figure
plot(SNRdB, log10(BER_ML))
hold on
plot(SNRdB, log10(BER_ZF))
plot(SNRdB, log10(BER_SIC))
title('Bit Error Rate of 2x2 MIMO i.i.d Rayleigh fading channel')
xlabel('SNR(dB)')
ylabel('log10(BER)')
legend('ML','ZF','ZF-SIC')

figure
plot(SNRdB,gd_ML)
hold on
plot(SNRdB,gd_ZF)
plot(SNRdB,gd_SIC)
title('Diversity Gain of 2x2 MIMO i.i.d Rayleigh fading channel')
xlabel('SNR(dB)')
ylabel('-log2(BER)/log2(SNR)')
ylim([0.5,2.5])
legend('ML','ZF','ZF-SIC')