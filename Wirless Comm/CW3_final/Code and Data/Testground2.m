clear all
H = [1 1 1 1; 1 1 1 1];
H = [1 -1 1 -1; -1 1 -1 1];
H = [1 -1 1 -1; j j j j];
SNR = 10;
% SNR = 10^22;
load Codebook.mat
R_1 = zeros(16,1);
R_2 = zeros(16,1);
for i=1:16
    R_1(i) = real(log2(det(eye(1,1)+SNR/1*W1(:,:,i)'*H'*H*W1(:,:,i))));
    R_2(i) = real(log2(det(eye(2,2)+SNR/2*W2(:,:,i)'*H'*H*W2(:,:,i))));
end

[M1,PMI_1] = max(R_1);
[M2,PMI_2] = max(R_2);

if M1 > M2
    RI = 1;
    PMI = PMI_1-1;% matlab counts 1-16, standard counst 0-15
else
    RI = 2;
    PMI = PMI_2-1;
end

