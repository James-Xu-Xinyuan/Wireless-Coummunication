clear all
U = zeros(4,16);
U(:,1) = [1;-1;-1;-1];
U(:,2) = [1;-j;1;j];
U(:,3) = [1;1;-1;1];
U(:,4) = [1;j;1;-j];
U(:,5) = [1;(-1-j)/sqrt(2);-j;(1-j)/sqrt(2)];
U(:,6) = [1;(1-j)/sqrt(2);j;(-1-j)/sqrt(2)];
U(:,7) = [1;(1+j)/sqrt(2);-j;(-1+j)/sqrt(2)];
U(:,8) = [1;(-1+j)/sqrt(2);j;(1+j)/sqrt(2)];
U(:,9) = [1;-1;1;1];
U(:,10) = [1;-j;-1;-j];
U(:,11) = [1;1;1;-1];
U(:,12) = [1;j;-1;j];
U(:,13) = [1;-1;-1;1];
U(:,14) = [1;-1;1;-1];
U(:,15) = [1;1;-1;-1];
U(:,16) = [1;1;1;1];

W = zeros(4,4,16);
W1 = zeros(4,1,16); % ne = 1 = NoLayer
for n = 1:16
    u = U(:,n);
    W(:,:,n) = eye(4,4)-2*u*u'/(u'*u);
    W1(:,1,n) =  W(:,1,n);
end

W2 = zeros(4,2,16); % ne = 2
W2(:,:,1)= [W(:,1,1) W(:,4,1)]/sqrt(2);
W2(:,:,2)= [W(:,1,2) W(:,2,2)]/sqrt(2);
W2(:,:,3)= [W(:,1,3) W(:,2,3)]/sqrt(2);
W2(:,:,4)= [W(:,1,4) W(:,2,4)]/sqrt(2);
W2(:,:,5)= [W(:,1,5) W(:,4,5)]/sqrt(2);
W2(:,:,6)= [W(:,1,6) W(:,4,6)]/sqrt(2);
W2(:,:,7)= [W(:,1,7) W(:,3,7)]/sqrt(2);
W2(:,:,8)= [W(:,1,8) W(:,3,8)]/sqrt(2);
W2(:,:,9)= [W(:,1,9) W(:,2,9)]/sqrt(2);
W2(:,:,10)= [W(:,1,10) W(:,4,10)]/sqrt(2);
W2(:,:,11)= [W(:,1,11) W(:,3,11)]/sqrt(2);
W2(:,:,12)= [W(:,1,12) W(:,3,12)]/sqrt(2);
W2(:,:,13)= [W(:,1,13) W(:,2,13)]/sqrt(2);
W2(:,:,14)= [W(:,1,14) W(:,3,14)]/sqrt(2);
W2(:,:,15)= [W(:,1,15) W(:,3,15)]/sqrt(2);
W2(:,:,16)= [W(:,1,16) W(:,2,16)]/sqrt(2);

save('Codebook','W1','W2')

% norm(W2(:,1,1))^2+norm(W2(:,2,1))^2 = 1
% norm(W1(:,:,1))^2 = 1
% norm(W2(:,:,1))^2 = 0.5
