% testing function: function CQI = Inst_User_Rate(H, Ai,Aj,Pn, Ptx)
clear all
nr = 2; % 1 or 2
nt = 4;
Ptx = 1*10^(-3)*10^(46/10); % 39.8107      
Pn = 1*10^(-3)*10^(-174/10); % 3.9811e-21 ? is it too small?     
K = 10; 
J =6;
t = 0.5; % space correlation
epsilon = 0.85;% time correlation
% rng('default')
[di,dj,Rt]=DropUser(0.5);
A0j_dB = 128.1+37.6*log10(dj); 
A0i_dB = 128.1+37.6*log10(di);
% channel initialization
H_tilt = sqrt(0.5)*(randn(nr,nt, K,J+1)+1i*randn(nr,nt, K,J+1));
H = zeros(nr,nt, K,J+1);

% shadowing
Sj_dB = randn(K,J)*8;   % std of 8dB
Si_dB = randn(K,1)*8; % user cell shadowing     % std of 8dB
% pass loss + shadowing
Aj = db2pow(A0j_dB+Sj_dB);
Ai = db2pow(A0i_dB+Si_dB); 

N = 1/sqrt(2)*(randn(nr,nt,K,J+1)+1i*randn(nr,nt,K,J+1));
H_tilt = epsilon*H_tilt + sqrt(1-epsilon^2)*N;

for k = 1:K
    for j = 1:J+1
        H(:,:,k,j) = H_tilt(:,:,k,j)*Rt(:,:,k,j);
    end
end

CQI_test = Inst_User_Rate(H(:,:,1,:), Ai(1),Aj(1,:),Pn, Ptx)
% CQI of user 1

