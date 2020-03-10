% Xinyuan Xu
% Wireless Commm 2020
% Task 3

clear all
load Codebook.mat
NoDrop = 30; % 100 is too much; took 8h to run with T=6000
T = 600; % number of time instances per drop
% tested in testground 4; long enough for this task
nr = 2; % receive antenna = 1 or 2
nt = 4; % transmit antenna
K = 10; % number of users in cell
J = 6; % number of interfering cells
Ptx = 1*10^(-3)*10^(46/10); % 39.8107      
Pn = 1*10^(-3)*10^(-174/10); % 3.9811e-21 ? is it too small?     
% SNR = Ptx/Pn;
t_space = 0.5; % space correlation
epsilon = 0.85;% time correlation
tc = 100;

Rate_average = zeros(K,NoDrop,2);
R_avg = zeros(K,T);
scheduled_user = zeros(T,1);
R_test = zeros(T+1,1) ;

for k =1:K
    R_avg(k,1) = 10^(-12); % initialized with a very small rate to avoid divide by 0
end

figure
for nr=1:2
    for drop = 1:NoDrop
        tic
        [di,dj,Rt]=DropUser(t_space);  
        % path loss
        A0j_dB = 128.1+37.6*log10(dj); 
        A0i_dB = 128.1+37.6*log10(di);
        % channel initialization
        H_tilt = sqrt(0.5)*(randn(nr,nt, K,J+1)+1i*randn(nr,nt, K,J+1));
        H = zeros(nr,nt, K,J+1);
    
        for t = 1:T
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
        
            for k = 1:K
                CQI(k) = Inst_User_Rate(H(:,:,k,:),Ai(k),Aj(k,:),Pn,Ptx);
                U(k) = CQI(k)/R_avg(k,t);
            end
    
            [M,q_star] = max(U); 
            scheduled_user(t) = q_star;

            for k = 1:K
                R_avg(k,t+1) = (1-1/tc)*R_avg(k,t);
            end
            R_avg(q_star,t+1) = R_avg(q_star,t+1) + 1/tc*CQI(q_star);
        end
    
        Rate_average(:,drop,nr) = R_avg(:,T);
        toc % for T=600, tc=100, Elapsed time is 22.348310 seconds.
    end
    
    DataPoints = reshape(Rate_average(:,:,nr),K*NoDrop,1);
    cdfplot(DataPoints)
    hold on
end


save('Task3','Rate_average')
title("CDF of user average rate, tc=100 T=600 drop=30")
xlabel("x = user average rate (b/s/Hz)")
ylabel("F(x) - CDF")
legend("nr=1","nr=2")

UAR_Nr1 = reshape(Rate_average(:,:,1),300,1);
UAR_Nr2 = reshape(Rate_average(:,:,2),300,1);
% value in comment: tc=100 T=600 drop=30
Mean_Nr1 = mean(UAR_Nr1) % 2.0257
Std_Nr1 = std(UAR_Nr1) % 0.3401
Mean_Nr2 = mean(UAR_Nr2) % 3.6403
Std_Nr2 = std(UAR_Nr2) % 0.6504