% Xinyuan Xu
% Wireless Commm 2020
% Task 5

clear all
load Codebook.mat
NoDrop = 30; 
T = 600; % number of time instances per drop
% tested in testground 4; long enough for this task
nr = 2; % receive antenna 
nt = 4; % transmit antenna
K = 10; % number of users in cell
J = 6; % number of interfering cells
Ptx = 1*10^(-3)*10^(46/10); % 39.8107      
Pn = 1*10^(-3)*10^(-174/10); % 3.9811e-21 ? is it too small?     
t_space = 0.5; % space correlation
tc = 100;

Epsilons = [0.02,0.18,0.45,0.72,0.98]; % time correlation

Rate_average = zeros(K,NoDrop,5);

figure
for counter=1:5
    tic
    epsilon = Epsilons(counter);
    
    R_avg = zeros(K,T);
    scheduled_user = zeros(T,1);
%     R_test = zeros(T+1,1) ;

    for k =1:K
        R_avg(k,1) = 10^(-12); 
    end
    
    for drop = 1:NoDrop
        [di,dj,Rt]=DropUser(t_space);  
        A0j_dB = 128.1+37.6*log10(dj); 
        A0i_dB = 128.1+37.6*log10(di);
        H_tilt = sqrt(0.5)*(randn(nr,nt, K,J+1)+1i*randn(nr,nt, K,J+1));
        H = zeros(nr,nt, K,J+1);
    
        for t = 1:T
            Sj_dB = randn(K,J)*8;   
            Si_dB = randn(K,1)*8; 
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
    
        Rate_average(:,drop,counter) = R_avg(:,T);
    end
    
    DataPoints = reshape(Rate_average(:,:,counter),K*NoDrop,1);
    cdfplot(DataPoints)
    hold on
    toc
end

save('Task5','Rate_average')
title("CDF of user average rate")
xlabel("x = user average rate (b/s/Hz)")
ylabel("F(x) - CDF")
legend("\epsilon=0.02","\epsilon=0.18","\epsilon=0.45","\epsilon=0.72","\epsilon=0.98")

% Elapsed time is 752.454472 seconds.
% Elapsed time is 679.568915 seconds.
% Elapsed time is 706.312207 seconds.
% Elapsed time is 708.883365 seconds.
% Elapsed time is 752.681832 seconds.

%  warning ('off','all');

