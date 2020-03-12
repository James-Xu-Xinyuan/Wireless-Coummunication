% Xinyuan Xu
% Wireless Commm 2020
% Task 4

clear all
load Codebook.mat
NoDrop = 30; 
nr = 2; % receive antenna 
nt = 4; % transmit antenna
K = 10; % number of users in cell
J = 6; % number of interfering cells
Ptx = 1*10^(-3)*10^(46/10); % 39.8107      
Pn = 1*10^(-3)*10^(-174/10); % 3.9811e-21 ? is it too small?     
t_space = 0.5; % space correlation
epsilon = 0.85;% time correlation

% tcs = [2,10,100,1000];
% Ts = [100,100,600,6000];
tcs = [2,10,100,1000];
Ts = [100,100,600,6000];

Rate_average = zeros(K,NoDrop,4);

figure
for counter=1:4
    tic
    tc = tcs(counter);
    T = Ts(counter);
    
    R_avg = zeros(K,T);
    scheduled_user = zeros(T,1);

    for k =1:K
        R_avg(k,1) = 10^(-12); 
    end
    
    for drop = 1:NoDrop
        [Ai,Aj,Rt]=DropUser(t_space);  
        H_tilt = sqrt(0.5)*(randn(nr,nt, K,J+1)+1i*randn(nr,nt, K,J+1));
        H = zeros(nr,nt, K,J+1);
    
        for t = 1:T
            N = 1/sqrt(2)*(randn(nr,nt,K,J+1)+1i*randn(nr,nt,K,J+1));
            H_tilt = epsilon*H_tilt + sqrt(1-epsilon^2)*N;
        
            for k = 1:K
                for j = 1:J+1
                    H(:,:,k,j) = H_tilt(:,:,k,j)*Rt(:,:,k,j)^(1/2);
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

save('Task4','Rate_average')
title("CDF of user average rate")
xlabel("x = user average rate (b/s/Hz)")
ylabel("F(x) - CDF")
legend("tc=2","tc=10","tc=100","tc=1000")

% Elapsed time is 111.107887 seconds.
% Elapsed time is 110.564017 seconds.
% Elapsed time is 220.675700 seconds.
% Elapsed time is 705.974525 seconds.
% Elapsed time is 1326.793080 seconds.