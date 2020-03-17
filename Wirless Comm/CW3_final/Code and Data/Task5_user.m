% Xinyuan Xu
% Wireless Commm 2020
% Task 5 - number of user part

clear all
load Codebook.mat
T = 300;     tc = 50;
nr = 2;      nt = 4; 
J = 6;         
Ptx = 1*10^(-3)*10^(46/10); % 39.8107      
Pn = 1*10^(-3)*10^(-174/10); % 3.9811e-21 
t_space = 0.5;      epsilon = 0.85;%

NoUser = [2,10,20,50]; 
Drops = [150,30,15,6];
% each cdf would have 300 points
Rate_average = zeros(50,300,4);

figure
for counter=1:4
    tic
    K = NoUser(counter);
    NoDrop = Drops(counter);
    
    R_avg = zeros(K,T);
    scheduled_user = zeros(T,1);

    for q =1:K %loop through user, initialize small non-zero rate
        R_avg(q,1) = 10^(-12); 
    end
    
    for drop = 1:NoDrop
        [Ai,Aj,Rt]=DropUser(t_space,K);  
        H_tilt = sqrt(0.5)*(randn(nr,nt, K,J+1)+1i*randn(nr,nt, K,J+1));
        H = zeros(nr,nt, K,J+1);
    
        for k = 1:T % increment of time index k
            N = 1/sqrt(2)*(randn(nr,nt,K,J+1)+1i*randn(nr,nt,K,J+1));
            H_tilt = epsilon*H_tilt + sqrt(1-epsilon^2)*N;
        
            for q = 1:K
                for j = 1:J+1
                    H(:,:,q,j) = H_tilt(:,:,q,j)*Rt(:,:,q,j)^(1/2);
                end
            end
        
            for q = 1:K
                CQI(q) = Inst_User_Rate(H(:,:,q,:),Ai(q),Aj(q,:),Pn,Ptx);
                U(q) = CQI(q)/R_avg(q,k);
            end
    
            [M,q_star] = max(U); 
            scheduled_user(k) = q_star;

            for q = 1:K
                R_avg(q,k+1) = (1-1/tc)*R_avg(q,k);
            end
            R_avg(q_star,k+1) = R_avg(q_star,k+1) + 1/tc*CQI(q_star);
        end
    
        Rate_average(1:K,drop,counter) = R_avg(:,T);
    end
    
    DataPoints = reshape(Rate_average(1:K,1:NoDrop,counter),300,1);
    cdfplot(DataPoints)
    hold on
    toc
end

save('Task5_user','Rate_average')
xlabel("x = user average rate (b/s/Hz)")
ylabel("F(x) - CDF")
legend("K=2","K=10","K=20","K=50")
% xlim([0,15])
xlim([0,5])

% Elapsed time is 631.633133 seconds.
% Elapsed time is 674.061636 seconds.
% Elapsed time is 605.744193 seconds.
% Elapsed time is 563.001382 seconds.