% testing proportional fair scheduling
% main objective: to find out how long the drop needs to be
clear all
nr = 2; % 1 or 2
nt = 4;
Ptx = 1*10^(-3)*10^(46/10); % 39.8107      
Pn = 1*10^(-3)*10^(-174/10);    
K = 10; 
J =6;
T = 100;  
tc = 2; 
t = 0.5; % space correlation
epsilon = 0.85;% time correlation

R_avg = zeros(K,T);
scheduled_user = zeros(T,1);
R_test = zeros(T+1,1) ;

for k =1:K
    R_avg(k,1) = 10^(-12); % initialized with a very small rate to avoid divide by 0
end

[Ai,Aj,Rt]=DropUser(t);  
H_tilt = sqrt(0.5)*(randn(nr,nt, K,J+1)+1i*randn(nr,nt, K,J+1));
H = zeros(nr,nt, K,J+1);

for t = 1:T
    CQI = zeros(K,1);
    U = zeros(K,1);% PF utility metric

    N = 1/sqrt(2)*(randn(nr,nt,K,J+1)+1i*randn(nr,nt,K,J+1));
    H_tilt = epsilon*H_tilt + sqrt(1-epsilon^2)*N;

    for k = 1:K
        for j = 1:J+1
            H(:,:,k,j) = H_tilt(:,:,k,j)*Rt(:,:,k,j)^(1/2);
        end
    end
    
    for k = 1:K
        % CQI_test = Inst_User_Rate(H(:,:,1,:), Ai(1),Aj(1,:),Pn, Ptx)
        % CQI of user 1
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
figure
for k = 1:K
    plot(R_avg(k,:))
    hold on
end
% Tested if tc = 50, nt = 1 or 2, T=300 is enough
% (tc, T): (10,60) (100,600) (500, 3000)

    