% Xinyuan Xu
% Wireless Commm 2020
% Task 3

clear all
NoDrop = 1000; % test 100 to 10k
T = 10^6; % number of time instances per drop
nr = 2; % receive antenna = 1 or 2
nt = 4; % transmit antenna
K = 10; % number of users in cell
J = 6; % number of interfering cells
Es = 1*10^(-3)*10^(46/10);      
Pn = 1*10^(-3)*10^(-174/10);    
SNR = Es/Pn;

for drop = 1:NoDrop
    [Ai,Aj]=DropUser(K,J);    
    
    for t = 1:T
        % specifically in this task: 
        % assume no time correlation and no space correlation
        H = 1/sqrt(2)*(randn(nr,nt,K)+1i*randn(nr,nt,K));
        
    end
      
end

