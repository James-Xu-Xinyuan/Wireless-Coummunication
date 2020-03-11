function [Ai,Aj,Rt]=DropUser(t,K)
% [Ai,Aj]=DropUser(10,6)
% pass loss + shadowing
% di: K x 1 vector
% dj: K x J matrix
% K : number of users in cells
% J : number of interfering cells
% t : % baseline spatial correlation magnitude
    if nargin < 2
        K = 10;
    end
    J = 6;
%     t = 0.5; % baseline spatial correlation magnitude
%     rng('default')
    
    % user coordinate
    theta_i = rand(K,1)*2*pi;
    di = (rand(K,1)*(250-35)+35)/1000;
    x = di.*cos(theta_i);
    y = di.*sin(theta_i);
    
    % interfering base station location
    theta_j = [0,1,2,3,4,5]*(2*pi/6);
    x_j = 500/1000*cos(theta_j);
    y_j = 500/1000*sin(theta_j);
    dj = zeros(K,J);
    
    for q = 1:K
        for j = 1:J
            dj(q,j) = sqrt((x_j(j)-x(q))^2+(y_j(j)-y(q))^2);
        end
    end
    
    % transmit correlation matrix
    Rt = zeros(4,4,K,J+1);
    Phi = 2*pi*rand(K,1);
    t_q0 = t*exp(1i*Phi);
    for k = 1:K
        % user and its cell
        tqi = t_q0(k);
        Rt(1,:,k,1) = [ 1 tqi tqi^2 tqi^3];
        Rt(2,:,k,1) = [conj(tqi) 1 tqi tqi^2];
        Rt(3,:,k,1) = [conj(tqi)^2 conj(tqi) 1 tqi];
        Rt(4,:,k,1) = [conj(tqi)^3 conj(tqi)^2 conj(tqi) 1];
        
        % assume no space transmit correlation 
        % between user and interfering cells
        for i = 2:J+1 % 2:7
            Rt(:,:,k,i) = eye(4,4); 
        end
    end
    
        
    % pretend A doesnt have the horizontal bar - path loss
    A0j_dB = 128.1+37.6*log10(dj); 
    Sj_dB = randn(K,J)*8;   % std of 8dB
    Aj = db2pow(A0j_dB+Sj_dB);% pass loss + shadowing
    
    A0i_dB = 128.1+37.6*log10(di);
    Si_dB = randn(K,1)*8; % user cell shadowing     % std of 8dB
    Ai = db2pow(A0i_dB+Si_dB); % pass loss + shadowing

