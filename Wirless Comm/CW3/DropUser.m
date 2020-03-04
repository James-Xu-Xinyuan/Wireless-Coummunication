function [Ai,Aj]=DropUser(K,J)
% [Ai,Aj]=DropUser(10,6)
% pass loss + shadowing
% Ai: K x 1 vector
% AJ: K x J matrix
% K : number of users in cells
% J : number of interfering cells
    rng('default')
    % user coordinate
    theta_i = rand(K,1)*2*pi;
    ri = (rand(K,1)*(250-35)+35)/1000;
    x = ri.*cos(theta_i);
    y = ri.*sin(theta_i);
    
    % interfering base station location
    theta_j = [0,1,2,3,4,5]*(2*pi/6);
    x_j = 500/1000*cos(theta_j);
    y_j = 500/1000*sin(theta_j);
    d = zeros(10,6);
    
    for q = 1:K
        for j = 1:J
            d(q,j) = sqrt((x_j(j)-x(q))^2+(y_j(j)-y(q))^2);
        end
    end
        
    % pretend A doesnt have the horizontal bar - path loss
    A0j_dB = 128.1+37.6*log10(d); 
    Sj_dB = randn(K,J)*8;   % std of 8dB
    Aj = db2pow(A0j_dB+Sj_dB);% pass loss + shadowing
    
    A0i_dB = 128.1+37.6*log10(ri);
    Si_dB = randn(K,1)*8; % user cell shadowing     % std of 8dB
    Ai = db2pow(A0i_dB+Si_dB); % pass loss + shadowing