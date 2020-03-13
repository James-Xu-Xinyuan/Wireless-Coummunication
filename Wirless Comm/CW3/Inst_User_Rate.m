function CQI = Inst_User_Rate(H, Ai,Aj,Pn, Ptx)
    % Ai: a scalar
    % Aj: a 6x1 vector
    % H(:,:,user_q,7) i+6j
    % calculate instantaenous rate of one user and feedback as CQI
    % coursework specific, no inter-user interference
    load Codebook.mat
    % [RI,PMI] = QuantizePrecoding(H, SNR);
    [nr,nt,K_user,J_user] = size(H);
%     rng('default')
    P_cell = zeros(4,2,7);
    for j = 2:7 % random interference precoding
        ri = round(rand+1);
        pmi = round(15*rand);    
        if ri == 1
            P_cell(:,1,j) = W1(:,:,pmi+1)*Ptx^(0.5);;
        else
            P_cell(:,:,j) = W2(:,:,pmi+1)*[Ptx 0 ; 0 Ptx]^(0.5);
            % precoder already divided by 1/sqrt(2) in codebook
        end
    end
    
    % P1 = W1(:,:,1)*Ptx^(0.5)
    % trace(P1'*P1)
    % P2 = W2(:,:,1)*[Ptx 0 ; 0 Ptx]^(0.5)
    % trace(P2'*P2) = 39.8107 = Ptx

    Rate = zeros(2,16);
    % try all the precoder
    for RI = 1:2
        for PMI = 0:15
            Hqi = H(:,:,1);
            Rni(:,:,1) = Pn*eye(nr,nr); % noise term
            Rni(:,:,2) = Pn*eye(nr,nr);
            if RI == 1  % no interstream interference: there is only 1 stream!
                P = W1(:,:,PMI+1)*Ptx^(0.5); 
                % matlab counts 1-16, standard counst 0-15
            else
                P = W2(:,:,PMI+1)*[Ptx 0 ; 0 Ptx]^(0.5); 
                % coursework specific, assume at max 2 streams
                % so only add 1 interference stream
                Rni(:,:,1) = Rni(:,:,1) + 1/Ai*Hqi*P(:,2)*(Hqi*P(:,2))';
                Rni(:,:,2) = Rni(:,:,2) + 1/Ai*Hqi*P(:,1)*(Hqi*P(:,1))';
            end

            % plus inter-cell interference term
            tmp = zeros(nr,nr);
            for j = 2:7
                tmp = tmp + 1/Aj(j-1).*H(:,:,j)*P_cell(:,:,j)*(H(:,:,j)*P_cell(:,:,j))';
            end
            Rni(:,:,1) = Rni(:,:,1) + tmp;
            Rni(:,:,2) = Rni(:,:,2) + tmp;
            % end of Rni

            % MMSE combiner 
            g(1,:) = Ai^(-1/2)*(Hqi*P(:,1))'*inv(Rni(:,:,1));
            if RI == 1
                g(2,:) = zeros(size(g(1,:)));
            else
                g(2,:) = Ai^(-1/2)*(Hqi*P(:,2))'*inv(Rni(:,:,2));
%                 Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  3.783505e-18. 
            end

            % SINR
            SINR(1) = 1/Ai*norm(g(1,:)*Hqi*P(:,1))^2/(g(1,:)*Rni(:,:,1)*g(1,:)');
            if RI ==1
                SINR(2) = 0;
            else
                SINR(2) = 1/Ai*norm(g(2,:)*Hqi*P(:,2))^2/(g(2,:)*Rni(:,:,2)*g(2,:)');
            end

            Rate(RI,PMI+1) = real(log2(1+SINR(1))+log2(1+SINR(2)));

        end
    end
    
    CQI = max(max(Rate));













