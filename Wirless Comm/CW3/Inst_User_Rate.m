function R = Inst_User_Rate(H, SNR,Ai,Aj,Pn)
% calculate instantaenous rate of one user
load Codebook.mat
% coursework specific, no inter-user interference
[RI,PMI] = QuantizePrecoding(H, SNR);

rng('default')
J = length(Aj);
P_cell = zeros(4,2,J);
for j = 1:J % random interference precoding
    ri = round(rand+1);
    pmi = round(15*rand);
    
    if ri == 1
        P_cell(j) = W1(:,:,pmi+1);
    else
        P_cell(j) = W2(:,:,pmi+1);
    end

end

% H_{q,j} : inter-cell link
[nr,nt]=size(H)
H_cell = 1/sqrt(2)*(randn(nr,nt,J)+1i*randn(nr,nt,J));

if RI == 1
    P = W1(:,:,PMI+1); % matlab counts 1-16, standard counst 0-15
    % no interstream interference: there is only 1 stream!
else
    P = W2(:,:,PMI+1); 
end

Rni = % to be continued