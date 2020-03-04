function answer=SINR_LT(Ai,Aj,Es, Pn)
% this function is coursework specific:
% assume tranmit power is the same for all cells (Es)
% Ai: path loss + shadowing of (central) cell i
% Aj: path loss + shadowing of interfering cells (vector)
% Pn = sigma^2 of noise, of user q

% Inter-cell interference term
Ic = 0;
for i = 1:length(Aj)
    Ic = Ic + Es/Aj(i);
end
answer = (Es/Ai)/(Pn + Ic);

