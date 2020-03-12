clear all
t_space = 0.5;
[Ai,Aj,Rt]=DropUser(t_space);  
% norm(Rt(:,:,1,1))
% 1     : 1
% 0.5   : 2.0856
% 0.9   : 3.5266
% 0.999 : 3.9950
Rt(:,:,1,1)=Rt(:,:,1,1)/norm(Rt(:,:,1,1));
norm(Rt(:,:,1,1))