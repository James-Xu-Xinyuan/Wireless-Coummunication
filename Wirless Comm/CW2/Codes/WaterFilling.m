% wirless comm coursework 2
% xx316, Xinyuan Xu, 2020 Feb

function s = WaterFilling(lemda, SNR)
%input: array of lemdas 
%output: array of optimal power allocation
lemda = sort(lemda, 'descend');
n = length(lemda);
s = zeros(1, n);
loop = 1;

i =1;
while(loop==1)
    loop = 0;
    sum=0;
    for k =1:n-i+1
        sum = sum +1/(SNR*lemda(k));
    end
    u = 1/(n-i+1)*(1+sum);
    for k = 1:n-i+1
        s(k)= u-1/(SNR*lemda(k));
        if s(k)<0
            loop =1;
            s(k)=0;
        end
    end
    
    if loop==1
        i=i+1;
    end
end