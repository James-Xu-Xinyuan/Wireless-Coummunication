clear all
theta_i = rand(1,1)*2*pi;

ri = (rand*(250-35)+35)/1000;

x = ri*cos(theta_i);
y = ri*sin(theta_i);

di = ri; %in km

theta_j = [0,1,2,3,4,5]*(2*pi/6);

x_j = 500/1000*cos(theta_j);
y_j = 500/1000*sin(theta_j);

dj = zeros(6,1);
for j = 1:6
    dj(j) = sqrt((x_j(j)-x)^2+(y_j(j)-y)^2);
end