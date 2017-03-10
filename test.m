
clear all; close all; clc;

VFOV= 30*pi;

pd= linspace(0.5,0.9999, 50);
lambda= linspace(1/100, 1/2000, 50);
[PD, LAMBDA]= meshgrid(pd,lambda);

out= LAMBDA .* (1-PD) ./ PD;

figure; hold on; grid on;
surf(PD,LAMBDA,out)
view(3)
xlabel('P_D');
ylabel('lambda')
zlabel('C')





























