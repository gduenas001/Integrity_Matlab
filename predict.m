function predict (Vn,Gn)
%function predict (v,g,Q,WB,dt)
%
% Inputs:
%   v, g - control inputs: velocity and gamma (steer angle)
%   Q - covariance matrix for velocity and gamma
%   WB - vehicle wheelbase
%   dt - timestep
%
% Outputs: 
%   XX, PX - predicted state and covariance (global variables)

global XX PX PARAMS

s= sin(Gn+XX(3)); c= cos(Gn+XX(3));
vts= Vn*PARAMS.dt*s; vtc= Vn*PARAMS.dt*c;

% jacobians   
Gv= [1 0 -vts;
     0 1  vtc;
     0 0 1];
Gu= [PARAMS.dt*c -vts;
     PARAMS.dt*s  vtc;
     PARAMS.dt*sin(Gn)/PARAMS.wheelbase Vn*PARAMS.dt*cos(Gn)/PARAMS.wheelbase];
  
% predict covariance
PX(1:3,1:3)= Gv*PX(1:3,1:3)*Gv' + Gu*PARAMS.Q*Gu';
if size(PX,1)>3
    PX(1:3,4:end)= Gv*PX(1:3,4:end);
    PX(4:end,1:3)= PX(1:3,4:end)';
end    

% predict state
XX(1:3)= [XX(1) + vtc; 
          XX(2) + vts;
         pi_to_pi(XX(3)+ Vn*PARAMS.dt*sin(Gn)/PARAMS.wheelbase)];
