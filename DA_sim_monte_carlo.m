
clear;
addpath('./utilities');

set(0,'DefaultFigureWindowStyle','docked')
dbstop if error

format compact
configfile; % ** USE THIS FILE TO CONFIGURE THE EKF-SLAM **

% Initialise states and other global variables
global XX PX lm IA

lm= [30, 30;
      0.5, -0.5];
%   lm= [30, 30, 30, 30;
%       1.5, -1.5, 3, -3];

wp= [0,80;
     0,  0];
xtrue= zeros(3,1);
XX= zeros(3,1);
PX= [0.3^2, 0, 0;
     0, 0.3^2, 0;
     0, 0, deg2rad(5)^2];

% Setup plots
if SWITCH_GRAPHICS
    fig=figure;
    plot(lm(1,:),lm(2,:),'b*')
    hold on, axis equal
    plot(wp(1,:),wp(2,:), 'g', wp(1,:),wp(2,:),'g.')
    xlabel('metres'), ylabel('metres')
    set(fig, 'name', 'EKF-SLAM Simulator')
    h= setup_animations;
    veh= [0 -WHEELBASE -WHEELBASE; 0 -2 2]; % vehicle animation
    plines=[]; % for laser line animation
    pcount=0;
end

% Initialise other variables and constants
dt= DT_CONTROLS;        % change in time between predicts
dtsum= 0;               % change in time since last observation
ftag= 1:size(lm,2);     % identifier for each landmark
da_table= zeros(1,size(lm,2)); % data association table 
iwp= 1;                 % index to first waypoint 
G= 0;                   % initial steer angle
QE= Q; RE= R; if SWITCH_INFLATE_NOISE, QE= 2*Q; RE= minR + 2*R; end % inflate estimated noises (ie, add stabilising noise)
if SWITCH_SEED_RANDOM, rand('state',SWITCH_SEED_RANDOM), randn('state',SWITCH_SEED_RANDOM), end
if SWITCH_PROFILE, profile on -detail builtin, end

Nloops= 12;
calcPCA_final= zeros(Nloops,4);    
realPCAt_final= zeros(Nloops,4);

lm_cell{1}= [30,30; 0.5,-0.5];
lm_cell{2}= [30,30; 0.6,-0.6];
lm_cell{3}= [30,30; 0.7,-0.7];
lm_cell{4}= [30,30; 0.8,-0.8];
lm_cell{5}= [30,30; 0.9,-0.9];
lm_cell{6}= [30,30; 1,-1];
lm_cell{7}= [30,30; 1.2,-1.2];
lm_cell{8}= [30,30; 1.4,-1.4];
lm_cell{9}= [30,30; 1.6,-1.6];
lm_cell{10}= [30,30; 1.8,-1.8];
lm_cell{11}= [30,30; 2,-2];
lm_cell{12}= [30,30; 2.4,-2.4];


for loop= 1:Nloops
disp(loop)
if SWITCH_SEED_RANDOM, rand('state',SWITCH_SEED_RANDOM), randn('state',SWITCH_SEED_RANDOM), end

lm= lm_cell{loop};
distance_lm(loop)= lm(2,1)*2;

step= 0;
PCA= zeros(5000,4);
PCAt= ones(5000,4);
realPCA= zeros(5000,4);
calcPCA= zeros(5000,4);
IA= zeros(1,4);

% Main loop 
while iwp ~= 0
    if step > NUMBER_STEPS, break, end
    step= step + 1;
    
    % sample from a  Gaussian
    XX= mvnrnd(xtrue,PX)';
    
    % Incorporate observation, (available every DT_OBSERVE seconds)
    dtsum= dtsum + dt;
    if dtsum >= DT_OBSERVE
        dtsum= 0;
        
        % get measurements
        [z,idft]= get_observations(xtrue, lm, ftag, MAX_RANGE); 
        z= add_observation_noise(z,R, SWITCH_SENSOR_NOISE);
        
        if ~isempty(z)
            % DA
            if SWITCH_ASSOCIATION == 1
                [zf,idf,zn, da_table]= data_associate_known(XX,z,idft, da_table);
            elseif SWITCH_ASSOCIATION == 0
                [zf,idf, zn]= data_associate(XX,PX,z,RE, GATE_REJECT, GATE_AUGMENT);
            elseif SWITCH_ASSOCIATION == 2
                [gamma,H,Y,idf,Noutliers,PCA(step),allOutliers]= data_associate_GD(XX,PX,z,RE,GATE);
            elseif SWITCH_ASSOCIATION == 3
                [PCA(step,:),PCAt(step,:)]= DA_monte_carlo(XX,PX,z,idft,RE,GATE);
            end
        else
            PCA(step)= 1;
        end
    end
    realPCA(step,:)= 1 - IA/step;
    calcPCA(step,:)= sum(PCA)/step;

    % Plots
    if SWITCH_GRAPHICS
        
        xt= transformtoglobal(veh, xtrue);
        set(h.xt, 'xdata', xt(1,:), 'ydata', xt(2,:))
        
        xv= transformtoglobal(veh, XX(1:3));
        pvcov= make_vehicle_covariance_ellipse(XX,PX);
        set(h.xv, 'xdata', xv(1,:), 'ydata', xv(2,:))
        set(h.vcov, 'xdata', pvcov(1,:), 'ydata', pvcov(2,:))
        drawnow
    end
    
    
end % end of main loop
if SWITCH_PROFILE, profile report, end


PCA(step+1:end,:)= [];
PCAt(step+1:end,:)= [];
realPCA(step+1:end,:)= [];
calcPCA(step+1:end,:)= [];

% save the last value
calcPCA_final(loop,:)= calcPCA(end,:);
realPCAt_final(loop,:)= realPCA(end,:);


end

% plots - PCA VS time
figure; hold on; grid on;
plot(PCA(:,1),'-b');
plot(PCAt(:,1),'or');
idx= find(PCAt(:,1) == 0);
for i= 1:length(idx)
    line([idx(i),idx(i)],[0,1],'color','red');
end
axis([0,step,0,1]);


% % plots - realPCA Vs calcPCA
% figure; hold on; grid on
% plot(realPCA(:,1),'k--','linewidth',4);
% 
% plot(calcPCA(:,1),'r-','linewidth',2);
% plot(calcPCA(:,2),'b-','linewidth',2);
% plot(calcPCA(:,3),'g-','linewidth',2);
% plot(calcPCA(:,4),'y-','linewidth',2);
% 
% legend('Real P(CA)','Areta','GD','BS','MJ','location','southeast');

% plots - finals
figure; hold on; grid on;
plot(distance_lm,realPCAt_final(:,1),'k--','linewidth',4);

plot(distance_lm,calcPCA_final(:,1),'r-','linewidth',2);
plot(distance_lm,calcPCA_final(:,2),'b-','linewidth',2);
plot(distance_lm,calcPCA_final(:,3),'g-','linewidth',2);
plot(distance_lm,calcPCA_final(:,4),'y-','linewidth',2);

legend('Real P(CA)','Areta','GD','BS','MJ','location','southeast');
xlabel('Distance between landmarks')
axis([distance_lm(1),distance_lm(end),0,1]);




