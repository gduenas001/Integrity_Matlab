%%% Configuration file
%%% Permits various adjustments to parameters of the SLAM algorithm.
%%% See ekfslam_sim.m for more information

addpath('./utilities');
set(0,'DefaultFigureWindowStyle','normal')
format compact

% Initialise states and other global variables
global  IA LB C_mesh_interp ny_mesh_interp PIA_terms
global XX PX LM DATA PARAMS SWITCH

% Load table for the lower bounds of the non-centrality parameter
load('LowerBound.mat');
load('PIA_terms_table');

% Load or create the map
% load('example_webmap');

LM= [30, 30;
    4, -4];
% LM= [30, 30, 30;
%     2, 0, -2];
% lm= [30, 30, 30, 30;
%     1.5, -1.5, 3, -3];

% Way points
wp= [0,5;
     0,  0];



% waypoint proximity
PARAMS.at_waypoint= 1.0; % metres, distance from current waypoint at which to switch to next waypoint
PARAMS.numLoops= 1; % number of loops through the waypoint list
PARAMS.numSteps= 100;


% control parameters
PARAMS.v= 0.5;
PARAMS.maxG= 30*pi/180;  % radians, maximum steering angle (-MAXG < g < MAXG)
PARAMS.maxRateG= 20*pi/180; % rad/s, maximum rate of change in steer angle
PARAMS.wheelbase= 4; % metres, vehicle wheel-base
PARAMS.dt_controls= 0.025; % seconds, time interval between control signals
PARAMS.dt= PARAMS.dt_controls;
PARAMS. veh= [0 -PARAMS.wheelbase -PARAMS.wheelbase;
                                        0 -2 2]; % vehicle animation
% control noises
PARAMS.sigmaV= 0.3; % m/s ( default= 0.3 )
PARAMS.sigmaG= deg2rad(3); % radians ( default= (3.0*pi/180) )
PARAMS.Q= [PARAMS.sigmaV^2, 0; 0, PARAMS.sigmaG^2];

% observation parameters
PARAMS.dz= 2; % d.o.f. of one measurement
PARAMS.maxRange= inf; % metres (default = 30)
PARAMS.V_FOV= pi*PARAMS.maxRange^2 / 2; % Total field of view (FOV) of the robot, constant over time
PARAMS.lambda= 1/1000; % density of new lm in the new FOV. Very small number until better modeled.
PARAMS.P_D= 0.95; % probability of detection of a lm.
PARAMS.dt_observe= 1*PARAMS.dt_controls; % seconds, time interval between observations
PARAMS.Const= (2*pi)^(PARAMS.dz/2) * PARAMS.lambda * (1- PARAMS.P_D) * PARAMS.P_D^(-1) ;

% observation noises
PARAMS.sigmaR= 0.1; % metres ( default 0.1 )
PARAMS.sigmaB= deg2rad(1); % radians ( default (1.0*pi/180) )
PARAMS.R= [PARAMS.sigmaR^2 0; 0 PARAMS.sigmaB^2];

     
% data association innovation gates (Mahalanobis distances)
PARAMS.gate= chi2inv(1-1e-15,PARAMS.dz); %0.9999999,dz);
PARAMS.gate= inf;
% For 2-D observation:
%   - common gates are: 1-sigma (1.0), 2-sigma (4.0), 3-sigma (9.0), 4-sigma (16.0)
%   - percent probability mass is: 1-sigma bounds 40%, 2-sigma 86%, 3-sigma 99%, 4-sigma 99.9%.


% switches
SWITCH.control_noise= 1; % if 0, velocity and gamma are perfect
SWITCH.sensor_noise= 1; % if 0, measurements are perfect
SWITCH.inflate_noise= 0; % if 1, the estimated Q and R are inflated (ie, add stabilising noise)
SWITCH.heading_known= 0; % if 1, the vehicle heading is observed directly at each iteration
SWITCH.batch_update= 1; % if 1, process scan in batch, if 0, process sequentially
SWITCH.seed_random= 1; % if not 0, seed the randn() with its value at beginning of simulation (for repeatability)
SWITCH.use_IEKF= 0; % if 1, use iterated EKF for updates, if 0, use normal EKF
SWITCH.profile= 0; % if 1, turn on MatLab profiling to measure time consumed by simulator functions
SWITCH.graphics= 1; % if 0, avoids plotting most animation data to maximise simulation speed
SWITCH.update_global= 0; % if 1, This alternative is a "global constraint" model devised by Jose Guivant, and may have better linearisation properties than the conventional range-bearing model.
SWITCH.association= 2; % if 0, associations are given, if 1, they are estimated using gates, if 2, more options

%% INITIALIZATIONS 
 
 % true & estiamted state
xtrue= zeros(3,1);
XX= zeros(3,1);

% initial pose covariance
PX= [0.0001^2, 0, 0;
     0, 0.0001^2, 0;
     0, 0, deg2rad(0.001)^2];


if SWITCH.seed_random, rand('state',SWITCH.seed_random), randn('state',SWITCH.seed_random), end
if SWITCH.profile, profile on -detail builtin, end

dtsum= 0;               % change in time since last observation
PARAMS.ftag= 1:size(LM,2);     % identifier for each landmark
da_table= zeros(1,size(LM,2)); % data association table 
iwp= 1;                 % index to first waypoint 
G= 0;                   % initial steer angle

% more initializations
step= 0; IA= 0;
DATA.PCA= zeros(5000,1);
DATA.PCA_MJ= zeros(5000,1);
DATA.PCA_MJ_outliers= zeros(5000,1);
DATA.PCAt= ones(5000,1);
DATA.realPCA= zeros(5000,1);
DATA.calcPCA= zeros(5000,1);
DATA.calcPCA_MJ= zeros(5000,1);
DATA.errorXX= zeros(5000,3);
DATA.stdXX= zeros(5000,3);

