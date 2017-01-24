%%% Configuration file
%%% Permits various adjustments to parameters of the SLAM algorithm.
%%% See ekfslam_sim.m for more information

addpath('./utilities');
set(0,'DefaultFigureWindowStyle','normal')
format compact

% control parameters
V= 0.5; % m/s
MAXG= 30*pi/180; % radians, maximum steering angle (-MAXG < g < MAXG)
RATEG= 20*pi/180; % rad/s, maximum rate of change in steer angle
WHEELBASE= 4; % metres, vehicle wheel-base
DT_CONTROLS= 0.025; % seconds, time interval between control signals

% control noises
sigmaV= 0.3; % m/s ( default= 0.3 )
sigmaG= deg2rad(3); % radians ( default= (3.0*pi/180) )
Q= [sigmaV^2, 0; 0, sigmaG^2];

% observation parameters
dz= 2; % d.o.f. of one measurement
MAX_RANGE= 50.0; % metres (default = 30)
DT_OBSERVE= 1*DT_CONTROLS; % seconds, time interval between observations

% observation noises
sigmaR= 0.1; % metres ( default 0.1 )
sigmaB= deg2rad(1); % radians ( default (1.0*pi/180) )
R= [sigmaR^2 0; 0 sigmaB^2];
% minR= [20^2, 0;
%          0,  deg2rad(5)^2];
     
% data association innovation gates (Mahalanobis distances)
GATE= chi2inv(0.99,dz);
GATE= inf;
GATE_REJECT= 4.0; % maximum distance for association
GATE_AUGMENT= 25.0; % minimum distance for creation of new feature
% For 2-D observation:
%   - common gates are: 1-sigma (1.0), 2-sigma (4.0), 3-sigma (9.0), 4-sigma (16.0)
%   - percent probability mass is: 1-sigma bounds 40%, 2-sigma 86%, 3-sigma 99%, 4-sigma 99.9%.

% waypoint proximity
AT_WAYPOINT= 1.0; % metres, distance from current waypoint at which to switch to next waypoint
NUMBER_LOOPS= 1; % number of loops through the waypoint list
NUMBER_STEPS= 1000;

% switches
SWITCH_CONTROL_NOISE= 1; % if 0, velocity and gamma are perfect
SWITCH_SENSOR_NOISE = 1; % if 0, measurements are perfect
SWITCH_INFLATE_NOISE= 0; % if 1, the estimated Q and R are inflated (ie, add stabilising noise)
SWITCH_HEADING_KNOWN= 0; % if 1, the vehicle heading is observed directly at each iteration
SWITCH_BATCH_UPDATE= 1; % if 1, process scan in batch, if 0, process sequentially
SWITCH_SEED_RANDOM= 1; % if not 0, seed the randn() with its value at beginning of simulation (for repeatability)
SWITCH_USE_IEKF= 0; % if 1, use iterated EKF for updates, if 0, use normal EKF
SWITCH_PROFILE= 0; % if 1, turn on MatLab profiling to measure time consumed by simulator functions
SWITCH_GRAPHICS= 1; % if 0, avoids plotting most animation data to maximise simulation speed
SWITCH_UPDATE_GLOBAL= 0; % if 1, This alternative is a "global constraint" model devised by Jose Guivant, and may have better linearisation properties than the conventional range-bearing model.
SWITCH_ASSOCIATION= 2; % if 0, associations are given, if 1, they are estimated using gates, if 2, more options


%% INITIALIZATIONS 

% Initialise states and other global variables
global XX PX lm IA step

lm= [30, 30;
      1.8, -1.8];
%   lm= [30, 30, 30, 30;
%       0.8, -0.8, 1.6, -1.6];

% Way points
wp= [0,80;
     0,  0];
 
 % true & estiamted state
xtrue= zeros(3,1);
XX= zeros(3,1);

% initial pose covariance
PX= [0.3^2, 0, 0;
     0, 0.3^2, 0;
     0, 0, deg2rad(5)^2];

dt= DT_CONTROLS;        % change in time between predicts
QE= Q; RE= R; if SWITCH_INFLATE_NOISE, QE= 2*Q; RE= minR + 2*R; end % inflate estimated noises (ie, add stabilising noise)
if SWITCH_SEED_RANDOM, rand('state',SWITCH_SEED_RANDOM), randn('state',SWITCH_SEED_RANDOM), end
if SWITCH_PROFILE, profile on -detail builtin, end

dtsum= 0;               % change in time since last observation
ftag= 1:size(lm,2);     % identifier for each landmark
da_table= zeros(1,size(lm,2)); % data association table 
iwp= 1;                 % index to first waypoint 
G= 0;                   % initial steer angle

% more initializations
step= 0;
PCA= zeros(5000,1);
PCAt= ones(5000,1);
realPCA= zeros(5000,1);
calcPCA= zeros(5000,1);
IA= 0;
errorXX= zeros(5000,3);
stdXX= zeros(5000,3);

