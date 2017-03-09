function [z,idf]= get_observations(x, idf)
% INPUTS:
%   x - vehicle pose [x;y;phi]
%   lm - set of all landmarks
%   idf - index tags for each landmark
%   rmax - maximum range of range-bearing sensor 
%
% OUTPUTS:
%   z - set of range-bearing observations
%   idf - landmark index tag for each observation



[lm_visible,idf]= get_visible_landmarks(x, idf);
z= compute_range_bearing(x,lm_visible);




function [lm_visible,idf]= get_visible_landmarks(x, idf)

global PARAMS LM

% Select set of landmarks that are visible within vehicle's semi-circular field-of-view
dx= LM(1,:) - x(1);
dy= LM(2,:) - x(2);
phi= x(3);

% incremental tests for bounding semi-circle
ii= find(abs(dx) < PARAMS.maxRange & abs(dy) < PARAMS.maxRange ... % bounding box
      & (dx*cos(phi) + dy*sin(phi)) > 0 ...  % bounding line
      & (dx.^2 + dy.^2) < PARAMS.maxRange^2);           % bounding circle
% Note: the bounding box test is unnecessary but illustrates a possible speedup technique
% as it quickly eliminates distant points. Ordering the landmark set would make this operation
% O(logN) rather that O(N).
  
% Eliminate the non-detected lm's using the probabitlity of detection
ind_detected= logical( binornd(1, PARAMS.P_D, [length(ii),1]) );

lm_visible= LM(:,ii);
lm_visible= LM(:,ind_detected');
idf= idf(ii);
idf= idf(ind_detected);

%
%

function z= compute_range_bearing(x,lm)
% Compute exact observation
dx= lm(1,:) - x(1);
dy= lm(2,:) - x(2);
phi= x(3);
z= [sqrt(dx.^2 + dy.^2);
    atan2(dy,dx) - phi];
    

