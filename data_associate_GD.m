

function [gamma_star,H_star,Y_star,idf,Noutliers,PCA,allOutliers]= data_associate_GD(x,P,z,R,GATE)

global T h H gamma ngamma y Y psi

allOutliers= 0;

Nxv= 3; % number of vehicle pose states
dz= 2;  % d.o.f. of each measurement
Nz= size(z,2); % number of measurements

gamma_star= [];
H_star= [];
Y_star= [];
idf= [];
PCA= 1;

% generate the table with associations
Noutliers= genereate_associations(x,P,z,R,GATE);

% If the only possible one is all outliers -> get out!
if Noutliers == Nz, allOutliers= 1; return, end

% generate the parameters for each association
generate_models(P,z,Nz,Noutliers,R);


if psi > 1
    % Candidate association
    [~,jstar]= min(ngamma);
    idf= T(jstar,1:Nz);
    
    % Eliminate the jstar association from ngamma
    ngamma(jstar)= [];
    min2= min(ngamma);
    
    % calculate PCA
    PCA= chi2cdf(min2, (Nz-Noutliers)*dz);
else
    jstar= 1;
    PCA= 1;
    
    idf= T(jstar,1:Nz);
end

% Save CA values to update KF
gamma_star= gamma{jstar};
H_star= H{jstar};
Y_star= Y{jstar};















