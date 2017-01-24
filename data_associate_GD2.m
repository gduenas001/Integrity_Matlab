
function [gamma_star,H_star,Y_star,idf,Noutliers,PCA,allOutliers]=...
    data_associate_GD2(x,P,z,R,GATE, option)

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
    switch option
        case 'GD2'
            % Candidate association
            [~,jstar]= min(ngamma);
            idf= T(jstar,1:Nz);
            
            y= cell(psi,1);
            ny= zeros(psi,1);
            % compute the centers
            for j= 1:psi
                y{j}= h{jstar} - h{j};
                ny(j)= y{j}'*(Y{j}\y{j}); % this is the non-centrality parameter
            end
            
            % Eliminate the jstar association from ny
            ny(jstar)= [];
            min2= min(ny);
            
            % calculate the PCA
            PCA= 1 - integral(@(x) PCAfun(x, (Nz-Noutliers)*dz, min2), 0, inf);
        case 'GD'
            % Candidate association
            [~,jstar]= min(ngamma);
            idf= T(jstar,1:Nz);
            
            % Eliminate the jstar association from ngamma
            ngamma(jstar)= [];
            min2= min(ngamma);
            
            % calculate PCA
            PCA= chi2cdf(min2, (Nz-Noutliers)*dz);
        case 'BS'
            
            
    end
else
    jstar= 1;
    PCA= 1;
    
    idf= T(jstar,1:Nz);
end
% Save CA values to update KF
gamma_star= gamma{jstar};
H_star= H{jstar};
Y_star= Y{jstar};


