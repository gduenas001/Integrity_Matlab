
function [PCA,PCAt]= DA_monte_carlo(x,P,z,idft,R,GATE)

global T h H gamma ngamma y Y psi IA 

Nxv= 3; % number of vehicle pose states
dz= 2;  % d.o.f. of each measurement
Nz= size(z,2); % number of measurements

PCA= 1;
PCAt= 1;

% generate the table with associations
Noutliers= genereate_associations(x,P,z,R,GATE);

% If the only possible one is all outliers -> get out!
if Noutliers == Nz
    PCA= [1,1,1,1];
    PCAt= [1,1,1,1];
    return
end

% generate the parameters for each association
generate_models(P,z,Nz,Noutliers,R);

if psi > 1
    
    % ------------------- Areta ------------------- %
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
    PCA(1)= 1 - integral(@(x) PCAfun(x, (Nz-Noutliers)*dz, min2), 0, inf);
    
    % check CA or IA
    boolCorrect= idf == idft;
    if sum(boolCorrect) ~= size(z,2)-Noutliers
        PCAt(1)= 0;
        IA(1) = IA(1) + 1;
    end
    
    
    % ------------------- GD ------------------- %
    % Candidate association
    [~,jstar]= min(ngamma);
    idf= T(jstar,1:Nz);
    
    % Eliminate the jstar association from ngamma
    ngamma2= ngamma;
    ngamma2(jstar)= [];
    min2= min(ngamma2);
    
    % calculate PCA
    PCA(2)= chi2cdf(min2, (Nz-Noutliers)*dz);
    
    % check CA or IA
    boolCorrect= idf == idft;
    if sum(boolCorrect) ~= size(z,2)-Noutliers
        PCAt(2)= 0;
        IA(2) = IA(2) + 1;
    end
    
    
    % ------------------- BS ------------------- %
    ptheta= zeros(psi,1);
    for j= 1:psi
        ptheta(j)= mvnpdf(gamma{j}, zeros((Nz- Noutliers)*dz, 1), Y{j});
    end
    [~,jstar]= max(ptheta);
    idf= T(jstar,1:Nz);
    
    % calculate PCA
    PCA(3)= ptheta(jstar) / sum(ptheta);
    
    % check CA or IA
    boolCorrect= idf == idft;
    if sum(boolCorrect) ~= size(z,2)-Noutliers
        PCAt(3)= 0;
        IA(3) = IA(3) + 1;
    end
    
    % ------------------- MJ ------------------- %
    % Candidate association
    [~,jstar]= min(ngamma);
    idf= T(jstar,1:Nz);
    
    y= cell(psi,1);
    ny= zeros(psi,1);
    % compute the centers
    for j= 1:psi
        y{j}= h{jstar} - h{j};
        ny(j)= y{j}'*(Y{j}\y{j});
    end
    
    % Eliminate the jstar association from ny
    ny(jstar)= [];
    min2= min(ny);
    
    % calculate PCA
    PCA(4)= chi2cdf(min2/4, (Nz-Noutliers)*dz + Nxv);
    
    % check CA or IA
    boolCorrect= idf == idft;
    if sum(boolCorrect) ~= size(z,2)-Noutliers
        PCAt(4)= 0;
        IA(4) = IA(4) + 1;
    end
    
% if there is only one possible association -> PCA= 1
else
    jstar= 1;
    PCA= [1,1,1,1];
    
    idf= T(jstar,1:Nz);
end




