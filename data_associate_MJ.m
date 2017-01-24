
function [gamma_star,H_star,Y_star,idf,Noutliers,PCA,allOutliers]= data_associate_MJ(x,P,z,R,GATE)

global lm

allOutliers= 0;

Nxv= 3; % number of vehicle pose states
dz= 2;  % d.o.f. of each measurement
Nlm= size(lm,2); % number of features already in map
Nz= size(z,2); % number of measurements

hlm= cell(Nlm,1);
Hlm= cell(Nlm,1);
Slm= cell(Nlm,1);
for l= 1:Nlm
    [hlm{l},Hlm{l}]= observe_model_localization(x,lm(:,l));
    Slm{l}= Hlm{l}*P*Hlm{l}' + R;
end

nis= zeros(Nz,Nlm);
for i= 1:Nz
    for l= 1:Nlm
        v= z(:,i) - hlm{l};
        nis(i,l)= v'*(Slm{l}\v);
    end
end

A= cell(Nz,1);
for i= 1:Nz
    A{i}= [0, find(nis(i,:) < GATE)];
end

% possible combinations
T= combvec(A{:})'; psi= size(T,1); 
T= [T,-ones(psi,1)]; T(1,end)= Nz; % Add a column with the number of outliers

% eliminate impossible associations
j= 2;
while j <= psi
    theta= sort(T(j,:));
    isoutlier= theta == 0;
    Noutliers= sum(isoutlier);
    if Noutliers
        ind= find(isoutlier,1,'last') + 1;
    else
        ind= 1;
    end
    
    isrep= ~all(diff(theta(ind:end)));    
    if isrep
        T(j,:)= [];
        psi= psi - 1;
    else
        T(j,Nz+1)= Noutliers;
        j= j + 1;
    end
end

% Select only the associations with less outliers
for Noutliers= 0:Nz
    ind= find(T(:,Nz+1) == Noutliers);
    if ~isempty(ind)
        T= T(ind,:);
        psi= length(ind);
        
        break
    end
end

% If the only possible one is all outliers -> get out!
if Noutliers == Nz
    gamma_star= [];
    H_star= [];
    Y_star= [];
    idf= [];
    PCA= 1;
    allOutliers= 1;
    
    return
end


% make R more dimensional
R= kron(eye(Nz-Noutliers),R);

% stuck the models in different orders
H= cell(psi,1);
h= cell(psi,1);
gamma= cell(psi,1);
ngamma= zeros(psi,1);
Y= cell(psi,1);
for j= 1:psi
    for i= 1:Nz
        if T(j,i) ~= 0
            H{j}= [H{j};Hlm{T(j,i)}];
            h{j}= [h{j};hlm{T(j,i)}];
        end
    end
    
    % If there is an outlier eliminate the measurement associated with it
    zj= z(:, T(j,1:Nz) ~= 0 );
    
    gamma{j}= zj(:) - h{j};
    Y{j}= R + H{j}*P*H{j}';
    ngamma(j)= gamma{j}'*(Y{j}\gamma{j}); % Compute weighted norms
end


if psi > 1
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
    PCA= chi2cdf(min2/4, (Nz-Noutliers)*dz + Nxv);
else
    jstar= 1;
    PCA= 1;
    
    idf= T(jstar,1:Nz);
end

% Save CA values to update KF
gamma_star= gamma{jstar};
H_star= H{jstar};
Y_star= Y{jstar};



