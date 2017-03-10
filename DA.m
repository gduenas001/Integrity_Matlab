
function [gamma_star, H_star, Y_star, R_star, PCA_GD, PCA_MJ]=...
    DA(z, R_2x2)

global T h hlm H Hlm gamma ngamma Y LB 
global PX PARAMS LM

Nz= size(z,2);
Nlm= size(LM,2);

% If no measurements -> return
if Nz == 0, gamma_star= []; H_star= []; R_out= []; PIA= 0; return, end

% Initialize lm models to be filled in the next function
hlm= cell(Nlm,1);
Hlm= cell(Nlm,1);

% Create T associations table
[psi,nis]= associations_table(z, Nz, Nlm, R_2x2, PARAMS.gate);

log_C2= log(PARAMS.Const^2);

beta= zeros(psi,1);
phi= zeros(1,psi);
H= cell(psi,1);
h= cell(psi,1);
gamma= cell(psi,1);
ngamma= zeros(psi,1);
Y= cell(psi,1);
for j= 1:psi
    
    % current possible association
    T_j= T(j,:);
    phi(j)= T_j(end);
    
    % If all outliers
    if phi(j) == Nz, beta(j)= (PARAMS.lambda* (1-PARAMS.P_D)/ PARAMS.P_D)^Nz; continue, end;
    
    % If there is an outlier eliminate the measurement associated with it
    zj= z(:, T_j(1:end-1) ~= 0 );
    
    % make R a block diagonal with higher dimension
    R= kron(eye(Nz-phi(j)),R_2x2);
    
    % create the models for each association
    for i= 1:Nz
        if T_j(i) ~= 0
            H{j}= [H{j};Hlm{T_j(i)}];
            h{j}= [h{j};hlm{T_j(i)}];
        end
    end
    gamma{j}= zj(:) - h{j};
    Y{j}= H{j}*PX*H{j}' + R;
    ngamma(j)= gamma{j}'*(Y{j}\gamma{j}); % Compute weighted norms
    
    % compute BETA 
    beta(j)= log(det(Y{j})) + ngamma(j) - phi(j)*log_C2;
end

[~,jstar]= min(beta);

gamma_star= gamma{jstar};
H_star= H{jstar};
Y_star= Y{jstar};
R_star= kron(eye(Nz-T(jstar,end)),R_2x2);

dof= PARAMS.dz .* (Nz - T(:,end));

PIA_term= zeros(psi,1);
ny= zeros(psi,1);
ny_extended= zeros(psi,1);
C= zeros(psi,1);
LHS= zeros(psi,1);

for j= 1:psi
    if j == jstar, continue, end;
%     if any(j == nonConflictAssoc), continue, end; % the corresponding PIA_term will be zero
    if phi(j) == Nz
        C(j)= ( PARAMS.lambda * (1 - PARAMS.P_D)/PARAMS.P_D )^Nz + phi(jstar)*log_C2...
            - log(det(Y{jstar}));
        C(j)= 
        PIA_term(j)= chi2cdf(C(j), dof(j), 'upper'); 
        ny(j)= 0;
    else
        if dof(j) == 0 % if all outliers -> ny = 0
            ny(j)= 0;
        else % Interpolate to lower bound the non-centrality parameters
            ny(j)= interp1( 2:1:200, squeeze(LB(1,dof(j)/2,:)), ngamma(j),'linear','extrap');
        end
        % Compute C_j
        C(j)= log( (det(Y{j})/det(Y{jstar})) * PARAMS.Const^(2*(phi(jstar)-phi(j))) );
        
        % Direct evaluation of P(IA) - integral
        fun= @(x) ncx2cdf(x-C(j),dof(j),ny(j)) .*  chi2pdf(x,dof(jstar));
        PIA_term(j) = integral( fun , 0, inf);
    end
    
    % Compute the LHS for the MJ P(IA)
    if C(j) >= 0 
        LHS(j)= ( sqrt(ny(j)) + sqrt(C(j)) / (1+sqrt(2)) )^2;
    else
        LHS(j)= 0.25* ( sqrt(ny(j)) - sqrt(-C(j)) )^2;
    end
end

PIA= sum(PIA_term);
PCA_GD= 1 - PIA;

% MJ approach with outliers - find minimum
LHS(jstar)= inf;
PCA_MJ= chi2cdf(min(LHS), Nz*PARAMS.dz);  
;

    
    
    
    


% Nxv= 3; % number of vehicle pose states
% Nz= size(z,2); % number of measurements
% 
% gamma_star= [];
% H_star= [];
% Y_star= [];
% idf= [];
% PCA= 1;
% PCAt= 1;
% 
% % generate the table with associations
% genereate_associations(z,R,GATE);
% 
% % generate the parameters for each association
% generate_models(z,dz,Nz,R,V_FOV,LAMBDA,P_D);
% 
% if psi > 1
%     switch option
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
%         case 'outliers'
%             [~,jstar]= min(beta); dof_0= dz*(Nz - T(jstar,end));
%             
%             % If the candidate association is all outliers -> notify
%             if T(jstar,end) == Nz, allOutliers= 1; end;
%             
%             nonConflictAssoc= nonConflictingAssociations(jstar);
%             nonConflictAssoc= []; % consider all conflicting
%             
%             dof= dz .* (Nz - T(:,end));
%             
%             C= zeros(1,psi);
%             y= cell(psi,1);
%             ny= zeros(psi,1);
%             PIA_term= zeros(psi,1);
%             for j= 1:psi
%                 if j == jstar, continue, end;
%                 if any(j == nonConflictAssoc), continue, end; % the corresponding PIA_term will be zero
%                 
%                 % Degrees of freedom
% %                 dof_j= dz*(Nz - T(j,end));
%                 
%                 if dof(j) == 0
%                     % if all outliers -> ny = 0
%                     ny(j)= 0;
%                 else               
%                     % Interpolate to lower bound the non-centrality parameters
%                     ny(j)= interp1( 2:1:200, squeeze(LB(5,dof(j)/2,:)), ngamma(j),'linear','extrap');
%                 end
%                 
%                 % Compute C_j
%                 C(j)= log( (det(Y{j})/det(Y{jstar})) * Const^(2*(phi(jstar)-phi(j))) );
% %                 if C(j) < 0 && T(jstar,end) > 0 
% %                     ;
% %                 end
%                 
% %                 % Direct evaluation of P(IA) - integral
% %                 fun= @(x) ncx2cdf(x-C(j),dof(j),ny(j)) .*  chi2pdf(x,dof_0);
% %                 PIA_term(j) = integral( fun , 0, inf);
%                 
%                 % Direct evaluation of P(IA) - table
%                 PIA_term(j)= interp2(C_mesh_interp,ny_mesh_interp,...
%                     squeeze(PIA_terms(:,dof(j)/2+1,dof_0/2+1,:))',C(j),ny(j),'spline',0);
%             end
%             
%             PIA= sum(PIA_term);
%             PCA= 1 - PIA;
%             idf= T(jstar,1:Nz);
%             
%             %% Comparizon with MJ approach - only compares with same #outliers associations
%             T_MJ= T;
%             Noutliers_MJ= T_MJ(jstar,end);
%             
%             % Keep only those associations with the same number of outliers
%             Tind_MJ= T_MJ(:,end) == Noutliers_MJ;
%             T_MJ= T_MJ(Tind_MJ, :);
%             psi_MJ= size(T_MJ,1);
%             
%             if psi_MJ == 1
%                 PCA_MJ= 1;
%             else
%                 ny_MJ= ny(Tind_MJ);
%                 [~,jstar_MJ]= min(ny_MJ);
%                 ny_MJ(jstar_MJ)= [];
%                 nymin2_MJ= min(ny_MJ);
%                 PCA_MJ= chi2cdf(nymin2_MJ/4, (Nz-Noutliers_MJ)*dz + Nxv);
%                 PI_MJ= 1- PCA_MJ;
%             end
%             %% Comparizon with MJ with outliers approach            
%             posInds= C > 0;
%             LHS= ( sqrt(ny(posInds)) + sqrt(C(posInds)') ).^2;
%             posdof= dof(posInds);
%             posPIA= zeros(sum(posInds),1);
%             for j= 1:sum(posInds)
%                 posPIA(j)= 1 - gamcdf(LHS(j), (dof_0+posdof(j))/2, 8);
%             end
%             
%             negInds= C < 0;
%             LHS= 0.5 .* ( sqrt(ny(negInds)) - sqrt(abs(C(negInds))') ).^2;
%             negdof= dof(negInds);
%             negPIA= zeros(sum(negInds),1);
%             for j= 1:sum(negInds)
%                 negPIA(j)= 1 - chi2cdf(LHS(j), dof_0+negdof(j));
%             end
%             PCA_MJ_outliers= 1 - max([posPIA; negPIA]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %            
%     end
% else
%     
%     
%     jstar= 1;
%     PCA= 1;
%     PCA_MJ= 1;
%     
%     if T(jstar,end) == Nz, allOutliers= 1; end;
%     idf= T(jstar,1:Nz);
% end
% 
% % Save CA values to update KF
% gamma_star= gamma{jstar};
% H_star= H{jstar};
% Y_star= Y{jstar};
% 
% % check CA or IA
% Noutlierst= 0; % Real number of outiers, always 0 with no outliers
% if T(jstar,end) < Noutlierst
%     IA= IA + 1; PCAt= 0;
% elseif T(jstar,end) == Noutlierst
%     if ~sum(idf == idft)
%         IA= IA + 1; PCAt= 0;
%     end
% elseif T(jstar,end) > Noutlierst
%     if any((idft - idf) .* idf)
%         IA= IA + 1; PCAt= 0;
%     end
% end
% 
% 




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [psi,nis]= associations_table(z, Nz, Nlm, R, gate)
global T


nis= ones(Nz,Nlm)*gate;
for i= 1:Nz
    for l= 1:Nlm
        [nis(i,l)]= lm_model_and_nis(z(:,i),R,l);        
    end
end

% Eliminate very improbable associations using individual gates
A= cell(Nz,1);
for i= 1:Nz
    A{i}= [0, find(nis(i,:) < gate)];
end

% Generate table with all possible combinations
T= allcomb(A{:}); psi= size(T,1); 
T= [T,-ones(psi,1)]; T(1,end)= Nz; % Add a column with the number of outliers

% Eliminate impossible associations, i.e. repeated landmarks
psi= eliminate_impossible_associations(psi, Nz);




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [nis]= lm_model_and_nis(z,R,lm_id)
global XX PX LM hlm Hlm

% return normalised innovation squared (ie, Mahalanobis distance) and normalised distance

% auxiliary values
dx= LM(1,lm_id)  -XX(1); 
dy= LM(2,lm_id) - XX(2);
d2= dx^2 + dy^2;
d= sqrt(d2);

xd= dx/d;
yd= dy/d;
xd2= dx/d2;
yd2= dy/d2;

% predict z
h= [d;
    atan2(dy,dx) - XX(3)];


v= z-h; v(2)= pi_to_pi(v(2));

% Quick check to discard associations
% if abs(v(1)) > 2 || abs(v(2)) > deg2rad(15), nis= inf; return, end;


% calculate H
H = [-xd -yd 0; 
      yd2 -xd2 -1];

S= H*PX*H' + R; % TODO: optimise this line -- H is sparse

nis= v'/S*v;

% Store the values in the global varibles
hlm{lm_id}= h;
Hlm{lm_id}= H;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [psi]= eliminate_impossible_associations(psi, Nz)
global T

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

