
function [gamma_star,H_star,Y_star,idf,PCA,PCA_MJ,PCAt,allOutliers]=...
    DA(x,P,z,idft,R,GATE,V_FOV,LAMBDA,P_D,option)

global T h H gamma ngamma y Y psi IA beta phi Const LB C_mesh_interp ny_mesh_interp PIA_terms step 

allOutliers= 0;

Nxv= 3; % number of vehicle pose states
dz= 2;  % d.o.f. of each measurement
Nz= size(z,2); % number of measurements

gamma_star= [];
H_star= [];
Y_star= [];
idf= [];
PCA= 1;
PCAt= 1;

% generate the table with associations
genereate_associations(x,P,z,R,GATE);

% generate the parameters for each association
generate_models(P,z,dz,Nz,R,V_FOV,LAMBDA,P_D);

if psi > 1
    switch option
        case 'direct'
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
            lambda2_mesh
            % Eliminate the jstar association from ny
            ny(jstar)= [];
            min2= min(ny);
            
            % calculate the PCA
            PCA= 1 - integral(@(x) PCAfun(x, (Nz-Noutliers)*dz, min2), 0, inf);
        case 'simplification'
            % Candidate association
            [~,jstar]= min(ngamma);
            idf= T(jstar,1:Nz);
            
            % Eliminate the jstar association from ngamma
            ngamma(jstar)= [];
            min2= min(ngamma);
            
            % calculate PCA
            PCA= chi2cdf(min2, (Nz-Noutliers)*dz);
        case 'BS'
            ptheta= zeros(psi,1);
            for j= 1:psi
                ptheta(j)= mvnpdf(gamma{j}, zeros((Nz- Noutliers)*dz, 1), Y{j});
            end
            [~,jstar]= max(ptheta);
            idf= T(jstar,1:Nz);
            
            % calculate PCA
            PCA= ptheta(jstar) / sum(ptheta);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
        case 'MJ'
            [~,jstar]= min(ngamma);
            if Noutliers == 0
                % Candidate association
                idf= T(jstar,1:Nz);
            else
                % find non-associated measurements
                nonAssociatedMeas= find( T(jstar,:) == 0 );
                T(:, nonAssociatedMeas)= 0;
                
                % Eliminate associations where the associated features in the
                % candidate are not associated
                i= 1;
                while i <= psi
                    if ~any( T(i,1:end-1) ) % if the row is all zeros -> eliminate
                        T(i,:)= []; ngamma(i)= []; h{i}= []; Y{i}= [];
                        h= h(~cellfun('isempty',h));Y= Y(~cellfun('isempty',Y));
                        psi= psi - 1;
                    else
                        i= i + 1;
                    end
                end
                [~,jstar]= min(ngamma);
                % Candidate association
                idf= T(jstar,1:Nz);
            end
            
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
        case 'outliers'
            [~,jstar]= min(beta);
            
            % If the candidate association is all outliers -> notify
            if T(jstar,end) == Nz, allOutliers= 1; end;
            
            nonConflictAssoc= nonConflictingAssociations(jstar);
            nonConflictAssoc= 1;
            
            C= zeros(1,psi);
            y= cell(psi,1);
            ny= zeros(psi,1);
            PIA_term= zeros(psi,1);
            for j= 1:psi
                if j == jstar, continue, end;
                if any(j == nonConflictAssoc), continue, end; % the corresponding PIA_term will be zero
                
                % Degrees of freedom
                dof= dz*(Nz - T(j,end));
                
                
                if dof == 0
                    % if all outliers -> ny = 0
                    ny(j)= 0;
                else               
                    % Interpolate to lower bound the non-centrality parameters
                    ny(j)= interp1( 2:1:200, squeeze(LB(5,dof/2,:)), ngamma(j),'linear','extrap');
                end
                
                % Compute C_j
                C(j)= log( (det(Y{j})/det(Y{jstar})) * Const^(2*(phi(jstar)-phi(j))) );
                
%                 % Direct evaluation of P(IA) - integral
%                 fun= @(x) ncx2cdf(x-C(j),dof,ny(j)) .*  chi2pdf(x,dof);
%                 PIA_term(j) = integral( fun , 0, inf);
                
                % Direct evaluation of P(IA) - table
                PIA_term(j)= ...
                    interp2(C_mesh_interp,ny_mesh_interp,squeeze(PIA_terms(:,dof/2,:))',C(j),ny(j),'linear',0);
            end
            
            PIA= sum(PIA_term);
            PCA= 1 - PIA;
            idf= T(jstar,1:Nz);
            
            %% Comparizon with MJ approach - only compares with same #outliers associations
            T_MJ= T;
            Noutliers_MJ= T_MJ(jstar,end);
            
            % Keep only those associations with the same number of outliers
            Tind_MJ= T_MJ(:,end) == Noutliers_MJ;
            T_MJ= T_MJ(Tind_MJ, :);
            psi_MJ= size(T_MJ,1);
            
            if psi_MJ == 1
                PCA_MJ= 1;
            else
                ny_MJ= ny(Tind_MJ);
                [~,jstar_MJ]= min(ny_MJ);
                ny_MJ(jstar_MJ)= [];
                nymin2_MJ= min(ny_MJ);
                PCA_MJ= chi2cdf(nymin2_MJ/4, (Nz-Noutliers_MJ)*dz + Nxv);
                PI_MJ= 1- PCA_MJ;
            end
    end
else
    
    
    jstar= 1;
    PCA= 1;
    PCA_MJ= 1;
    
    if T(jstar,end) == Nz, allOutliers= 1; end;
    idf= T(jstar,1:Nz);
end

% Save CA values to update KF
gamma_star= gamma{jstar};
H_star= H{jstar};
Y_star= Y{jstar};

% check CA or IA
Noutlierst= 0; % Real number of outiers, always 0 with no outliers
if T(jstar,end) < Noutlierst
    IA= IA + 1; PCAt= 0;
elseif T(jstar,end) == Noutlierst
    if ~sum(idf == idft)
        IA= IA + 1; PCAt= 0;
    end
elseif T(jstar,end) > Noutlierst
    if any((idft - idf) .* idf)
        IA= IA + 1; PCAt= 0;
    end
end




