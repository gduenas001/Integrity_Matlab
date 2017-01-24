
function [gamma_star,H_star,Y_star,idf,Noutliers,PCA,PCAt,allOutliers]=...
    DA(x,P,z,idft,R,GATE,option)

global T h H gamma ngamma y Y psi IA step

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
Noutliers= generate_models(P,z,Nz,R);

% If the only possible one is all outliers -> get out!
if Noutliers == Nz, allOutliers= 1; return, end


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
            ptheta= zeros(psi,1);
            for j= 1:psi
                ptheta(j)= mvnpdf(gamma{j}, zeros((Nz- Noutliers)*dz, 1), Y{j});
            end
            [~,jstar]= max(ptheta);
            idf= T(jstar,1:Nz);
            
            % calculate PCA
            PCA= ptheta(jstar) / sum(ptheta);
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

% check CA or IA
if ~allOutliers
    % check CA or IA
    boolCorrect= idf == idft;
    if sum(boolCorrect) ~= Nz-Noutliers % this is not a good check
        PCAt= 0;
        IA = IA + 1;
    end
else
    PCAt= 1;
end





