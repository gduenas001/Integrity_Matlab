

function  Noutliers= generate_models(P,z,Nz,R_2x2)

global T h hlm H Hlm gamma ngamma Y psi

% possible number of outliers- array from 0 to the association with less outliers
Noutliers_possible= unique(T(:,end));

% Generate the variables: h, H, gamma, ngamma
for Noutliers_ind= 1:length(Noutliers_possible)
    Noutliers= Noutliers_possible(Noutliers_ind);
    if Noutliers == Nz, disp('all outliers'); return, end;
    
    ind= find(T(:,end) == Noutliers);
    T_temporal= T(ind,:);
    psi= length(ind);
    
    % make R a block diagonal with higher dimension
    R= kron(eye(Nz-Noutliers),R_2x2);

    % stuck the models in different orders
    H= cell(psi,1);
    h= cell(psi,1);
    gamma= cell(psi,1);
    ngamma= zeros(psi,1);
    Y= cell(psi,1);
    for j= 1:psi
        for i= 1:Nz
            if T_temporal(j,i) ~= 0
                H{j}= [H{j};Hlm{T_temporal(j,i)}];
                h{j}= [h{j};hlm{T_temporal(j,i)}];
            end
        end
        
        % If there is an outlier eliminate the measurement associated with it
        zj= z(:, T_temporal(j,1:Nz) ~= 0 );
        
        gamma{j}= zj(:) - h{j};
        Y{j}= R + H{j}*P*H{j}';
        
        ngamma(j)= gamma{j}'*(Y{j}\gamma{j}); % Compute weighted norms      
    end
    
    if any( ngamma < chi2inv(0.99, (Nz-Noutliers)*2) )
        T= T_temporal;
        break;
    end;
end
    
    
%     
% % Select only the associations with less outliers
% for Noutliers= 0:Nz
%     ind= find(T(:,Nz+1) == Noutliers);
%     
%     if ~isempty(ind)
%         T= T(ind,:);
%         psi= length(ind);
%         
%         break
%     end
% end
% 