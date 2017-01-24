
function [zf, idf, zn, jstar, PCA]= data_associate_paper2(x,P,z,R)

zf= []; zn= [];
idf= []; PCA= 1;
jstar= 0;

Nxv= 3; % number of vehicle pose states
Nlm= (length(x) - Nxv)/2; % number of features already in map

% if both lm are not in the map, put them
if Nlm ~= 2
    zn= z; return
else
    zf= z;
end

% make R more dimensional
R= kron(eye(Nlm),R);

% get models
[h_1,Hz_1,Hx_1]= observe_model2(x,z,1);
[h_2,Hz_2,Hx_2]= observe_model2(x,z,2);

z= z(:); % convert to vector of measurements

% permutations
permuts= perms(1:2);

% stuck the models in different orders
h{1}= [h_2;h_1];
h{2}= [h_1;h_2];
Hx{1}= blkdiag(Hx_2,Hx_1);
Hx{2}= blkdiag(Hx_1,Hx_2);
Hz{1}= blkdiag(Hz_2,Hz_1);
Hz{2}= blkdiag(Hz_1,Hz_2);

% inoovation
gamma{1}= z - h{1};
gamma{2}= z - h{2};

% normalization parameter and its chol
Lambda{1}= Hz{1}*R*Hz{1}' + Hx{1}*kron(eye(2),P(1:3,1:3))*Hx{1}';
Lambda{2}= Hz{2}*R*Hz{2}' + Hx{2}*kron(eye(2),P(1:3,1:3))*Hx{2}';
L{1}= chol(Lambda{1},'lower');
L{2}= chol(Lambda{2},'lower');

% normalize gamma
gamma_u{1}= L{1}\gamma{1};
gamma_u{2}= L{2}\gamma{2};

% compute norms
C(1)= norm(gamma_u{1});
C(2)= norm(gamma_u{2});

% DA - minimun norm
[~,jstar]= min(C);
idf= permuts(jstar,:);

disp(jstar);































