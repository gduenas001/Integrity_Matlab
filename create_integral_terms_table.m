
clear; clc; close all;

syms 

C= -10:5:100;
p= 2:2:10;
lambda2= 0:1:100;

lenC= length(C);
lenp= length(p);
lenlambda2= length(lambda2);

PIA_terms= zeros(lenC,lenp,lenlambda2);
tic
parfor i= 1:lenC
    for j= 1: lenp
        disp(['C: ', num2str(C(i)), '      p: ', num2str(p(j))]);
        for k= 1:lenlambda2
            
            fun= @(x) ncx2cdf(x-C(i),p(j),lambda2(k)) .*  chi2pdf(x,p(j));
            PIA_terms(i,j,k)= integral( fun, 0, inf );
            
        end
    end
end
toc

[C_mesh_interp,ny_mesh_interp]= meshgrid(C,lambda2);
interp2(C_mesh_interp,ny_mesh_interp,squeeze(PIA_terms(:,4,:))',5,7.9)
%%

% [P,LAMBDA2]= meshgrid(p,lambda2);
% figure; hold on; grid on;
% xlabel('degrees of freedom');
% ylabel('non-centrality parameter');
% zlabel('P(IA)');
% surf(P, LAMBDA2, squeeze(PIA_terms(15,:,:))');





