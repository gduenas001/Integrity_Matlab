
clear; clc; close all;

syms 

C= -20:5:100;
dof_j= 2:2:10;
dof_0= 2:2:10;
lambda2= [[0:0.1:10],[10.5:0.5:20],[21:1:40],[42:2:50],[55:5:100]];
% lambda2= [[0:0.5:10],[11:1:15],[20:5:40],[50:10:100]];

lenC= length(C);
lendof_j= length(dof_j);
lendof_0= length(dof_0);
lenlambda2= length(lambda2);

PIA_terms= zeros(lenC,lendof_j,lendof_0,lenlambda2);
tic
for i= 1:lenC
    for j= 1: lendof_j
        disp(['C: ', num2str(C(i)), '      dof_j: ', num2str(dof_j(j))]);
        for k= 1:lendof_0
            for l= 1:lenlambda2
                
                fun= @(x) ncx2cdf(x-C(i),dof_j(j),lambda2(l)) .*  chi2pdf(x,dof_0(k));
                PIA_terms(i,j,k,l)= integral( fun, 0, inf );
                
            end
        end
    end
end
toc

[C_mesh_interp,ny_mesh_interp]= meshgrid(C,lambda2);

C_test= 3;
dof_j_test= 4;
dof_0_test= 2;
lambda2_test= 5;

PIA_test_table= interp2...
    (C_mesh_interp,ny_mesh_interp,squeeze(PIA_terms(:,dof_j_test,dof_0_test,:))',C_test,lambda2_test,...
    'spline')

fun= @(x) ncx2cdf(x-C_test,dof_j_test,lambda2_test) .*  chi2pdf(x,dof_0_test);
PIA_test_int= integral( fun, 0, inf )


%%

% [P,LAMBDA2]= meshgrid(p,lambda2);
% figure; hold on; grid on;
% xlabel('degrees of freedom');
% ylabel('non-centrality parameter');
% zlabel('P(IA)');
% surf(P, LAMBDA2, squeeze(PIA_terms(15,:,:))');





