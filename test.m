
clear; close all; clc;

Ntrials= 10000;
dof= 6;

psi= 3;
gamma= zeros(psi-1,1);

% noncetrality parameters
ncx(1)= 40;
ncx(2)= 60;

IA= zeros(Ntrials, 1);
for trial= 1:Ntrials
    disp(trial)
    % generate samples
    gamma0= chi2rnd(dof);
    for i= 1:psi-1
        gamma(i)= ncx2rnd(dof,ncx(i));
    end
    
    % check association
    IA(trial)= any(gamma < gamma0);
end

% P(IA) direct evaluation
fun= @(x) ( ncx2cdf(x,dof,ncx(1)) + ncx2cdf(x,dof,ncx(2)) ) .* chi2pdf(x,dof);
direct_evaluation_PIA = integral( fun , 0, inf)

% P(IA) MJ upper bound
MJ_PIA = chi2cdf( min(ncx)/4 , dof + 3, 'upper')

% actual P(IA)
P_IA= sum(IA) / Ntrials






















