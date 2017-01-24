
function [out]= PCAfun(x,dof, noncentr)

out= ncx2cdf(x, dof, noncentr) .* chi2pdf(x, dof);