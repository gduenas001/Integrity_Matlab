function p= make_vehicle_covariance_ellipse(x,P)
% compute ellipses for plotting vehicle covariances
N= 10;
inc= 2*pi/N;
phi= 0:inc:2*pi;
circ= 3*[cos(phi); sin(phi)];

p= make_ellipse(x(1:2), P(1:2,1:2), circ);


