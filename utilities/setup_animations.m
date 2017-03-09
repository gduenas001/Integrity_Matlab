function h= setup_animations(wp)

global LM

% Create lm names
str= {};
for l= 1:size(LM,2)
    str= [str, strcat('lm', num2str(l))];
end

% Create figure
fig=figure; hold on; axis equal;
xlabel('metres'), ylabel('metres')
set(fig, 'name', 'EKF-localization Simulator')

% Plot lm and waypoints
plot(LM(1,:),LM(2,:),'b*');
text(LM(1,:)-1,LM(2,:)-2,str,'FontSize',5,'color','b');
plot(wp(1,:),wp(2,:), 'g', wp(1,:),wp(2,:),'g.')

% Initialize dynamic plots
h.xt= patch(0,0,'b'); % vehicle true
h.xv= patch(0,0,'r'); % vehicle estimate
h.pth= plot(0,0,'k.','markersize',2); % vehicle path estimate
h.obs= plot(0,0,'y'); % observations
h.xf= plot(0,0,'r+'); % estimated features
h.vcov= plot(0,0,'r'); % vehicle covariance ellipses
h.fcov= plot(0,0,'r'); % feature covariance ellipses
