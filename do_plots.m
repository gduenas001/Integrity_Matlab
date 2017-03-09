
function do_plots(xtrue, z, h, dtsum)

global XX PARAMS SWITCH

% Plots
if SWITCH.graphics
    
    xt= transformtoglobal(PARAMS.veh, xtrue);
    set(h.xt, 'xdata', xt(1,:), 'ydata', xt(2,:))
    
    xv= transformtoglobal(PARAMS.veh, XX(1:3));
    pvcov= make_vehicle_covariance_ellipse();
    set(h.xv, 'xdata', xv(1,:), 'ydata', xv(2,:))
    set(h.vcov, 'xdata', pvcov(1,:), 'ydata', pvcov(2,:))
    
    %         pcount= pcount+1;
    %         if pcount == 120 % plot path infrequently
    %             pcount=0;
    %             set(h.pth, 'xdata', DATA.path(1,1:DATA.i), 'ydata', DATA.path(2,1:DATA.i))
    %         end
    
    if dtsum==0 && ~isempty(z) % plots related to observations
        set(h.xf, 'xdata', XX(4:2:end), 'ydata', XX(5:2:end))
        plines= make_laser_lines (z,XX(1:3));
        set(h.obs, 'xdata', plines(1,:), 'ydata', plines(2,:))
    end
    drawnow
end












































