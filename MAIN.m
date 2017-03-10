
% This script runs a localization simulator with DA and P(IA)

dbstop if error
dbclear if error


clear all; close all;
configfile; % ** USE THIS FILE TO CONFIGURE THE EKF-SLAM **

h= setup_animations(wp);


tic
% *****************    MAIN LOOP    *****************
while iwp ~= 0
    if step > PARAMS.numSteps, break, end;
    if iwp==0 && PARAMS.numLoops > 1, iwp=1; PARAMS.numLoops= PARAMS.numLoops-1; end % perform loops: if final waypoint reached, go back to first

    step= step + 1;
    
    % Compute true data & add noise
    [G,iwp]= compute_steering(xtrue, wp, iwp,  G);
    xtrue= vehicle_model(xtrue, G);
    [Vn,Gn]= add_control_noise(G);
    
    % EKF predict step
    predict (Vn,Gn);
        
    % Incorporate observation, (available every DT_OBSERVE seconds)
    dtsum= dtsum + PARAMS.dt;
    if dtsum >= PARAMS.dt_observe
        dtsum= 0;
        
        % get measurements
        [z,idft]= get_observations(xtrue, PARAMS.ftag); 
        z= add_observation_noise(z);
        
        % DA
        if ~isempty(z)
            if SWITCH.association == 0
                [zf,idf,zn, da_table]= data_associate_known(XX,z,idft, da_table);
            elseif SWITCH.association == 1
                [zf,idf, zn]= data_associate_localNN(XX,PX,z,RE, GATE_REJECT, GATE_AUGMENT);
            elseif SWITCH.association == 2
                [gamma, H, Y, R, DATA.PCA(step), DATA.PCA_MJ(step)]= DA(z, PARAMS.R);
            end
            
            % update the state
            if ~isempty(gamma)
                K= PX*H'/Y;
                XX= XX + K*gamma;
                PX= PX - K*H*PX; % or PX= (eye(3) - K*H)*PX
            end
        else
            DATA.PCA(step)= 1;
            DATA.PCA_MJ(step)= 1;
        end
    end
    DATA.errorXX(step,:)= abs(xtrue - XX)';
    DATA.stdXX(step,:)= 3*sqrt(diag(PX))';
%     DATA.realPCA(step)= 1 - IA/step;
    DATA.calcPCA(step)= sum(DATA.PCA)/step;
    DATA.calcPCA_MJ(step)= sum(DATA.PCA_MJ)/step;
    
    % Plots    
    do_plots(xtrue, z, h, dtsum);
end 
% *****************  END OF MAIN LOOP    *****************
toc

% Post-processing

if SWITCH.profile, profile report, end

post_processing_and_plots(step)




