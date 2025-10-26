function [event_value,isterminal,direction] = HalfCar_EventsFcn(t_hc,Z_hc,b,delta,tau,traj)
% clc;
b.verbose = false;
[~,~,fv] = HalfCarBekker_F3(t_hc,Z_hc,b,delta,tau,traj);

b.verbose = true;
event_value = [Z_hc(4) - 0.1;                 % event 1: dx < 0.1 m/s
    b.R*(Z_hc(12)+Z_hc(11)) - 0.1;      % event 2: w_r + w_f < 8 rad/s
    fv(end,3) - 20;                    % event 3: Nf > 0
    fv(end,7) - 20;                   % event 4: Nr > 0
    double(isreal(Z_hc))];
isterminal = [1;1;1;1;1];  % Halt integration
direction = [-1;-1;-1;-1;0];   % The zero can be approached from either direction
% event_value
end
