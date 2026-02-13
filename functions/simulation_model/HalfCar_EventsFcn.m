function [event_value,isterminal,direction] = HalfCar_EventsFcn(t_hc,Z_hc,b,delta,tau,traj)
% clc;
b.verbose = false;
[~,~,fv] = HalfCarBekker_F3(t_hc,Z_hc,b,delta,tau,traj);
% clc;
% fprintf('slip = %2.2f',fv(end,21));
% b.verbose = true;
event_value = [Z_hc(4) - 0.1;                 % event 1: dx < 0.1 m/s from pos to neg
    b.R*(Z_hc(12)+Z_hc(11)) - 0.1;      % event 2: w_r + w_f < 8 rad/s from pos to neg
    fv(end,3) - 20;                    % event 3: Nf > 0  from pos to neg
    fv(end,7) - 20;                   % event 4: Nr > 0   from pos to neg
    fv(end,21)-0.95;                    % event 4: slip_r > .95  from neg to pos
    double(isreal(Z_hc))];
isterminal = [1;1;1;1;1;1];  % Halt integration
direction = [-1;-1;-1;-1;1;0];   
% event_value
end
