%% ================================================================
%  KOOPMAN-BASED NONLINEAR MPC WITH REAL-TIME REFRESH (10 Hz)
%  - MPC solve at 10 Hz
%  - Internal control loop at 100 Hz
%  - Koopman model refreshed each second using GP mapping
%  - Prediction horizon = 2 sec   (20 steps @ dt=0.1s)
%  - Control horizon = 1 sec      (10 steps)
% ================================================================

import casadi.*

%% USER INPUTS (plug in your Koopman matrices)
K = ...;     % r x r
B = ...;     % r x m
C = ...;     % p x r   -> must output [x;y;psi;v]

r = size(K,1);
m = size(B,2);

dt   = 0.1;      % sampling time = 10Hz
Np   = 20;       % prediction horizon (2 sec)
Nc   = 10;       % control horizon (1 sec)

% Control bounds
u_min = [-0.5; -1.0];
u_max = [ 0.5;  1.0];

goal = [5;5];    % example goal (x_goal,y_goal)

obstacles = [2 2 0.5;  3 4 0.6];  
% each row: [xo, yo, radius]

%% COST WEIGHTS
Q_goal = 50;     % weight for distance-to-goal
R_u    = diag([1 1]);   % control effort
R_du   = diag([10 10]); % control rate penalty
W_obs  = 3;             % obstacle repulsion weight
eps_obs = 0.2;          % small softening (avoid divide-by-zero)

%% BUILD CASADI OPTIMIZATION
opti = Opti();

% Decision variables
Z  = opti.variable(r, Np+1);    % lifted state trajectory
U  = opti.variable(m, Np);      % control sequence

%% Cost
J = 0;

for k = 1:Np
    % Koopman output
    yk = C * Z(:,k);
    xk = yk(1);
    ykpos = yk(2);

    % Goal cost
    J = J + Q_goal * ((xk - goal(1))^2 + (ykpos - goal(2))^2);

    % Control effort
    if k <= Nc
        J = J + U(:,k)' * R_u * U(:,k);
        if k > 1
            J = J + (U(:,k)-U(:,k-1))' * R_du * (U(:,k)-U(:,k-1));
        end
    end

    % Obstacle avoidance cost (soft barrier)
    for i = 1:size(obstacles,1)
        xo = obstacles(i,1);
        yo = obstacles(i,2);
        ro = obstacles(i,3);

        d2 = (xk - xo)^2 + (ykpos - yo)^2;
        J = J + W_obs * 1/( (d2 - ro^2 + eps_obs) );
    end
end

% assign cost
opti.minimize(J);

%% DYNAMICS CONSTRAINTS
for k = 1:Np
    uk = U(:,k);

    % Only first Nc controls are independent; later controls = last one
    if k > Nc
        uk = U(:,Nc);
    end

    % Koopman linear propagation
    opti.subject_to(Z(:,k+1) == K*Z(:,k) + B*uk);
end

%% CONTROL BOUNDS
opti.subject_to(u_min <= U(:,1:Nc) <= u_max);

%% Solver
opts.ipopt.print_level = 0;
opts.print_time = 0;
opti.solver('ipopt',opts);

%% REAL-TIME LOOP (10 Hz MPC, 100 Hz internal sim)

% INITIAL REAL STATE
x_real = [x0; y0; psi0; v0];

for outer_loop = 1:TotalSteps   % << runs at 10 Hz
    
    % --- (1) GP lifting of real state to Koopman state ----------------
    z0 = GP_lift(x_real);   % Your GP returns r×1 mean vector

    opti.set_initial(Z(:,1), z0);
    opti.subject_to(Z(:,1) == z0);

    % --- (2) Solve MPC problem ---------------------------------------
    sol = opti.solve();

    % extract FIRST control for this second's worth of execution
    u_star = sol.value(U(:,1));

    % --- (3) Apply control in a fast 100 Hz loop ----------------------
    for fast_iter = 1:10      % 10×100Hz = 1 control horizon
        x_real = nonlinear_vehicle_update(x_real, u_star, 0.01);
    end

    % --- (4) Next iteration will refresh Koopman state with GP -------
end
