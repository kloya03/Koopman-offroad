%% ================================================================
%  KOOPMAN-BASED NONLINEAR MPC WITH REAL-TIME REFRESH (10 Hz)
%  - MPC solve at 10 Hz
%  - Internal control loop at 100 Hz
%  - Koopman model refreshed each second using GP mapping
%  - Prediction horizon = 2 sec   (20 steps @ dt=0.1s)
%  - Control horizon = 1 sec      (10 steps)
% ================================================================

clc;
clear;
addpath('/home/kloya/Documents/casadi-3.7.2-linux64-matlab2018b/')
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
    files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
    filename = fullfile(folder, files(123).name);
    model = load(filename,"MDL_fitr","A","B","Bc1","C",...
        "Cc1","K_obs","mean_std_out","rr","b","valData"); % load the file
%% USER INPUTS (plug in your Koopman matrices)
K = model.A;     % r x r
B = model.B;     % r x m
C = model.C;     % p x r   -> must output [x;y;psi;v]
Cc1 = model.Cc1;
Bc1 = model.Bc1;

r = size(K,1);
m = size(B,2);

dt   = 0.01;      % sampling time = 10Hz
Np   = 200;       % prediction horizon (2 sec)
Nc   = 50;       % control horizon (1 sec)

% Control bounds
u_min = [-0.35; 0];
u_max = [ 0.35;  135.0];

goal = [50;50];    % example goal (x_goal,y_goal)

obstacles = [20 20 5;  30 20 2];  
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
Xcart = opti.variable(3,Np+1);
%% Cost
J = 0;

for k = 1:Np
    % Koopman output
    yk = C * Z(:,k) + Cc1;
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
    opti.subject_to(Z(:,k+1) == K*Z(:,k) + B*uk + Bc1);
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
