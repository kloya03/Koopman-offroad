%% ================================================================
%  KOOPMAN-BASED MPC FOR TIME-VARYING TRAJECTORY TRACKING
%  -------------------------------------------------------
%  - Uses Koopman model: z_{k+1} = K z_k + B u_k
%  - Output mapping:      y_k = C z_k   (e.g., [x;y;psi;v])
%  - Tracks a reference trajectory y_ref(:,k)
%  - Works with GP lift: z0 = GP_lift(x_real)
%  ---------------------------------------------------------------
%  Author: Kartik Loya (Koopman Online Learning Framework)
% ================================================================
% KOOPMAN MPC WITH KINEMATICS + NORMALIZATION BIAS
%        -------------------------------------------------
%  Koopman is trained at 100Hz
%  MPC runs at 10Hz  → compress 10 Koopman steps per MPC step
%
%  State:
%     Z      = Koopman lifted state (dimension r)
%     Xcart  = [X; Y; psi]  (kinematic states)
%
%  Output:
%     yK = C*Z + Cc1   (body-frame u,v,r)
%     yFull = [X;Y;psi; u;v;psi_dot]
%
%  Dynamics:
%     Z⁺ = A10 Z + B u + C10     (100→10 Hz compression)
%     X⁺ = X + f_kin(uB,vB,rB,psi) dt
%
% ================================================================
clc;
clear;

addpath('/home/kloya/Documents/casadi-3.7.2-linux64-matlab2018b/')
addpath('../../functions/utility/')
addpath('../../functions/simulation_model')
load("../../datasets/Offroad_InputsSignals.mat")
% load('sine_traj_clay.mat')
model ="clay";
if model == "sandy"
    folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
    files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
    filename = fullfile(folder, files(124).name);
    model = load(filename,"MDL_fitr","A","B","Bc1","C",...
        "Cc1","K_obs","mean_std_out","rr","b","valData"); % load the file
elseif model =="clay"
    folder = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
    files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
    filename = fullfile(folder, files(123).name);
    model = load(filename,"MDL_fitr","A","B","Bc1","C",...
        "Cc1","K_obs","mean_std_out","rr","b","valData"); % load the file
end


%% === TIMING PARAMETERS ===
dt_mpc   = 0.01;      % MPC step (10Hz)
b.dt = 0.01;  % model at 100 Hz
b.tstart = 0;
b.tstop = 20;
b.tspan = b.tstart:dt_mpc:b.tstop;
b.nt = size(b.tspan);
%% Trajectory Simulation for tracking

% b.elev = 0;
% b.terrain="_sandyloam";
% run('loadCommonParams.m');
% b.verbose = true;
% b.hr = nan;
% b.hf = nan;
% traj = randi(1000)  %102
% traj_name(:,traj)
% delta_inp = 0.35*sin(2*b.tspan); % DELTA(:,traj);
% tau_inp = 80 + 30*sin(2*b.tspan); %TAU(:,traj);
% 
% delta = @(t) zoh1d_scalar(t, 0, b.dt, delta_inp, 2001);
% tau   = @(t) zoh1d_scalar(t, 0, b.dt, tau_inp, 2001);
% xd = 5;
% Zhc0 = [0;            % X        [1]
%     0;              % Y        [2]
%     1.1335;%-pi + 2*pi*rand;% psi      [3]
%     xd;             % dx       [4]
%     0;              % dy       [5]
%     0;              % dpsi     [6]
%     0;              % z        [7]
%     0;              % dz       [8]
%     0.0037;%0.01*rand;      % theta    [9]
%     0;              % dtheta  [10]
%     (1 + 0.1*rand) * xd / b.R;% omega_f [11]
%     (1 + 0.1*rand) * xd / b.R];% omega_r [12]
% wrappedEventFcn = @(t, x) HalfCar_EventsFcn(t, x, b,delta,tau,traj);
% opts = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);
% [t_hc,Z_hc,te,xe,ie] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,delta,tau,traj),b.tspan,Zhc0,opts);
% [~,~,fv,~] = cellfun(@(t,Z) HalfCarBekker_F3(t,Z.',b,delta,tau,traj),...
%     num2cell(t_hc),num2cell(Z_hc,2),'uni',0);
% fv = cell2mat(fv);
% events = [te,xe,ie];
% Data = iddata([Z_hc(:,1:6)],[fv(:,12), fv(:,11)],dt_mpc);
% Data.OutputName = {'x','y','$\psi$','$u$', '$v$', '$\dot{\psi}$'};
% Data.InputName = {'$\delta$', '$\tau$'};
figure(4)
[ypred,yout]  = K_RSSID_prediction(getexp(Data,1),...
    model.MDL_fitr,model.A,model.B,model.Bc1,...
    model.C,model.Cc1,model.K_obs,model.mean_std_out,250);
% figure(1)
hold on;
% plot(yout(:,1),yout(:,2),'LineWidth',4);hold on;
plot(ypred(:,1),ypred(:,2),'LineWidth',4); grid on; 
clc;
% save('sine_traj_clay.mat')
%%
Np   = 200;        % prediction horizon
Nc   = 200;        % control horizon
block = 50;
yref_10Hz = Z_hc(:,1:6);   % @ 10hz
x0 = Zhc0(1:6,:).';
% yref_10Hz = yref(1:10:end, :);   % size ~201 x 6
TotalSteps = size(yref_10Hz, 1);            % number of 10Hz samples

%% === USER SUPPLIED MATRICES ===
K = model.A;     % r x r
B = model.B;     % r x m
C = model.C;     % p x r    (must output [X;Y;psi;u:v:psi_dot] or whatever you track)
Cc1 = model.Cc1;   % unnormalized term to get real output
Bc1 = model.Bc1;    % unnormalized term to get real input effect

r = size(K,1);
m = size(B,2);
p = size(C,1);

%% === BUILD OPTIMIZER ===
opti = casadi.Opti();
Z     = opti.variable(r, Np+1);    % Koopman states
Xcart = opti.variable(3, Np+1);    % [X;Y;psi]
U     = opti.variable(m, Np);      % control inputs
Yref  = opti.parameter(6, Np);     % full reference: [X;Y;psi; u;v;r]
Z0par = opti.parameter(r,1);
X0par = opti.parameter(3,1);
U0prev = opti.parameter(m,1);
% Uref = opti.parameter(2,Np);
% === COST WEIGHTS ===
Q   = diag([1 1 1 1 15 15]);   % tracking cost for [X,Y,psi,u,v,r]
Ru  = diag([1e-2 1e-6]);             % input effort
Rdu = diag([1000 1]);
% R0du = diag([1000 1]);
% === COST FUNCTION ===
J_track = 0;   % tracking error
J_u     = 0;   % input magnitude
J_du    = 0;   % input rate
J_du0 = 0;
for k = 1:Np
    % Koopman output (body-frame)
    yK = C*Z(:,k) + Cc1;   % gives [u;v;r]
    % Full predicted output
    yFull = [ Xcart(:,k) ; yK ];
    % Tracking error cost
    e = yFull - Yref(:,k);
    J_track = J_track + e' * Q * e;
    % Input cost
    if k <= Nc
        % Control effort
        J_u = J_u + U(:,k)' * Ru * U(:,k);
        % du = U(:,k) - Uref(:,k);
        % Smooth input (Δu)
        
        if k > 1
            du = U(:,k) -U(:,k-1);
            J_du = J_du + du' * Rdu * du;
        else
            du = U(:,1) - U0prev;
            J_du = J_du + du' * Rdu * du;
        end

    end
end
J = J_track + J_u + J_du;% + J_du0;
opti.minimize(J);

% === DYNAMICS CONSTRAINTS ===
opti.subject_to(Z(:,1)     == Z0par);
opti.subject_to(Xcart(:,1) == X0par);
for k = 1:Np
    uk = U(:, min(k, Nc));  % hold last input after control horizon
    % ---- Koopman multirate propagation ----
    opti.subject_to( Z(:,k+1) == K * Z(:,k) + B*uk + Bc1 );
    % ---- Koopman body velocities ----
    yK = C*Z(:,k+1) + Cc1;
    uB = yK(1);
    vB = yK(2);
    rB = yK(3);
    psi_k = Xcart(3,k);
    % ---- Kinematic Integration (10Hz) ----
    opti.subject_to( Xcart(1,k+1) == Xcart(1,k) + (uB*cos(psi_k) - vB*sin(psi_k)) * dt_mpc );
    opti.subject_to( Xcart(2,k+1) == Xcart(2,k) + (uB*sin(psi_k) + vB*cos(psi_k)) * dt_mpc );
    opti.subject_to( Xcart(3,k+1) == Xcart(3,k) + rB * dt_mpc );
end

% === CONTROL CONSTRAINTS
u_min = [-0.35; 0];
u_max = [ 0.35;  130.0];
opti.subject_to(u_min <= U(:,1:Nc) <= u_max);

% === SOLVER ===
solver_opts.ipopt.print_level = 0;
solver_opts.print_time = 0;
solve_opts = struct();   % empty
opti.solver('ipopt', solver_opts, solve_opts);

%% ================================================================
%                  REAL-TIME EXECUTION LOOP
% ================================================================
% Initial real nonlinear state

% (1) Lift to Koopman state via GP
z0 = zeros(model.rr,1);
for i = 1:model.rr
    X0n = (x0(1,model.K_obs) - model.mean_std_out(1,:))./model.mean_std_out(2,:);  % normalize back for GP
    Zi_mean = predict(model.MDL_fitr(i).gprMDL,X0n);
    z0(i,:) = Zi_mean;
end
ytest(:,i) = x0(1,model.K_obs).' - C*z0 - Cc1;
x_real = Zhc0;    % real nonlinear state [X,Y,psi,u,v,r]
xsol = x_real.';
b.verbose = false;
usol = []; u_last = [0;0];
tic; et_realsim = [];
opts = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);
wrappedEventFcn = @(t, x) HalfCar_EventsFcn(t, x, b,delta_mpc,tau_mpc,traj);
refresh = 1; t_mpc_start = 0;
for MPC_iter = 1:block:TotalSteps
    % MPC_iter
    rows_idx = MPC_iter + (0:Np);
    per_rows_idx = MPC_iter + (0:Np);
    % clamp all indices to T10 so we don't go out of bounds
    rows_idx(rows_idx > TotalSteps) = TotalSteps;

    % (2) Provide the reference trajectory for next Np steps and initial
    % condition
    Yref_val = yref_10Hz(rows_idx(1:end-1),:).';
    opti.set_value(Yref, Yref_val);
    opti.set_value(Z0par, z0);
    opti.set_value(X0par, x_real(1:3));
    opti.set_value(U0prev, u_last);
    % opti.set_value(Uref, [fv(rows_idx(1:end-1),12).';fv(rows_idx(1:end-1),11).']);
    opti.set_initial(U, [fv(rows_idx(1:end-1),12).';fv(rows_idx(1:end-1),11).']);
    % opti.set_initial(Xcart, Z_hc(rows_idx,1:3).');
    % (3) Solve MPC
    sol = opti.solve();   % actual solve
    % ---- evaluate cost at solution ----
    J_val       = sol.value(J);
    J_track_val = sol.value(J_track);
    J_u_val     = sol.value(J_u);
    J_du_val    = sol.value(J_du);
    J_du0_val    = sol.value(J_du0);
    stats = opti.stats();
    % Extract u_nc
    u_star = sol.value(U(:,1:block));
    u_last = u_star(:,end);
    X_pos = sol.value(Xcart);
    yfull = [X_pos; C*sol.value(Z)+Cc1].';

    % real system feedback
    et_1 = toc;
    % (4) Apply control in high-rate real dynamics
    t_mpc_start = (MPC_iter-1)*dt_mpc;
    delta_mpc = @(t) zoh1d_scalar(t, t_mpc_start, dt_mpc, sol.value(U(1,:)), rows_idx(end));
    tau_mpc   = @(t) zoh1d_scalar(t, t_mpc_start, dt_mpc, sol.value(U(2,:)), rows_idx(end));
    t_mpc =  t_mpc_start + (0:block) * dt_mpc;

    [t_hc,X_real,te,xe,ie] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,delta_mpc,tau_mpc,traj),t_mpc,x_real,opts);
    events = [te,xe,ie];
    x_real = X_real(end,:);    % real nonlinear state [X,Y,psi,u,v,r]
    x0 = X_real(end,1:6);
    % (1) Lift to Koopman state via GP
    z0 = zeros(model.rr,1);
    for i = 1:model.rr
        X0n = (x0(1,model.K_obs) - model.mean_std_out(1,:))./model.mean_std_out(2,:);  % normalize back for GP
        Zi_mean = predict(model.MDL_fitr(i).gprMDL,X0n);
        z0(i,:) = Zi_mean;
    end
    ytest = [ytest,x0(1,model.K_obs).' - C*z0 - Cc1];
    et_2 = toc;
    et_realsim = [et_realsim, et_2-et_1];


    usol = [usol, u_star];
    xsol = [xsol; X_real(2:end,:)];

    % ---- print one-line summary ----
    fprintf(['MPC iter %3d | J = %.2f (track = %.2f, u = %.2f, du = %.2f, du0 = %.2f) ', ...
        '| iters = %d | status = %s\n'], ...
        MPC_iter, J_val, J_track_val, J_u_val, J_du_val, J_du0_val,...
        stats.iter_count, stats.return_status);

    % reconstruction prediction plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    clf(1)
    lw = 5;
    
    subplot(2,2,1:2)
    plot(Z_hc(:,1),Z_hc(:,2),'LineWidth',lw); hold on;
    plot(yfull(:,1),yfull(:,2),'LineWidth',lw); hold on;
    plot(xsol(:,1),xsol(:,2),'LineWidth',lw); hold on;
    plot(x0(:,1),x0(:,2),'ok','LineWidth',lw); grid on;
    xlabel('X')
    ylabel('Y')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    grid on;
    hold off;

    subplot(2,2,3)
    plot(b.tspan,fv(:,11),'LineWidth',lw); hold on;
    plot(rows_idx(1:end-1)*dt_mpc,sol.value(U(2,:)),'LineWidth',lw); hold on;
    plot((1:size(usol,2))*dt_mpc,usol(2,:),'LineWidth',lw); grid on;
    xlabel('Time [s]')
    ylabel('$\tau$','Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    grid on;

    subplot(2,2,4)
    plot(b.tspan,fv(:,12),'LineWidth',lw); hold on;
    plot(rows_idx(1:end-1)*dt_mpc,sol.value(U(1,:)),'LineWidth',lw); hold on;
    plot((1:size(usol,2))*dt_mpc,usol(1,:),'LineWidth',lw); grid on;
    xlabel('Time [s]')
    ylabel('$\delta$','Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    grid on;
    % for jj=1:3
    %     subplot(2,2,jj+1)
    %     plot(b.tspan,Z_hc(:,model.K_obs(jj)),'LineWidth',lw); hold on;
    %     plot(b.tspan,xsol(:,model.K_obs(jj)),'LineWidth',lw); grid on;
    %     xlabel('Time [s]')
    %     ylabel(model.valData.OutputName{model.K_obs(jj)},'Interpreter','latex')
    %     box on;
    %     set(gca, 'LineWidth', 1.5)
    %     ax = gca;   % Get the current axes handle
    %     ax.FontSize = 30; % Set the font size to 14 points
    %     grid on;
    %     hold off;
    % end
end
et_mpc1 = toc;
et_mpc = et_mpc1 - sum(et_realsim);
%% Plotting
cc = 'b';
lw = 5;
figure(2)
subplot(2,2,1)
hold on;
plot(Z_hc(:,1),Z_hc(:,2),'k','LineWidth',lw);  hold on;
plot(xsol(1:2001,1),xsol(1:2001,2),cc,'LineWidth',lw-1); grid on;
xlabel('X')
ylabel('Y')
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 30; % Set the font size to 14 points
legend('$y_{ref}$','K-MPC','interpreter','latex','fontsize',30)
grid on;
% hold off;

for jj=1:3
    subplot(2,2,jj+1)
    hold on;
    plot(b.tspan,Z_hc(:,model.K_obs(jj)),'k','LineWidth',lw); hold on;
    plot(b.tspan,xsol(1:2001,model.K_obs(jj)),cc,'LineWidth',lw); grid on;
    xlabel('Time [s]')
    ylabel(model.valData.OutputName{model.K_obs(jj)},'Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    grid on;
    % hold off;
end

figure(3)
for jj=1:2
    subplot(1,2,jj)
    % hold on;
    % plot(0,0,'LineWidth',lw-3);
    plot(b.tspan,usol(jj,1:2001),'r','LineWidth',lw); grid on;
    xlabel('Time [s]')
    ylabel(model.valData.InputName{jj},'Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    grid on;
    hold off;
end

save('sine_traj_sandy_on_clay.mat')
%%
%
% model ='sandy';
% testset = 'sandy';
% if model == 'sandy'
%     folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
%     files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
%     filename = fullfile(folder, files(123).name);
%     model = load(filename,"MDL_fitr","A","B","Bc1","C",...
%         "Cc1","K_obs","mean_std_out","rr","b"); % load the file
% elseif model =='clay'
%     folder = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
%     files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
%     filename = fullfile(folder, files(124).name);
%     model = load(filename,"MDL_fitr","A","B","Bc1","C",...
%         "Cc1","K_obs","mean_std_out","rr","b"); % load the file
% end
% if testset=='sandy'
%     load('../../datasets/sandyloam_100hz_no_elev_experiment_1579.mat',...
%         'valData','exp_x','allindices');
% elseif testset=='clay'
%     load('../../datasets/sandyloam_100hz_no_elev_experiment_1579.mat',...
%         'valData','exp_x','allindices');
% end
% %
% % dts = 14,5
% dts = 14%randi(50)
% data =  getexp(valData,dts);
% yref = data.OutputData;
% x0 = yref(1,:);
% yref_10Hz = yref(1:10:end, :);   % size ~201 x 6
% TotalSteps = size(yref_10Hz, 1);            % number of 10Hz samples
% % plot(yref(:,1),yref(:,2))
% % axis equal
% % hold off
% Ts = dt_mpc;               % your Koopman dt
% Q = diag([5 1 1]);
% Qx = (model.C.')*Q*model.C;
% R = diag([1e-4 1e-6]);
% D = zeros(size(model.C,1), size(model.B,2));
% sys = ss(model.A,model.B,model.C,D,Ts);
% [K,S,E] = lqr(sys, Qx, R);   % discrete-time LQR
%
% u = -Kx