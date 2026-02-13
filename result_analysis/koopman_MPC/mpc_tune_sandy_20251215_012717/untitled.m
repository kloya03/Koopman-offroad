clc;
clear;
addpath('/home/kloya/Documents/casadi-3.7.2-linux64-matlab2018b/')
addpath('../../../functions/utility/')
addpath('../../../functions/simulation_model')
vid = 1;
%% === Test Trajectory ===
% b.dt = 0.01;  % model at 100 Hz
% b.tstart = 0;
% b.tstop = 20;
% b.tspan = b.tstart:b.dt:b.tstop;
% b.nt = size(b.tspan);
% b.verbose = true;
% b.hr = nan;
% b.hf = nan;
% b.elev = 0;
% b.traj = 102;
% b.terrain="_sandyloam";
% run('loadCommonParams.m');

% [Data,Zhc0,events] = find_traj(b);
% save('sandy_data_102.mat');
% clc;
% clear;
% load("sandy_102_v2.mat")
load("clay_data_high_u.mat")
%% Model prediction

folder = '../../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(124).name);
models{1} = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr"); % load the file

folder = '../../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(123).name);
models{2} = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr"); % load the file
clearvars filename files folder
models{1}.filename = "F_Ksandy_on_sandy_presen.avi";
models{2}.filename = "Kclay_on_sandy_present.avi";


% for i=1:length(models)
%     ig1 = figure(1);
%     ig1.WindowState = 'maximized';   % (newer MATLAB versions)
%     [ypred,yout]  = K_RSSID_prediction(getexp(Data,1),...
%         models{i}.MDL_fitr,models{i}.A,models{i}.B,models{i}.Bc1,...
%         models{i}.C,models{i}.Cc1,models{i}.K_obs,models{i}.mean_std_out,200);
%     plot(ypred(:,1),ypred(:,2),'LineWidth',4); hold on;
%     clearvars model
% end
% plot(yout(:,1),yout(:,2),'k','LineWidth',4); grid on;
% legend('K-SSID sandy','K-SSID clay','True')
%% MPC parameters
dt_mpc   = 0.1;      % MPC step (10Hz)
Np   = 20;        % prediction horizon
Nc   = 20;        % control horizon
block = 1;
inner = dt_mpc / b.dt; % =10 Koopman steps per MPC step  100 Hz model -->10 Hz MPC
yref = Data.OutputData;   % @ 10hz
x0 = yref(1:6,:).';
yref_nHz = yref(1:inner:end, :);   % size ~201 x 6
TotalSteps = size(yref_nHz, 1);            % number of 10Hz samples
tspan = b.tstart:dt_mpc:Data.SamplingInstants(end,1);

for model_iter = 2:length(models)
    x_real = []; usol = []; xsol = [];
    %% === USER SUPPLIED MATRICES ===
    K = models{model_iter}.A;     % r x r
    B = models{model_iter}.B;     % r x m
    C = models{model_iter}.C;     % p x r    (must output [X;Y;psi;u:v:psi_dot] or whatever you track)
    Cc1 = models{model_iter}.Cc1;   % unnormalized term to get real output
    Bc1 = models{model_iter}.Bc1;    % unnormalized term to get real input effect
    r = size(K,1);
    m = size(B,2);
    p = size(C,1);

    %% Model matrices w.r.t time step
    A10 = K^inner;
    S = zeros(size(K));   % sum of K^i
    Ki = eye(size(K));
    for i = 1:inner
        S = S + Ki;      % accumulate K^i
        Ki = Ki*K;
    end
    B10  = S * B;        % effective B for 10 steps
    BC10  = S * Bc1;      % effective constant term

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
    Q   = diag([15 15 15 1 15 15]);   % tracking cost for [X,Y,psi,u,v,r]
    Ru  = diag([1e-2 1e-6]);             % input effort
    Rdu = diag([100 1]);
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

    %% === DYNAMICS CONSTRAINTS ===
    opti.subject_to(Z(:,1)     == Z0par);
    opti.subject_to(Xcart(:,1) == X0par);
    for k = 1:Np
        uk = U(:, min(k, Nc));  % hold last input after control horizon
        % ---- Koopman multirate propagation ----
        opti.subject_to( Z(:,k+1) == A10 * Z(:,k) + B10*uk + BC10 );
        % ---- Koopman body velocities ----
        yK = C*Z(:,k+1) + Cc1;
        uB = yK(1);
        vB = yK(2);
        rB = yK(3);
        psi_k = Xcart(3,k);
        % ---- Kinematic Integration (10Hz) ----
        opti.subject_to( Xcart(1,k+1) == Xcart(1,k) +...
            (uB*cos(psi_k) - vB*sin(psi_k)) * dt_mpc );
        opti.subject_to( Xcart(2,k+1) == Xcart(2,k) + ...
            (uB*sin(psi_k) + vB*cos(psi_k)) * dt_mpc );
        opti.subject_to( Xcart(3,k+1) == Xcart(3,k) + rB * dt_mpc );
    end

    %% === CONTROL CONSTRAINTS
    u_min = [-0.35; 0];
    u_max = [ 0.35;  130.0];
    opti.subject_to(u_min <= U(:,1:Nc) <= u_max);

    % === SOLVER ===
    solver_opts.ipopt.print_level = 0;
    solver_opts.print_time = 0;
    % Hard limits
    solver_opts.ipopt.max_iter     = 200;     % cap iterations
    solver_opts.ipopt.max_cpu_time = 1;    % seconds (pick your real-time budget)
    solve_opts = struct();   % empty
    opti.solver('ipopt', solver_opts, solve_opts);

    %% ================================================================
    %                  REAL-TIME EXECUTION LOOP
    % ================================================================
    % Initial real nonlinear state
    if vid==1
        outputfilename = models{model_iter}.filename;
        v = VideoWriter(outputfilename,'Motion JPEG AVI'); %
        v.FrameRate = round(length(1:block:TotalSteps)/20); % Define the frame rate (FPS)
        open(v);
    end


    x_real = Zhc0;    % real nonlinear state [X,Y,psi,u,v,r]
    xsol = x_real.';
    b.verbose = false;
    usol = []; slip_r=[];
    u_last = Data.InputData(1,:).';
    Jval = 0;
    tic; et_realsim = [];sol_time=[];

    refresh = 1; t_mpc_start = 0;
    for MPC_iter = 1:block:TotalSteps
        % MPC_iter
        rows_idx = MPC_iter + (0:Np);
        per_rows_idx = MPC_iter + (0:Np);
        % clamp all indices to T10 so we don't go out of bounds
        rows_idx(rows_idx > TotalSteps) = TotalSteps;
        
        sol1 =toc;
        % (1) Lift to Koopman state via GP
        z0 = zeros(models{model_iter}.rr,1);
        for i = 1:models{model_iter}.rr
            X0n = (x0(1,models{model_iter}.K_obs) - ...
                models{model_iter}.mean_std_out(1,:))./models{model_iter}.mean_std_out(2,:);  % normalize back for GP
            Zi_mean = predict(models{model_iter}.MDL_fitr(i).gprMDL,X0n);
            z0(i,:) = Zi_mean;
        end

        % (2) Provide the reference trajectory for next Np steps and initial condition
        Yref_val = yref_nHz(rows_idx(1:end-1),:).';
        opti.set_value(Yref, Yref_val);
        opti.set_value(Z0par, z0);
        opti.set_value(X0par, x_real(1:3));
        opti.set_value(U0prev, u_last);
        % opti.set_value(Uref, [fv(rows_idx(1:end-1),12).';fv(rows_idx(1:end-1),11).']);
        opti.set_initial(U, [Data.InputData(rows_idx(1:end-1),1).';...
            Data.InputData(rows_idx(1:end-1),2).']);
        % opti.set_initial(Xcart, Z_hc(rows_idx,1:3).');
        % (3) Solve MPC
        
        sol = opti.solve();   % actual solve
        sol2 = toc;
        sol_time(MPC_iter) = sol2-sol1;
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
        delta_mpc = @(t) zoh1d_scalar(t, t_mpc_start, dt_mpc, ...
            sol.value(U(1,:)), rows_idx(end));
        tau_mpc   = @(t) zoh1d_scalar(t, t_mpc_start, dt_mpc, ...
            sol.value(U(2,:)), rows_idx(end));
        t_mpc =  t_mpc_start + (0:block) * dt_mpc;

        wrappedEventFcn = @(t, x) HalfCar_EventsFcn(t, x, b,delta_mpc,tau_mpc,0);
        opts = odeset('RelTol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);
        sol_ode = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,...
            delta_mpc,tau_mpc,0),t_mpc,x_real,opts);
        events = [sol_ode.xe,sol_ode.ye.',sol_ode.ie];
        if ~isempty(events)
            TotalSteps_break = length(usol);
            break
        end
        X_real = deval(sol_ode, t_mpc).';   % (length(t_mpc) x nStates)
        t_hc   = t_mpc(:);
        [~,~,fv,~] = cellfun(@(t,Z) HalfCarBekker_F3(t,Z.',b,delta_mpc,tau_mpc,0),...
            num2cell(t_hc),num2cell(X_real,2),'uni',0);
        fv = cell2mat(fv);
        slip_r = [slip_r, fv(1:end-1,21).'];

        x_real = X_real(end,:);    % real nonlinear state [X,Y,psi,u,v,psi_dot,z,dz,theta,dtheta,w_f,w_r]
        x0 = X_real(end,1:6);

        et_2 = toc;
        et_realsim = [et_realsim, et_2-et_1];


        usol = [usol, u_star];
        xsol = [xsol; X_real(2:end,:)];
        Jval = Jval +J_val;

        % ---- print one-line summary ----
        fprintf(['MPC iter %3d | J = %.2f (track = %.2f, u = %.2f, du = %.2f, du0 = %.2f) ', ...
            '| iters = %d | status = %s| Terrain: %s on %d \n'], ...
            MPC_iter, J_val, J_track_val, J_u_val, J_du_val, J_du0_val,...
            stats.iter_count, stats.return_status, b.terrain,models{model_iter}.rr);

        % reconstruction prediction plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ig = figure(model_iter+1); clf
        set(ig, 'Units','pixels', 'Position',[100 100 1074 647]);  % pick a fixed size
        set(ig, 'Resize','off');
        % ig.WindowState = 'maximized';   % (newer MATLAB versions)
        lw = 5;
        axTraj = subplot(2,2,[1 3]);   % big trajectory plot (left)
        axU1   = subplot(2,2,2);       % control 1 (top-right)
        axU2   = subplot(2,2,4);       % control 2 (bottom-right)

        % (Optional) set styles once
        grid(axU1,'on'); grid(axU2,'on');
        plotControlsMPC(axU1, axU2, ...
            b.tspan, Data.InputData, ...
            rows_idx, dt_mpc, sol.value(U), ...
            usol, lw);
        % --- inside your loop ---
        plotTrajWithCar(axTraj, Data.OutputData, yfull, x0, xsol, 6, 3);

        if vid==1
            frame = getframe(gcf); % gcf gets the current figure
            writeVideo(v, frame);
        end

    end
    if vid==1
        close(v);
    end

    et_mpc1 = toc;
    et_mpc = et_mpc1 - sum(et_realsim);
    
    if ~isempty(events)
        models{model_iter}.tspan = b.tstart:dt_mpc:events(1,1);
        models{model_iter}.usol = usol;
        models{model_iter}.xsol = xsol;
        models{model_iter}.Jval = Jval./MPC_iter;
        models{model_iter}.slip_r = slip_r;
    else
        models{model_iter}.tspan = tspan;
        models{model_iter}.usol = usol(:,1:TotalSteps);
        models{model_iter}.xsol = xsol(1:TotalSteps,:);
        models{model_iter}.Jval = Jval./TotalSteps;
        models{model_iter}.slip_r = slip_r(:,1:TotalSteps);
    end
    SolTime = sum(sol_time);

    models{model_iter}.sol_time = [SolTime, mean(sol_time),et_mpc, et_mpc1];
    % clearvars usol xsol et_mpc sol_time J_val
end

% save('MPC_onsandy_20_1_v2.mat','-v7.3')

%% Plotting
cc = {'b','r'};
lw = 5;
figure(4)
subplot(2,2,1)
hold on;
plot(Data.OutputData(:,1),Data.OutputData(:,2),'k','LineWidth',lw);   hold on;
plot(models{2}.xsol(:,1),models{2}.xsol(:,2),cc{2},'LineWidth',lw-1); hold on;
plot(models{1}.xsol(:,1),models{1}.xsol(:,2),cc{1},'LineWidth',lw-1); 
grid on;
xlabel('X')
ylabel('Y')
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 25; % Set the font size to 14 points
legend('$y_{ref}$ clay','K-MPC clay','K-MPC sandy loam','interpreter','latex','fontsize',25)
% axis equal 
axis([-15 15 0 40 ])
for jj=1:3
    subplot(2,2,jj+1)
    nk = models{1}.K_obs(jj);
    hold on;
    plot(b.tspan,Data.OutputData(:,nk),'k','LineWidth',lw); hold on;
    plot(models{2}.tspan,models{2}.xsol(:,nk),cc{2},'LineWidth',lw); hold on;
    plot(models{1}.tspan,models{1}.xsol(:,nk),cc{1},'LineWidth',lw); 
    grid on;
    xlabel('Time [s]')
    ylabel(Data.OutputName{models{model_iter}.K_obs(jj)},'Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 25; % Set the font size to 14 points
    grid on;
    % hold off;
end

figure(5)
for jj=1:2
    subplot(2,2,jj)
    % hold on;
    plot(0,0,'LineWidth',lw-3);
    plot(models{2}.tspan,models{2}.usol(jj,:),cc{2},'LineWidth',lw); hold on;
    plot(models{1}.tspan(2:end),models{1}.usol(jj,:),cc{1},'LineWidth',lw);
    grid on;
    xlabel('Time [s]')
    ylabel(Data.InputName{jj},'Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 25; % Set the font size to 14 points
    hold on;
end

% plot(0,0,'LineWidth',lw-3);
figure(6)
plot(models{2}.tspan,models{2}.slip_r,cc{2},'LineWidth',lw); hold on;
plot(models{1}.tspan(2:end),models{1}.slip_r,cc{1},'LineWidth',lw); grid on;
xlabel('Time [s]')
ylabel('$s_r$','Interpreter','latex')
box on;
set(gca, 'LineWidth', 1.5)
legend('K-MPC clay','K-MPC sandy loam','interpreter','latex','fontsize',35)
ax = gca;   % Get the current axes handle
ax.FontSize = 50; % Set the font size to 14 points
hold on;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Data,Zhc0,fv,events] = find_traj(b)
% Trajectory Simulation for tracking
load("../../../datasets/Offroad_InputsSignals.mat")
delta_inp =  DELTA(:,b.traj);
tau_inp = TAU(:,b.traj);
% delta_inp = 0.35*sin(2*b.tspan); % DELTA(:,traj);
% tau_inp = 100 + 20*sin(-b.tspan); %TAU(:,traj);
delta = @(t) zoh1d_scalar(t, 0, b.dt, delta_inp, 2001);
tau   = @(t) zoh1d_scalar(t, 0, b.dt, tau_inp, 2001);
xd = 7;
Zhc0 = [0;            % X        [1]
    0;              % Y        [2]
    1.1335;%-pi + 2*pi*rand;% psi      [3]
    xd;             % dx       [4]
    0;              % dy       [5]
    0;              % dpsi     [6]
    0;              % z        [7]
    0;              % dz       [8]
    0.0037;%0.01*rand;      % theta    [9]
    0;              % dtheta  [10]
    (1 + 0.1) * xd / b.R;% omega_f [11]
    (1 + 0.1) * xd / b.R];% omega_r [12]
wrappedEventFcn = @(t, x) HalfCar_EventsFcn(t, x, b,delta,tau,b.traj);
opts = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);
[t_hc,Z_hc,te,xe,ie] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,delta,tau,b.traj),...
    b.tspan,Zhc0,opts);
[~,~,fv,~] = cellfun(@(t,Z) HalfCarBekker_F3(t,Z.',b,delta,tau,b.traj),...
    num2cell(t_hc),num2cell(Z_hc,2),'uni',0);
fv = cell2mat(fv);
% if ~isempty(events)
%     events = [te,xe,ie];
% end

Data = iddata([Z_hc(:,1:6)],[fv(:,12), fv(:,11)],b.dt);
Data.OutputName = {'x','y','$\psi$','$u$', '$v$', '$\dot{\psi}$'};%,'$Z$',...
% '$\dot{Z}$','$\theta$','$\dot{\theta}$','$\omega_f$','$\omega_r$'};
Data.InputName = {'$\delta$', '$\tau$'};
end


function plotTrajWithCar(ax, outXY, yfull, x0, xsol, L, W)
% Call every loop. Creates graphics once per-axes, then updates.

lw = 5;

% Store handles per-axes (so it works with subplots)
ud = ax.UserData;

needCreate = isempty(ud) || ~isstruct(ud) || ...
    ~isfield(ud,'p_out') || ~isgraphics(ud.p_out);

if needCreate
    cla(ax); hold(ax,'on')

    ud.p_out   = plot(ax, outXY(:,1), outXY(:,2),'k', 'LineWidth', lw);
    ud.p_xsol  = plot(ax, nan, nan, 'LineWidth', lw);
    ud.p_yfull = plot(ax, yfull(:,1), yfull(:,2), 'LineWidth', lw);
    ud.p_x0    = plot(ax, x0(:,1),  x0(:,2), 'k-o', 'LineWidth', lw-2);

    ud.carPatch = patch(ax, nan, nan, 'r', 'FaceAlpha', 0.2, ...
        'EdgeColor', 'r', 'LineWidth', 2);
    ud.carHead  = plot(ax, nan, nan, 'r-', 'LineWidth', 2);

    grid(ax,'on'); box(ax,'on'); axis(ax,'equal')
    xlabel(ax,'X'); ylabel(ax,'Y');
    set(ax, 'LineWidth', 1.5, 'FontSize', 30);

    % (Fix legend label count to match handles)
    ud.leg = legend(ax, [ud.p_out, ud.p_yfull, ud.p_xsol, ud.p_x0, ud.carPatch], ...
        {'$y_{ref}$','K-mpc pred','K-mpc','x0','car'}, ...
        'Interpreter','latex','Location','best');
    ud.leg.FontSize = 15;

    ax.UserData = ud;
end

% --- Update every call ---
set(ud.p_yfull, 'XData', yfull(:,1), 'YData', yfull(:,2), 'LineWidth', lw-2);
set(ud.p_xsol , 'XData', xsol(:,1), 'YData', xsol(:,2), 'LineWidth', lw-2);
set(ud.p_x0   , 'XData', x0(1,1),    'YData', x0(1,2),    'LineWidth', lw-2);

x   = x0(1,1);
y   = x0(1,2);
yaw = x0(1,3);   % radians

corn = [ L/2,  W/2;
    L/2, -W/2;
    -L/2, -W/2;
    -L/2,  W/2 ]';

R  = [cos(yaw) -sin(yaw);
    sin(yaw)  cos(yaw)];

cw = R*corn + [x; y];
px = [cw(1,:) cw(1,1)];
py = [cw(2,:) cw(2,1)];

set(ud.carPatch, 'XData', px, 'YData', py);

nose = R*[L/2; 0] + [x; y];
set(ud.carHead, 'XData', [x nose(1)], 'YData', [y nose(2)]);

drawnow limitrate
end


%%
function plotControlsMPC(axU1, axU2, b_tspan, InputData, rows_idx, dt_mpc, U, usol, lw)

% --- tau subplot (InputData(:,2), U(2,:), usol(2,:)) ---
updateCtrlAxes(axU1, ...
    b_tspan, InputData(:,2), ...
    rows_idx, dt_mpc, U(2,:), ...
    usol(2,:), ...
    lw, '$\tau$');

% --- delta subplot (InputData(:,1), U(1,:), usol(1,:)) ---
updateCtrlAxes(axU2, ...
    b_tspan, InputData(:,1), ...
    rows_idx, dt_mpc, U(1,:), ...
    usol(1,:), ...
    lw, '$\delta$');

drawnow limitrate
end

%%
function updateCtrlAxes(ax, t_meas, u_meas, rows_idx, dt_mpc, u_pred, u_model, lw, ylab)

ud = ax.UserData;

needCreate = isempty(ud) || ~isstruct(ud) || ~isfield(ud,'l_meas') || ~isgraphics(ud.l_meas);

if needCreate
    cla(ax); hold(ax,'on');

    % ud.l_meas  = plot(ax, nan, nan, 'LineWidth', lw);
    ud.l_pred  = plot(ax, nan, nan, 'LineWidth', lw);
    ud.l_model = plot(ax, nan, nan, 'LineWidth', lw);

    grid(ax,'on'); box(ax,'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, ylab, 'Interpreter','latex');

    set(ax, 'LineWidth', 1.5, 'FontSize', 30);

    % ud.leg = legend(ax, [ud.l_meas, ud.l_pred, ud.l_model], ...
    %     {'True','K-mpc pred','K-mpc'}, 'Location','best');
    ud.leg = legend(ax, [ ud.l_pred, ud.l_model], ...
        {'K-mpc pred','K-mpc'}, 'Location','best');
    ud.leg.FontSize = 15;
    ax.UserData = ud;
end

% Build time vectors (match your original intent)
t_pred  = rows_idx(1:end-1) * dt_mpc;
t_model = (1:numel(u_model)) * dt_mpc;

% Update line data (NO hold-on accumulation)
% set(ud.l_meas , 'XData', t_meas , 'YData', u_meas );
set(ud.l_pred , 'XData', t_pred , 'YData', u_pred );
set(ud.l_model, 'XData', t_model, 'YData', u_model);

% Optional: keep x-limits sensible
tmax = max([t_pred(:); t_model(:)]);
tmin = min([t_pred(:); t_model(:)]);
% tmax = max([t_meas(:); t_pred(:); t_model(:)]);
% tmin = min([t_meas(:); t_pred(:); t_model(:)]);
if isfinite(tmin) && isfinite(tmax) && tmax > tmin
    xlim(ax, [tmin tmax]);
end

end


% function plotTrajWithCar(fig_id, outXY, yfull, x0, xsol, L, W)
% % One-call plotting function.
% % Call it every loop; it will create graphics once and then just update.
%
% persistent h
% lw = 5;
% if isempty(h) || ~isvalid(h.fig) || h.fig.Number ~= fig_id
%     % --- First time (create) ---
%     h.fig = figure(fig_id); clf
%     h.ax  = axes('Parent', h.fig); hold(h.ax,'on')
%     h.p_out   = plot(h.ax, outXY(:,1), outXY(:,2), 'LineWidth', lw);
%     h.p_yfull = plot(h.ax, yfull(:,1), yfull(:,2), 'LineWidth', lw);
%     h.p_xsol   = plot(h.ax, nan, nan, 'LineWidth', lw);
%     h.p_x0    = plot(h.ax, x0(:,1),  x0(:,2), 'LineWidth', lw);
%     h.carPatch = patch(h.ax, nan, nan, 'r', 'FaceAlpha', 0.2, ...
%                        'EdgeColor', 'r', 'LineWidth', 2);
%     h.carHead  = plot(h.ax, nan, nan, 'r-', 'LineWidth', 2);
%
%     grid(h.ax,'on'); box(h.ax,'on'); axis(h.ax,'equal')
%     xlabel(h.ax,'X'); ylabel(h.ax,'Y');
%     set(h.ax, 'LineWidth', 1.5, 'FontSize', 30);
%
%     h.leg = legend(h.ax, ...
%     [h.p_out, h.p_yfull, h.p_xsol, h.p_x0, h.carPatch], ...
%     {'$y_{ref}$','K-mpc pred','K-mpc'},'Interpreter','latex','Location','best');
%     h.leg.FontSize = 22;   % adjust if you want
% end
%
% % --- Update every call ---
% set(h.p_xsol, 'XData', xsol(:,1), 'YData', xsol(:,2),'linewidth',lw-2);
% set(h.p_yfull,'XData', yfull(:,1), 'YData', yfull(:,2),'linewidth',lw-2);
% set(h.p_x0 , 'XData', x0(:,1), 'YData', x0(:,2),'Marker', 'o', 'Color', 'k', 'LineWidth', lw-3);
% x   = x0(:,1);%xsol(end,1);
% y   = x0(:,2);%xsol(end,2);
% yaw = x0(:,3);%xsol(end,3);   % radians
%
% corn = [ L/2,  W/2;
%          L/2, -W/2;
%         -L/2, -W/2;
%         -L/2,  W/2 ]';
%
% R  = [cos(yaw) -sin(yaw);
%       sin(yaw)  cos(yaw)];
%
% cw = R*corn + [x; y];
% px = [cw(1,:) cw(1,1)];
% py = [cw(2,:) cw(2,1)];
%
% set(h.carPatch, 'XData', px, 'YData', py);
%
% nose = R*[L/2; 0] + [x; y];
% set(h.carHead, 'XData', [x nose(1)], 'YData', [y nose(2)]);
%
% drawnow limitrate
% end


% figure(model_iter+3)
% lw = 5;
% subplot(2,2,1:2)
% plot(Data.OutputData(:,1),Data.OutputData(:,2),'LineWidth',lw); hold on;
% plot(yfull(:,1),yfull(:,2),'LineWidth',lw); hold on;
% plot(models{model_iter}.xsol(:,1),models{model_iter}.xsol(:,2),...
%     'LineWidth',lw); hold on;
% plot(x0(:,1),x0(:,2),'ok','LineWidth',lw); grid on;
% xlabel('X')
% ylabel('Y')
% box on;
% set(gca, 'LineWidth', 1.5)
% ax = gca;   % Get the current axes handle
% ax.FontSize = 30; % Set the font size to 14 points
% grid on;
% hold off;

% subplot(2,2,2)
% plot(b.tspan,Data.InputData(:,2),'LineWidth',lw); hold on;
% plot(rows_idx(1:end-1)*dt_mpc,sol.value(U(2,:)),'LineWidth',lw); hold on;
% plot((1:size(models{model_iter}.usol,2))*dt_mpc,...
%     models{model_iter}.usol(2,:),'LineWidth',lw); grid on;
% xlabel('Time [s]')
% ylabel('$\tau$','Interpreter','latex')
% box on;
% set(gca, 'LineWidth', 1.5)
% ax = gca;   % Get the current axes handle
% ax.FontSize = 30; % Set the font size to 14 points
% grid on;
% %
% subplot(2,2,4)
% plot(b.tspan,Data.InputData(:,1),'LineWidth',lw); hold on;
% plot(rows_idx(1:end-1)*dt_mpc,sol.value(U(1,:)),'LineWidth',lw); hold on;
% plot((1:size(models{model_iter}.usol,2))*dt_mpc,...
%     models{model_iter}.usol(1,:),'LineWidth',lw); grid on;
% xlabel('Time [s]')
% ylabel('$\delta$','Interpreter','latex')
% box on;
% set(gca, 'LineWidth', 1.5)
% ax = gca;   % Get the current axes handle
% ax.FontSize = 30; % Set the font size to 14 points
% grid on;

% for jj=1:3
%     subplot(2,2,jj+1)
%     plot(b.tspan,Z_hc(:,models{model_iter}.K_obs(jj)),'LineWidth',lw); hold on;
%     plot(b.tspan,models{model_iter}.xsol(:,models{model_iter}.K_obs(jj)),'LineWidth',lw); grid on;
%     xlabel('Time [s]')
%     ylabel(Data.OutputName{models{model_iter}.K_obs(jj)},'Interpreter','latex')
%     box on;
%     set(gca, 'LineWidth', 1.5)
%     ax = gca;   % Get the current axes handle
%     ax.FontSize = 30; % Set the font size to 14 points
%     grid on;
%     hold off;
% end