clc;
clear;
addpath('/home/kloya/Documents/casadi-3.7.2-linux64-matlab2018b/')
addpath('../../../functions/utility/')
addpath('../../../functions/simulation_model')

%% === Test Trajectory ===
b.dt = 0.01;  % model at 100 Hz
b.tstart = 0;
b.tstop = 20;
b.tspan = b.tstart:b.dt:b.tstop;
b.nt = size(b.tspan);
b.verbose = true;
b.hr = nan;
b.hf = nan;
b.elev = 0;
b.terrain="_sandyloam";
run('loadCommonParams.m');

% [Data,Zhc0] = find_traj(b);
% save('test_data_1.mat','Data','Zhc0');
load("test_data_1.mat")
%% Model prediction

folder = '../../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(124).name);
models{1} = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","b"); % load the file

folder = '../../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(123).name);
models{2} = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","b"); % load the file
clearvars filename files folder
models{1}.filename = "sandy_sandy_MPC_1.avi";
models{2}.filename = "sandy_clay_MPC_1.avi";


% for i=1:length(models)
%     figure(1)
%     [ypred,yout]  = K_RSSID_prediction(getexp(Data,1),...
%         models{i}.MDL_fitr,models{i}.A,models{i}.B,models{i}.Bc1,...
%         models{i}.C,models{i}.Cc1,models{i}.K_obs,models{i}.mean_std_out,250);
%     plot(ypred(:,1),ypred(:,2),'LineWidth',4); hold on;
%     clearvars model
% end
% plot(yout(:,1),yout(:,2),'k','LineWidth',4); grid on;

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

for model_iter = 1:length(models)

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
    
    outputfilename = models{model_iter}.filename;
    v = VideoWriter(outputfilename,'Motion JPEG AVI'); %
    v.FrameRate = round(length(1:block:TotalSteps)/20); % Define the frame rate (FPS)
    open(v); 

    % (1) Lift to Koopman state via GP
    z0 = zeros(models{model_iter}.rr,1);
    for i = 1:models{model_iter}.rr
        X0n = (x0(1,models{model_iter}.K_obs) - models{model_iter}.mean_std_out(1,:))./models{model_iter}.mean_std_out(2,:);  % normalize back for GP
        Zi_mean = predict(models{model_iter}.MDL_fitr(i).gprMDL,X0n);
        z0(i,:) = Zi_mean;
    end
    ytest(:,i) = x0(1,models{model_iter}.K_obs).' - C*z0 - Cc1;
    x_real = Zhc0;    % real nonlinear state [X,Y,psi,u,v,r]
    xsol = x_real.';
    b.verbose = false;
    usol = []; u_last = [50;0];Jval = 0;
    tic; et_realsim = [];sol_time=0;

    refresh = 1; t_mpc_start = 0;
    for MPC_iter = 1:block:TotalSteps
        % MPC_iter
        rows_idx = MPC_iter + (0:Np);
        per_rows_idx = MPC_iter + (0:Np);
        % clamp all indices to T10 so we don't go out of bounds
        rows_idx(rows_idx > TotalSteps) = TotalSteps;

        % (2) Provide the reference trajectory for next Np steps and initial
        % condition
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
        sol1 =toc;
        sol = opti.solve();   % actual solve
        sol2 = toc;
        sol_time = sol_time + sol2-sol1;
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
        opts = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);
        [t_hc,X_real,te,xe,ie] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,...
            delta_mpc,tau_mpc,0),t_mpc,x_real,opts);
        events = [te,xe,ie];
        x_real = X_real(end,:);    % real nonlinear state [X,Y,psi,u,v,psi_dot,z,dz,theta,dtheta,w_f,w_r]
        x0 = X_real(end,1:6);

        % (1) Lift to Koopman state via GP
        z0 = zeros(models{model_iter}.rr,1);
        for i = 1:models{model_iter}.rr
            X0n = (x0(1,models{model_iter}.K_obs) - ...
                models{model_iter}.mean_std_out(1,:))./models{model_iter}.mean_std_out(2,:);  % normalize back for GP
            Zi_mean = predict(models{model_iter}.MDL_fitr(i).gprMDL,X0n);
            z0(i,:) = Zi_mean;
        end

        et_2 = toc;
        et_realsim = [et_realsim, et_2-et_1];


        usol = [usol, u_star];
        xsol = [xsol; X_real(2:end,:)];
        Jval = Jval +J_val;

        % ---- print one-line summary ----
        fprintf(['MPC iter %3d | J = %.2f (track = %.2f, u = %.2f, du = %.2f, du0 = %.2f) ', ...
            '| iters = %d | status = %s\n'], ...
            MPC_iter, J_val, J_track_val, J_u_val, J_du_val, J_du0_val,...
            stats.iter_count, stats.return_status);

        % reconstruction prediction plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plotTrajWithCar(model_iter+1, Data.OutputData, yfull,...
            x0, xsol, 9, 4);


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
        % 
        % subplot(2,2,3)
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
        % 
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
        frame = getframe(gcf); % gcf gets the current figure
        writeVideo(v, frame);
    end
    et_mpc1 = toc;
    et_mpc = et_mpc1 - sum(et_realsim);
    close(v);
    models{model_iter}.usol = usol(:,1:TotalSteps);
    models{model_iter}.xsol = xsol(1:TotalSteps,:);
    models{model_iter}.Jval = Jval./TotalSteps;
    models{model_iter}.sol_time = [sol_time, et_mpc, et_mpc1];
    % clearvars usol xsol et_mpc sol_time J_val
end

save('MPC_onSL_1.mat')

    %% Plotting
cc = {'b','r'};
lw = 5;
figure(4)
subplot(2,2,1)
hold on;
plot(Data.OutputData(:,1),Data.OutputData(:,2),'k','LineWidth',lw);   hold on;
plot(models{1}.xsol(:,1),models{1}.xsol(:,2),cc{1},'LineWidth',lw-1); hold on;
plot(models{2}.xsol(:,1),models{2}.xsol(:,2),cc{2},'LineWidth',lw-1); grid on;
xlabel('X')
ylabel('Y')
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 30; % Set the font size to 14 points
legend('$y_{ref}$','K-MPC sandy loam','K-MPC clay','interpreter','latex','fontsize',20)

for jj=1:3
    subplot(2,2,jj+1)
    hold on;
    plot(b.tspan,Data.OutputData(:,models{1}.K_obs(jj)),'k','LineWidth',lw); hold on;
    plot(b.tspan,models{1}.xsol(:,models{1}.K_obs(jj)),cc{1},'LineWidth',lw); hold on;
    plot(b.tspan,models{2}.xsol(:,models{2}.K_obs(jj)),cc{2},'LineWidth',lw); grid on;
    xlabel('Time [s]')
    ylabel(Data.OutputName{models{model_iter}.K_obs(jj)},'Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    grid on;
    % hold off;
end

figure(5)
for jj=1:2
    subplot(1,2,jj)
    % hold on;
    % plot(0,0,'LineWidth',lw-3);
    plot(b.tspan,models{1}.usol(jj,:),cc{1},'LineWidth',lw); hold on;
    plot(b.tspan,models{2}.usol(jj,:),cc{2},'LineWidth',lw); grid on;
    xlabel('Time [s]')
    ylabel(Data.InputName{jj},'Interpreter','latex')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    hold on;
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Data,Zhc0,events] = find_traj(b)
% Trajectory Simulation for tracking
% load("../../../datasets/Offroad_InputsSignals.mat")
delta_inp = 0.35*sin(2*b.tspan); % DELTA(:,traj);
tau_inp = 80 + 30*sin(2*b.tspan); %TAU(:,traj);
delta = @(t) zoh1d_scalar(t, 0, b.dt, delta_inp, 2001);
tau   = @(t) zoh1d_scalar(t, 0, b.dt, tau_inp, 2001);
xd = 5;
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
wrappedEventFcn = @(t, x) HalfCar_EventsFcn(t, x, b,delta,tau,0);
opts = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);
[t_hc,Z_hc,te,xe,ie] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,delta,tau,0),...
    b.tspan,Zhc0,opts);
[~,~,fv,~] = cellfun(@(t,Z) HalfCarBekker_F3(t,Z.',b,delta,tau,0),...
    num2cell(t_hc),num2cell(Z_hc,2),'uni',0);
fv = cell2mat(fv);
events = [te,xe,ie];
Data = iddata([Z_hc(:,1:6)],[fv(:,12), fv(:,11)],b.dt);
Data.OutputName = {'x','y','$\psi$','$u$', '$v$', '$\dot{\psi}$'};%,'$Z$',...
% '$\dot{Z}$','$\theta$','$\dot{\theta}$','$\omega_f$','$\omega_r$'};
Data.InputName = {'$\delta$', '$\tau$'};
end


function plotTrajWithCar(fig_id, outXY, yfull, x0, xsol, L, W)
% One-call plotting function.
% Call it every loop; it will create graphics once and then just update.

persistent h
lw = 5;
if isempty(h) || ~isvalid(h.fig) || h.fig.Number ~= fig_id
    % --- First time (create) ---
    h.fig = figure(fig_id); clf
    h.ax  = axes('Parent', h.fig); hold(h.ax,'on')
    h.p_out   = plot(h.ax, outXY(:,1), outXY(:,2), 'LineWidth', lw);
    h.p_yfull = plot(h.ax, yfull(:,1), yfull(:,2), 'LineWidth', lw);
    h.p_xsol   = plot(h.ax, nan, nan, 'LineWidth', lw);
    h.p_x0    = plot(h.ax, x0(:,1),  x0(:,2), 'LineWidth', lw);
    h.carPatch = patch(h.ax, nan, nan, 'r', 'FaceAlpha', 0.2, ...
                       'EdgeColor', 'r', 'LineWidth', 2);
    h.carHead  = plot(h.ax, nan, nan, 'r-', 'LineWidth', 2);

    grid(h.ax,'on'); box(h.ax,'on'); axis(h.ax,'equal')
    xlabel(h.ax,'X'); ylabel(h.ax,'Y');
    set(h.ax, 'LineWidth', 1.5, 'FontSize', 30);

    h.leg = legend(h.ax, ...
    [h.p_out, h.p_yfull, h.p_xsol, h.p_x0, h.carPatch], ...
    {'$y_{ref}$','K-mpc pred','K-mpc'},'Interpreter','latex','Location','best');
    h.leg.FontSize = 22;   % adjust if you want
end

% --- Update every call ---
set(h.p_xsol, 'XData', xsol(:,1), 'YData', xsol(:,2),'linewidth',lw-2);
set(h.p_yfull,'XData', yfull(:,1), 'YData', yfull(:,2),'linewidth',lw-2);
set(h.p_x0 , 'XData', x0(:,1), 'YData', x0(:,2),'Marker', 'o', 'Color', 'k', 'LineWidth', lw-3);
x   = x0(:,1);%xsol(end,1);
y   = x0(:,2);%xsol(end,2);
yaw = x0(:,3);%xsol(end,3);   % radians

corn = [ L/2,  W/2;
         L/2, -W/2;
        -L/2, -W/2;
        -L/2,  W/2 ]';

R  = [cos(yaw) -sin(yaw);
      sin(yaw)  cos(yaw)];

cw = R*corn + [x; y];
px = [cw(1,:) cw(1,1)];
py = [cw(2,:) cw(2,1)];

set(h.carPatch, 'XData', px, 'YData', py);

nose = R*[L/2; 0] + [x; y];
set(h.carHead, 'XData', [x nose(1)], 'YData', [y nose(2)]);

drawnow limitrate
end
