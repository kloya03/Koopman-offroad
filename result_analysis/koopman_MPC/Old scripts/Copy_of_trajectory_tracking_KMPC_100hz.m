%% ================================================================
%  Koopman MPC Weight Tuning (No Uref cost) + Save Video/Results
%  ------------------------------------------------------------
%  - Builds ONE CasADi Opti problem with weight PARAMETERS
%  - Loops over candidate Q,R,Rdu sets
%  - Normalizes input penalties using known ranges:
%      delta in [-0.35, 0.35]
%      tau   in [20, 135]
%  - Saves each run to a new folder
% ================================================================
clc; clear;

addpath('/home/kloya/Documents/casadi-3.7.2-linux64-matlab2018b/')
addpath('../../functions/utility/')
addpath('../../functions/simulation_model')

load("../../datasets/Offroad_InputsSignals.mat")

%% ================= SIMULATE ONE REFERENCE TRAJ =================
b = model.b;
b.dt = dt_mpc;
b.tstart = 0;
b.tstop  = 20;
b.tspan  = b.tstart:dt_mpc:b.tstop;
b.verbose = 1;
b.elev = 0;
if soil=="sandy"; b.terrain="_sandyloam"; else; b.terrain="_clay"; end
run('loadCommonParams.m');

traj = randi(1000);
delta_inp = DELTA(:,traj);
tau_inp   = TAU(:,traj);

delta = @(t) zoh1d_scalar(t, 0, b.dt, delta_inp, 2001);
tau   = @(t) zoh1d_scalar(t, 0, b.dt, tau_inp, 2001);

xd = 5;
Zhc0 = [0; 0; -pi + 2*pi*rand; xd; 0; 0; 0; 0; 0.01*rand; 0; ...
        (1 + 0.1*rand) * xd / b.R; (1 + 0.1*rand) * xd / b.R];

wrappedEventFcn0 = @(t, x) HalfCar_EventsFcn(t, x, b, delta, tau, traj);
opts0 = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn0);

[t_hc,Z_hc] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,delta,tau,traj), b.tspan, Zhc0, opts0);

[~,~,fv,~] = cellfun(@(t,Z) HalfCarBekker_F3(t,Z.',b,delta,tau,traj), ...
    num2cell(t_hc), num2cell(Z_hc,2), 'uni',0);
fv = cell2mat(fv);

% Reference: [X Y psi u v r]
yref = Z_hc(:,1:6);
TotalSteps = size(yref,1);
save('tets')


%% ================= USER SETTINGS =================
soil = "sandy";                 % "sandy" or "clay"
model_idx = 124;                % sandy default
if soil=="clay"; model_idx=123; end

dt_mpc = 0.01;                  % keep as your script (0.01 s)
Np = 300;
Nc = 300;
block = 100;

MAKE_VIDEO = true;              % saves mpc_run.mp4 for each case (can be large!)
VIDEO_FPS  = 20;

% Input bounds
u_min = [-0.35; 0];
u_max = [ 0.35; 135];

% Input scales for cost normalization (so delta/tau comparable)
u_scale  = [0.35; 135];         % magnitude scaling
du_scale = [0.70; 115];         % range length scaling (0.35-(-0.35)=0.7, 135-20=115)

%% ============== OUTPUT ROOT FOLDER ===============
run_tag = datestr(now,'yyyymmdd_HHMMSS');
rootdir = fullfile(pwd, ['mpc_tune_' char(soil) '_' datestr(now,'yyyymmdd_HHMMSS')]);
rootdir = char(rootdir);   % <-- force to character vector for exist/mkdir
if ~exist(rootdir,'dir')
    mkdir(rootdir);
end

%% ================= LOAD KOOPMAN MODEL =================
if soil == "sandy"
    folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
else
    folder = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
end
files = dir(fullfile(folder,'*.mat'));
filename = fullfile(folder, files(model_idx).name);

model = load(filename,"MDL_fitr","A","B","Bc1","C","Cc1","K_obs","mean_std_out","rr","b","valData");

K   = model.A;
B   = model.B;
Bc1 = model.Bc1;
C   = model.C;
Cc1 = model.Cc1;

r = size(K,1);
m = size(B,2);


%% ================= BUILD OPTI ONCE =================
% load("tets.mat")
import casadi.*
opti = casadi.Opti();

Z     = opti.variable(r, Np+1);
Xcart = opti.variable(3, Np+1);
U     = opti.variable(m, Np);

Yref  = opti.parameter(6, Np);
Z0par = opti.parameter(r,1);
X0par = opti.parameter(3,1);

% --- weight parameters (diagonal entries) ---
wQ  = opti.parameter(6,1);   % [X Y psi u v r]
wRu = opti.parameter(2,1);   % [delta tau]
wRdu= opti.parameter(2,1);   % [d_delta d_tau]
wQf = opti.parameter(6,1);   % terminal weights

% initial constraints
opti.subject_to(Z(:,1)     == Z0par);
opti.subject_to(Xcart(:,1) == X0par);

% dynamics + cost
J_track = 0; J_u = 0; J_du = 0;

for k = 1:Np
    uk = U(:, min(k,Nc));

    % Koopman propagation
    opti.subject_to( Z(:,k+1) == K*Z(:,k) + B*uk + Bc1 );

    % Koopman output (body-frame)
    yK = C*Z(:,k) + Cc1;        % [u; v; r]
    uB = yK(1); vB = yK(2); rB = yK(3);

    % Kinematics
    psi_k = Xcart(3,k);
    opti.subject_to( Xcart(1,k+1) == Xcart(1,k) + (uB*cos(psi_k) - vB*sin(psi_k))*dt_mpc );
    opti.subject_to( Xcart(2,k+1) == Xcart(2,k) + (uB*sin(psi_k) + vB*cos(psi_k))*dt_mpc );
    opti.subject_to( Xcart(3,k+1) == Xcart(3,k) + rB*dt_mpc );

    % tracking error
    yFull = [Xcart(:,k); yK];
    e = yFull - Yref(:,k);

    % wrap heading error
    e(3) = atan2(sin(e(3)), cos(e(3)));

    % tracking cost (diagonal weights)
    J_track = J_track + sum(wQ .* (e.^2));

    % input effort cost (normalized by range scale)
    un = uk ./ u_scale;
    J_u = J_u + sum(wRu .* (un.^2));

    % input rate cost (normalized by range scale)
    if k > 1
        du  = uk - U(:,k-1);
        dun = du ./ du_scale;
        J_du = J_du + sum(wRdu .* (dun.^2));
    end
end

% terminal cost
yK_N = C*Z(:,Np+1) + Cc1;
yF_N = [Xcart(:,Np+1); yK_N];
eN   = yF_N - Yref(:,Np);
eN(3)= atan2(sin(eN(3)), cos(eN(3)));
J_term = sum(wQf .* (eN.^2));

J = J_track + J_u + J_du + J_term;
opti.minimize(J);

% input constraints
opti.subject_to(u_min <= U(:,1:Nc) <= u_max);

% solver
solver_opts.ipopt.print_level = 0;
solver_opts.print_time = 0;
opti.solver('ipopt', solver_opts);

%% ================= CANDIDATE WEIGHT SETS =================
% Q weights: [X Y psi u v r]
Q_list = [
    1   1   2   0.5  8   2
    2   2   4   0.5  10  2
    5   5   6   1    12  3
    10  10  8   1    15  4
    5   5   10  0.5  20  5
];

% Ru weights act on normalized inputs [delta/0.35, tau/135]
Ru_list = [
    1e-3  1e-3
    5e-3  1e-3
    1e-2  5e-3
];

% Rdu weights act on normalized rates [d_delta/0.7, d_tau/115]
Rdu_list = [
    1   1
    5   5
    10  10
];

term_mult = [5 10];   % terminal multiplier

%% ================= SUMMARY STORAGE =================
summary = [];
summary_hdr = ["case_id","rmseXY","rmseX","rmseY","rmsePsi","rmseU","rmseV","rmseR", ...
               "Q_X","Q_Y","Q_psi","Q_u","Q_v","Q_r","Ru_d","Ru_tau","Rdu_d","Rdu_tau","termMult"];

case_id = 0;

%% ================================================================
%                     GRID SEARCH LOOP
% ================================================================
for iq = 1:size(Q_list,1)
for iru = 1:size(Ru_list,1)
for irdu = 1:size(Rdu_list,1)
for itf = 1:numel(term_mult)

    case_id = case_id + 1;
    casedir = fullfile(rootdir, sprintf('case_%03d',case_id));
    if ~exist(casedir,'dir'); mkdir(casedir); end

    % set weights
    Qw  = Q_list(iq,:).';
    Ruw = Ru_list(iru,:).';
    Rduw= Rdu_list(irdu,:).';
    Qfw = term_mult(itf) * Qw;

    opti.set_value(wQ,  Qw);
    opti.set_value(wRu, Ruw);
    opti.set_value(wRdu,Rduw);
    opti.set_value(wQf, Qfw);

    % reset state
    x_real = Zhc0(:).';     % nonlinear full state
    x0 = x_real(1:6);

    % initial GP lift
    z0 = zeros(model.rr,1);
    for i = 1:model.rr
        X0n = (x0(model.K_obs) - model.mean_std_out(1,:))./model.mean_std_out(2,:);
        z0(i) = predict(model.MDL_fitr(i).gprMDL, X0n);
    end

    xsol = x_real;
    usol = [];
    cost_log = [];
    iter_log = [];
    status_log = strings(0);

    % video
    if MAKE_VIDEO
        fig_mpc = figure(200); clf(fig_mpc);
        set(fig_mpc,'Color','w','Position',[50 50 1600 900]);
        vid = VideoWriter(fullfile(casedir,'mpc_run.avi'), 'Motion JPEG AVI');
        vid.FrameRate = VIDEO_FPS;
        open(vid);
    end

    % ===== MPC loop =====
    for MPC_iter = 1:block:(TotalSteps-1)

        rows_idx = MPC_iter + (0:Np);
        rows_idx(rows_idx > TotalSteps) = TotalSteps;

        Yref_val = yref(rows_idx(1:end-1), :).';

        opti.set_value(Yref, Yref_val);
        opti.set_value(Z0par, z0);
        opti.set_value(X0par, x_real(1:3).');

        % warm start (use data inputs)
        Uinit = [ fv(rows_idx(1:end-1),12).'; fv(rows_idx(1:end-1),11).' ];
        opti.set_initial(U, Uinit);
        opti.set_initial(Xcart, yref(rows_idx,1:3).');

        sol = opti.solve();
        stats = opti.stats();

        Jval = sol.value(J);
        cost_log = [cost_log; Jval];
        iter_log = [iter_log; stats.iter_count];
        status_log(end+1,1) = string(stats.return_status);

        % apply first block controls
        u_star = sol.value(U(:,1:block));
        usol   = [usol, u_star];

        t_mpc_start = (MPC_iter-1)*dt_mpc;
        delta_mpc = @(t) zoh1d_scalar(t, t_mpc_start, dt_mpc, u_star(1,:), block+1);
        tau_mpc   = @(t) zoh1d_scalar(t, t_mpc_start, dt_mpc, u_star(2,:), block+1);

        t_mpc = t_mpc_start + (0:block)*dt_mpc;

        wrappedEventFcn = @(t, x) HalfCar_EventsFcn(t, x, b, delta_mpc, tau_mpc, traj);
        opts = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);

        [~,X_real] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,delta_mpc,tau_mpc,traj), t_mpc, x_real, opts);

        x_real = X_real(end,:);
        xsol   = [xsol; X_real(2:end,:)];

        % re-lift
        x0 = x_real(1:6);
        z0 = zeros(model.rr,1);
        for i = 1:model.rr
            X0n = (x0(model.K_obs) - model.mean_std_out(1,:))./model.mean_std_out(2,:);
            z0(i) = predict(model.MDL_fitr(i).gprMDL, X0n);
        end

        % plot & video frame
        if MAKE_VIDEO
            clf(fig_mpc); lw=3;
            subplot(2,2,1:2)
            plot(yref(:,1),yref(:,2),'LineWidth',lw); hold on;
            plot(xsol(:,1),xsol(:,2),'LineWidth',lw);
            plot(x_real(1),x_real(2),'ok','LineWidth',2);
            grid on; box on; xlabel('X'); ylabel('Y');
            title(sprintf('Case %03d | t = %.2f s', case_id, t_mpc_start));
            legend('Ref','Executed','Now','Location','best');

            subplot(2,2,3)
            tt=(0:size(usol,2)-1)*dt_mpc;
            plot(tt,usol(2,:),'LineWidth',lw); grid on; box on;
            xlabel('Time [s]'); ylabel('\tau');

            subplot(2,2,4)
            plot(tt,usol(1,:),'LineWidth',lw); grid on; box on;
            xlabel('Time [s]'); ylabel('\delta');

            drawnow;
            writeVideo(vid, getframe(fig_mpc));
        end
    end

    if MAKE_VIDEO; close(vid); end

    % ===== METRICS =====
    T = min(size(yref,1), size(xsol,1));
    err = xsol(1:T,1:6) - yref(1:T,1:6);
    err(:,3) = atan2(sin(err(:,3)), cos(err(:,3)));

    rmse = sqrt(mean(err.^2,1));              % [X Y psi u v r]
    rmseXY = sqrt(mean(sum(err(:,1:2).^2,2)));

    % ===== SAVE FINAL FIGURES =====
    lw=3;

    fig2 = figure(301); clf(fig2); set(fig2,'Color','w','Position',[50 50 1600 900]);
    subplot(2,2,1)
    plot(yref(:,1),yref(:,2),'LineWidth',lw); hold on;
    plot(xsol(1:T,1),xsol(1:T,2),'LineWidth',lw);
    grid on; box on; xlabel('X'); ylabel('Y'); legend('Ref','Executed','Location','best');
    title(sprintf('Case %03d | RMSE_{XY}=%.3f', case_id, rmseXY));

    for jj=1:3
        subplot(2,2,jj+1)
        plot(b.tspan(1:T), yref(1:T, model.K_obs(jj)), 'LineWidth', lw); hold on;
        plot(b.tspan(1:T), xsol(1:T, model.K_obs(jj)), 'LineWidth', lw);
        grid on; box on;
        xlabel('Time [s]');
        ylabel(model.valData.OutputName{model.K_obs(jj)},'Interpreter','latex');
        legend('Ref','Executed','Location','best');
    end
    saveas(fig2, fullfile(casedir,'final_states.png'));
    savefig(fig2, fullfile(casedir,'final_states.fig'));

    fig3 = figure(302); clf(fig3); set(fig3,'Color','w','Position',[50 50 1200 800]);
    tt=(0:size(usol,2)-1)*dt_mpc;
    subplot(2,1,1)
    plot(tt, usol(1,:), 'LineWidth', lw); grid on; box on;
    xlabel('Time [s]'); ylabel('\delta');
    subplot(2,1,2)
    plot(tt, usol(2,:), 'LineWidth', lw); grid on; box on;
    xlabel('Time [s]'); ylabel('\tau');
    saveas(fig3, fullfile(casedir,'final_inputs.png'));
    savefig(fig3, fullfile(casedir,'final_inputs.fig'));

    % ===== SAVE DATA =====
    save(fullfile(casedir,'mpc_results.mat'), ...
        'soil','traj','dt_mpc','Np','Nc','block', ...
        'Qw','Ruw','Rduw','Qfw','term_mult', ...
        'rmse','rmseXY','xsol','usol','cost_log','iter_log','status_log', ...
        'filename');

    % ===== ADD TO SUMMARY =====
    row = [case_id, rmseXY, rmse, ...
           Qw.', Ruw.', Rduw.', term_mult(itf)];
    summary = [summary; row]; %#ok<AGROW>

    fprintf('Case %03d done | RMSE_XY=%.4f | RMSE=[%.3f %.3f %.3f %.3f %.3f %.3f]\n', ...
        case_id, rmseXY, rmse(1),rmse(2),rmse(3),rmse(4),rmse(5),rmse(6));

end
end
end
end

%% ================= SAVE SUMMARY + PRINT BEST =================
save(fullfile(rootdir,'summary.mat'), 'summary', 'summary_hdr', 'Q_list','Ru_list','Rdu_list','term_mult','soil','traj');

% write CSV
csvpath = fullfile(rootdir,'summary.csv');
fid = fopen(csvpath,'w');
fprintf(fid, "%s", summary_hdr(1));
for i=2:numel(summary_hdr); fprintf(fid,",%s", summary_hdr(i)); end
fprintf(fid,"\n");
fclose(fid);
dlmwrite(csvpath, summary, '-append');

% best by rmseXY
[~,bestIdx] = min(summary(:,2));
bestRow = summary(bestIdx,:);
fprintf('\nBEST CASE = %03d | RMSE_XY=%.4f\n', bestRow(1), bestRow(2));
disp('Saved all results under:');
disp(rootdir);
