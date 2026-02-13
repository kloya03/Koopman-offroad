addpath("../../functions/simulation_model")
foldername = "";
b.terrain = "_clay";
b.elev = 0;
run('loadCommonParams.m');
b.verbose = true;
b.hr = nan;
b.hf = nan;

%% Time step
b.dt = 0.01;
b.tstart = 0;
b.tstop = 20;
b.tspan = b.tstart:b.dt:b.tstop;
b.nt = size(b.tspan);
ngrid = 1600;
%% Input values for delta and Fu
load("../../datasets/Offroad_InputsSignals.mat")
rng(42)

for traj=10
    delta_inp = DELTA(:,traj);
    tau_inp = TAU(:,traj);
    N = numel(b.tspan);
    delta = @(t) foh1d_scalar(t, b.tstart, b.dt, delta_inp, N);
    tau   = @(t) foh1d_scalar(t, b.tstart, b.dt, tau_inp, N);

    try
        
        % clc;
        % traj
        % Z = [X, Y, psi, dx, dy, dpsi, z, dz, theta, dtheta, omega_f, omega_r]
        %      1, 2,   3,  4,  5,    6, 7,  8,     9,     10
        xd = 10 * rand;  % ensure scalar
        Rf = b.R;        % ensure scalar or assign explicitly
        Z0 = [0;            % X        [1]
            0;              % Y        [2]
            -pi + 2*pi*rand;% psi      [3]
            xd;             % dx       [4]
            0;              % dy       [5]
            0;              % dpsi     [6]
            0;              % z        [7]
            0;              % dz       [8]
            0.01*rand;      % theta    [9]
            0;              % dtheta  [10]
            (1 + 0.1*rand) * xd / Rf;% omega_f [11]
            (1 + 0.1*rand) * xd / Rf];% omega_r [12]

        L = length(Z0);

        wrappedEventFcn = @(t, x) HalfCar_EventsFcn(t, x, b,delta,tau,traj);
        opts = odeset('Reltol',1e-5,'AbsTol',1e-6,'Events', wrappedEventFcn);

        [t_hc,Z_hc,te,xe,ie] = ode113(@(t,Z) HalfCarBekker_F3(t,Z,b,delta,tau,traj),b.tspan,Z0,opts);

        [dZdt,bekk_h,fv,sig_tau] = cellfun(@(t,Z) HalfCarBekker_F3(t,Z.',b,delta,tau,traj),...
            num2cell(t_hc),num2cell(Z_hc,2),'uni',0);
        fv = cell2mat(fv);
        bekk_h = cell2mat(bekk_h);
        sig_tau1 =  cat(3,sig_tau{:});
        dZdt = cell2mat(dZdt.').';

        % fv = [Flf, Fcf, Nf, Frr_f, Flr, Fcr, Nr, Frr_r, f_adx, f_ady, b.tau_t, delta,...
        %                         vlf, vcf, omega_f, b.sf, vlr, vcr, omega_r, b.sr]
        % bekk_h = [-int_sig_sin_f, int_tau_cos_f, thf_f, thr_f, thm_f, hf_f,
        %                 -int_sig_sin_r, int_tau_cos_r, thf_r, thr_r, thm_r, hf_r]
        % sig_tau = [sig_f, tau_xf, sig_r, tau_xr]

        m=matfile(sprintf('%s/data_%d.mat',foldername,traj),'writable',true);
        m.dZdt = dZdt;
        m.Z_hc = Z_hc;
        m.t_hc = t_hc;
        m.Z0 = Z0;
        m.INN = INN;
        m.bekk_h = bekk_h;
        m.sig_tau = sig_tau1;
        m.suff = suff;
        m.fv = fv;
        m.b = b;
        m.events = [te,xe,ie];
        m.input = [delta_inp,tau_inp];
        fprintf('Finished traj %d\n', traj);

    catch ME
        warning('Error in traj %d: %s', traj, ME.message);
        % Optional: log error to file or display full stack trace
        % disp(getReport(ME));
    end
end