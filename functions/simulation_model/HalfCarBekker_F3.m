function [dZdt, out1, out2, out3] = HalfCarBekker_F3(t,Z,b,delta,tau,traj)
% clc;
% t
% Z = [X, Y, psi, dx, dy, dpsi, z, dz, theta, dtheta]
%      1, 2,   3,  4,  5,    6, 7,  8,     9,     10
X = Z(1);
Y = Z(2);
psi = Z(3);
dx = Z(4);
dy = Z(5);
dpsi = Z(6);
z = Z(7);
dz = Z(8);
th = Z(9);
dth = Z(10);
omega_f = Z(11);
omega_r = Z(12);

eps = 0.05;

% max torque until base speed, then linearly drop
tau_t = min(max(tau(t),0),b.tau_curve(omega_r));

delta_t = delta(t);

% Velocities in wheel frame
vcf = (dy + b.Lf*dpsi)*cos(delta_t) - dx*sin(delta_t);          % cornering velocity, front wheel
vlf = (dy + b.Lf*dpsi)*sin(delta_t) + dx*cos(delta_t);          % longitudinal velocity, front wheel
vcr = dy - b.Lr*dpsi;                                           % cornering velocity, rear wheel
vlr = dx;                                                       % longitudinal velocity, rear wheel

% Position / Velocity Kinematics
Xf = X + b.Lf*cos(psi);
Yf = Y + b.Lf*sin(psi);
Xr = X - b.Lr*cos(psi);
Yr = Y - b.Lr*sin(psi);

dX = dx*cos(psi) - dy*sin(psi);
dY = dx*sin(psi) + dy*cos(psi);
dXf = dx*cos(psi) - (dy + b.Lf*dpsi)*sin(psi);
dYf = dx*sin(psi) + (dy + b.Lf*dpsi)*cos(psi);
dXr = dx*cos(psi) - (dy - b.Lr*dpsi)*sin(psi);
dYr = dx*sin(psi) + (dy - b.Lr*dpsi)*cos(psi);

% Vertical dynamics
z_r = z - b.Lr*sin(th);
z_f = z + b.Lf*sin(th);
z_rg = b.h(Xr,Yr);
z_fg = b.h(Xf,Yf);

dz_r = dz - b.Lr*dth*cos(th);
dz_f = dz + b.Lf*dth*cos(th);
dz_rg = b.dh_x(Xr,Yr)*dXr + b.dh_y(Xr,Yr)*dYr;
dz_fg = b.dh_x(Xf,Yf)*dXf + b.dh_y(Xf,Yf)*dYf;

Nfi = b.m*b.g/2 - b.kf*(z_f - z_fg) - b.cf*(dz - dz_fg);
Nri = b.m*b.g/2 - b.kr*(z_r - z_rg) - b.cr*(dz - dz_rg);

maxErr = 1e4;
errTol = 1e-4;
ct = 0;
    
while maxErr > errTol && ct < 25
    % displayFormula('max_error = maxErr')
    % Front wheel force computation 
    Wf = Nfi;
    beta_f = atan2(round(vcf,2),round(vlf,2));    % using round for not getting large value when both are small
    if  abs(b.R*omega_f) < eps && abs(vlf) < eps
        sf = 0;
    elseif abs(omega_f*b.R) >= abs(vlf)
        sf = (b.R*omega_f-vlf)/(b.R*omega_f);
    else
         sf = (b.R*omega_f-vlf)/(vlf);
    end
    % b.sf = b.s;
    [Force,hf,int_f,vec_force_r] = bekkerforce_front(b,Wf,beta_f,sf);
    b.hf = hf;
    Flf = Force(1);        % longitudinal force, front wheel
    Fcf = Force(2);        % cornering force, front wheel
    Frr_f = b.P_rr^(b.alpha_rr)*Wf^(b.beta_rr)*(b.A_rr + b.B_rr*vlf + b.C_rr*vlf^2);
    
    % Rear wheel force computation
    Wr = Nri;
    beta_r = atan2(round(vcr,2),round(vlr,2));
    if abs(b.R*omega_f) < eps && abs(vlf) < eps
        sr =0;
    elseif abs(omega_r*b.R) >= abs(vlr)
        sr = (b.R*omega_r-vlr)/(b.R*omega_r);
    else
         sr = (b.R*omega_r-vlr)/(vlr);
    end
    % b.sr = b.s;
    [Force,hr,int_r,vec_force_f] = bekkerforce_rear(b,Wr,beta_r,sr);
    b.hr = hr;
    Flr = Force(1);        % longitudinal force, rear wheel
    Fcr = Force(2);        % cornering force, rear wheel 
    Frr_r = b.P_rr^(b.alpha_rr)*Wr^(b.beta_rr)*(b.A_rr + b.B_rr*vlr + b.C_rr*vlr^2);

    f_adx = 0.5*b.rho_air*b.Afx*b.Cd*dx^2;
    f_ady = 0.5*b.rho_air*b.Afy*b.Cd*dy^2;

    % Wheel angular accelerations
    domega_r = (-b.R*(Flr + Frr_r) + tau_t)./b.Iw;
    domega_f = (-b.R*(Flf + Frr_f))./b.Iw;

    % Accelerations
    ddx = (b.m*dy*dpsi*cos(th) + b.m*dz*dth +...
        (Flf*cos(delta_t) - Fcf*sin(delta_t) + Flr) - f_adx) / b.m;
    ddy = (-b.m*dx*dpsi*cos(th) + b.m*dz*dpsi*sin(th) + Flf*sin(delta_t) + Fcf*cos(delta_t) + Fcr- f_ady) / b.m;
    ddpsi = ((Flf*sin(delta_t) + Fcf*cos(delta_t))*b.Lf - Fcr*b.Lr)/b.Iz;

    % Wheel accelerations
    ddXf = (ddx - dy*dpsi - b.Lf*dpsi^2)*cos(psi) - (ddy + dx*dpsi + b.Lf*ddpsi)*sin(psi);
    ddYf = (ddx - dy*dpsi - b.Lf*dpsi^2)*sin(psi) + (ddy + dx*dpsi + b.Lf*ddpsi)*cos(psi);
    ddXr = (ddx - dy*dpsi + b.Lr*dpsi^2)*cos(psi) - (ddy + dx*dpsi - b.Lr*ddpsi)*sin(psi);
    ddYr = (ddx - dy*dpsi + b.Lr*dpsi^2)*sin(psi) + (ddy + dx*dpsi - b.Lr*ddpsi)*cos(psi);
    
    % Vertical displacement of wheels
    z_rg = b.h(Xr,Yr) - hr;
    z_fg = b.h(Xf,Yf) - hf;

    % Normal forces
    Nf = b.m*b.g/2 - b.kf*(z_f - z_fg) - b.cf*(dz - dz_fg) + b.mw*(...
         b.ddh_xx(Xf,Yf)*dXf^2 + b.dh_x(Xf,Yf)*ddXf + b.ddh_yy(Xf,Yf)*dYf^2 + b.dh_y(Xf,Yf)*ddYf + 2*dXf*dYf*b.ddh_xy(Xf,Yf));

    Nr = b.m*b.g/2 - b.kr*(z_r - z_rg) - b.cr*(dz - dz_rg) + b.mw*(...
         b.ddh_xx(Xr,Yr)*dXr^2 + b.dh_x(Xr,Yr)*ddXr + b.ddh_yy(Xr,Yr)*dYr^2 + b.dh_y(Xr,Yr)*ddYr + 2*dXr*dYr*b.ddh_xy(Xr,Yr));

    % Check normal force errors
    errF = abs(Nf - Nfi);
    errR = abs(Nr - Nri);
    maxErr = max(errF,errR);
    
    % if (Nf < 0) 
    %     % fprintf('Front wheel Lifted')
    % elseif (Nr < 0)
    % %     fprintf('Rear wheel Lifted')
    % else
    Nfi = Nf;
    Nri = Nr;
    % end
    
    ct = ct + 1;
    % fprintf('%0.0f\t %0.0f\t %0.0f\t %0.5f\t %0.5f\n',ct, Nf, Nr, hf, hr)
end
% fprintf('ct = %0.0f\n',ct)

ddz = (-b.m*dy*dpsi*sin(th) - b.m*dx*dth - b.kr*(z_r - z_rg) - b.cr*(dz_r - dz_rg) - b.kf*(z_f - z_fg) - b.cf*(dz_f - dz_fg)) / b.m;
ddth = (b.kr*(z_r - z_rg)*b.Lr*cos(th) + b.cr*(dz_r - dz_rg)*b.Lr*cos(th) - b.kf*(z_f - z_fg)*b.Lf*cos(th) - b.cf*(dz_f - dz_fg)*b.Lf*cos(th)) / b.Iy;

% Z = [X, Y, psi, dx, dy, dpsi, z, dz, theta, dtheta, d_omega_f, d_omega_r]
%      1, 2,   3,  4,  5,    6, 7,  8,     9,     10,   11,         12

dZdt = [dX;
        dY;
        dpsi;
        ddx;
        ddy;
        ddpsi;
        dz;
        ddz;
        dth;
        ddth;
        domega_f;
        domega_r];

forces = [Flf, Fcf, Nf, Frr_f, Flr, Fcr, Nr, Frr_r, f_adx, f_ady, tau_t, delta_t];
vel = [vlf, vcf, omega_f, sf, beta_f, vlr, vcr, omega_r, sr, beta_r];
out1 = [int_f,int_r];  
out2 = [forces, vel];
out3 = [vec_force_f.',vec_force_r.'];


% clc

if b.verbose == true  
    fprintf('Tr = %f==== Time: %4.10f\t ========================================\n',traj,t)
    fprintf('sf: %3.2f\t, Vf:%3.2f\t, omega_f: %3.2f\t, Vcf: %3.2f\t, beta_f %3.2f\t \n',...
        sf,vlf,omega_f,vcf,atan2(vcf,vlf));
    fprintf('sr: %3.2f\t, Vr:%3.2f\t, omega_r: %3.2f\t, Vcr: %3.2f\t, beta_r %3.2f\t \n',...
        sr,vlr,omega_r,vcr,atan2(vcr,vlr));
    fprintf('Flf: %4.2f\t, Fcf: %4.2f\t, Nf: %7.2f\t, Flr: %4.2f\t, Fcr: %4.2f\t, Nr: %5.2f\t \n',...
        Flf,Fcf,Nf,Flr,Fcr,Nr);
    fprintf('delta: %3.2f\t, tau_r: %3.2f\t, int_tau_r: %3.2f\t, hr: %3.2f\t, int_tau_f: %3.2f\t, hf: %3.2f\t\n',...
        delta_t,tau_t,int_r(2),int_r(end),int_f(2),int_f(end));

% fprintf('%4.3f\t ddx: %7.2f \n',t,ddx)
end
% fprintf('%0.3f\n',dXdt)
% fprintf('----------------\n')
end    


% disp([t,Nf,Nr,Flf,Fcf,Flr,Fcr])
% fprintf('%4.2f\t %4.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t\n',t,ddy, -b.m*dx*dpsi, Flf*sin(delta), Fcf*cos(delta),Fcr)
% fprintf('%4.2f\t %4.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f \t %6.2f\t\n',t,ddx, b.m*dy*dpsi, Flf*cos(delta), - Fcf*sin(delta),Flr,Fx)
% disp([t,ddx,Fe,Fx_y,Fu,Fd])
% fprintf('%4.2f\t %4.4f\t %6.4f\t %6.2f\t %6.2f\t %6.2f\t\n',t,ddx, b.m*dy*dpsi*cos(th) + b.m*dz*dth ,...
%     1.25*(Flf*cos(delta_t) - Fcf*sin(delta_t) + Flr),b.Fx_t,f_ad)