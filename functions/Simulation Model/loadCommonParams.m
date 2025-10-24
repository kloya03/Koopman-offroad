
%% Vehicle Parameters
% https://military.polaris.com/en-us/mrzr-d4-military-tan/specs/

%% Wheel parameters
b.R = 0.33;                          % Tire radius [m]
b.b = 0.2286;              % Wheel width [m]
b.mw = 30;                           % Unsprung mass
b.mwi = 5;
b.Iw = 0.5 * b.mwi * b.R^2;          % Inertia wheel
b.kf = 0.5*1e4;                      % Stiffness (front)
b.kr = 0.5*1e4;                      % Stiffness (rear)
b.cf = 3e2;                          % Damping  (front)
b.cr = 3e2;                          % Damping  (rear)

%% Vehicle body parameters
b.m = 952.54/2;                      % Sprung Mass
L = 3.556;
W = 1.511;
H = 1.875;
wb = 2.719;                          % wheelbase
b.Lf = wb/2;                         % length (front)
b.Lr = wb/2;                         % length (rear)
b.Iy = 1/12*b.m*(L^2 + H^2);
b.Iz = 1/12*b.m*(L^2 + W^2);
b.g = 9.81;                          % Gravity

%% Air drag coeffs
% f_adx = 0.5*b.rho_air*b.Afx*b.Cd*dx^2;
% f_ady = 0.5*b.rho_air*b.Afy*b.Cd*dy^2;
b.Afx = 2.5;                          % front Area
b.Afy = 4.5;
b.rho_air = 1.2;                     % air density Kg/m^3
b.Cd = 0.5;                         % Air drag coeffs

%% Rolling resistance parameters
% F_rr = b.P_rr^(b.alpha_rr)*F_z^(b.beta_rr)*(b.A_rr + b.B_rr*vlf + b.C_rr*vlf^2)
b.alpha_rr = -0.4;
b.beta_rr = 0.8;
b.P_rr = 0.5*(22.0627/33);
b.A_rr = 84e-4;
b.B_rr = 6.2e-4;
b.C_rr = 1.6e-4;

%% soft soil Parameters for Bekker Model
if b.terrain == 1
    % % Terrain parameters [Clay]
    b.kc = 13200;          % Cohesive modulus [N/m^{n+1}]
    b.kphi = 692200;       % Frictional modulus [N/m^{n+2}]c
    b.n = 0.5;
    b.kx = 0.01;          % Shear deformation modulus [m]
    b.ky = b.kx;           %                           ?????
    b.co = 4140;           % Cohesion [Pa]
    b.phi = 0.2269;       % Angle of internal friction [rad]
    b.a0 = 0.43;          % Coefficients for theta_m
    b.a1 = 0.32;          % ''  ''  ''  ''
    suff = "_clay"

elseif b.terrain == 2
    % % Terrain parameters [Sandy Loam]
    b.kc = 5300;          % Cohesive modulus [N/m^{n+1}]
    b.kphi = 1515000;     % Frictional modulus [N/m^{n+2}]
    b.n = 0.7;
    b.kx = 0.025;          % Shear deformation modulus [m]
    b.ky = b.kx;           %                           ?????
    b.co = 1700;           % Cohesion [Pa]
    b.phi = 0.5061;       % Angle of internal friction [rad]
    b.a0 = 0.18;          % Coefficients for theta_m
    b.a1 = 0.32;          % ''  ''  ''  ''
    suff = "_sandy_loam"
else
    % % Terrain parameters [Sand]
    b.kc = 1000;          % Cohesive modulus [N/m^{n+1}]
    b.kphi = 1528600;     % Frictional modulus [N/m^{n+2}]
    b.n = 1.08;         % Sinkage exponent [-]
    b.n = 1.08;
    b.kx = 0.024;          % Shear deformation modulus [m]
    b.ky = b.kx;           %                           ?????
    b.co = 1000;           % Cohesion [Pa]  200
    b.phi = 0.4712;       % Angle of internal friction [rad]
    b.a0 = 0.18;          % Coefficients for theta_m ????
    b.a1 = 0.32;          % ''  ''  ''  ''           ????
    suff = "_sand"

end
% b.b0 = -0.2;          % Coefficients for theta_r
% b.b1 = 0;             % ''  ''  ''  ''
b.lam = 0.08;
b.ret_xy = 1;          % logical for whether to return x,y forces

% b.beta = 0.3;         % side slip angle
% b.r = 0.8;            % Wheel radius [m]
% b.W = 2500;           % Normal load [N]
% b.W = (b.m1 + b.m2)*b.g;
% b.s = 0.2;            % Longitudinal slip [-]
% b.v = 5.5;            % Longitudinal speed [-]

%% engine torque
% b.omega_max = 1050;
% Define RPM values (in revolutions per minute)
rpm_vals = [50, 200, 300, 350, 400, 450, 550, 650, 700];
% Corresponding torque values (in Newton-meters)
tau_vals = [90, 105, 110, 120, 125, 130, 135, 130,  120];
rad_per_sec = 2*pi*rpm_vals/60;
speeds = rad_per_sec*0.33;
engine_speeds = linspace(0, rad_per_sec(end)+20, 500);
b.tau_curve = @(engine_speeds) max(interp1(rad_per_sec, tau_vals, engine_speeds , 'spline', 'extrap'),0);
estimated_torque = b.tau_curve(engine_speeds);
% plot(engine_speeds*0.33,estimated_torque); hold on;
% plot(speeds,tau_vals,'or')

%% Terrain elevation definition
% Parameters
if b.elev == 0
    b.H = 0
else
    b.H = 0.1         % Maximum height in (m)
end
b.Lx = 50;          % Domain size in x (m)
b.Ly = 50;          % Domain size in y (m)
b.N = 4;            % Number of Fourier modes per axis (low = smoother)
% Random coefficients for cosine and sine terms
rng(42);  % for reproducibility
A = randn(b.N,b.N);
B = randn(b.N,b.N);
C = randn(b.N,b.N);
D = randn(b.N,b.N);
b.wf = 2;
kx = (1:b.N)' * pi * b.wf / b.Lx;
ky = (1:b.N)  * pi * b.wf / b.Ly;

% Terrain height function: normalized to max |height| = H
raw_terrain = @(x, y) arrayfun(@(xi, yi) ...
    sum(sum( ...
        A .* cos(kx * xi) .* cos(ky * yi) + ...
        B .* sin(kx * xi) .* cos(ky * yi) + ...
        C .* cos(kx * xi) .* sin(ky * yi) + ...
        D .* sin(kx * xi) .* sin(ky * yi) ...
    )), x, y);

% Evaluate raw terrain on a reference grid to find its max
[xref, yref] = meshgrid(linspace(0, b.Lx, 100), linspace(0, b.Ly, 100));
zref = raw_terrain(xref, yref);
zmax = max(abs(zref(:)));   % Maximum absolute height
% Final height function scaled to max height H
b.h = @(x, y) b.H * raw_terrain(x, y) / zmax;

b.dh_x = @(x, y) b.H / zmax * arrayfun(@(xi, yi) ...
    sum(sum( ...
        -A .* (kx .* sin(kx * xi)) .* cos(ky * yi) + ...
         B .* (kx .* cos(kx * xi)) .* cos(ky * yi) - ...
         C .* (kx .* sin(kx * xi)) .* sin(ky * yi) + ...
         D .* (kx .* cos(kx * xi)) .* sin(ky * yi) ...
    )), x, y);

b.dh_y = @(x, y) b.H/ zmax * arrayfun(@(xi, yi) ...
    sum(sum( ...
        -A .* cos(kx * xi) .* (ky .* sin(ky * yi)) - ...
         B .* sin(kx * xi) .* (ky .* sin(ky * yi)) + ...
         C .* cos(kx * xi) .* (ky .* cos(ky * yi)) + ...
         D .* sin(kx * xi) .* (ky .* cos(ky * yi)) ...
    )), x, y);

b.ddh_xx = @(x, y) b.H/ zmax * arrayfun(@(xi, yi) ...
    sum(sum( ...
        -A .* (kx.^2 .* cos(kx * xi)) .* cos(ky * yi) - ...
         B .* (kx.^2 .* sin(kx * xi)) .* cos(ky * yi) - ...
         C .* (kx.^2 .* cos(kx * xi)) .* sin(ky * yi) - ...
         D .* (kx.^2 .* sin(kx * xi)) .* sin(ky * yi) ...
    )), x, y);

b.ddh_yy = @(x, y) b.H/ zmax * arrayfun(@(xi, yi) ...
    sum(sum( ...
        -A .* cos(kx * xi) .* (ky.^2 .* cos(ky * yi)) - ...
         B .* sin(kx * xi) .* (ky.^2 .* cos(ky * yi)) - ...
         C .* cos(kx * xi) .* (ky.^2 .* sin(ky * yi)) - ...
         D .* sin(kx * xi) .* (ky.^2 .* sin(ky * yi)) ...
    )), x, y);

b.ddh_xy = @(x, y) b.H/ zmax * arrayfun(@(xi, yi) ...
    sum(sum( ...
         A .* (kx .* sin(kx * xi)) .* (ky .* sin(ky * yi)) - ...
         B .* (kx .* cos(kx * xi)) .* (ky .* sin(ky * yi)) + ...
         C .* (kx .* sin(kx * xi)) .* (ky .* cos(ky * yi)) - ...
         D .* (kx .* cos(kx * xi)) .* (ky .* cos(ky * yi)) ...
    )), x, y);


% [xq, yq] = meshgrid(linspace(-1, b.Lx, 200), linspace(-1, b.Ly, 200));
% zq = arrayfun(b.h, xq, yq);
% surf(xq, yq, zq); shading interp; colormap jet;
% xlabel('x (m)'); ylabel('y (m)'); zlabel('Height (m)');
% title('Smooth Random Terrain via Low-Frequency Fourier Series');
% axis equal



% % H = 0;
% % w1 = 0.5;
% % w2 = 1.5;
% H = 0.15;
% w1 = 0.5;
% w2 = 0.2;
% b.h = @(x,y) H*sin(w1*x).^2.*cos(w2.*y);
% b.dh_x = @(x,y) 2*H*w1*cos(w1*x)*cos(w2*y)*sin(w1*x);
% b.dh_y = @(x,y) -H*w2*sin(w1*x)^2*sin(w2*y);
% b.ddh_xx = @(x,y) 2*H*w1^2*cos(2*w1*x)*cos(w2*y);
% b.ddh_yy = @(x,y) -H*w2^2*cos(w2*y)*sin(w1*x)^2;
% b.ddh_xy = @(x,y) -H*w1*w2*sin(2*w1*x)*sin(w2*y);
% %
% % z = terrain_elev(0.1);
% % b.h = z.h; %@(x,y) H*sin(w1*x).^2.*cos(w2.*y);
% % b.dh_x = z.hx; %@(x,y) 2*H*w1*cos(w1*x)*cos(w2*y)*sin(w1*x);
% % b.dh_y = z.hy; %@(x,y) -H*w2*sin(w1*x)^2*sin(w2*y);
% % b.ddh_xx = z.hxx; %@(x,y) 2*H*w1^2*cos(2*w1*x)*cos(w2*y);
% % b.ddh_yy = z.hyy; %@(x,y) -H*w2^2*cos(w2*y)*sin(w1*x)^2;
% % b.ddh_xy = z.hxy; %@(x,y) -H*w1*w2*sin(2*w1*x)*sin(w2*y);



