function [F,h,int, out2] = bekkerforce_front(b,W,beta,s)
if ~isnan(b.hf)
    h0 = b.hf;
else
    h0 = (3*W / (b.b*(3-b.n)*(b.kc/b.b + b.kphi)*sqrt(2*b.R)))^(2/(2*b.n+1));  % initial guess sinkage
end

    delh = 1e-6 * h0;
tol = 1e-7;     % [N]
temp = b.ret_xy;
b.ret_xy = 0;
h = newtraph(h0,delh,tol,b,W,beta,s);
b.ret_xy = temp;
[F, int, out2] =  bekker(h,b,W,beta,s);
end


function [out,out1,out2] = bekker(hf,b,W,beta,s)

th_f = acos(1-hf/b.R);                  % front contact angle [rad]
th_m = (b.a0 + b.a1*s)*th_f;          % max stress angle [rad]
%  th_r = (b.b0 + b.b1*s)*th_f;          
th_r = -acos(1-b.lam*hf/b.R);           % rear contact angle [rad]
Nth = 100;                             % theta discretization
th = linspace(th_r, th_f, Nth);
th_e = th_f - (th - th_r)*(th_f - th_m)/(th_m - th_r);  % equivalent rear angle [rad]

h = zeros(1,Nth);    % sinkage [m]
hG = b.R*(cos(th) - cos(th_f));
hL = b.R*(cos(th_e) - cos(th_f));
h(th >= th_m) = hG(th >= th_m);
h(th <  th_m) = hL(th <  th_m);

sig = (b.kc/b.b + b.kphi)*(h.^b.n);         %*(1/((b.b).^(b.n-1)));       % normal stress [Pa]
jx = b.R*((th_f - th) - (1-s)*(sin(th_f) - sin(th)));
tau_max = (b.co + sig*tan(b.phi));     
tau_x = tau_max.*(1-exp(-(jx)/b.kx));
int_Fz = b.R*b.b*(tau_x.*sin(th) + sig.*cos(th));
Fz = trapz(th,int_Fz);                  % reaction force [N]

int_sig = trapz(th,b.R*b.b*(-sig.*sin(th)));
int_tau = trapz(th,b.R*b.b*(tau_x.*cos(th)));
vec_force = [sig; tau_x];

if b.ret_xy
    jy = b.R*(1-s)*(th_f-th)*tan(beta);
    tau_y = tau_max.*(1-exp(-(jy)/b.ky));
    Fy = -trapz(th,b.R*b.b*tau_y);
    int_Fx = b.R*b.b*(-sig.*sin(th) + tau_x.*cos(th));
    Fx = trapz(th,int_Fx);
    out = [Fx, Fy, Fz];
    out1 = [int_sig, int_tau, th_f, th_r, th_m, hf];
    out2 = vec_force;
else
    out = Fz;
end
end


function h = newtraph(h0,delh,tol,b,W,beta,s)

h1 = h0;
h0 = 0.5;      % arbitrary large number
ct = 0;

while abs(h1 - h0) > tol && ct<100
    h0 = h1;
    Fh0 = bekker(h0,b,W,beta,s) - W;
    Fh0_p = bekker(h0 + delh,b,W,beta,s) - W;
    Fh0_m = bekker(h0 - delh,b,W,beta,s) - W;

    dFh0 = (Fh0_p - Fh0_m)/(2*delh);
    h1 = h0 - Fh0/dFh0;
    ct = ct + 1;
end
h=h1;
end






