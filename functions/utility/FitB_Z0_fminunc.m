function [Br, Zlift0, fval, out] = FitB_Z0_fminunc(Phi_B, Phi_Z, BZ_init, Yhr, nB, opts)
if nargin < 6, opts = struct; end
if ~isfield(opts,'MaxIter'),   opts.MaxIter   = 5000;    end
if ~isfield(opts,'Display'),   opts.Display   = 'iter'; end
if ~isfield(opts,'lambda_B'),  opts.lambda_B  = 0;      end
if ~isfield(opts,'lambda_Z'),  opts.lambda_Z  = 0;      end

ntr = size(Phi_B,3);
rr = size(Phi_Z,2);
nu = size(Phi_B,2)./rr;
ny = size(Phi_B,1)./nB;

Zlift0   = reshape(BZ_init(nu*rr+1:end,:), rr, ntr);
vecB0 = sum(BZ_init(1:nu*rr,:), 2) / ntr;
x0     = [vecB0; Zlift0(:)];
lambda_B = opts.lambda_B;  lambda_Z = opts.lambda_Z;

options = optimoptions('fminunc','SpecifyObjectiveGradient',true, ...
    'Algorithm','quasi-newton','MaxIterations',opts.MaxIter, ...
    'Display',opts.Display,'PlotFcn',{'optimplotfval'});

fun = @(x) BZ_cost_func(x, Phi_B, Phi_Z, Yhr, nB);
[x_opt, fval, exitflag, output] = fminunc(fun, x0, options);

vecB   = x_opt(1:nu*rr);
Zlift0 = reshape(x_opt(nu*rr+1:end), rr, ntr);
Br     = reshape(vecB, rr, nu);

out.fval = fval; out.exitflag = exitflag; out.output = output;
out.rr = rr; out.nu = nu; out.ntr = ntr; out.ny = ny;
end
