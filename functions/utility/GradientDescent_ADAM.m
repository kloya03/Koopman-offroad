%% 6. Gradient Descent (Adam)
function [Br,Zlift0,del_cost,total_cost] = GradientDescent_ADAM(Phi_B,Phi_Z,BZ,Yhr,nB,opts)

ntr = size(BZ,2);
rr  = size(Phi_Z,2);
nu  = size(Phi_B,2)./rr;

% ----- Initial Guess -----
Zlift0 = reshape(BZ(nu*rr+1:end,:), rr, ntr);   % parameter 2 (lifted initial states)
vec_Br = sum(BZ(1:nu*rr,:),2)/ntr;             % parameter 1 (B matrix, vectorized)
x      = [vec_Br; Zlift0(:)];                  % stacked parameter vector

% ----- Adam Hyperparameters -----
if isfield(opts,'alpha0')
    alpha0 = opts.alpha0;      % base learning rate
else
    alpha0 = 0.01;             % good starting point if data are normalized
end

if isfield(opts,'beta1')
    beta1 = opts.beta1;
else
    beta1 = 0.9;               % default Adam β1
end

if isfield(opts,'beta2')
    beta2 = opts.beta2;
else
    beta2 = 0.999;             % default Adam β2
end

if isfield(opts,'eps_adam')
    eps_adam = opts.eps_adam;
else
    eps_adam = 1e-8;           % numerical stability
end

num_iter = opts.max_iter;
cost_tol = opts.del_cost_tol;

% Adam moment vectors (same size as x)
m = zeros(size(x));    % first moment (mean of gradients)
v = zeros(size(x));    % second moment (mean of squared gradients)

total_cost = nan(1,num_iter+1);
del_cost   = nan(1,num_iter);

for iter = 1:num_iter

    % ----- Cost and gradient wrt x -----
    [f,g] = BZ_cost_func(x,Phi_B,Phi_Z,Yhr,nB);
    total_cost(1,iter) = f;

    % ----- Adam moments update -----
    m = beta1*m + (1 - beta1)*g;
    v = beta2*v + (1 - beta2)*(g.^2);

    % Bias-corrected moments
    m_hat = m ./ (1 - beta1^iter);
    v_hat = v ./ (1 - beta2^iter);

    % Adam parameter update
    step = alpha0 * m_hat ./ (sqrt(v_hat) + eps_adam);
    x    = x - step;

    % ----- Cost change and early stopping -----
    if iter > 2
        del_cost(iter-1) = total_cost(1,iter) - total_cost(1,iter-1);

        if mod(iter,500) == 0
            fprintf(['Adam GD:: Iter %4d | Cost: %e | ΔJ: %e | ' ...
                'mean(|step|): %e | max(|step|): %e\n'], ...
                iter, total_cost(iter), del_cost(iter-1), ...
                mean(abs(step)), max(abs(step)));
        end
    else
        fprintf('Adam GD:: Iter %4d | Cost: %e\n',iter, total_cost(iter));

        % if abs(del_cost(iter-1)) < cost_tol
        %     fprintf('Adam GD:: Stopping early at iter %d (cost plateau: |ΔJ| = %.3e)\n', ...
        %         iter, abs(del_cost(iter-1)));
        %     break;
        % end
    end
end

% Unpack x back into vec_Br and Zlift0
vec_Br = x(1:nu*rr,:);
Zlift0 = reshape(x(nu*rr+1:end,:), rr, ntr);

% Final cost at converged x
total_cost(1,iter+1) = BZ_cost_func(x,Phi_B,Phi_Z,Yhr,nB);

% Reshape B back to (rr x nu)
Br = reshape(vec_Br, rr, nu);

end
