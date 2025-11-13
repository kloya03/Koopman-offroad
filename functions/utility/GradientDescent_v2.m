%% 6. Gradient Descent
function [Br,Zlift0,del_cost,total_cost] = GradientDescent_v2(Phi_B,Phi_Z,BZ,Yhr,nB,opts)

ntr = size(BZ,2);
rr = size(Phi_Z,2);
nu = size(Phi_B,2)./rr;

% Initial Guess
Zlift0 = reshape(BZ(nu*rr+1:end,:),rr,ntr);  % Zlift0 = zeros(rr,ntr);  % parameter 2
vec_Br = sum(BZ(1:nu*rr,:),2)/ntr;   % Br = zeros(rr,1);  % parameter 1 
x = [vec_Br;Zlift0(:)];  % Combine parameters for cost function
alpha_Z = 0.99;          % Learning rate
alpha_B = 0.01;
num_iter = opts.max_iter;        % Number of iterations
cost_tol  = opts.del_cost_tol;        % stop if cost change is tiny

for iter=1:num_iter
    % Compute the gradient with respect to the j-th feature and total cost
    [f,g] = BZ_cost_func(x,Phi_B,Phi_Z,Yhr,nB);
    gradient_B = g(1:nu*rr,:);
    gradient_Zj = reshape(g(nu*rr+1:end,:),rr,ntr);
    total_cost(1,iter) = f;

    % Update the Zlift0 j-th  and B parameter
    Zlift0 = Zlift0 - alpha_Z * gradient_Zj;
    vec_Br = vec_Br - alpha_B * gradient_B;
    x = [vec_Br;Zlift0(:)];  % Combine parameters for cost function

    % Compute the del_cost
    if iter > 2
        del_cost(iter-1) = total_cost(1,iter) - total_cost(1,iter-1);

        if mod(iter,200)==0
            fprintf('GD:: Iteration %d | Cost: %f| Cost delta: %f\n', iter, total_cost(iter), del_cost(end));
        end
        
        if abs(del_cost(iter-1)) < cost_tol
            fprintf('GD:: Stopping early at iter %d (cost plateau: |ΔJ| = %.3e)\n', ...
                iter, abs(del_cost(iter-1)));
            break;
        end
    end
    
    
end
total_cost(1,iter+1) = BZ_cost_func(x,Phi_B,Phi_Z,Yhr,nB);
Br = reshape(vec_Br,rr,nu);

end
