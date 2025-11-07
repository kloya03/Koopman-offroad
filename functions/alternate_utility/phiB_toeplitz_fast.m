function [Phi_B_N, BZ] = phiB_toeplitz_fast(Ex_obs, U_tr, Y_tr, nB, Phi_Z)

lambda_ridge = 1e-5;
[ny, rr] = size(Ex_obs);
[nu_nB, ntr] = size(U_tr); 
nu = nu_nB / nB;

Phi_B_N = zeros(ny*nB, nu*rr, ntr, 'like', U_tr);
want_BZ = nargout > 1;
if want_BZ
    BZ = zeros(nu*rr + rr, ntr); else, BZ = []; 
end

% Pre-extract E{l}
E = cell(1, nB); for l = 1:nB, E{l} = Ex_obs(:,:,l); end

% Ensure a pool exists (cheap if already open)
try gcp('nocreate'); catch, end

parfor i = 1:ntr
    UU = reshape(U_tr(:,i), nu, nB);   % [nu x nB]
    YY = Y_tr(:,i);                    % [ny*nB x 1]

    phi_B = zeros(ny*nB, nu*rr);
    for ii = 1:nB
        % accumulate lagged contributions; specialize nu=2 if you like
        S = zeros(ny, nu*rr);
        col = 1;
        for j = 1:nu
            Sj = zeros(ny, rr);
            for l = 1:ii
                Sj = Sj + UU(j, ii+1-l) * E{l};  % no kron
            end
            S(:, col:col+rr-1) = Sj; col = col + rr;
        end
        r0 = (ii-1)*ny + 1;
        phi_B(r0:r0+ny-1, :) = S;
    end

    % Normalize columns with guard
    sig = std(phi_B, 0, 1);  sig = sig + eps(class(sig));
    Phi_B_i = (phi_B) ./ sig;
    Phi_B_N(:,:,i) = Phi_B_i;

    if want_BZ
        phi = [Phi_B_i, Phi_Z];
        G = phi.'*phi;
        if lambda_ridge > 0
            G(1:size(G,1)+1:end) = G(1:size(G,1)+1:end) + lambda_ridge;
        end
        rhs = phi.'*YY;
        [R,p] = chol(G);
        if p > 0
            jitter = max(1e-12, 1e-9 * trace(G)/size(G,1));
            R = chol(G + jitter*eye(size(G), 'like', G));
        end
        BZ(:,i) = R \ (R.' \ rhs);
    end
end
end
