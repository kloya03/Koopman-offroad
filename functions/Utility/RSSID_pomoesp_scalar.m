%% 4.  RSSID pomoesp scalar update          %%%%%%
function [Xi_N2,S1] = RSSID_pomoesp_scalar(Y_N,U_N,Phi_N,Xi_N,S)

% Single-column (rank-1) recursive PO-MOESP update (Algorithm 4.1)
%
% Inputs (new Hankel-block column):
%   u_nu : r*nu x 1   (stacked inputs u_{k-ν+1:k})
%   y_nu : m*nu x 1   (stacked outputs y_{k-ν+1:k})
%   phi  : (r+m)*nu x 1 (regressor built from past input+output)
%
% State (fields of S) at step N:
%   P        : (UN*UN')^{-1}            [r*nu x r*nu]
%   Psi      : (PhiN * PiU * PhiN')^{-1}[(r+m)*nu x (r+m)*nu]
%   YU       : YN * UN'                 [m*nu x r*nu]
%   PhiU     : PhiN * UN'               [(r+m)*nu x r*nu]
%   Xi       : compressed I/O matrix    [m*nu x m*nu]
%   YPiPhiT  : YN * PiU * PhiN'         [m*nu x (r+m)*nu]  (recommended to cache)
%
% Output:
%   S (all fields updated to step N+1)
%
% Notes:
% - Uses the paper's direct rank-1 Xi update: Xi_{N+1} = Xi_N + α e e' − β f f'
%   with f = (Y Π_U^⊥ Φ^T) Ψ q + e, computed using *old* YPiPhiT and Psi.
% - If you prefer, you can recompute Xi via: Xi = (Y Π_U^⊥ Φ^T) Ψ (Y Π_U^⊥ Φ^T)'.
% - For numerical stability, avoid explicit inverses; this code only inverts scalars.
Nc = size(Y_N,2);
% Unpack
P    = S.P;
Psi  = S.Psi;
YU   = S.YU;
PhiU = S.PhiU;
YPiPhi = S.YPiPhi;
Xi_N2 = Xi_N;
for iter = 1:Nc
    u_nu = U_N(:,iter);
    y_nu = Y_N(:,iter);
    phi_nu = Phi_N(:,iter);

    % ===== (1) Gains & innovations (Alg. 4.1: (64)-(67)) =====
    Pu   = P * u_nu;
    denA = 1.0 + u_nu.' * Pu;
    q     = (PhiU * Pu) - phi_nu;       % q_{N+1}
    e     = y_nu - YU * Pu;             % e_{N+1}
    Psiq = Psi * q;
    denB = denA + q.' * Psiq;

    alpha = 1.0 / denA;                 % α_{N+1}
    beta  = 1.0 / denB;                 % β_{N+1}

    % ===== (2) Xi update (Remark 4.3: (74)-(75)) =====
    f = (YPiPhi * Psiq) + e;
    Xi_N2 = Xi_N2 + alpha * (e * e.') - beta * (f * f.');

    % ===== (3) Carry recursions (Alg. 4.1: (71)-(73)) =====
    P     = P - alpha * (Pu * Pu.');           % P_{N+1}
    YU    = YU + y_nu * u_nu.';                              % (Y U^T)_{N+1}
    PhiU  = PhiU + phi_nu  * u_nu.';                            % (Phi U^T)_{N+1}

    % ===== (4) Update YΠU⊥Φ^T and Ψ (Alg. 4.1: (68)-(69)) =====
    % (Y Π_U^⊥ Φ^T)_{N+1} = (Y Π_U^⊥ Φ^T)_N - α e q'
    YPiPhi = YPiPhi - alpha * (e * q.');         % rank-1 downdate
    % Ψ_{N+1} = Ψ_N - β Ψ q q' Ψ
    Psi = Psi - beta * (Psiq * Psiq.');   % symmetric rank-1 update

end
    
% ===== (5) Save back =====
S1.P    = P;
S1.YU   = YU;
S1.PhiU = PhiU;
S1.YPiPhi = YPiPhi;
S1.Psi  = Psi;

end
