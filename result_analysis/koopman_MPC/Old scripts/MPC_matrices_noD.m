function [H,f,Phi,Gamma] = MPC_matrices_noD(Z0,A,B,C,yref,ri,rNp,qi,qNp,Np)
% Relevant lengths
r = size(A,1);        % number of states
nu = size(B,2);        % number of inputs
ny = size(C,1);        % number of outputs

% Full System Matrices
Phi = zeros(r*(Np),r);
Gamma = zeros(r*(Np),nu*(Np));
Pbar = zeros(r*(Np),nu*(Np));
Qbar = zeros(r*(Np),r*(Np));
Rbar = zeros(nu*(Np),nu*(Np));

Gamrow = []; Gbar = []; Lam_Z = [];
for i = 1:Np
    Phi((i-1)*r+1:r*(i),:) = A^(i);
    Gamrow = [(A^(i-1))*B, Gamrow];
    Gamma(r*(i-1)+1:r*(i),1:size(Gamrow,2)) = Gamrow;
    if i < Np
        Lam_Z = [Lam_Z, -2*(yref(:,i).')*qi*C];
    end
end
Lam_Z=[Lam_Z,-2*(yref(:,end).')*qNp*C];

Qbar(1:r*(Np-1),1:r*(Np-1)) = kron(eye(Np-1),(C.')*qi*C);
Qbar(r*(Np-1)+(1:r),r*(Np-1)+(1:r)) = (C.')*qNp*C;

Rbar(1:nu*(Np-1),1:nu*(Np-1)) = eye(nu)*ri;
Rbar(nu*(Np-1)+(1:nu),nu*(Np-1)+(1:nu)) =eye(nu)*rNp;

H = 2*(Rbar + Gamma.'*Qbar*Gamma);
f = 2*Gamma.'*Qbar.'*Phi*Z0 + Gamma.'*Lam_Z';

end
