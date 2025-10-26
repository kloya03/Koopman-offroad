clc;
clear;
% ny : number of outputs
% nu : number of inputs
% nl : time delay--length of rows in Hankel Matrix
% sy, su: Past time delay--length O/I   < nl
% nr, nr_i : (nx*nl)  No. of rows in Hankel Matrix
% Nts : No. of time steps in a trajectory
% nc : No. of columns in Hankel Matrix
% nB : No of rows for B matrix computation
filename ='Datasets/sandyloam_100hz_no_elev_experiment_1.mat';
load(filename,"b","trainData","valData","testData","numTest","numVal",...
    "numTrain","t_hc");
global ny nu nl sy su Ntr mean_std_inp mean_std_out idx_data
%% normalize data
for i=1:2
    inp = cell2mat(trainData(:,:,i).InputData);
    mean_std_inp(:,i) = [mean(inp(:));std(inp(:))];
end
for i=1:3
    out = cell2mat(trainData(:,i+2).OutputData);
    mean_std_out(:,i) = [mean(out(:));std(out(:))];
end
clearvars out inp

%% Parameter selection
tic
ny = size(trainData(:,3:5),2);     % number of outputs
nu = size(trainData,3);       % number of inputs
nl = 600;                     % time delay--length of Hankel Matrix   *************
sy = 200;
su = sy;
nB = nl-sy;
N4horizon = [nl,sy,su];
Ntr = randi(size(trainData,4));
n_stride = 5;
idx_data = 1:n_stride;
prev_GrassDist = [];
ct = [];
cut_off = 0.95;
%% Initialize with 5 trajectory data
[~,~,~,Xi_N1,SN1] = initialize_RSSID(trainData(:,3:5,:,idx_data),...
    nl,sy,mean_std_inp,mean_std_out);

[Gam_Xi_R,rr] = find_ExObs(Xi_N1,cut_off);

%% Recursive SSID

for iter =1+n_stride:n_stride:numTrain

    %%%%% Check Subspace distance for new data%%%%%
    traj = iter:iter+n_stride;
    [Y_N,U_N,Phi_N,Xi_i] = initialize_RSSID(trainData(:,3:5,:,traj),...
        nl,sy,mean_std_inp,mean_std_out);
    [Gam_Xi_i,ri] = find_ExObs(Xi_i,cut_off);
    GrDR_N = subspace(Gam_Xi_i,Gam_Xi_R);
    prev_GrassDist = [prev_GrassDist; GrDR_N];

    %%%% Updating the Subspace %%%%%
    if GrDR_N > 0.01
        ct = [ct;iter GrDR_N];
        idx_data = [idx_data, traj];
        %% Recursive subspace Identification
        [Xi_N1,SN1] = RSSID_pomoesp_scalar(Y_N,U_N,Phi_N,Xi_N1,SN1);

        %% Find reduced order subspace
        [Gam_Xi_R,rr] = find_ExObs(Xi_N1,cut_off);

        %% Find Koopman Matrices and realizations of latent initial values
        [A,C,B,XGpr,ZGpr, ytest,del_cost,total_cost] = find_KoopmanMtrices(...
            trainData(:,3:5,:,idx_data),Gam_Xi_R,ny,rr,nu,nB,mean_std_inp,mean_std_out);

        %% Fit GP basis functions
        for i =1:rr
            opts = statset('fitrgp');
            opts.TolFun = 1e-08; opts.MaxIter = maxiter;
            MDL_fitr(i).gprMDl = fitrgp(XGpr,ZGpr(:,i),'verbose',0,...
                'FitMethod','exact','PredictMethod','exact',...
                'KernelFunction','ardsquaredexponential',...
                'optimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions',struct('UseParallel',true,...
                'ShowPlots',0));
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. create Hankel Matrix    %%%%%%
function [H1,H2,H3] = createHankelMatrix(data,nl,sy,mean_std)
% H1: past hankel matrix
% H2: future hankel matrix
% H3: Hankel matrix with no repetition
% data: {1, No. of trajectory}(numPoints, numFeatures)
% mean: [mean of fetures, 1;
%        std of same features];
% nl: row size of hankel
% sy: divide into past=sy and future=nl-sy hankel
% Construct Hankel matrix using vectorized indexing
H = [];
H3 = [];
if iscell(data)
    for i=1:size(data,2)
        D = data{i};
        D_norm = ( D - mean_std(1,:) )./mean_std(2,:);
        [numPoints, numFeatures] = size(D);
        nc = numPoints-nl+1;
        nr = numFeatures*nl;
        idx = (0:nl-1)'; % Create index offsets for time steps
        HH = reshape(D_norm(idx + (1:nc), :).', nr, nc);
        H = [H, HH];
        if nargout > 2
            H3 = [H3, HH(:,1:nl:end)];
        end
    end
else
    [numPoints, numFeatures] = size(data);
    nc = numPoints-nl+1;
    nr = numFeatures*nl;
    idx = (0:nl-1)'; % Create index offsets for time steps
    data_norm = ( data - mean_std(1,:) )./mean_std(2,:);
    H = reshape(data_norm(idx + (1:nc), :).', nr, nc);
end
H1 = H(1:numFeatures*sy,:);
H2 = H(numFeatures*sy+1:end,:);
end


%% 2. initialize RSSID      %%%%%%
function [Y_N,U_N,Phi_N,Xi_N,S] = initialize_RSSID(Train_Data,...
    nl,sy,mean_std_inp,mean_std_out)
% Using inverse (QR factorization) to calculate terms Pi_orth_U and Psi_N
% PiOrt_U = eye(size(U_N,2))-(U_N.')*(pinv(U_N*U_N.'))*U_N = I - Q1*Q1.';
% Psi_N = (Phi_N*Pi_Uc*(Phi_N.'))^(-1) = (R22*R22.')^(-1);
% P_N = (U_N*U_N.')^(-1) = (R11*R11.')^(-1)

% Build Hankel blocks on (already) normalized data
[Yp,Y_N] = createHankelMatrix(Train_Data.OutputData,...
    nl,sy,mean_std_out );
[Up,U_N] = createHankelMatrix(Train_Data.InputData,...
    nl,sy,mean_std_inp);
Phi_N = [Up;Yp];
% Sizes
lU   = size(U_N,1);          % r*nu
lPhi = size(Phi_N,1);        % (r+m)*nu
T1    = size(U_N,2);
lY = size(Y_N,1);
% ========= RQ via QR(S.'): correct stack & correct partitions =========
[Q,R]  = qr([U_N;Up;Yp;Y_N].');
R = R.';
Q1 = Q(:,1:lU);
% Q2 = Q(:, lU+(1:lPhi));
R32  = R(lU+lPhi+(1:lY), lU+(1:lPhi));            % lY x lPhi
Xi_N = R32 *R32.';

if nargout > 4
    R11 = R(1:lU,1:lU);
    R22 = R(lU+(1:lPhi),lU+(1:lPhi));
    S.PiOrt_U = eye(T1) - Q1*Q1.';
    S.YPiPhi = Y_N*S.PiOrt_U*Phi_N.';
    S.YU = Y_N*U_N.';
    S.PhiU = Phi_N*U_N.';

    Iu   = eye(size(R11,1));
    Iphi = eye(size(R22,1));
    if istril(R11)                      % LOWER triangular
        S.P  = R11' \ (R11 \ Iu);       % == (R11*R11.').^(-1) stably
    else                                % UPPER triangular
        S.P  = R11  \ (R11' \ Iu);
    end

    if istril(R22)                      % LOWER triangular
        S.Psi = R22' \ (R22 \ Iphi);    % == (R22*R22.').^(-1) stably
    else
        S.Psi = R22  \ (R22' \ Iphi);
    end
end

% Gphi = R22*R22.';                    % == Φ Π_U^⊥ Φ^T (by construction)
% Gphi = (Gphi+Gphi.')/2;
% % Cholesky (or tiny ridge if needed)
% [Rc,p] = chol(Gphi);
% if p==0
%     Psi_N = Rc \ (Rc' \ eye(lPhi));  % == (Φ Π_U^⊥ Φ^T)^{-1}
% else
%     lam = 1e-8 * trace(Gphi)/lPhi;
%     Rc  = chol(Gphi + lam*eye(lPhi));
%     Psi_N = Rc \ (Rc' \ eye(lPhi));
% end
% Xi_N = Y_N*PiOrt_U*(Phi_N.')*Psi_N*Phi_N*PiOrt_U;
% PiOrth_U = eye(T1)-(U_N.')*(pinv(U_N*U_N.'))*U_N;
% P_N = (R11*R11.')^(-1);
% Psi_N = R22\((R22.')\eye(lPhi));  % (R22*R22.')^(-1);

%%%%% error analysis %%%%%
% rcond(Gphi)
% rcond(P_N)
% rcond(Psi_N)
% % rcond(PiOrth_U)
% % rcond(PiOrt_U)
% rcond(U_N*U_N.')
% rank(U_N)==size(U_N,1)
% % Verify identities (after the fix)
% PiU_perp_direct = eye(T1) - Q1*Q1.';
% E1 = norm(U_N*U_N.' - R11*R11.','fro')/norm(U_N*U_N.','fro');
% E2 = norm(Phi_N*PiU_perp_direct*Phi_N.' - R22*R22.','fro')/norm(R22*R22.','fro');
% E3 = norm(Xi_N - R32Q2T,'fro') / max(1,norm(R32Q2T,'fro'));
% fprintf('rel.err UU^T: %.2e,   rel.err ΦΠΦ^T: %.2e\n,  rel.err Xi Xi^T: %.2e\n', E1, E2, E3);
end

%% 3.  find Extended subspace      %%%%%%

% Observability Matrix and system rank
function [Gam_Xi_R,rr] = find_ExObs(Xi_N,tau)
% tau = 0.95;
[U,S,~] = svd(Xi_N,'econ');     % S is diag of singular values
sig = diag(S);                  % for PSD, sig == eigenvalues
Cu  = cumsum(sig) / sum(sig);   % NO sqrt here

rr   = find(Cu >= tau, 1, 'first');
Gam_Xi_R = U(:,1:rr)*sqrt(S(1:rr,1:rr));

end

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
clearvars S
Xi_N2 = Xi_N;
for iter = 1:Nc
    u_nu = U_N(:,iter);
    y_nu = Y_N(:,iter);
    phi_nu = Phi_N(:,iter);
    % ===== (1) Gains & innovations (Alg. 4.1: (64)-(67)) =====
    alpha = 1.0 / (1.0 + u_nu.' * P * u_nu);                 % α_{N+1}
    e     = y_nu - YU * (P * u_nu);                         % e_{N+1}
    q     = (PhiU * (P * u_nu)) - phi_nu;                      % q_{N+1}
    beta  = 1.0 / ( (1.0/alpha) + q.' * Psi * q );           % β_{N+1}

    % ===== (2) Xi update (Remark 4.3: (74)-(75)) =====
    f = (YPiPhi * (Psi * q)) + e;
    Xi_N2 = Xi_N2 + alpha * (e * e.') - beta * (f * f.');

    % ===== (3) Carry recursions (Alg. 4.1: (71)-(73)) =====
    P     = P - alpha * (P * u_nu) * (u_nu.' * P);           % P_{N+1}
    YU    = YU + y_nu * u_nu.';                              % (Y U^T)_{N+1}
    PhiU  = PhiU + phi  * u_nu.';                            % (Phi U^T)_{N+1}

    % ===== (4) Update YΠU⊥Φ^T and Ψ (Alg. 4.1: (68)-(69)) =====
    % (Y Π_U^⊥ Φ^T)_{N+1} = (Y Π_U^⊥ Φ^T)_N - α e q'
    YPiPhi = YPiPhi - alpha * (e * q.');         % rank-1 downdate
    % Ψ_{N+1} = Ψ_N - β Ψ q q' Ψ
    Psi = Psi - beta * (Psi * q) * (q.' * Psi);   % symmetric rank-1 update

end

% ===== (5) Save back =====
S1.P    = P;
S1.YU   = YU;
S1.PhiU = PhiU;
S1.YPiPhi = YPiPhi;
S1.Psi  = Psi;

end


%% 5. Find Koopman matrices and latent state realization

function [A,C,B,XGpr,ZGpr,ytest,del_cost,total_cost] = find_KoopmanMtrices(Train_Data,...
    Gamma_Xi_R,ny,rr,nu,nB,mean_std_inp,mean_std_out)

% A,C computation for Recursive SSID
C = Gamma_Xi_R(1:ny,:);
A = pinv(Gamma_Xi_R(1:end-ny,:))*Gamma_Xi_R(ny+1:end,:);

[~,~,U_tr] = createHankelMatrix(Train_Data.InputData,nB,0,mean_std_inp);
[~,~,Y_tr] = createHankelMatrix(Train_Data.OutputData,nB,0,mean_std_out);

[Phi_Z_N, Ex_obs] = phiZ(C,A,nB,ny,rr); 
[Phi_B_N, BZ]  = phiB(Ex_obs,U_tr,Y_tr,nB,ny,nu,rr,Phi_Z_N);

[B,Zlift,del_cost,total_cost] = GradientDescent(Phi_B_N,Phi_Z_N,BZ,Y_tr,nu,rr,nB);
clearvars Phi_B_N phi_Z_N
[~,~,~,min_D] = kmeans(Y_tr(1:nx,:).', 800);
[~,idx_GP] = min(min_D);
XGp = Y_tr(:,idx_GP);
ZGp = Zlift(:,idx_GP); 
ytest(:,:)= XGp - C*ZGp;
XGpr = XGp.';
ZGpr = ZGp.';

end

%% 6. Gradient Descent
function [Br,Zlift0,del_cost,total_cost] = GradientDescent(Phi_B,Phi_Z,BZ,Yhr,nu,rr,nB)

ntr = size(BZ,2);
% Initial Guess
Zlift0 = reshape(BZ(nu*rr+1:end,:),rr,ntr);  % Zlift0 = zeros(rr,ntr);  % parameter 2
vec_Br = sum(BZ(1:nu*rr,:),2)/ntr;   % Br = zeros(rr,1);  % parameter 1 
alpha_Z = 0.99;          % Learning rate
alpha_B = 0.01;
num_iter = 5000;        % Number of iterations
% Compute the predicted output , - phi_B*B_rep - Phi_Z_norm*Z0 prediction is (nx*nB,ntr) 
prediction = squeeze(pagemtimes(Phi_B,repmat(vec_Br,1,1,ntr)))+Phi_Z*Zlift0;
residual = Yhr - prediction;
total_cost = [(1/2*nB*ntr)*sum(diag((residual.')*residual)), zeros(1,num_iter)];
del_cost = total_cost(1,1);
% fprintf('Iteration %d | Cost: %f\n', 0, total_cost(1));

for iter=1:num_iter
    % Compute the gradient with respect to the j-th feature
    gradient_Zj = (-1/(nB*ntr))*(Phi_Z.')*residual;
    gradient_Bs = (-1/(nB*ntr))*pagemtimes(pagetranspose(Phi_B),reshape(residual,ny*nB,1,ntr)); 
    gradient_B = sum(gradient_Bs,3);

    % Update the Zlift0 j-th  and B parameter
    Zlift0 = Zlift0 - alpha_Z * gradient_Zj;
    vec_Br = vec_Br - alpha_B * gradient_B;

    % Compute the predicted output
    % - phi_B_norm*B_rep - Phi_Z_norm*Z0; prediction is (nx*nB,ntr) 
    prediction = squeeze(pagemtimes(Phi_B,repmat(vec_Br,1,1,ntr)))+Phi_Z*Zlift0; 
    residual = Yhr - prediction; % size is (nx*nB,ntr)

    % Compute the total cost
    total_cost(iter+1) = sum(diag((1/(2*ntr*nB))*(residual.')*residual));
    del_cost(iter) = total_cost(1,iter+1) - total_cost(1,iter);
    % fprintf('Iteration %d | Cost: %f\n', iter, total_cost(iter+1));
end
Br = reshape(vec_Br,rr,nu);

end

%% 7. PhiB

function [Phi_B_N,BZ] = phiB(Ex_obs,U_tr,Y_tr,nr,ny,nu,rr,Phi_Z)        %%%%%%

ntr = size(U_tr,2);
for i=1:ntr
    UU = U_tr(:,i);
    YY = Y_tr(:,i);
    phi_B = zeros(ny*nr,nu*rr);
    for ii = 1:nr
        PHIB = zeros(ny,nu*rr);
        for kk=1:ii
            phib = kron(UU((kk-1)*nu+1:nu*kk).',Ex_obs(:,:,ii+1-kk));
            PHIB = PHIB + phib;
        end
        phi_B(ny*(ii-1)+1:ny*ii,:) =PHIB;
    end

    % Phi_B(:,:,(iter-1)*ntr+i) = phi_B;
    Phi_B_N(:,:,i) = (phi_B-...
        mean(phi_B, 1))./std(phi_B,0,1);    % normalize
    if nargout > 1
        phi = [phi_B Phi_Z];
        BZ(:,i) = pinv((phi.')*phi)*(phi.')*YY(:,i);
    end
end

end

%% 8. PhiZ
function [phi_Z_N,Ex_obs] = phiZ(Cr,Kr,nr,nx,rr)      %%%%%%
Ex_obs = zeros(nx,rr,nr+1);
phi_Z = zeros(nx*nr,rr);
for ii = 1:nr
    Ex_obs(:,:,ii+1) =  Cr*(Kr^(ii-1));
    phi_Z(nx*(ii-1)+1:nx*ii,:) = Cr*(Kr^(ii-1));
end
phi_Z_N = (Phi_Z - mean(Phi_Z, 1))./std(Phi_Z, 0, 1); % normalize
end


%% 9. Validation Function