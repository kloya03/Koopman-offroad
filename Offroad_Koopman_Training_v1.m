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
filename ='Datasets/sandyloam_100hz_no_elev_experiment.mat';
load(filename,"b","trainData","valData","testData","numTest","numVal",...
    "numTrain","t_hc");
addpath("functions/Utility")
global ny nu nl sy su Ntr mean_std_inp mean_std_out idx_data K_obs
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
tic

%% Parameter selection
tic
ny = size(trainData(:,K_obs),2);     % number of outputs
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
cut_off = 5;
K_obs = 4:6;
%% Initialize with 5 trajectory data
[~,~,~,Xi_N1,SN1] = initialize_RSSID(trainData(:,K_obs,:,idx_data),...
    nl,sy,mean_std_inp,mean_std_out);

[Gam_Xi_R,rr] = find_ExObs(Xi_N1,cut_off);

%% Recursive SSID
toc
for iter =1+n_stride:n_stride:numTrain

    %%%%% Check Subspace distance for new data%%%%%
    traj = iter:iter+n_stride;
    [Y_N,U_N,Phi_N,Xi_i] = initialize_RSSID(trainData(:,K_obs,:,traj),...
        nl,sy,mean_std_inp,mean_std_out);
    [Gam_Xi_i,ri] = find_ExObs(Xi_i,cut_off);
    GrDR_N = subspace(Gam_Xi_i,Gam_Xi_R);
    Gam_Xi_rold = Gam_Xi_R;
    prev_GrassDist = [prev_GrassDist; GrDR_N];

    %%%% Updating the Subspace %%%%%
    if GrDR_N > 0.01
        ct = [ct;iter GrDR_N];
        idx_data = [idx_data, traj];

        % Recursive subspace Identification
        [Xi_N1,SN1] = RSSID_pomoesp_scalar(Y_N,U_N,Phi_N,Xi_N1,SN1);

        % Find reduced order subspace
        [Gam_Xi_R,rr] = find_ExObs(Xi_N1,cut_off);
        check_sub = subspace(Gam_Xi_i,Gam_Xi_R);
        % check_sub1 = subspace(Gam_Xi_rold,Gam_Xi_R);
        % fprintf('Iteration %d | Counter: %f | Time: %d \n', iter, ct(end,:), toc);
        fprintf('Iteration %d-%d | sytem order: %d | Gr Dist: %.2f | check Dist: %.2f  \n', iter,iter+n_stride,rr,ct(end,end),check_sub);
    end
end


%% Find Koopman Matrices and realizations of latent initial values
[A,C,B,XGpr,ZGpr, ytest,del_cost,total_cost] = find_KoopmanMtrices(...
    trainData(:,K_obs,:,idx_data),Gam_Xi_R,nB,mean_std_inp,mean_std_out);

%% Fit GP basis functions
for i =1:rr
    opts = statset('fitrgp');
    opts.TolFun = 1e-08; opts.MaxIter = maxiter;
    MDL_fitr(i).gprMDL = fitrgp(XGpr,ZGpr(:,i),'verbose',0,...
        'FitMethod','exact','PredictMethod','exact',...
        'KernelFunction','ardsquaredexponential',...
        'optimizeHyperparameters','auto',...
        'HyperparameterOptimizationOptions',struct('UseParallel',true,...
        'ShowPlots',0));
end

%% Validate the model using the validation dataset
errors = Koopman_validation(valData,MDL_fitr.gprMDL,K,B,C,K_obs);