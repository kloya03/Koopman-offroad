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
filename ='../../datasets/sandyloam_100hz_no_elev_experiment_1579.mat';
load(filename,"b","trainData","valData","testData","numTest","numVal",...
    "numTrain","t_hc");
addpath("../../functions/utility")
% global ny nu nl sy su Ntr mean_std_inp mean_std_out idx_data K_obs
K_obs = 4:6;  % Only velocities as the observable

%% normalize data
for i=1:2
    inp = cell2mat(trainData(:,:,i).InputData);
    mean_std_inp(:,i) = [mean(inp(:));std(inp(:))];
end
for i=K_obs
    out = cell2mat(trainData(:,i).OutputData);
    mean_std_out(:,i+1-K_obs(1,1)) = [mean(out(:));std(out(:))];
end
clearvars out inp
%% Parameters
tic
% parameters selection 
nl = 600;                     % time delay--length of Hankel Matrix   *************
sy = 300;
cut_off = 7;
nB = 300;
save_filename = "Koopman_model_"+nl+"_"+sy+"_"+cut_off+"_"+nB;%d_%d_%d',nl,sy,cut_off,nB

ny = size(trainData(:,K_obs),2);     % number of outputs
nu = size(trainData,3);       % number of inputs
su = sy;
N4horizon = [nl,sy,su];
Ntr = randi(size(trainData,4));
n_stride = 5;
idx_data = 1:n_stride;
prev_GrassDist = [];
ct = [];

%% Initialize with n_stride trajectory data
[~,~,~,Xi_N1,SN1] = initialize_RSSID(trainData(:,K_obs,:,idx_data),...
    nl,sy,mean_std_inp,mean_std_out);

[Gam_Xi_R,rr] = find_ExObs(Xi_N1,cut_off);
%% Recursive SSID
toc
for iter =1+n_stride:n_stride:numTrain
    
    %%%%% Check Subspace distance for new data%%%%%
    traj = iter:min(iter+n_stride-1,numTrain);
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
        % [~,~,~,Xi_NN] = initialize_RSSID(trainData(:,K_obs,:,1:traj(:,end)),nl,sy,mean_std_inp,mean_std_out);
        
        % Find reduced order subspace
        [Gam_Xi_R,rr] = find_ExObs(Xi_N1,cut_off);
        check_sub = subspace(Gam_Xi_i,Gam_Xi_R);
        % check_sub1 = subspace(Gam_Xi_rold,Gam_Xi_R);
        % fprintf('Iteration %d | Counter: %f | Time: %d \n', iter, ct(end,:), toc);
        fprintf('RSSID:: Iteration %d-%d | sytem order: %d | Gr Dist: %.2f | check Dist: %.2f  \n', iter,iter+n_stride-1,rr,ct(end,end),check_sub);
    end
end
clc;
fprintf('RSSID:: Iteration %d-%d | sytem order: %d | Gr Dist: %.2f | check Dist: %.2f  \n', iter,iter+n_stride-1,rr,ct(end,end),check_sub);
et1 = toc;
%% Find Koopman Matrices and realizations of latent initial values

opts.maxiter = 10000;
[A,Cn,Bn,XGprn,ZGpr, ytest,del_cost,total_cost] = find_KoopmanMatrices(...
    trainData(:,K_obs,:,idx_data),Gam_Xi_R,nB,mean_std_inp,mean_std_out,opts);
et2 = toc;

% un-normalize from the mean and std of I/O data
B = Bn./mean_std_inp(2,:);
C = Cn.*mean_std_out(2,:).';
Bc1 = -B*mean_std_inp(1,:).';
Cc1 = mean_std_out(1,:).';
% XGpr = (XGprn - mean_std_out(1,:))./mean_std_out(2,:);
%% Fit GP basis functions
maxiter = 1000;
for i =1:rr
    opts = statset('fitrgp');
    opts.TolFun = 1e-08; opts.MaxIter = maxiter;
    MDL_fitr(i).gprMDL = fitrgp(XGprn,ZGpr(:,i),'verbose',0,...
        'FitMethod','exact','PredictMethod','exact',...
        'KernelFunction','ardsquaredexponential',...
        'optimizeHyperparameters','auto',...
        'HyperparameterOptimizationOptions',struct('UseParallel',true,...
        'ShowPlots',0));
end
et3 = toc;
%% Validate the model using the validation dataset

refresh = [25,50,75,100,125,150,175,200,225,250];
test_ntr = size(testData,4);

for jj = 1:size(refresh,2)
    time_error = zeros(size(testData(:,:,:,1).OutputData)).';
        total_rmse = zeros(size(testData(:,:,:,1).OutputData,2),1);
    for i=1:test_ntr
        [Y_pred,y_out,error,Y_pred_95]  = K_RSSID_prediction(testData(:,:,:,i),...
            MDL_fitr,A,B,Bc1,C,Cc1,K_obs,mean_std_out,refresh(jj));

        time_error = time_error + error.withTime;
        total_rmse = total_rmse + error.overallRMSE;  % make sure time stamp...
                                              %  for each traj is same otherwise ...
                                              % rmse = (n1*rmse1+ n2*rmse2)/(n1+n2)
    end
    error_with_time(:,:,jj) = time_error./test_ntr;
    overall_error(:,jj) = total_rmse./test_ntr;

end

save(save_filename,'-v7.3')