%% Find Koopman Matrices and realizations of latent initial values
clc;
clear;
addpath("../../functions/utility")
tic
load("train_test_2.mat")
opts.max_iter = 10000;
[A,C,B,XGpr,ZGpr, ytest,del_cost,total_cost] = find_KoopmanMatrices(...
    trainData(:,K_obs,:,idx_data),Gam_Xi_R,nB,mean_std_inp,mean_std_out,opts);
et2 = toc

%% Fit GP basis functions
maxiter = 1000;
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
save('train_test_3','-v7.3')