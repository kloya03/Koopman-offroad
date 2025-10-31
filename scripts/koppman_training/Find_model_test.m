%% Find Koopman Matrices and realizations of latent initial values
clc;
clear;
tic
load("train_test_2.mat")
[A,C,B,XGpr,ZGpr, ytest,del_cost,total_cost] = find_KoopmanMatrices(...
    trainData(:,K_obs,:,idx_data),Gam_Xi_R,nB,mean_std_inp,mean_std_out);
et2 = toc
save('train_test_3','-v7.3')