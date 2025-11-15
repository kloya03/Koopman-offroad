%% 5. Find Koopman matrices and latent state realization

function [An,Cn,Bn,XGprn,ZGpr,ytest,del_cost,total_cost] = find_KoopmanMatrices(Train_Data,...
    Gamma_Xi_R,nB,mean_std_inp,mean_std_out,opts)

ny = size(mean_std_out,2);
% A,C computation for Recursive SSID
Cn = Gamma_Xi_R(1:ny,:);
An = pinv(Gamma_Xi_R(1:end-ny,:))*Gamma_Xi_R(ny+1:end,:);

[~,~,U_tr] = createHankelMatrix(Train_Data.InputData,nB,0,mean_std_inp);   %% time consuming
[~,~,Y_tr] = createHankelMatrix(Train_Data.OutputData,nB,0,mean_std_out);  %% time consuming

[Phi_Z_N,D_PhiZ,Ex_obs] = phiZ(Cn,An,nB); 
[Phi_B_N,D_PhiB,BZ]  = phiB_parallel(Ex_obs,U_tr,Y_tr,nB,Phi_Z_N);   %% time consuming

[BN,ZliftN,del_cost,total_cost] = GradientDescent_ADAM(Phi_B_N,Phi_Z_N,BZ,Y_tr,nB,opts);
% [BN, ZliftN, total_cost, del_cost] = FitB_Z0_fminunc(Phi_B_N, Phi_Z_N, BZ, Y_tr, nB);
clearvars Phi_B_N phi_Z_N
rr = size(BN, 1); 
nu = size(BN, 2);

% un-normalize from the std of feature marix phi_B and phi_Z
Bn = reshape(diag(D_PhiB)\BN(:),rr,nu);   %% un-normalize the B matrix 
Zlift = diag(D_PhiZ)\ZliftN;             %% un-normalize the Z matrix 

Yini = Y_tr(1:ny,:); %*diag(mean_std_out(2,:))).' + mean_std_out(1,:).';
[~,~,~,min_D] = kmeans(Yini.', min(800,size(Yini,2)) );
[~,idx_GP] = min(min_D);
XGp = Yini(:,idx_GP);
ZGp = Zlift(:,idx_GP); 
ytest(:,:)= XGp - (C)*ZGp;
XGprn = XGp.';
ZGpr = ZGp.';
end