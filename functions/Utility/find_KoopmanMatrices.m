%% 5. Find Koopman matrices and latent state realization

function [A,C,B,XGpr,ZGpr,ytest,del_cost,total_cost] = find_KoopmanMtrices(Train_Data,...
    Gamma_Xi_R,nB,mean_std_inp,mean_std_out)

ny = size(mean_std_out,2);
% A,C computation for Recursive SSID
C = Gamma_Xi_R(1:ny,:);
A = pinv(Gamma_Xi_R(1:end-ny,:))*Gamma_Xi_R(ny+1:end,:);

[~,~,U_tr] = createHankelMatrix(Train_Data.InputData,nB,0,mean_std_inp);
[~,~,Y_tr] = createHankelMatrix(Train_Data.OutputData,nB,0,mean_std_out);

[Phi_Z_N, Ex_obs] = phiZ(C,A,nB); 
[Phi_B_N, BZ]  = phiB(Ex_obs,U_tr,Y_tr,nB,Phi_Z_N);

[B,Zlift,del_cost,total_cost] = GradientDescent(Phi_B_N,Phi_Z_N,BZ,Y_tr,nB);
clearvars Phi_B_N phi_Z_N

Yini = ((Y_tr(1:ny,:).')*diag(mean_std_out(2,:))).' + mean_std_out(1,:).';
[~,~,~,min_D] = kmeans(Yini.', min(800,size(Yini,2)) );
[~,idx_GP] = min(min_D);
XGp = Yini(:,idx_GP);
ZGp = Zlift(:,idx_GP); 
ytest(:,:)= XGp - C*ZGp;
XGpr = XGp.';
ZGpr = ZGp.';

end