%% 10. Koopman Prediction

function [Y_pred,y_out,error,Y_pred_95] = K_RSSID_prediction(data,MDL_fitr,...
    A,B,Bc1,C,Cc1,K_obs,mean_std_out,refresh)
% data : 1 experiment with n sample time points and input output

rr = size(A,1);
y_out = data.OutputData;
u_in = (data.InputData).';
ts = data.Ts;
tspan = data.SamplingInstants - ts;
Y_pred = zeros(size(y_out));
Y_cov = zeros(size(K_obs,2),size(y_out,1));
% GP initial condition
Z_Gpr = zeros(rr,2);
for i = 1:rr
    X0n = (y_out(1,K_obs) - mean_std_out(1,:))./mean_std_out(2,:);  % normalize back for GP 
    [Zi_mean,Zi_sd,Zi_95] = predict(MDL_fitr(i).gprMDL,X0n);
    Z_Gpr(i,:) = [Zi_mean Zi_sd.^2];
    Z_95(i,:) = Zi_95;
end
pos0 = (y_out(1,1:3)).';
Ymean = C*Z_Gpr(:,1)+Cc1;            % added normalized terms for unnormalized
Y_pred(1,:) = [pos0.', Ymean.'];
Y_cov(:,1) = diag(C*diag(Z_Gpr(:,2))*C.');
for i = 2:size(tspan,1) 

    if mod(i,refresh)==0
        X0_corrn = ((y_out(i,K_obs) - mean_std_out(1,:))./mean_std_out(2,:)).'; % normalize back for GP 
        Ymean = y_out(i,K_obs).';
        pos_cur = y_out(i,1:3).';
        for j = 1:rr
            [Zi_mean,Zi_sd,Zi_95] = predict(MDL_fitr(j).gprMDL,X0_corrn.');
            Z_Gpr(j,:) = [Zi_mean Zi_sd.^2];
            Z_95(j,:) = Zi_95;
        end
        
    else
        % dx = Ymean(1,1).*cos(pos0(3,1)) - Ymean(2,1).*sin(pos0(3,1));
        % dy = Ymean(1,1).*sin(pos0(3,1)) + Ymean(2,1).*cos(pos0(3,1));
        % pos_cur = pos0+[dx;dy;Ymean(3,end)].*ts;
        Z_Gpr(:,1) = A*Z_Gpr(:,1) + B*u_in(:,i-1) + Bc1;  % added normalized terms for unnormalized
        Z_Gpr(:,2) = diag(A*diag(Z_Gpr(:,2))*A.');
        Y_cov(:,i) = diag(C*diag(Z_Gpr(:,2))*C.');
        Ymean = C*Z_Gpr(:,1)+Cc1;                 % added normalized terms for unnormalized
        dx = Ymean(1,1).*cos(pos0(3,1)) - Ymean(2,1).*sin(pos0(3,1));
        dy = Ymean(1,1).*sin(pos0(3,1)) + Ymean(2,1).*cos(pos0(3,1));
        pos_cur = pos0+[dx;dy;Ymean(3,end)].*ts;
        
    end
    pos0 = pos_cur;
    Y_pred(i,:) =  [pos0.', Ymean.'];  
    Y_cov(:,i) = diag(C*diag(Z_Gpr(:,2))*C.');
end

if nargout > 3
    % Confidence Interval
    Ystd = sqrt(Y_cov);
    CI95 = tinv([0.025 0.975],inf);
    for i=1:size(K_obs,2)
        YCI_95(:,:,i) = bsxfun(@times, Ystd(i,:), CI95(:));
        Y_pred_95(:,:,i) = (Y_pred(:,K_obs(i))).' + YCI_95(:,:,i);
    end
end

error.withTime = (Y_pred - y_out).^2;
error.overallRMSE = rmse(Y_pred,y_out,1);

end