%% 10. Koopman Prediction

function [Ymean,Ycov,pos,error_vel,error_pos] = K_RSSID_prediction(data,MDL_fitr,K,B,C,K_obs)
% data : 1 experiment with n sample time points and input output

rr = size(K,1);
y_out = data.OutputData;
u_in = (data.InputData).';
ts = data.Ts;
tspan = data.SamplingInstants - ts;

% GP initial condition
Z_Gpr = zeros(rr,2);
for i = 1:rr
    [Zi_mean,Zi_sd,Zi_95] = predict(MDL_fitr(i).gprMDL,y_out(1,K_obs));
    Z_Gpr(i,:) = [Zi_mean Zi_sd.^2];
    Z_95(i,:) = Zi_95;
end
pos0 = (y_out(1,1:3)).';
pos = pos0;
Ymean = [(y_out(1,1:3)).';C*Z_Gpr(:,1)];
Ycov = diag(C*Z_Gpr(:,2)*C.');
for i = 2:size(tspan,1) 

    if mod(i,1000)==0
        Ymean_corr = y_out(i,K_obs).';
        Ymean = [Ymean Ymean_corr];
        for j = 1:rr
            [Zi_mean,Zi_sd,Zi_95] = predict(MDL_fitr(j).gprMDL,Ymean_corr.');
            Z_Gpr(j,:) = [Zi_mean Zi_sd.^2];
            Z_95(j,:) = Zi_95;
        end
    end

    Z_Gpr(:,1) = K*Z_Gpr(:,1) + B*u_in(:,i-1);
    Z_Gpr(:,2) = diag(K*diag(Z_Gpr(:,2))*K.');
    Ymean = [Ymean C*Z_Gpr(:,1)];
    Ycov = [Ycov C*diag(Z_Gpr(:,2))*C.'];
    dx = Ymean(1,i).*cos(pos0(3,1)) - Ymean(2,i).*sin(pos0(3,1));
    dy = Ymean(1,i).*sin(pos0(3,1)) + Ymean(2,i).*cos(pos0(3,1));
    pos_cur = pos0+[dx;dy;Ymean(3,i)].*ts;
    pos =  [pos pos_cur];
    pos0 = pos_cur;
end

error_vel = rmse(Ymean,y_out(:,K_obs).',2);
error_pos = rmse(pos,y_out(:,1:3).',2);

end