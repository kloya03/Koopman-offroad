%% 9. Validation Function

function [err_val] = Koopman_validation(Data,MDL_fitr, K,B,C,K_obs)
% error_val : No. of observation x no. of data samples (trajectories)

ntr = size(Data,4);

for i=1:ntr
    [Ymean,pos,y_out,Ycov,error_vel,error_pos] = K_RSSID_prediction(Data(:,:,:,i),MDL_fitr,K,B,C,K_obs);
    
    err_val(:,i) = [error_vel;error_pos];
end

end