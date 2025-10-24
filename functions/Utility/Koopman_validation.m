%% 9. Validation Function

function [err_val] = Koopman_validation(valData,MDL_fitr, K,B,C,K_obs)

ntr = size(valData,4);

for i=1:ntr
    [Ymean,Ycov,pos,error_vel,error_pos] = K_RSSID_prediction(valData(:,:,:,i),MDL_fitr,K,B,C,K_obs);
    
    err_val(:,i) = [error_vel;error_pos];
end

end