clc;
clear;
addpath('../../functions/utility/')
folder = '/scratch/kloya/Koopman-offroad/scripts/koopman_training/results/sandyloam_noelev_models/models/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
tic
for k = randi(length(files),1)
    k
    filename = fullfile(folder, files(k).name);
    disp(['Loading: ' filename]);

    load(filename);   % load the file

    refresh = [25,50,75,100,125,150,175,200,225,250];
    test_ntr = size(testData,4);
    traj = randi(test_ntr);
    
    for jj = 1:size(refresh,2)

        time_error = zeros(size(testData(:,:,:,1).OutputData));
        total_rmse = zeros(size(testData(:,:,:,1).OutputData,2),1);
        for i=traj  %1:test_ntr
            clc;
            [k,jj,i]
            [Y_pred,y_out,error,Y_pred_95]  = K_RSSID_prediction(testData(:,:,:,i),...
                MDL_fitr,A,B,Bc1,C,Cc1,K_obs,mean_std_out,refresh(jj));
    
            time_error = time_error + error.withTime;
            total_rmse = total_rmse + error.overallRMSE.';  % make sure time stamp...
                                                  %  for each traj is same otherwise ...
                                                  % rmse = (n1*rmse1+ n2*rmse2)/(n1+n2)
        end
        error_with_time(:,:,jj) = sqrt(time_error./test_ntr);
        overall_error(:,jj) = total_rmse./test_ntr;
    end
end
display(overall_error)
et_val = toc
