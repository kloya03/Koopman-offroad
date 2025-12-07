clc;
clear;
addpath('../../functions/utility/')
% folder = '/scratch/kloya/Koopman-offroad/scripts/koopman_training/results/sandyloam_noelev_models/models/';
folder = 'results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(1).name);
load(filename,"testData");
% testData = merge(valData,testData);
test_ntr = size(testData,4);
refresh = 250; %[25,50,75,100,125,150,175,200,225,250];
tic
parfor k = 1:length(files)
    allErr = []; allY = [];
    k
    filename = fullfile(folder, files(k).name);
    % disp(['Loading: ' filename]);

    model = load(filename,"MDL_fitr","A","B","Bc1","C",...
        "Cc1","K_obs","mean_std_out","rr","base_name");   % load the file
    model_complexity(k,:) =[model.rr];
    model_name{k} = model.base_name;
    for i=1:test_ntr
        clc;
        [k,i]
        [ypred,yout]  = K_RSSID_prediction(getexp(testData,i),...
            model.MDL_fitr,model.A,model.B,model.Bc1,...
            model.C,model.Cc1,model.K_obs,model.mean_std_out,refresh);
        err = yout - ypred;   % (2001 x 6)
        allErr = [allErr; err];       % grows to (2001*100) x 6
        allY = [allY; yout];        % (2001*100)x6
        % overallRMSE = sqrt(mean(err.^2, 1));    % 1 x 6
        % total_rmse = total_rmse + overallRMSE;  % make sure time stamp...
        %  for each traj is same otherwise ...
        % rmse = (n1*rmse1+ n2*rmse2)/(n1+n2)
    end
    range_y = max(allY) - min(allY);   % 1x6
    % ---- Metrics across ALL experiments ----
    RMSE(k,:) = sqrt(mean(allErr.^2, 1));        % 1x6 RMSE per state
    err_var(k,:) = var(allErr, 1);               % error variance per state (1x6)
    err_mean(k,:) = mean(allErr, 1);             % bias per state (1x6)
    
    % ---- Normalized error and variance ----
    NRMSE(k,:) = RMSE(k,:) ./ range_y;           % normalized RMSE
    Nerr_Var(k,:) = err_var(k,:) ./ (range_y.^2);

    % ---- Total RMSE across everything ----
    total_rmse(k,:) = sqrt(mean(allErr(:).^2));
    Total_Nrmse_var(k,:) = [mean(NRMSE(k,:)),...
        mean(RMSE(k,:)), mean(err_var(k,:)), mean(Nerr_Var(k,:))];
end

Total_Nrmse_var = [model_complexity,total_rmse,Total_Nrmse_var];
et_val = toc
%%

[rmse_val, rmse_ind] = sortrows(Total_Nrmse_var,2);
[NRMSE_val, NRMSE_ind] = sortrows([model_complexity,NRMSE],5);
[RMSE_val, RMSE_ind] = sortrows([model_complexity,RMSE],5);

save('sandy_loam_errors','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choosing exp 76 as the best for sandy loam
% 76: parameters = rr = 52,  600 200 200 8

%% RMSE VS Model Complexity
% kk = 55:3:80;
% mc = model_complexity(kk,1);
% Nerr = NRMSE(kk,:);
% lw = 2;
% for jj=1:size(NRMSE,2)
%     plot(mc,Nerr(:,jj),'-o','linewidth',lw)
%     hold on;
% end
% xlabel('Model Complexity');
% ylabel('Normalized RMSE');
% title('Normalized RMSE vs Model Complexity');
% legend(data.OutputName);
% hold off;

%% RMSE VS refresh rate
% kk = 76;
% filename = fullfile(folder, files(kk).name);
% 
% model = load(filename,"error_with_time");
% mc = model_complexity(kk,1);
% Nerr = NRMSE(kk,:);
% lw = 2;
% for jj=1:size(NRMSE,2)
%     plot(mc,Nerr(:,jj),'-o','linewidth',lw)
%     hold on;
% end
% xlabel('Model Complexity');
% ylabel('Normalized RMSE');
% title('Normalized RMSE vs Model Complexity');
% legend(data.OutputName);
% hold off;

%% RMSE VS Model Complexity
% kk = 76;
% filename = fullfile(folder, files(kk).name);
% 
% model = load(filename,"error_with_time");
% mc = model_complexity(kk,1);
% Nerr = NRMSE(kk,:);
% lw = 2;
% for jj=1:size(NRMSE,2)
%     plot(mc,Nerr(:,jj),'-o','linewidth',lw)
%     hold on;
% end
% xlabel('Model Complexity');
% ylabel('Normalized RMSE');
% title('Normalized RMSE vs Model Complexity');
% legend(data.OutputName);
% hold off;