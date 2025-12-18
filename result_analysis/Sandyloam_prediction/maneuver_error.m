clc;
clear;
addpath('../../functions/utility/')
load('../expt_name.mat')
load('../../datasets/sandyloam_100hz_no_elev_experiment_1579.mat','allindices')
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(123).name);
model = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","task_id","base_name","trainData"); 
refresh = 250; %[25,50,75,100,125,150,175,200,225,250];
nf = length(model.trainData.OutputName);
maneuvers = {'straight','circle','multisine','slalom','fishhook'};
ny = length(maneuvers);
range_y = zeros(ny,nf); % 10x6x5
RMSE = zeros(ny,nf); % --"--
err_var = zeros(ny,nf); % --"--
err_mean = zeros(ny,nf); % --"--
NRMSE = zeros(ny,nf); % --"--
Nerr_Var = zeros(ny,nf); % --"--
total_rmse = zeros(1,nf); % 10x1x135
Total_Nrmse_var = zeros(ny,4); % 10x6x135
tic
allErr = []; allY = [];
expname = cell2mat(expt_name(1,:));
for k=1:ny
    maneuvers{k}
    test_ntr = find(expname(1,1:1423)==maneuvers{k});
    for i=1:length(test_ntr)
        traj = find(test_ntr(i)==allindices);
        [ypred,yout]  = K_RSSID_prediction(getexp(model.trainData,test_ntr(i)),...
            model.MDL_fitr,model.A,model.B,model.Bc1,...
            model.C,model.Cc1,model.K_obs,model.mean_std_out,refresh);
        err = yout - ypred;   % (2001 x 6)
        allErr = [allErr; err];       % grows to (2001*100) x 6
        allY = [allY; yout];        % (2001*100)x6
    end
    range_y(k,:) = max(allY) - min(allY);   % 1x6
    % ---- Metrics across ALL experiments ----
    RMSE(k,:) = sqrt(mean(allErr.^2, 1));        % 1x6 RMSE per state
    err_var(k,:) = var(allErr, 1);               % error variance per state (1x6)
    err_mean(k,:) = mean(allErr, 1);             % bias per state (1x6)
    
    % ---- Normalized error and variance ----
    NRMSE(k,:) = RMSE(k,:) ./ range_y(k,:);           % normalized RMSE
    Nerr_Var(k,:) = err_var(k,:) ./ (range_y(k,:).^2);
    
    % ---- Total RMSE across everything ----
    total_rmse(k,:) = sqrt(mean(allErr(:).^2));
    Total_Nrmse_var(k,:) = [mean(NRMSE(k,:)),...
        mean(RMSE(k,:)), mean(err_var(k,:)), mean(Nerr_Var(k,:))];
end
save('RMSE_maneuvers_sl_1','-v7.3')


%%
clc;
clear;
addpath('../../functions/utility/')
load('../expt_name.mat')
load('../../datasets/clay_100hz_no_elev_experiment_1472.mat','allindices')
folder = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(124).name);
model = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","task_id","base_name","trainData"); 
refresh = 250; %[25,50,75,100,125,150,175,200,225,250];
nf = length(model.trainData.OutputName);
maneuvers = {'straight','circle','multisine','slalom','fishhook'};
ny = length(maneuvers);
range_y = zeros(ny,nf); % 10x6x5
RMSE = zeros(ny,nf); % --"--
err_var = zeros(ny,nf); % --"--
err_mean = zeros(ny,nf); % --"--
NRMSE = zeros(ny,nf); % --"--
Nerr_Var = zeros(ny,nf); % --"--
total_rmse = zeros(1,nf); % 10x1x135
Total_Nrmse_var = zeros(ny,4); % 10x6x135
tic
allErr = []; allY = [];
expname = cell2mat(expt_name(1,:));
for k=1:ny
    maneuvers{k}
    test_ntr = find(expname(1,1:1143)==maneuvers{k});
    for i=1:length(test_ntr)
        traj = find(test_ntr(i)==allindices);
        [ypred,yout]  = K_RSSID_prediction(getexp(model.trainData,test_ntr(i)),...
            model.MDL_fitr,model.A,model.B,model.Bc1,...
            model.C,model.Cc1,model.K_obs,model.mean_std_out,refresh);
        err = yout - ypred;   % (2001 x 6)
        allErr = [allErr; err];       % grows to (2001*100) x 6
        allY = [allY; yout];        % (2001*100)x6
    end
    range_y(k,:) = max(allY) - min(allY);   % 1x6
    % ---- Metrics across ALL experiments ----
    RMSE(k,:) = sqrt(mean(allErr.^2, 1));        % 1x6 RMSE per state
    err_var(k,:) = var(allErr, 1);               % error variance per state (1x6)
    err_mean(k,:) = mean(allErr, 1);             % bias per state (1x6)
    
    % ---- Normalized error and variance ----
    NRMSE(k,:) = RMSE(k,:) ./ range_y(k,:);           % normalized RMSE
    Nerr_Var(k,:) = err_var(k,:) ./ (range_y(k,:).^2);
    
    % ---- Total RMSE across everything ----
    total_rmse(k,:) = sqrt(mean(allErr(:).^2));
    Total_Nrmse_var(k,:) = [mean(NRMSE(k,:)),...
        mean(RMSE(k,:)), mean(err_var(k,:)), mean(Nerr_Var(k,:))];
end

save('RMSE_maneuvers_clay_1','-v7.3')