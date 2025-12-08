clc;
clear;


addpath('../../functions/utility/')
% folder = '/scratch/kloya/Koopman-offroad/scripts/koopman_training/results/clay_noelev_models/models/';
folder = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(1).name);
load(filename,"testData");
test_ntr = size(testData,4);
refresh = [25,50,75,100,125,150,175,200,225,250];
ny = length(testData.OutputName);
nf = 1;%length(files);
% Initialize
range_y = zeros(size(refresh,2),ny,nf); % 10x6x135
RMSE = zeros(size(refresh,2),ny,nf); % --"--
err_var = zeros(size(refresh,2),ny,nf); % --"--
err_mean = zeros(size(refresh,2),ny,nf); % --"--
NRMSE = zeros(size(refresh,2),ny,nf); % --"--
Nerr_Var = zeros(size(refresh,2),ny,nf); % --"--
total_rmse = zeros(size(refresh,2),1,nf); % 10x1x135
Total_Nrmse_var = zeros(size(refresh,2),6,nf); % 10x6x135
model_complexity = zeros(size(refresh,2),1,nf); %10x1x135

tic
k=41;
parfor jj = 1:length(refresh)
    clc;
    k
    % for jj = 1%1:length(files)

        [k,jj]
        allErr = []; allY = [];
        filename = fullfile(folder, files(k).name);
        model = load(filename,"MDL_fitr","A","B","Bc1","C",...
            "Cc1","K_obs","mean_std_out","rr","base_name");   % load the file
        model_complexity(jj,:,k) =model.rr;
        for i=1:test_ntr
            [ypred,yout]  = K_RSSID_prediction(getexp(testData,i),...
                model.MDL_fitr,model.A,model.B,model.Bc1,...
                model.C,model.Cc1,model.K_obs,model.mean_std_out,refresh(jj));
            err = yout - ypred;   % (2001 x 6)
            allErr = [allErr; err];       % grows to (2001*100) x 6
            allY = [allY; yout];        % (2001*100)x6
        end
        range_y(jj,:,k) = max(allY) - min(allY);   % 1x6
        % ---- Metrics across ALL experiments ----
        RMSE(jj,:,k) = sqrt(mean(allErr.^2, 1));        % 1x6 RMSE per state
        err_var(jj,:,k) = var(allErr, 1);               % error variance per state (1x6)
        err_mean(jj,:,k) = mean(allErr, 1);             % bias per state (1x6)
    
        % ---- Normalized error and variance ----
        NRMSE(jj,:,k) = RMSE(jj,:,k) ./ range_y(jj,:,k);           % normalized RMSE
        Nerr_Var(jj,:,k) = err_var(jj,:,k) ./ (range_y(jj,:,k).^2);
    
        % ---- Total RMSE across everything ----
        total_rmse(jj,:,k) = sqrt(mean(allErr(:).^2));
        Total_Nrmse_var(jj,:,k) = [model_complexity(jj,:,k),...
                                    total_rmse(jj,:,k),...
                                   mean(NRMSE(jj,:,k)), ...
                                   mean(RMSE(jj,:,k)), ...
                                   mean(err_var(jj,:,k)), ...
                                   mean(Nerr_Var(jj,:,k))];
    % end
end
et_val = toc

[rmse_val, rmse_ind] = sortrows(Total_Nrmse_var(:,:,end),2);
[NRMSE_val, NRMSE_ind] = sortrows([model_complexity(:,:,end),NRMSE(:,:,end)],5);
[RMSE_val, RMSE_ind] = sortrows([model_complexity(:,:,end),RMSE(:,:,end)],5);

save('clay_errors_41','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
load("clay_errors.mat")
% Choosing exp 76 as the best for sandy loam
% 76: parameters = rr = 52,

%% RMSE VS Model Complexity
% mm = {'*clay_nl400_sy200_nB200_c*.mat',...
%     '*clay_nl400_sy200_nB300_c*.mat',...
%     '*clay_nl400_sy200_nB400_c*.mat',...
%     '*clay_nl400_sy300_nB200_c*.mat',...
%     '*clay_nl400_sy300_nB300_c*.mat',...
%     '*clay_nl400_sy300_nB400_c*.mat'};

% mm = {'*clay_nl600_sy200_nB200_c*.mat',...
%     '*clay_nl600_sy200_nB300_c*.mat',...
%     '*clay_nl600_sy200_nB400_c*.mat',...
%     '*clay_nl600_sy300_nB200_c*.mat',...
%     '*clay_nl600_sy300_nB300_c*.mat',...
%     '*clay_nl600_sy300_nB400_c*.mat'...
%     '*clay_nl600_sy400_nB200_c*.mat',...
%     '*clay_nl600_sy400_nB300_c*.mat',...
%     '*clay_nl600_sy400_nB400_c*.mat'};
mm = {'*clay_nl600_sy400_nB200_c*.mat'}%...
    % '*clay_nl400_sy200_nB300_c*.mat'};
folders = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';

for i=1:size(mm,2)
    clc;
    files_int = dir(fullfile(folders, mm{i}));
    [nameMatch] = ismember({files.name}, {files_int.name});
    idx = find(nameMatch==1);
    mc = [idx.',model_complexity(idx,1)];
    idxx = [sortrows(mc,2,'ascend')];
    Nerr = [Total_Nrmse_var(idxx(:,1),2)];
    Nrmse = [NRMSE(idxx(:,1),:)];
    idxx1 = []
    if nnz(Nrmse(:,4:6)>2)>0
        continue;
    else
        figure(i+6)
        lw = 3;
        % plot(idxx(:,2),Nerr(:,1),'-o','linewidth',lw); hold on;
        for jj=4:6
            plot(idxx(:,2),Nrmse(:,jj),'-o','linewidth',lw);
            hold on;
        end
        xlabel('Model Complexity');
        ylabel('RMSE');
        grid on;
        legend(testData.OutputName(4:6),'Interpreter','latex');
        ax = gca;   % Get the current axes handle
        ax.FontSize = 25; % Set the font size to 14 points
    end

end


%% RMSE VS refresh rate
clc;
clear;
kk = 24;
addpath('../../functions/utility/')
folder = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","base_name","testData");   % load the file
base_name
test_ntr = size(testData,4);
refresh = [25,50,75,100,125,150,175,200,225,250];
parfor k=1:size(refresh,2)        
    allErr = []; allY = [];
    for i=1:test_ntr
        clc;
        [k,i]
        [ypred,yout]  = K_RSSID_prediction(getexp(testData,i),...
            MDL_fitr,A,B,Bc1,...
            C,Cc1,K_obs,mean_std_out,refresh(k));
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
    % ---- Normalized error and variance ----
    NRMSE(k,:) = RMSE(k,:) ./ range_y;           % normalized RMSE
end
save('clay_errors_refresh_41','-v7.3')
%%
clc;
clear;
load('clay_errors_refresh_41.mat')
lw = 3;
figure(2)
for jj=K_obs
    plot(0.01*refresh,RMSE(:,jj),'-o','linewidth',lw)
    hold on; grid on;
end
xlabel('Refreshing time [s]');
ylabel('RMSE');
% title('Normalized RMSE vs Refresh rate');
legend(testData.OutputName(K_obs),'Interpreter','latex');
ax = gca;   % Get the current axes handle
ax.FontSize = 20; % Set the font size to 14 points

hold off;

%% error with time
kk = 76;
filename = fullfile(folder, files(kk).name);
load(filename,"error_with_time","t_hc","K_obs");
N_err_wTime = error_with_time(:,:,end);
lw = 2;
for jj=K_obs
    plot(t_hc(1:249),N_err_wTime(1:249,jj),'-','linewidth',lw)
    hold on; grid on;
end
xlabel('Time [s]');
ylabel('RMSE');
axis([0 2 0 1])
% title('Normalized RMSE vs time');
legend(testData.OutputName(K_obs),'Interpreter','latex');
ax = gca;   % Get the current axes handle
ax.FontSize = 20; % Set the font size to 14 points

hold off;

%% Random trajectories

clc;
clear;
kk = 76;
folder = 'results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","base_name","testData","t_hc");   % load the file
test_ntr = size(testData,4);
refresh = [25,50,75,100,125,150,175,200,225,250];
k=10; tstep = 0+(1:249); lw=2;
trajj = [2,3,5,6,8,25,26,27,49,67,73,76,79,80,84,88,94,97,98,99];
for i=10%randi(size(trajj,2))%2:2:size(trajj,2)
    trajj(i)
    allErr = []; allY = [];
    [ypred,yout]  = K_RSSID_prediction(getexp(testData,trajj(i)),...
        MDL_fitr,A,B,Bc1,...
        C,Cc1,K_obs,mean_std_out,refresh(k));
    figure(i)
    subplot(2,2,1)
    plot(yout(tstep,1),yout(tstep,2),'b','linewidth',lw); hold on;
    plot(ypred(tstep,1),ypred(tstep,2),'--r','linewidth',lw);grid on;
    legend('True','K-SSID')
    xlabel('X')
    ylabel('Y')
    ax = gca;   % Get the current axes handle
    ax.FontSize = 15; % Set the font size to 14 points
    axis equal

    for jj=K_obs
        subplot(2,2,jj-2)
        plot(t_hc(tstep),yout(tstep,jj),'b','linewidth',lw); hold on;
        plot(t_hc(tstep),ypred(tstep,jj),'--r','linewidth',lw);grid on;
        ylabel(testData.OutputName(jj),'Interpreter','latex','FontSize',15)
        xlabel('Time [s]')
        xlim([0 tstep(end)*0.01])
        ax = gca;   % Get the current axes handle
        ax.FontSize = 15; % Set the font size to 14 points
    end
   

end


%% 