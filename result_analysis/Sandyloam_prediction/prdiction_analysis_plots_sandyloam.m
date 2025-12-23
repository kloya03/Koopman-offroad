clc;
clear;
addpath('../../functions/utility/')
% folder = '/scratch/kloya/Koopman-offroad/scripts/koopman_training/results/sandyloam_noelev_models/models/';
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(1).name);
load(filename,"testData");
% testData = merge(valData,testData);
test_ntr = size(testData,4);
refresh = 250; %[25,50,75,100,125,150,175,200,225,250];
ny = length(testData.OutputName);
nf = length(files);
range_y = zeros(ny,nf); % 10x6x135
RMSE = zeros(ny,nf); % --"--
err_var = zeros(ny,nf); % --"--
err_mean = zeros(ny,nf); % --"--
NRMSE = zeros(ny,nf); % --"--
Nerr_Var = zeros(ny,nf); % --"--
total_rmse = zeros(1,nf); % 10x1x135
Total_Nrmse_var = zeros(6,nf); % 10x6x135
model_complexity = zeros(1,nf); %10x1x135
tic
for k = 1:length(files)
    allErr = []; allY = [];
    k
    filename = fullfile(folder, files(k).name);
    % disp(['Loading: ' filename]);

    model = load(filename,"MDL_fitr","A","B","Bc1","C",...
        "Cc1","K_obs","mean_std_out","rr","task_id","base_name");   % load the file
    model_complexity(k,:) =[model.rr];
    taskNo(k,:) = model.task_id;
    % model_name{k} = model.base_name;
    for i=1:test_ntr
        [ypred,yout]  = K_RSSID_prediction(getexp(testData,i),...
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
    Total_Nrmse_var(k,:) = [taskNo(k,:),model_complexity(k,:),mean(NRMSE(k,:)),...
        mean(RMSE(k,:)), mean(err_var(k,:)), mean(Nerr_Var(k,:))];

end
et_val = toc

[rmse_val, rmse_ind] = sortrows(Total_Nrmse_var,2);
[NRMSE_val, NRMSE_ind] = sortrows([taskNo,model_complexity,NRMSE],6);
[RMSE_val, RMSE_ind] = sortrows([taskNo,model_complexity,RMSE],6);

save('sandy_loam_errors_all','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
load("sandy_loam_errors.mat")
% Choosing exp 76 as the best for sandy loam
% 76: parameters = rr = 52,

%% RMSE VS Model Complexity
% mm = {'*sandyloam_nl400_sy200_nB200_c*.mat',...1
% '*sandyloam_nl600_sy300_nB400_c*.mat'...12
% '*sandyloam_nl400_sy300_nB200_c*.mat',...4
% '*sandyloam_nl400_sy300_nB300_c*.mat',...5
% '*sandyloam_nl400_sy300_nB400_c*.mat'...6
% '*sandyloam_nl600_sy400_nB300_c*.mat',...14
% '*sandyloam_nl600_sy400_nB400_c*.mat'}; %15
mm = {'*sandyloam_nl400_sy200_nB300_c*.mat',...2
    '*sandyloam_nl400_sy200_nB400_c*.mat',...3
    '*sandyloam_nl600_sy200_nB200_c*.mat',...7
    '*sandyloam_nl600_sy200_nB300_c*.mat',...8
    '*sandyloam_nl600_sy200_nB400_c*.mat',...9
    '*sandyloam_nl600_sy300_nB200_c*.mat',...10
    '*sandyloam_nl600_sy300_nB300_c*.mat',...11
    '*sandyloam_nl600_sy400_nB200_c*.mat'};%,...13
% mm = {'*sandyloam_nl400_sy300_nB200_c*.mat'};% 4
addpath('../../functions/utility/')
% folder = '/scratch/kloya/Koopman-offroad/scripts/koopman_training/results/clay_noelev_models/models/';
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.

for i=6%1:size(mm,2)
    clc;
    files_int = dir(fullfile(folder, mm{i}));
    [nameMatch] = ismember({files.name}, {files_int.name});
    idx = find(nameMatch==1);
    mc = [idx.',model_complexity(idx,1)];
    idxx = [sortrows(mc,2,'ascend')];
    % Nerr = [Total_Nrmse_var(idxx(:,1),2)];
    % Nerr1 = [Nerr(1:5,1);0.212;Nerr(6:end,1)];
    idxx1 = [idxx(1:5,:);[76,52];idxx(6:end,:)];
    Nrmse = [NRMSE(idxx(:,1),:)];
    Nrmse = Nrmse./max(Nrmse);
    MNrmse = mean(Nrmse,2);
    % Nrmse(7,4) = 0.051;
    % idxx(3:4,2) = idxx(3:4,2)-18;
    % if nnz(Nrmse(:,4:6)>5)<0
    %     continue;
    % else
    figure(1)
    lw = 5;
    
    % plot(idxx(:,2),Nerr(:,1),'-o','linewidth',lw); hold on;
    for jj=4:6
        plot(idxx(:,2)-4,Nrmse(:,jj),'-o','linewidth',lw);
        hold on;
    end
    plot(idxx(:,2)-4,MNrmse,'--k','linewidth',lw); hold on;
    xlabel('Model Order (r)' );
    ylabel('N-RMSE');
    grid on; box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    % labels = [testData.OutputName{4:6},{'mean'}]
    legend(testData.OutputName{4:6},'Interpreter','latex','FontSize',35);

    % end

end


%% RMSE VS refresh rate
clc;
clear;
kk = 123;
addpath('../../functions/utility/')
folder = '../../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","base_name");   % load the file
load('../../datasets/sandyloam_100hz_elev_experiment_1572.mat','testData','b');
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
save('clay_on_sandyloam_errors_noelev_refresh_123','-v7.3')
%%
% clc;
% clear;
% models = {'sandyloam_errors_refresh_124.mat', 'sandyloam_errors_elev_refresh_124.mat'}
for jk=1:2
    load(models{jk})
    lw = 5;
    figure(5)
    for jj=K_obs
        plot(0.01*refresh,RMSE(:,jj),'-o','linewidth',lw)
        hold on; grid on;
    end
    xlabel('Refreshing time [s]');
    ylabel('RMSE');
    % title('Normalized RMSE vs Refresh rate');
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 35; % Set the font size to 14 points
    legend(testData.OutputName(K_obs),'Interpreter','latex','FontSize',35);
    axis([0 3 0 0.8])
    hold off;
end

%%
clc;
clear;
models = {'sandyloam_errors_refresh_124.mat', 'sandyloam_errors_elev_refresh_124.mat'}
noelev = load(models{1},'NRMSE','refresh','mean_std_out');
elev = load(models{2},'NRMSE');
perc = 100*(elev.NRMSE-noelev.NRMSE)./noelev.NRMSE;
AA = perc(:,4:6)
mean(AA)
%% error with time
kk = 124;
addpath('../../functions/utility/')
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
load(filename,"error_with_time","t_hc","K_obs");
N_err_wTime = error_with_time(:,:,end);
lw = 5;
figure(kk)
for jj=K_obs
    plot(t_hc(1:249),N_err_wTime(1:249,jj),'-','linewidth',lw)
    hold on; grid on;
end
xlabel('Time [s]');
ylabel('RMSE (t) ');
axis([0 2.5 0 0.6])
% title('Normalized RMSE vs time');
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 35; % Set the font size to 14 points
legend(testData.OutputName(K_obs),'Interpreter','latex','FontSize',35);
hold off;

%% Random trajectories

clc;
clear;
kk = 124;
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","base_name","testData","t_hc");   % load the file
test_ntr = size(testData,4);
refresh = [25,50,75,100,125,150,175,200,225,250];
k=6; tstep = 0+(1:2001); lw=4;
trajj = [2,3,5,6,8,25,26,27,49,67,73,76,79,80,84,88,94,97,98,99];
for i=4;%randi(size(trajj,2))%2:2:size(trajj,2)
    trajj(i)
    allErr = []; allY = [];
    [ypred,yout,~,ypred_95]  = K_RSSID_prediction(getexp(testData,trajj(i)),...
        MDL_fitr,A,B,Bc1,...
        C,Cc1,K_obs,mean_std_out,refresh(k));
    figure(1)
    subplot(2,2,1)
    plot(yout(tstep,1),yout(tstep,2),'b','linewidth',lw); hold on;
    plot(ypred(tstep,1),ypred(tstep,2),'r','linewidth',lw-1);grid on;
    legend('True','K-SSID')
    xlabel('X')
    ylabel('Y')
    box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 35; % Set the font size to 14 points
    hold off;

    for jj=K_obs
        subplot(2,2,jj-2)
        plot(t_hc(tstep),yout(tstep,jj),'b','linewidth',lw); hold on;
        plot(t_hc(tstep),ypred(tstep,jj),'r','linewidth',lw);hold on;
        % plot(t_hc(tstep),squeeze(ypred_95(1,tstep,jj-3)),'-k','linewidth',lw);hold on;
        % plot(t_hc(tstep),squeeze(ypred_95(2,tstep,jj-3)),'-k','linewidth',lw);grid on;
        ylabel(testData.OutputName(jj),'Interpreter','latex','FontSize',35)
        xlabel('Time [s]')
        xlim([0 tstep(end)*0.01])
        box on;
        set(gca, 'LineWidth', 1.5)
        ax = gca;   % Get the current axes handle
        ax.FontSize = 35; % Set the font size to 14 points
        grid on;
        hold off;
        if jj==4
            ylim([min(yout(tstep,jj))-2 2+max(yout(tstep,jj))])
        end
    end


end
%%
figure(1)
for i=randi(20)
    i
    yout = getexp(testData,i);
    plot(yout.OutputData(:,1),yout.OutputData(:,2),'b','linewidth',lw); hold off;
end


%% eigenvalue plot
% clc;
% clear;
kk =124
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
load(filename,"A")
eigvals = eig(A);     % typically you get them from a matrix
% Compute magnitudes
mag = abs(eigvals);
figure; hold on; axis equal;
% --- Plot unit circle ---
theta = linspace(0, 2*pi, 400);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.5);  % unit circle
% --- Scatter eigenvalues with color = magnitude ---
scatter(real(eigvals), imag(eigvals), 80, mag, 'filled');
% --- Colorbar ---
cb = colorbar('Ticks',[0.95, 0.96, 0.97, 0.98, 0.99, 1]);
ylabel(cb, 'Eigenvalue Magnitude');
% --- Formatting ---
xlabel('Real');
ylabel('Imaginary');
% title('Eigenvalues in Complex Plane with Magnitude Colorbar');
grid on;
box on;
set(gca, 'LineWidth', 1.5)
cb.Ticks = [0.95, 0.96, 0.97, 0.98, 0.99, 1];
xlim([-1.05 1.05]);
ylim([-1.05 1.05]);
ax = gca;   % Get the current axes handle
ax.FontSize = 20; % Set the font size to 14 points