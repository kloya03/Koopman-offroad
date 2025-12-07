clc;
clear;
addpath('../../../functions/utility/')
folder = 'sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, 'task46*.mat'));   % or *.txt, *.csv, etc.
tic
for k = randi(length(files),1)
    k
    filename = fullfile(folder, files(k).name);
    disp(['Loading: ' filename]);

    load(filename);   % load the file

    refresh = [25,50,75,100,125,150,175,200,225,250];
    test_ntr = size(testData,4);
    traj = randi(test_ntr);

    for jj = 10

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

%%
labl = {'x', 'y', '$\psi$','u','v','$\dot{\psi}$'};
numVars = 6;
figure;
tl = 1:249;

for i = 4:numVars
    subplot(2, 2, i-3);  % creates a 3-by-2 grid of subplots
    plot(t_hc(tl),Y_pred(tl, i), 'LineWidth', 1.5);
    hold on;
    plot(t_hc(tl),y_out(tl, i), 'LineWidth', 1.5);

    % title(['Variable ' num2str(i)]);
    xlabel('time');
    ylabel(labl{i},'Interpreter','latex');
    legend('Prediction', 'Output');
    grid on;
end

hold on;
subplot(2, 2, 4);  % creates a 3-by-2 grid of subplots
plot(Y_pred(1:249, 1),Y_pred(1:249, 2), 'LineWidth', 1.5);
hold on;
plot(y_out(1:249, 1),y_out(1:249, 2), 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('trajectory');
axis equal
