%% Data Gather and plotting
clc;
clear;
close all;

startIndex = 1;  % Starting index
endIndex = 1600; % Ending index
missingFiles=[];
foldername = '../../datasets/sandyloam_noElev_dataset_100hz';
load('../../datasets/Offroad_InputsSignals','traj_name')
fails=[]; incomp =[]; comp = [];
exp_u_incomp = [];exp_y_incomp = {};
exp_x_incomp = {};exp_t_incomp = {};
exp_bekk_h_incomp = {};exp_dZdt_incomp = {};
exp_fv_incomp = {};exp_sig_tau_incomp = {};
ct = 0; ctf = 0;
for iter=startIndex:endIndex
%     clc;
    iter

    % Generate the filename
    fileName = sprintf('data_%d.mat', iter);
    fileName = sprintf('%s/%s',foldername,fileName);

    % Check if the file exists in the current directory
    if ~isfile(fileName)
        % Add the missing file number to the list
        missingFiles = [missingFiles, iter, 0];
        continue
    end
    load(fileName)
    % Ht = b.h(Z_hc(:,1),Z_hc(:,2));
    if  abs(max(Z_hc(:,4))) > 20
        fails = [fails, iter]
    elseif size(t_hc,1) == 2001
        i=1;
        ct = ct+1;
        plotFigures(i,Z_hc,t_hc,fv,"complete")
        Data_states = Z_hc(:,1:12);
        % Data_observables = [dZdt(:,4:5), dZdt(:,8), dZdt(:,3), dZdt(:,9)];
        Data_observables = [Z_hc(:,1:6)];%,Z_hc(:,4:5), Z_hc(:,6)];
        Data_inputs = [fv(:,12), fv(:,11)];
        exp_u{1,ct} = Data_inputs;
        exp_y{1,ct} = Data_observables;
        exp_x{1,ct} = Data_states;
        exp_bekk_h{1,ct} = bekk_h;
        exp_dZdt{1,ct} = dZdt;
        exp_fv{1,ct} = fv;
        exp_sig_tau{1,ct} = sig_tau;
        comp = [comp, iter];
        exp_name{1,ct} = traj_name(:,iter);

    else
        clc;
        i=10;
        ctf = ctf +1
        % plotFigures(i,Z_hc,t_hc,fv,"incomplete")
        incomp = [incomp, iter];
        Data_states = Z_hc(:,1:12);
        % Data_observables = [dZdt(:,4:5), dZdt(:,8), dZdt(:,3), dZdt(:,9)];
        Data_observables = [Z_hc(:,1:6)];%,Z_hc(:,4:5), Z_hc(:,6)];
        Data_inputs = [fv(:,12), fv(:,11)];
        exp_u_incomp{1,ctf} = Data_inputs;
        exp_y_incomp{1,ctf} = Data_observables;
        exp_x_incomp{1,ctf} = Data_states;
        exp_t_incomp{1,ctf} = t_hc;
        exp_bekk_h_incomp{1,ctf} = bekk_h;
        exp_dZdt_incomp{1,ctf} = dZdt;
        exp_fv_incomp{1,ctf} = fv;
        exp_sig_tau_incomp{1,ctf} = sig_tau;
    end

end
ind_ref{1} = "X, Y, psi, u, v, dpsi, z, dz, theta, dtheta, wf, wr";
ind_ref{2} = "fv = [Flf, Fcf, Nf, Frr_f, Flr, Fcr, Nr, Frr_r, f_adx, f_ady,"+...
    "b.tau_t, delta, vlf, vcf, omega_f, b.sf, beta_f vlr, vcr, omega_r, b.sr, beta_r]";
ind_ref{3} = "bekk_h = [-int_sig_sin_f, int_tau_cos_f, thf_f, thr_f, thm_f, hf_f,"+...
    "-int_sig_sin_r, int_tau_cos_r, thf_r,thr_r, thm_r, hf_r]";
ind_ref{4} = "sig_tau = [sig_f, tau_xf, sig_r, tau_xr]";


%% Divide into Test, Training and validation set

numTraj = size(exp_y,2);
numTest = 100;
numVal = 50;
numTrain = numTraj-numVal-numTest;
allindices = randperm(numTraj);
expt_name = exp_name(:,allindices);
% training Data
trainData = iddata(exp_y(:, allindices(1:numTrain)),...
    exp_u(:, allindices(1:numTrain)), b.dt);
trainData.OutputName = {'x','y','$\psi$','$u$', '$v$', '$\dot{\psi}$'};
trainData.InputName = {'$\delta$', '$\tau$'};

% validation data
valData = iddata(exp_y(:, allindices(numTrain+1:numTrain+numVal)),...
    exp_u(:, allindices(numTrain+1:numTrain+numVal)), b.dt);
valData.OutputName = {'x','y','$\psi$','$u$','$v$','$\dot{\psi}$'};
valData.InputName = {'$\delta$', '$\tau$'};

% test Data
testData = iddata(exp_y(:, allindices(numTrain+numVal+1:end)),...
    exp_u(:, allindices(numTrain+numVal+1:end)), b.dt);
testData.OutputName = {'x','y','$\psi$','$u$','$v$','$\dot{\psi}$'};
testData.InputName = {'$\delta$','$\tau$'};
% % 
% figure(5)
% idplot(valData)
% figure(6)
% idplot(testData)
% % 
test_xy = reshape(cell2mat(testData(:,[1,2]).OutputData),...
    size(testData.OutputData{1},1),2,size(testData.OutputData,2));
val_xy = reshape(cell2mat(valData(:,[1,2]).OutputData),...
    size(valData.OutputData{1},1),2,size(valData.OutputData,2));
train_xy = reshape(cell2mat(trainData(:,[1,2]).OutputData),...
    size(trainData.OutputData{1},1),2,size(trainData.OutputData,2));

%% Save all
% save('../../datasets/clay_100hz_no_elev_experiment_1572.mat',...
%     "exp_x_incomp","exp_t_incomp","exp_y_incomp","exp_u_incomp","incomp",...
%     "exp_x","exp_y","exp_u","missingFiles","b","t_hc","exp_sig_tau_incomp",...
%     "exp_fv_incomp","fails","exp_dZdt_incomp","exp_bekk_h_incomp",...
%     "exp_sig_tau","exp_fv","exp_dZdt","exp_bekk_h","ind_ref","numTraj",...
%     "numVal","numTrain","numTest","comp","allindices","trainData",...
%     "testData","valData","test_xy","train_xy","val_xy",'-v7.3');


%% plot training validation and testing data
% for i = 1:size(trainData.OutputData,2)
% 
%     figure(7)
% 
%     subplot(3,1,1)
%     hold on;
%     grid on;
%     plot(train_xy(:,1,i),train_xy(:,2,i), 'linewidth', 1.5);
%     xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
%     ylabel('$Y$', 'fontsize', 20, 'interpreter', 'latex');
%     title('training dataset traj')
%     axis equal;
% 
%     if i <= size(testData.OutputData,2)
%         subplot(3,1,2)
%         hold on;
%         grid on;
%         plot(val_xy(:,1,i),val_xy(:,2,i), 'linewidth', 1.5);
%         xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
%         ylabel('$Y$', 'fontsize', 20, 'interpreter', 'latex');
%         title('validation dataset traj')
%         axis equal;
% 
%         subplot(3,1,3)
%         hold on;
%         grid on;
%         plot(test_xy(:,1,i),test_xy(:,2,i), 'linewidth', 1.5);
%         xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
%         ylabel('$Y$', 'fontsize', 20, 'interpreter', 'latex');
%         title('testing dataset traj')
%         axis equal;
%     end
% end

%% save figures
% figs = findall(0, 'Type', 'figure');
% for i = 1:length(figs)
%     savefig(figs(i), sprintf('figure_%d.fig', i));
%     exportgraphics(figs(i), sprintf('figure_%d.png', i));
% end

%% Functions

function plotFigures(i,Z_hc, t_hc, fv, plot_title)

% XY trajectory
figure(i);
hold on;
grid on;
plot(Z_hc(:,1), Z_hc(:,2), 'linewidth', 1.5);
plot(Z_hc(1,1), Z_hc(1,2), 'ok', 'linewidth', 3);
xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$Y$', 'fontsize', 20, 'interpreter', 'latex');
title(plot_title)
axis equal;

% Velocities
% figure(i+1);
% ylabs = {'$\psi$', '$u$', '$v$', '$\dot{\psi}$'};
% for jj = 3:6
%     subplot(2, 2, jj-2);
%     hold on;
%     grid on;
%     plot(t_hc, Z_hc(:,jj), 'linewidth', 1.5);
%     ylabel(ylabs{jj-2}, 'fontsize', 20, 'interpreter', 'latex');
%     grid on;
% end
% title(plot_title)

% Inputs
% figure(i+2);
% subplot(2, 1, 1);
% stairs(t_hc, fv(:,12), 'linewidth', 1.5);
% hold on;
% ylabel('$\delta$', 'fontsize', 20, 'Interpreter', 'latex');
% grid on;
% subplot(2, 1, 2);
% stairs(t_hc, fv(:,11), 'linewidth', 1.5);
% hold on;
% ylabel('$\tau$', 'fontsize', 20, 'Interpreter', 'latex');
% grid on;
% title(plot_title)

% % Observables
% figure(i+3)
% hold on;
% hold on; grid on;
% scatter3(Z_hc(:,4),Z_hc(:,5),Z_hc(:,6),'b','filled','MarkerFaceAlpha',0.1,'MarkerEdgeColor','none')
% xlabel('$u$','fontsize',20,'interpreter','latex')
% ylabel('$v$','fontsize',20,'interpreter','latex')
% zlabel('$\dot{\psi}$','fontsize',20,'interpreter','latex')
% title(plot_title)

end


%% 3) Quick PE diagnostics on a random subset (Hankel rank on 5 s windows)
% nu = 700;                              % past window (tune to your ID plan)
% kchk = 1:1582;%randperm(n_traj, min(50,n_traj));
% ct = 0;cc = 0;
% for k = kchk
%     clc;
%     k
%     U = createHankelMatrix(exp_u{1,k}, nu, nu);
%     rU = rank(U);
%     if rU < 2*nu
%         ct = ct+1;
%         Ut = createHankelMatrix(exp_u{1,k}(:,2), nu,nu);rUt = rank(Ut);
%         Ud = createHankelMatrix(exp_u{1,k}(:,1), nu,nu);rUd = rank(Ud);
%         fail_tr(ct,:) = [k,rUt,rUd,rU];
%         fprintf('traj%4d:rank(U)=%d(target~%d),rank(U_t,D)=%d,%d,(target~%d)\n',...
%             fail_tr(ct), rU, nu, rUt, rUd, nu);
%     else
%         cc = cc+1;
%         pass_tr(cc,:) = [k,rU];
%     end
% end

% function [H1,H2] = createHankelMatrix(data,nl,sy)
% % Construct Hankel matrix using vectorized indexing
%     H = [];
%     if iscell(data)
%         for i=1:size(data,2)
%             D = data{i};
%             [numPoints, numFeatures] = size(D);
%             nc = numPoints-nl+1;
%             nr = numFeatures*nl;
%             idx = (0:nl-1)'; % Create index offsets for time steps
%             H = [H, reshape(D(idx + (1:nc), :).', nr, nc)];
%         end
%     else
%         [numPoints, numFeatures] = size(data);
%         nc = numPoints-nl+1;
%         nr = numFeatures*nl;
%         idx = (0:nl-1)'; % Create index offsets for time steps
%         H = reshape(data(idx + (1:nc), :).', nr, nc);
%     end
%     H1 = H(1:numFeatures*sy,:);
%     H2 = H(numFeatures*sy+1:end,:);
% end