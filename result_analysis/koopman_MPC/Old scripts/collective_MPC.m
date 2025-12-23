%%%%%%%%% Main function for MPC %%%%%%% %%
% function [U] = MPC_koopman(K,B,C,x0,yref,[tsim,tpred,tcont],str,sys)
addpath('../../functions/utility/')
folder = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(124).name);
model = load(filename,"MDL_fitr","A","B","Bc1","C",...
        "Cc1","K_obs","mean_std_out","rr","valData","b"); % load the file
expNo = randi(50)
data = getexp(model.valData,expNo);
dt = data.Ts{1};
yref = model.valData.OutputData
tsim =10; % simulation time
tpred =2; % prediction horizon
tcont = 1; % control horizon
T = [tsim,tpred,tcont,dt];
t1 = 0:dt:tsim;
yref = valData;
x0 = [-0.7160;-0.9789]; %
cost_param = [5,10,0.0001;5,10,0]; % Q Qn R

%% SSID Koopman MPC with 35 basis functions
tic
[U_GPK,X_GPK,cost_GPK,ff_GPK] = MPC_koopman(model,x0,yref,T,str,cost_param);
et_GPK = toc;

%% Nonlinear MPC
% tic
% TNL = T; 
% TNL(:,4) = 0.1;
% [U_NL,X_NL,cost_NL,ff_NL,t_NL] = NMPC_duf(x0,TNL,yref,CV);
% et_NL = toc;

%% [co,C1,C2,Cu] = CostFcn(X_NL.',U_NL.',CV,yref);
% [co,C1,C2,Cu] = CostFcn(X_GPK.',U_GPK,CV,yref);
% [co,C1,C2,Cu] = CostFcn(X_K6.',U_K6,CV,yref);

%% Plots
lw=3; ls=30;
figure
plot(t1,X_GPK(1,:),'-','LineWidth',lw); hold on
plot(t1,X_K6(1,:),':','LineWidth',lw); hold on
plot(t_NL,X_NL(1,:),'--','LineWidth',lw); hold on
plot(t1(1:80:end),yref(1,1:80:end),'-squarek','LineWidth',lw-2);
xlabel('Time [s]')
ylabel('$x_1$','Interpreter','latex')
grid on
set(gca,'fontsize',ls)
xlim([0 8])
hleg = legend('GP-SSID','K-EDMD',...
    'Nonlinear Model','Target','fontsize',15,'location','best');
set(hleg,'color','none');
ax1 = gca;
% 
figure
plot(t1,X_GPK(2,:),'-','LineWidth',lw); hold on
plot(t1,X_K6(2,:),':','LineWidth',lw); hold on
plot(t_NL,X_NL(2,:),'--','LineWidth',lw); hold on
plot(t1(1:80:end),yref(1,1:80:end),'-squarek','LineWidth',lw-2);
xlabel('Time [s]')
ylabel('$x_2$','Interpreter','latex')
grid on
set(gca,'fontsize',ls)
xlim([0 8])
ax2 = gca;

figure
stairs(t1(1:end-1),U_GPK,'-','LineWidth',lw); hold on
stairs(t1(1:end-1),U_K6,':','LineWidth',lw); hold on
stairs(t_NL,U_NL,'--','LineWidth',lw);
xlabel('Time [s]')
ylabel('$u$','Interpreter','latex')
grid on
set(gca,'fontsize',ls)
xlim([0 8])
axu=gca;

% X_NL = zeros(size(X_GPK));
figure
plot(X_GPK(1,:),X_GPK(2,:),'-','LineWidth',lw); hold on
plot(X_K6(1,:),X_K6(2,:),':','LineWidth',lw); hold on
plot(X_NL(1,:),X_NL(2,:),'--','LineWidth',lw); hold on
plot(yref(1,end),yref(2,end),'-squarek','LineWidth',lw+1);hold on
plot(x0(1,:),x0(2,:),'o','LineWidth',lw+1); hold on
grid on
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
xlim([-2.2 2.2])
ylim([-2.2 2.2])
set(gca,'fontsize',ls)
legend('','','','Target','Start','fontsize',20,'location','best')
axf = gca;

%% Single trajectory data
tic
K = Sys.K; B = Sys.B; C = Sys.C;
sys = 1;
[U_GPKs,X_GPKs,cost_GPKs,ff_GPKs] = MPC_koopman(K,B,C,x0,yref,T,str,CV,sys);
et_GPK = toc;

%% Plots
% X_GPK = X_GPKs; U_GPK = U_GPKs;
% X_K6 = zeros(size(X_GPKs)); X_NL = X_K6;
lw=3; ls=30;
figure
plot(t1,X_GPK(1,:),'-','LineWidth',lw); hold on
plot(t1,X_K6(1,:),':','LineWidth',lw); hold on
plot(t_NL,X_NL(1,:),'--','LineWidth',lw); hold on
plot(t1(1:80:end),yref(1,1:80:end),'-squarek','LineWidth',lw-2);
xlabel('Time [s]')
ylabel('$x_1$','Interpreter','latex')
grid on
set(gca,'fontsize',ls)
xlim([0 8])
ax1 = gca;

figure
plot(t1,X_GPK(2,:),'-','LineWidth',lw); hold on
plot(t1,X_K6(2,:),':','LineWidth',lw); hold on
plot(t_NL,X_NL(2,:),'--','LineWidth',lw); hold on
plot(t1(1:80:end),yref(1,1:80:end),'-squarek','LineWidth',lw-2);
xlabel('Time [s]')
ylabel('$x_2$','Interpreter','latex')
grid on
set(gca,'fontsize',ls)
xlim([0 8])
ax2 = gca;

figure
stairs(t1(1:end-1),U_GPK,'-','LineWidth',lw); hold on
stairs(t1(1:end-1),U_K6,':','LineWidth',lw); hold on
stairs(t_NL,U_NL,'--','LineWidth',lw);
xlabel('Time [s]')
ylabel('$u$','Interpreter','latex')
grid on
set(gca,'fontsize',ls)
xlim([0 8])
ylim([-3 3])
axu=gca;

figure
plot(X_GPK(1,:),X_GPK(2,:),'-','LineWidth',lw); hold on
plot(X_K6(1,:),X_K6(2,:),':','LineWidth',lw); hold on
plot(X_NL(1,:),X_NL(2,:),'--','LineWidth',lw); hold on
plot(yref(1,end),yref(2,end),'squarek','LineWidth',lw-1);hold on
plot(x0(1,:),x0(2,:),'o','LineWidth',lw+1); hold on
grid on
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
xlim([-2.2 2.2])
ylim([-2.2 2.2])
set(gca,'fontsize',ls)
legend('','','','Target','Start','fontsize',20,'location','best')
axf = gca;

% figure
% plot(X_GPKs(1,:),X_GPKs(2,:),'k','LineWidth',lw); hold on
% plot(yref(1,end),yref(2,end),'squarek','LineWidth',lw);
% plot(x0(1,:),x0(2,:),'og','LineWidth',lw); hold on
% grid on
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% % xlim([-2.2 2.2])
% % ylim([-2.2 2.2])
% set(gca,'fontsize',ls)
% legend('','Start','Target','fontsize',20,'location','best')
% axf = gca;
%%
exportgraphics(ax1,'Fig6a.eps')
exportgraphics(ax2,'Fig6b.eps')
exportgraphics(axu,'Fig6c.eps')
exportgraphics(axf,'Fig6d.eps')

%%
function dx = DuffingStateFcn(x, u)
        % State equations of the Duffing oscillator.
        % States:
        %   x(1)  x inertial coordinate of center of mass
        %   x(2)  dx velocity
        % Inputs:
        %   u(1)
        % Parameters
        alpha = 1; beta = -4; delta = 0.02;
        % State equations
        dx(1,:) = x(2,:);
        dx(2,:) = -beta*x(1,:) -delta*x(2,:)- alpha*x(1,:).^3 + u;
end