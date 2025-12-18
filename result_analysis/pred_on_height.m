clc;
clear;

clay_elev = load("../datasets/clay_100hz_elev_experiment_1524.mat",'valData','b');

i = randi(50)  % 76
experiment_clay = getexp(clay_elev.valData,i);

kk = 123;
folder = '../scripts/koopman_training/results/clay_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
model_clay = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","base_name","t_hc","b");   % load the file
refresh = [25,50,75,100,125,150,175,200,225,250];
k=10; tstep = 0+(1:2001); lw=4;
[ypred,yout,~,ypred_95]  = K_RSSID_prediction(experiment_clay,...
    model_clay.MDL_fitr,model_clay.A,model_clay.B,model_clay.Bc1,...
    model_clay.C,model_clay.Cc1,model_clay.K_obs,model_clay.mean_std_out,...
    refresh(k));
% terrain height
xy = [ypred(:,1:2);yout(:,1:2)];
[xq, yq] = meshgrid(linspace(min(xy(:,1)), max(xy(:,1)), 200),...
    linspace(min(xy(:,2)), max(xy(:,2)), 200));
zq = arrayfun(clay_elev.b.h, xq, yq);
% plotting
figure(2)
subplot(2,2,1)
surf(xq, yq, zq); shading interp; colormap parula; colorbar; hold on;
z_out  = interp2(xq, yq, zq, yout(:,1),  yout(:,2),  'linear');   % or 'nearest'
z_pred = interp2(xq, yq, zq, ypred(:,1), ypred(:,2), 'linear');
plot3(yout(tstep,1),yout(tstep,2),z_out,'k','linewidth',lw); hold on;
plot3(ypred(tstep,1),ypred(tstep,2),z_pred,'r','linewidth',lw);grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('H (m)');
axis([min(xy(:,1)), max(xy(:,1)), min(xy(:,2)), max(xy(:,2))])
% axis equal;
ax = gca;   % Get the current axes handle
ax.FontSize = 25; % Set the font size to 14 points
hold off;
for jj=model_clay.K_obs
    subplot(2,2,jj-2)
    plot(model_clay.t_hc(tstep),yout(tstep,jj),'k','linewidth',lw); hold on;
    plot(model_clay.t_hc(tstep),ypred(tstep,jj),'r','linewidth',lw);hold on;
    % plot(t_hc(tstep),squeeze(ypred_95(1,tstep,jj-3)),'-k','linewidth',lw);hold on;
    % plot(t_hc(tstep),squeeze(ypred_95(2,tstep,jj-3)),'-k','linewidth',lw);grid on;
    ylabel(experiment_clay.OutputName(jj),'Interpreter','latex','FontSize',15)
    xlabel('Time [s]')
    xlim([0 tstep(end)*0.01]);
    grid on; box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    
    hold off;
    if jj==4
        ylim([min(yout(tstep,jj))-2 2+max(yout(tstep,jj))])
    end
end
% sgtitle("Clay soil")
legend('True','K-SSID')

%% %%%%%%%%%%
% clc;
clear;
% SL_elev = load("../datasets/sandyloam_100hz_no_elev_experiment_1579.mat",'testData','b');
SL_elev = load("../datasets/sandyloam_100hz_elev_experiment_1572.mat",'valData','b');
i = randi(50)   %28
experiment_SL = getexp(SL_elev.valData,i);
kk = 124;
folder = '../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files = dir(fullfile(folder, '*.mat'));   % or *.txt, *.csv, etc.
filename = fullfile(folder, files(kk).name);
model_sandyloam = load(filename,"MDL_fitr","A","B","Bc1","C",...
    "Cc1","K_obs","mean_std_out","rr","base_name","t_hc");   % load the file
refresh = [25,50,75,100,125,150,175,200,225,250];
k=10; tstep = 0+(1:2001); lw=4;
allErr = []; allY = [];
[ypred,yout,~,ypred_95]  = K_RSSID_prediction(experiment_SL,...
    model_sandyloam.MDL_fitr,model_sandyloam.A,model_sandyloam.B,model_sandyloam.Bc1,...
    model_sandyloam.C,model_sandyloam.Cc1,model_sandyloam.K_obs,...
    model_sandyloam.mean_std_out,refresh(k));
% terrain height
xy = [ypred(:,1:2);yout(:,1:2)];
[xq, yq] = meshgrid(linspace(min(xy(:,1)), max(xy(:,1)), 200),...
    linspace(min(xy(:,2)), max(xy(:,2)), 200));
zq = arrayfun(SL_elev.b.h, xq, yq);
% plotting
figure(1)
subplot(2,2,1)
surf(xq, yq, zq); shading interp; colormap parula; colorbar; hold on;
z_out  = interp2(xq, yq, zq, yout(:,1),  yout(:,2),  'linear');   % or 'nearest'
z_pred = interp2(xq, yq, zq, ypred(:,1), ypred(:,2), 'linear');
plot3(yout(tstep,1),yout(tstep,2),z_out,'k','linewidth',lw); hold on;
plot3(ypred(tstep,1),ypred(tstep,2),z_pred,'r','linewidth',lw);grid on;
% legend('True','K-SSID')
xlabel('X (m)'); ylabel('Y (m)'); zlabel('H (m)');
axis([min(xy(:,1)), max(xy(:,1)), min(xy(:,2)), max(xy(:,2))])
ax = gca;   % Get the current axes handle
ax.FontSize = 25; % Set the font size to 14 points
hold off;
for jj=model_sandyloam.K_obs
    subplot(2,2,jj-2)
    plot(model_sandyloam.t_hc(tstep),yout(tstep,jj),'k','linewidth',lw); hold on;
    plot(model_sandyloam.t_hc(tstep),ypred(tstep,jj),'r','linewidth',lw);hold on;
    % plot(t_hc(tstep),squeeze(ypred_95(1,tstep,jj-3)),'-k','linewidth',lw);hold on;
    % plot(t_hc(tstep),squeeze(ypred_95(2,tstep,jj-3)),'-k','linewidth',lw);grid on;
    ylabel(experiment_SL.OutputName(jj),'Interpreter','latex','FontSize',15)
    xlabel('Time [s]')
    xlim([0 tstep(end)*0.01]);
    grid on; box on;
    set(gca, 'LineWidth', 1.5)
    ax = gca;   % Get the current axes handle
    ax.FontSize = 30; % Set the font size to 14 points
    hold off;
    if jj==4
        ylim([min(yout(tstep,jj))-2 2+max(yout(tstep,jj))])
    end
end
legend('True','K-SSID')
% sgtitle("Sandy loam soil")