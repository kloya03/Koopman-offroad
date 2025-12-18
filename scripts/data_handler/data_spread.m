
%%
clc;
clear;
load('../../datasets/sandyloam_100hz_no_elev_experiment_1579.mat','exp_fv','exp_bekk_h')


beta = []; sr = [];v=[];hf=[];
for i=1:size(exp_fv,2)
    exp_1 = exp_fv{i};
    exp_2 = exp_bekk_h{i};
    beta = [beta;exp_1(:,17);exp_1(:,22)];
    sr = [sr;exp_1(:,16);exp_1(:,21)];
    v = [v;exp_1(:,14);exp_1(:,19)];
    hf = [hf;exp_2(:,6);exp_2(:,12)];
end
% histogram(beta);
% xlabel('slip angle');
% ylabel('Frequency');
% % title('Histogram of all sinkage values');

figure(1)
subplot(2,2,1)
histogram(sr);
xlabel('slip ratio');
ylabel('Frequency');
grid on;
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 30; % Set the font size to 14 points

subplot(2,2,2)
histogram(v);
xlabel('cornering velocity (m/s)');
ylabel('Frequency');
grid on;
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 30; % Set the font size to 14 points

% histogram(hf);
% xlabel('sinkage (m)');
% ylabel('Frequency');
% % title('Histogram of all sinkage values');


%%
clc;
clear;
load('../../datasets/clay_100hz_no_elev_experiment_1472.mat','exp_fv','exp_bekk_h')

beta = []; sr = [];v=[];hf=[];
for i=1:size(exp_fv,2)
    exp_1 = exp_fv{i};
    exp_2 = exp_bekk_h{i};
    beta = [beta;exp_1(:,17);exp_1(:,22)];
    sr = [sr;exp_1(:,16);exp_1(:,21)];
    v = [v;exp_1(:,14);exp_1(:,19)];
    hf = [hf;exp_2(:,6);exp_2(:,12)];
end
% histogram(beta);
% xlabel('slip angle');
% ylabel('Frequency');
% % title('Histogram of all sinkage values');
figure(1)
subplot(2,2,3)
histogram(sr);
xlabel('slip ratio');
ylabel('Frequency');
grid on;
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 30; % Set the font size to 14 points


subplot(2,2,4)
histogram(v);
xlabel('cornering velocity (m/s)');
ylabel('Frequency');
grid on;
box on;
set(gca, 'LineWidth', 1.5)
ax = gca;   % Get the current axes handle
ax.FontSize = 30; % Set the font size to 14 points


% histogram(hf);
% xlabel('sinkage (m)');
% ylabel('Frequency');
% % title('Histogram of all sinkage values');