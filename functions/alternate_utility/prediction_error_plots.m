clc;
clear all;


% Folder containing .mat files
folder = 'sandyloam_noelev_models/models_with_error';     % <-- change to your folder

% Pattern to match only desired files
files = dir(fullfile(folder, 'task46_sandyloam_*_error.mat'));

% Preallocate
numFiles = length(files);

for k = 1:numFiles
    fname = fullfile(folder, files(k).name);

    % Load only specific variables (faster, cleaner)
    load(fname, 'overall_error','task_id','rr','nl','sy','nB','ytest','cut_off',...
        'A');
    overall(k,:) = [k, mean(overall_error,"all"), ...
        mean(overall_error(:,end),"all"), mean(abs(ytest),"all"),...
        str2num(task_id),rr,nl,sy,nB,cut_off];
    eigenvals{k} = eig(A);
end
B = overall(:,2) + overall(:,3);
aa = find(min(B)==B);
overall(:,end+1) =B;
BA = sortrows(overall, 5);   % sort by column 2 (ascending)


%%
% overall_error is a 6x10 matrix
[numRows, numCols] = size(overall_error);

figure;
hold on;
referesh = [25,50,75,100,125,150,175,200,225,250];
for i = 1:numRows
    plot(referesh, overall_error(i, :), 'LineWidth',3);
end

xlabel('refresh rate (number of time step)');
ylabel('Error');
title('Overall Error forvariable with refresh rate');
legend('x','y','$\psi$','u','v','$\dot{\psi}$','interpreter','latex');
grid on;
hold off;
%%
clc;
clf;

% jj = [1:3:27; 2:3:27; 3:3:27;...
%     28:3:54; 29:3:54; 30:3:54;...
%     55:3:81; 56:3:81; 57:3:81;...
%     82:3:108; 83:3:108; 84:3:108;
%     109:3:135; 110:3:135; 111:3:135];
jj=1

for kk=1,%randi(15)
    BD = sortrows(BA(jj(kk,:),:),10);
    figure(1)
    plot(BD(:,6),BD(:,11)); hold on;
    % axis([0 300 0.3 0.5])
    display(BD)
end

%% 

th = linspace(0, 2*pi, 400);
figure(1)


for i=46
    i,overall(:)
    eigvals = eigenvals{i};
    eig_mag = abs(eigvals);
    figure(1)
    plot(cos(th), sin(th), 'k--', 'LineWidth', 1.2);  
    axis equal;
    hold on;
    scatter(real(eigvals), imag(eigvals), 80, eig_mag, 'filled');
    axis equal;
    grid on;
    xlabel('Real Part');
    ylabel('Imaginary Part');
    title('Eigenvalues with Magnitude Colorbar');
    col = colorbar;
    ylabel(col, '|\lambda|');
    
    % Limits for nicer view
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    clf
end






