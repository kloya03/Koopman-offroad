%% Eigenvalue plot + CENTER inset zoom (cluster near 1+0i)
kk = 124;
folder   = '../../scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/';
files    = dir(fullfile(folder,'*.mat'));
filename = fullfile(folder, files(kk).name);

load(filename,"A");
eigvals = eig(A);
mag     = abs(eigvals);

% ----- MAIN PLOT -----
figure; hold on; axis equal; grid on; box on;

% Unit circle
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1.5);

% Eigenvalues (color = magnitude)
sc = scatter(real(eigvals), imag(eigvals), 80, mag, 'filled');

cb = colorbar('Ticks',[0.95 0.96 0.97 0.98 0.99 1.00]);
ylabel(cb,'Eigenvalue Magnitude');

xlabel('Real'); ylabel('Imaginary');
set(gca,'LineWidth',1.5,'FontSize',30);
xlim([-1.05 1.05]); ylim([-1.05 1.05]);

axMain = gca;  % keep handle to main axes

% ----- DEFINE ZOOM REGION: CLUSTER NEAR 1+0i -----
axMain = gca; fig = gcf;

c = 1 + 0i;          % target cluster near eigenvalue 1
r0 = 0.2;           % radius around 1 (increase if needed)
imagMax = 0.25;      % keep near real-axis cluster (set [] to disable)

dist = abs(eigvals - c);
mask = (dist <= r0);
if ~isempty(imagMax)
    mask = mask & (abs(imag(eigvals)) <= imagMax);
end

% Auto-relax if too few points selected
if nnz(mask) < 8
    r0 = 0.20;
    mask = (abs(eigvals - c) <= r0);
    if ~isempty(imagMax)
        mask = mask & (abs(imag(eigvals)) <= imagMax);
    end
end

xZ = real(eigvals(mask));
yZ = imag(eigvals(mask));

% Fallback
if isempty(xZ)
    xZ = real(eigvals); yZ = imag(eigvals);
end

% Bounds + padding
pad = 0.13;   % increase to enlarge the small rectangle
x1 = min(xZ) - pad;  x2 = max(xZ) + pad;
y1 = min(yZ) - pad;  y2 = max(yZ) + pad;

% Make it SQUARE (looks better with axis equal and avoids skinny inset)
cx = 0.5*(x1+x2);  cy = 0.5*(y1+y2);
half = 0.55*max((x2-x1),(y2-y1));   % 0.55 adds a bit more margin
x1 = cx - half*0.5;  x2 = cx + half;
y1 = cy - half;  y2 = cy + half;

% Draw zoom rectangle on MAIN axes
hRect = rectangle(axMain, 'Position', [x1 y1 (x2-x1) (y2-y1)], ...
    'EdgeColor','k','LineWidth',1.8);

% ----- INSET AXES: BIGGER and CENTERED -----
posM = axMain.Position;            % [l b w h] of main axes in normalized units
wI = 0.15 * posM(3);               % inset width relative to main axes
hI = 0.6 * posM(4);               % inset height relative to main axes
leftI   = posM(1) + posM(3)/2 - wI/2;
bottomI = posM(2) + posM(4)/2 - hI/2;

axInset = axes('Position',[leftI bottomI wI hI]);
hold(axInset,'on'); box(axInset,'on'); grid(axInset,'on'); %axis(axInset,'equal');
set(axInset,'Color','w','LineWidth',1.6,'FontSize',18);

% plot in inset
th = linspace(0,2*pi,400);
plot(axInset, cos(th), sin(th), 'k--', 'LineWidth', 1.2);
scatter(axInset, real(eigvals), imag(eigvals), 50, mag, 'filled');

xlim(axInset,[x1 x2-0.1]);
ylim(axInset,[y1 y2]);
colormap(axInset, colormap(axMain));

% ----- CONNECTOR LINES: TOUCH BOTH BOXES CLEANLY -----
% Since the zoom box is near x~1 (right side), connect from its LEFT edge
% to the RIGHT edge of the inset (looks clean and avoids weird crossings).

ptMain_TL = data2norm(axMain, [x1 y2]);   % left-top corner of zoom rect (data->norm)
ptMain_BL = data2norm(axMain, [x1 y1]);   % left-bottom corner of zoom rect

posI = axInset.Position;                  % inset in normalized figure coords
ptInset_TR = [posI(1)+posI(3), posI(2)+posI(4)];  % inset top-right
ptInset_BR = [posI(1)+posI(3), posI(2)];          % inset bottom-right

annotation(fig,'line', [ptMain_TL(1)-0.095 ptInset_TR(1)], [ptMain_TL(2)+0.005 ptInset_TR(2)], ...
    'Color','k','LineWidth',1.4);
annotation(fig,'line', [ptMain_BL(1)-0.095 ptInset_BR(1)], [ptMain_BL(2)+0.01 ptInset_BR(2)], ...
    'Color','k','LineWidth',1.4);

% helper: data coords -> normalized figure coords
function p = data2norm(ax, xy)
    axPos = ax.Position; xl = ax.XLim; yl = ax.YLim;
    xn = (xy(1)-xl(1))/(xl(2)-xl(1));
    yn = (xy(2)-yl(1))/(yl(2)-yl(1));
    p = [axPos(1)+xn*axPos(3), axPos(2)+yn*axPos(4)];
end
