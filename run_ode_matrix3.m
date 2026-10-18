%% 
% clear all
clc
load('simulink_rotx_1e3.mat');

%% femodel generation
stlFile = 'cad/blanket_obr2.STL';
E = 2e11*1e-6;          % Young's modulus in MPa
nu = 0.3;               % Poisson's ratio (nondimensional)
rho = 3625*1e-12;      % Mass density in t/mm^3
origins = [0 0 1567/2];
numFrames = size(origins,1);

femodel = createpde('structural','modal-solid');
importGeometry(femodel,stlFile);
structuralProperties(femodel, ...
    'YoungsModulus',E, ...
    'PoissonsRatio',nu, ...
    'MassDensity',rho);
generateMesh(femodel, 'GeometricOrder', 'quadratic', 'Hmax', 1e3, 'Hmin', 8e2, 'Hgrad', 2);


faceIDs = [10];    % List in the same order as the interface frame origins

structuralBC(femodel, ...
    'Face',faceIDs(1), ...
    'Constraint','multipoint', ...
    'Reference',origins(1,:));

rom = reduce(femodel,'FrequencyRange',[0 8]*(2*pi));
TransformationMatrix=...
    rom.MPCTransformationMatrix*rom.TransformationMatrix;

volumes = zeros(1,size(femodel.Mesh.Elements,2));
for i= 1:size(femodel.Mesh.Elements,2)
    volumes(i) = volume(femodel.Mesh,i)/volume(femodel.Mesh);
end

weights= zeros(1,size(femodel.Mesh.Nodes,2));
for i= 1:size(femodel.Mesh.Nodes,2)
    weights(i)=sum(volumes(sum(femodel.Mesh.Elements==i)==1))/10;
end
m_i = volume(femodel.Mesh)*femodel.MaterialProperties.MaterialAssignments.MassDensity*1e3;



AttachmentDoF = rom.RetainedDoF;
CondensedDOF = true(3*size(femodel.Mesh.Nodes,2), 1);
CondensedDOF(AttachmentDoF) = false;
CondensedDoF = find(CondensedDOF);
DOF = [AttachmentDoF, CondensedDoF'];

Psi=TransformationMatrix(:,end);
% --- MAIN SCRIPT ---
z = [zeros(6,1); 0; zeros(6,1); 0];
% u_input = zeros(6,1);
u_input = @(t) [0; 0; 0; 1e3; 0; 0];

% Setup ODE options
% 1. Pass the Mass matrix function (it must accept t and Z as arguments)
% 2. Declare 'strong' state dependence so the solver computes the proper Jacobians
% opts = odeset('OutputFcn', @printProgress, ...
%     'Stats', 'on', ...
%     'Mass', @(t, Z) assembleMassMatrix(femodel.Mesh.Nodes, origins, Psi, DOF, m_i, weights, Z), ...
%     'MStateDependence', 'strong');


% 1. Define the Jacobian Pattern (14x14)
% Z = [translations (3); rotations (3); modal (1); trans_vel (3); rot_vel (3); modal_vel (1)]
S = sparse(14, 14);
% Top 7 states (kinematics) depend only on the bottom 7 states (velocities)
S(1:7, 8:14) = speye(7); 
% Bottom 7 states (dynamics) depend on all states (positions and velocities)
S(8:14, 1:14) = 1; 

% 2. Setup ODE options with JPattern and your MEX files (if compiled)
opts = odeset('Mass', @(t, Z) assembleMassMatrix(femodel.Mesh.Nodes, origins, Psi, DOF, m_i, weights, Z), ...
    'MStateDependence', 'strong', ...
    'JPattern', S);
% Solve the system using the updated RHS function
[time, states] = ode15s(@(t,Z) cb_dynamics_rhs(t, Z, u_input(t), rom.K, femodel.Mesh.Nodes, origins, Psi, DOF, m_i, weights), [0, 1], z, opts);
elapsedTime = toc;


% save solution.mat time states
fprintf('The script took %.4f seconds to run.\n', elapsedTime);

fig1 = figure(1);
defaultLatex(fig1);
fig1.PaperPosition = [0, 0, 16, 12]; % Adjusted size for 6 subplots
tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'tight');

% 1. Translation X
ax1 = nexttile;
plot(time, states(:, 1), 'LineWidth', 1.5, 'Color', [0, 0.5, 0]); hold on;
plot(out.rigid.Time, squeeze(out.rigid.Data(1,:,:)), 'LineWidth', 1, 'Color', 'k','LineStyle', '--');
title('Rigid Translation X'); ylabel('Meters'); 
legend('Proposed method', 'Simulink', 'Location', 'best'); grid on; ylim([-1,1]);

% 2. Translation Y
ax2 = nexttile;
plot(time, states(:, 2), 'LineWidth', 1.5, 'Color', [0, 0.5, 0]); hold on;
plot(out.rigid.Time, squeeze(out.rigid.Data(2,:,:)), 'LineWidth', 1, 'Color', 'k','LineStyle', '--');
title('Rigid Translation Y'); ylabel('Meters'); grid on;

% 3. Translation Z
ax3 = nexttile;
plot(time, states(:, 3), 'LineWidth', 1.5, 'Color', [0, 0.5, 0]); hold on;
plot(out.rigid.Time, squeeze(out.rigid.Data(3,:,:)), 'LineWidth', 1, 'Color', 'k','LineStyle', '--');
title('Rigid Translation Z'); ylabel('Meters'); grid on;

% 4. Rotation X
ax4 = nexttile;
plot(time, states(:, 4), 'LineWidth', 1.5, 'Color', [0, 0.5, 0]); hold on;
plot(out.rigid.Time, squeeze(out.rigid.Data(6,:,:)), 'LineWidth', 1, 'Color', 'k','LineStyle', '--');
title('Rigid Rotation X'); ylabel('Rad'); grid on; ylim([-1,1]);

% 5. Rotation Y
ax5 = nexttile;
plot(time, states(:, 5), 'LineWidth', 1.5, 'Color', [0, 0.5, 0]); hold on;
plot(out.rigid.Time, squeeze(out.rigid.Data(5,:,:)), 'LineWidth', 1, 'Color', 'k','LineStyle', '--');
title('Rigid Rotation Y'); ylabel('Rad'); grid on; ylim([-1,1]);

% 6. Rotation Z
ax6 = nexttile;
plot(time, states(:, 6), 'LineWidth', 1.5, 'Color', [0, 0.5, 0]); hold on;
plot(out.rigid.Time, squeeze(out.rigid.Data(4,:,:)), 'LineWidth', 1, 'Color', 'k','LineStyle', '--');
title('Rigid Rotation Z'); ylabel('Rad'); grid on; ylim([-1,1]);

fontSetting(fig1);
saveFigPdf(fig1, 'rigid_motions.pdf', [16, 12], '-r1000');


% ---------------------------------------------------------
% Figure 2: Modal Deflection with Inset Zoom
% ---------------------------------------------------------
fig2 = figure(2);
defaultLatex(fig2);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'tight'); 
fig2.PaperPosition = [0, 0, 8, 6];

ax7 = nexttile;

hold on;
plot(out.modal.Time, 1e-3*squeeze(out.modal.Data(1,:,:)), 'Color', [0, 0, 0], 'linewidth', .1);
plot(time, states(:, 7), 'Color', [0, 0.5, 0], 'linewidth', .1);

xlabel('Time [s]');
ylabel('Modal Deflection [m]');
xlim([time(1), time(end)]);
grid on;

lgd = legend({'Simulink ODE15s','Proposed method'}, 'Location', 'northeast');
lgd.FontSize = 6;

% --- Inset Zoom Setup ---
% Define zoom time window (Adjust tmin/tmax based on where the interesting dynamics occur)
tmin = 0.0; 
tmax = 0.05; 

% Automatically find Y limits for the zoom box
idx_zoom = (time >= tmin) & (time <= tmax);
if any(idx_zoom)
    ymin = min(states(idx_zoom, 7));
    ymax = max(states(idx_zoom, 7));
else
    yl = ylim(ax7); ymin = yl(1); ymax = yl(2);
end

% Draw dashed rectangle on main plot
rectangle(ax7, 'Position', [tmin, ymin, (tmax - tmin), (ymax - ymin)], ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.8, 'LineStyle','--');

axes('position', [.25 .2 .3 .2]); 
box on; grid on; hold on;

plot(time(idx_zoom), states(idx_zoom, 7), 'Color', [0, 0.5, 0], 'linewidth', .5); 
idx_zoom_sim = (out.modal.Time >= tmin) & (out.modal.Time <= tmax);

sim_modal_data = 1e-3 * squeeze(out.modal.Data(1,:,:));
plot(out.modal.Time(idx_zoom_sim), sim_modal_data(idx_zoom_sim), 'Color', [0, 0, 0], 'linewidth', .5);
idx_zoom_sim = (out.modal.Time >= tmin) & (out.modal.Time <= tmax);
sim_modal_data = 1e-3 * squeeze(out.modal.Data(1,:,:));


fontSetting(fig2);
saveFigPdf(fig2, 'modal_deflection_inset.pdf', [8, 6], '-r1000');

% --- LOCAL FUNCTIONS ---

function rhs = cb_dynamics_rhs(t, Z, u, K, nodes, origins, Psi, DOF, m_i, weights)

F_in = [eye(6,6); zeros(1,6)];
gbar = [zeros(7,6); F_in];

A_c = computeRotationMatrix(Z);
rhs = computeForce(K,nodes, origins, Psi, DOF, m_i, weights, Z) + gbar * u;
end

function Q_aug = computeForce(K,nodes, origin_ref, Psi, DOF, m_i, weights, Z)
    n = 7;
    Q_c= computeCoriolisForce(nodes, origin_ref, Psi, DOF, m_i, weights, Z);

    Q_st = [zeros(n-1,1);-K(end,end) * 1e3 * Z(n)];
    Q = Q_c+Q_st;
    euler_rates = Z(11:13);
    Q_aug = [Z(8:10); euler_rates; Z(14); Q];
end

function saveFigPdf(fig, filename, figsize,resolution)
if ~exist('figsize', 'var') || isempty(figsize)
    figsize = [5, 7];
end

if ~exist('resolution', 'var') || isempty(resolution)
    resolution = '-r300';
end

fig.PaperUnits = 'centimeters';
fig.PaperPositionMode = 'manual';
fig.PaperSize = figsize;
fig.PaperPosition = [0, 0, figsize(1), figsize(2)];
print(fig, '-dpdf', filename,resolution);
end

function defaultLatex(fig)
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(fig,'defaultLegendAutoUpdate','off')
end
function fontSetting(fig)
allaxes = findall(fig, 'type', 'axes');
set(allaxes,'TickLabelInterpreter','latex');
set(allaxes,'FontUnits','points');
set(allaxes,'FontName','times');
set(allaxes,'FontWeight','normal');
set(allaxes,'FontSize',8);
end