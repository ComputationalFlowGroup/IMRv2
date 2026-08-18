%% basic simulation input
clear all
clc
close all

addpath src\common\
addpath ../cmap/
% % addpath ../Anisotropic_material_IMR/IMRv2/
% %%
% % load data and process
% load("../data/high_stretch/R1R2_data.mat")
% R1 = table2array(R1_Lamb(:,2));
% t1 = table2array(R1_Lamb(:,1));
% R2 = table2array(R2_Lamb(:,2));
% t2 = table2array(R2_Lamb(:,1));
% 
% 
% 
% % load("../data/stiff/R1_and_R_2.mat")
% % R1 = table2array(R1_anastas(:,2));
% % t1 = table2array(R1_anastas(:,1));
% % R2 = table2array(R2_anastas(:,2));
% % t2 = table2array(R2_anastas(:,1));
% 
% % load("../data/R1_and_R2_soft.mat")
% % R1 = table2array(R1_soft(:,2));
% % t1 = table2array(R1_soft(:,1));
% % R2 = table2array(R2_soft(:,2));
% % t2 = table2array(R2_soft(:,1));
% 
% % Remove duplicate time points, keeping the first occurrence
% [t1_unique, idx1] = unique(t1, 'stable');
% R1_unique = R1(idx1);
% 
% [t2_unique, idx2] = unique(t2, 'stable');
% R2_unique = R2(idx2);
% 
% % Optional: overwrite originals
% t1 = t1_unique;
% R1 = R1_unique;
% 
% t2 = t2_unique;
% R2 = R2_unique;
% 
% % tshare = (t1+t2)./2;
% 
% tshare = t2;
% 
% R1interp = interp1(t1, R1, tshare);
% R2interp = interp1(t2, R2, tshare);
% 
% theta = [0, pi/2]; Y20 = sqrt(5/(16*pi))*(3*cos(theta).^2 - 1);
% 
% M = [1 Y20(1); 1 Y20(2)];
% for i = 1:length(R1interp)
%     b = [R1interp(i); R2interp(i)];
%     x = M \ b;
%     Rbar(i) = x(1); ep2(i) = x(2)./Rbar(i);
% end
% 
% Rmax = max(Rbar).*1e-6; 
% 
% R1interp = R1interp./max(Rbar);
% R2interp = R2interp./max(Rbar);
% Rbar = Rbar./max(Rbar);
% 
% tc =  Rmax*sqrt(1048/101325);
% tshare = tshare.*1e-6./tc;
% 
% figure
% % plot(tshare, R1interp, '^--')
% hold on
% % plot(tshare, R2interp, '^--')
% plot(tshare, Rbar, 'o')
% plot(tshare, ep2, 'o')
 

% load("../data/Sims_Brown_Surya/FEM_Equil_analysis.mat")

% Req = amp_extractf(1,:).*1e-6;
% epnmeq =  0.*amp_extractf(3:end,:).*1e-6./Req;
% epnmeq =  0.*amp_extractf(3:end,:).*1e-6;



load("../data/Sims_Brown_Surya/aniso_sim_FEA_sphequil.mat")
% load("../data/Sims_Brown_Surya/iso_sim_FEA_slight.mat")
% load("../data/Sims_Brown_Surya/aniso_sim_FEA_new_props.mat")
% load("../data/Sims_Brown_Surya/aniso_sim_FEA_reanalyzed.mat")
Req =  amp_extractf(1,end).*1e-6;
Rexp = amp_extractf(1,:).*1e-6;
epnmeq =  amp_extractf(3:end,end).*1e-6./Req;

Rmax = Rexp(1);
k = 0;
idxs = [3 5];%size(amp_extractf,1);
for i = idxs
    k = k+1;
    amp(k,:) = amp_extractf(i,:)./amp_extractf(1,:);
end
texp = tStep;
[maxamp, maxampidx] = max(amp(:,1))

%%
addpath src/common/
tic
% -------- Radial Solver ----------------------------------------%
% Rmax = 150e-6;
% Rmax = 50e-6;
% Req = Rmax;
% Rmax = 100e-6;
% Req = Rmax;
% mu =  0.2625;
% G = 200e3;
% alph = 0.0;
% ani = [2.5 0];
Rmax = 1.1*Req;


mu =  0.001;
G = 105e3;
alph = 0.0;
ani = [5 0];

% mu =  0.2625;
% G = 105e3;
% alph = 0.0;
% ani = [5 0];

sig = 0.00;
p_a = -1.15*101325; f_a = 50e3;
rho = 1000;
p8 = 101325;
tcLIC = Rmax*sqrt(rho/p8);
pertmod = 0;
tf_nd = 10;
tsteps = 3000; ultra = false;

% Optional rotational full-model comparison from ../../IMR_nonspherical_dynamics.
% This is only run when ani = [0 0]. Set Enabled=false to skip it.
fullModel = struct();
fullModel.Enabled = true;
fullModel.Root = ""; % Empty uses auto-detection for ../../IMR_nonspherical_dynamics.
fullModel.xN = 256;
fullModel.L = 5;
fullModel.MaxSteps = 1000;
fullModel.TimeSteppingMethod = 2;
fullModel.ForcedEp = 'F';
fullModel.Model = "me";
fullModel.Verbose = true;
fullModel.RunOnlyWhenIsotropic = true;
fullModel.IsotropicTolerance = 100*eps;


% -------- perturbation solver initial conditions ---------------%
% Mode numbers
n = mode_extractf(idxs,1)';
m = zeros(size(n));
N = n;
ep0 = amp(:,1);
ep0(1) = 0.1;
epd0 = zeros(size(ep0));
epeq = epnmeq(idxs-2);
% ep0 = epeq.*0;
% Rmax = Req.*1.5;



t = linspace(0, tf_nd, tsteps);
[t, R, epnm] = f_call_IMRv2(Rmax, Req, ep0, epd0, epeq, n, m, ...
    mu, G, alph, ani, sig, p_a, f_a, tf_nd, tsteps, ultra, ...
    'pertmod', pertmod);

hasDistinctIsotropicModel = any(ani ~= 0);
if hasDistinctIsotropicModel
    [tiso, Riso, epnmiso] = f_call_IMRv2(Rmax, Req, ep0, epd0, ...
        epeq, n, m, mu, G, alph, [0 0], sig, p_a, f_a, ...
        tf_nd, tsteps, ultra, 'pertmod', pertmod);
else
    tiso = [];
    Riso = [];
    epnmiso = [];
end

fullModelSolution = runFullNonsphericalModelForBasicSimulation( ...
    fullModel, Rmax, Req, ep0, epd0, epeq, n, mu, G, alph, sig, ...
    p_a, f_a, rho, p8, tf_nd, tsteps, ultra, ani);



% Rsiminterp = interp1(t, R, tshare(tshare < max(t)), 'linear');
% epnmsiminterp = interp1(t, epnm, tshare(tshare < max(t)), 'linear');
% 
% Rbar(isnan(Rbar)) = 0;
% ep2(isnan(ep2)) = 0;
% 
% rmseR = rmse(Rsiminterp,Rbar(tshare < max(t))')
% rmseep = rmse(epnmsiminterp,ep2(tshare < max(t))')

%%
% load("../data/Sims_Brown_Surya/aniso_model_FEA_rad_aniso.mat")
ms.AxesFontSize = 14;
ms.LabelFontSize = 16;
ms.LegendFontSize = 11;
ms.LineWidth = 2.0;
ms.LineWidthAlt = 1.4;

nModes = min(length(n), 10);%size(epnm, 2);
nCols = 3;
nModeRows = ceil(nModes / nCols);
cmap = viridis(nModes + 2);
cmap = cmap(2:end-1, :);
tExpNd = texp ./ tcLIC;
if fullModelSolution.success
    xLimits = [0, max([t(:); tExpNd(:); fullModelSolution.t(:)])];
else
    xLimits = [0, max([t(:); tExpNd(:)])];
end

if hasDistinctIsotropicModel
    modelLabel = 'IMRv2 anisotropic model';
else
    modelLabel = 'IMRv2 isotropic model';
end

figComparison = figure('Name', 'IMR and FEM mode comparison', ...
    'Color', 'w', 'Units', 'pixels', 'Position', [100 40 950 1450]);
tl = tiledlayout(figComparison, nModeRows + 1, nCols, ...
    'TileSpacing', 'compact', 'Padding', 'loose');
comparisonAxes = gobjects(nModes + 1, 1);

axR = nexttile(tl, [1 nCols]);
comparisonAxes(1) = axR;
hold(axR, 'on')
box(axR, 'on')
grid(axR, 'on')
radialColor = [0.10 0.42 0.58];
hModel = plot(axR, t, R, '-', 'Color', radialColor, ...
    'LineWidth', ms.LineWidth);
if hasDistinctIsotropicModel
    hIsotropic = plot(axR, tiso, Riso, '--', 'Color', radialColor, ...
        'LineWidth', ms.LineWidthAlt);
else
    hIsotropic = gobjects(0);
end
if fullModelSolution.success
    hFull = plot(axR, fullModelSolution.t, fullModelSolution.R, ':', ...
        'Color', [0.08 0.08 0.08], 'LineWidth', ms.LineWidth);
else
    hFull = gobjects(0);
end
hData = scatter(axR, tExpNd, Rexp ./ Rexp(1), 32, ...
    'MarkerFaceColor', [0.45 0.25 0.55], ...
    'MarkerEdgeColor', [0.30 0.15 0.40], ...
    'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55);
xlim(axR, xLimits)
radialValues = [R(:); Riso(:); fullModelSolution.R(:); ...
    Rexp(:) ./ Rexp(1)];
radialValues = radialValues(isfinite(radialValues));
radialPadding = 0.05 * max(eps, ...
    max(radialValues) - min(radialValues));
ylim(axR, [max(0, min(radialValues) - radialPadding), ...
    max(radialValues) + radialPadding])
xlabel(axR, '$t^*$', 'Interpreter', 'latex', ...
    'FontSize', ms.LabelFontSize)
ylabel(axR, '$R/R_{\max}$', 'Interpreter', 'latex', ...
    'FontSize', ms.LabelFontSize)
radialLegendHandles = hModel;
radialLegendLabels = {modelLabel};
if hasDistinctIsotropicModel
    radialLegendHandles = [radialLegendHandles hIsotropic];
    radialLegendLabels{end + 1} = 'IMRv2 isotropic model';
end
if fullModelSolution.success
    radialLegendHandles = [radialLegendHandles hFull];
    radialLegendLabels{end + 1} = 'Full nonspherical model';
end
radialLegendHandles = [radialLegendHandles hData];
radialLegendLabels{end + 1} = 'FEM data';
legend(axR, radialLegendHandles, radialLegendLabels, ...
    'Interpreter', 'latex', 'FontSize', ms.LegendFontSize, ...
    'Location', 'best', 'NumColumns', ...
    min(4, numel(radialLegendLabels)))

for i = 1:nModes
    ax = nexttile(tl);
    comparisonAxes(i + 1) = ax;
    hold(ax, 'on')
    box(ax, 'on')
    grid(ax, 'on')
    col = cmap(i, :);

    scatter(ax, tExpNd, amp(i, :), 32, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 0.65 .* col, ...
        'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55);
    hModeData = ax.Children(1);
    if hasDistinctIsotropicModel
        hModeIso = plot(ax, tiso, epnmiso(:, i), '--', 'Color', col, ...
            'LineWidth', ms.LineWidthAlt);
    else
        hModeIso = gobjects(0);
    end
    hModeModel = plot(ax, t, epnm(:, i), '-', ...
        'Color', 0.55 .* col, 'LineWidth', ms.LineWidth);
    [fullModeValues, hModeFull] = plotFullModelModeIfAvailable(ax, ...
        fullModelSolution, n(i), ms);

    if hasDistinctIsotropicModel
        modeValues = [amp(i, :).'; epnm(:, i); ...
            epnmiso(:, i); fullModeValues(:)];
    else
        modeValues = [amp(i, :).'; epnm(:, i); ...
            fullModeValues(:)];
    end
    modeValues = modeValues(isfinite(modeValues));
    if isempty(modeValues)
        modeLimits = [-1 1];
    else
        modeMin = min(modeValues);
        modeMax = max(modeValues);
        modeSpan = modeMax - modeMin;
        modeScale = max(abs(modeValues));
        if modeSpan <= 100 * eps(max(modeScale, eps))
            modePadding = 0.05 * modeScale;
            if modePadding == 0
                modePadding = 1e-12;
            end
        else
            modePadding = 0.06 * modeSpan;
        end
        modeLimits = [modeMin - modePadding, modeMax + modePadding];
    end
    xlim(ax, xLimits)
    ylim(ax, modeLimits)
    xlabel(ax, '$t^*$', 'Interpreter', 'latex', ...
        'FontSize', ms.LabelFontSize)
    ylabel(ax, sprintf('$\\epsilon_{%.0f}$', n(i)), ...
        'Interpreter', 'latex', 'FontSize', ms.LabelFontSize)
    if i == 1
        modeLegendHandles = [hModeModel hModeData];
        modeLegendLabels = {modelLabel, 'FEM data'};
        if hasDistinctIsotropicModel
            modeLegendHandles = [hModeModel hModeIso hModeData];
            modeLegendLabels = {modelLabel, 'IMRv2 isotropic model', ...
                'FEM data'};
        end
        if ~isempty(hModeFull) && isgraphics(hModeFull)
            modeLegendHandles = [modeLegendHandles hModeFull];
            modeLegendLabels{end + 1} = 'Full nonspherical model';
        end
        legend(ax, modeLegendHandles, modeLegendLabels, ...
            'Interpreter', 'latex', 'FontSize', ms.LegendFontSize, ...
            'Location', 'best')
    end
end

for i = 1:numel(comparisonAxes)
    ax = comparisonAxes(i);
    ax.FontSize = ms.AxesFontSize - (i > 1);
    ax.TickLabelInterpreter = 'latex';
    ax.Box = 'on';
    ax.Layer = 'top';
    ax.GridAlpha = 0.18;
end

% xlim([0 1])

% Stop here for the main basic-simulation workflow. The sections below are
% older exploratory plotting snippets that load different datasets.
return

%%
clear all
close all
addpath ../cmap/

cmap = viridis(5);

clc
load("../data/soft/data_and_sim.mat")
figure
plot(tshare, Rbar, 'o', 'MarkerFaceColor',cmap(1,:), 'Color',cmap(1,:))
hold on
plot(tshare, ep2, '^', 'MarkerFaceColor',cmap(2,:), 'Color',cmap(3,:))
plot(t,R, '-','LineWidth',2, 'Color',cmap(3,:))
plot(t, epnm(:,1), '-', 'LineWidth',2, 'Color',cmap(4,:))


% load("../data/soft/fullmodel_noaniso.mat")
% plot(t./max(R), R./max(R), '--')
% hold on
% plot(t./max(R), ep, 'r--')
% hold on
% plot(t./max(R), epirr(1, 1:length(t)), '-.')
ylim([-0.25 1])
xlim([0 1])
ax = gca;

ylabel('Interface Quantity', 'Interpreter','latex', 'FontSize',16)
xlabel('Time', 'Interpreter','latex', 'FontSize',16)
ax.TickLabelInterpreter = 'latex';
legend('$\overline{R}$ Tzoumaka 2022', '$\epsilon_2$ Tzoumaka 2022', ...
    '$\overline{R}$ Current model', '$\epsilon_2$ Current model', 'Interpreter', ...
    'latex', 'FontSize', 14)
ax = gca;
ax.FontSize = 14;
ax.LabelFontSizeMultiplier = 1.5;

%%
clear all
clc
close all

load("../data/high_stretch/data_and_sim.mat")

R1sim = R;
R2sim = R;
tf = max(t);

x = cos(theta);

for i = 1:length(n)
    Pnt = legendre(n(i),x);
    Rmod = Pnt(1,:);
    R1sim = R1sim + Rmod(1).*R.*epnm(:,i).*sqrt((2*n(i)+1)/(4*pi));
    R2sim = R2sim + Rmod(2).*R.*epnm(:,i).*sqrt((2*n(i)+1)/(4*pi));
end

figure(1)
plot(tshare, R1interp, 'o')
hold on
plot(tshare, R2interp, 'o')
plot(t, R1sim)
plot(t, R2sim)

R1simint = interp1(t, R1sim, tshare);
R2simint = interp1(t, R2sim, tshare);


figure(2)
semilogy(tshare,abs(R1interp-R1simint)./R1interp, 'LineWidth',2)
hold on
figure(3)
semilogy(tshare,abs(R2interp-R2simint)./R2interp, 'LineWidth',2)
hold on

avg_relerr_pa1 = mean((abs(R1interp(tshare < tf)-R1simint(tshare < tf))./R1interp(tshare < tf)));
avg_relerr_pa2 = mean((abs(R2interp(tshare < tf)-R2simint(tshare < tf))./R2interp(tshare < tf)));

0.5*(avg_relerr_pa1+avg_relerr_pa2)


load("../data/high_stretch/fullmodel_noaniso.mat")
t = t./max(R);
R = R./max(R);


R1sim = R;
R2sim = R;

x = cos(theta);

for i = 1:length(n)
    Pnt = legendre(n(i),x);
    Rmod = Pnt(1,:);
    R1sim = R1sim + Rmod(1).*R.*ep(:).*sqrt((2*n(i)+1)/(4*pi));
    R2sim = R2sim + Rmod(2).*R.*ep(:).*sqrt((2*n(i)+1)/(4*pi));
end
figure(1)
plot(t, R1sim, '--')
hold on
plot(t, R2sim, '--')
ylim([-1 1.5])


R1simint = interp1(t, R1sim, tshare);
R2simint = interp1(t, R2sim, tshare);



figure(2)
semilogy(tshare,abs(R1interp-R1simint)./R1interp, '--')
hold on
figure(3)
semilogy(tshare,abs(R2interp-R2simint)./R2interp,  '--')
hold on

avg_relerr_f1 = mean((abs(R1interp(tshare < tf)-R1simint(tshare < tf))./R1interp(tshare < tf)));
avg_relerr_f2 = mean((abs(R2interp(tshare < tf)-R2simint(tshare < tf))./R2interp(tshare < tf)));

0.5*(avg_relerr_f1+avg_relerr_f2)

load("../data/high_stretch/fullmodel_noaniso.mat")
t = t./max(R);
R = R./max(R);

R1sim = R;
R2sim = R;

x = cos(theta);

for i = 1:length(n)
    Pnt = legendre(n(i),x);
    Rmod = Pnt(1,:);
    R1sim = R1sim + Rmod(1).*R.*epirr(:,i).*sqrt((2*n(i)+1)/(4*pi));
    R2sim = R2sim + Rmod(2).*R.*epirr(:,i).*sqrt((2*n(i)+1)/(4*pi));
end
figure(1)
plot(t, R1sim, '-.')
plot(t, R2sim, '-.')
ylim([-1 1.5])


R1simint = interp1(t, R1sim, tshare);
R2simint = interp1(t, R2sim, tshare);



figure(2)
semilogy(tshare,abs(R1interp-R1simint)./R1interp, 'k-.')
hold on
figure(3)
semilogy(tshare,abs(R2interp-R2simint)./R2interp, 'k-.')
hold on


avg_relerr_p1 = mean((abs(R1interp(tshare < tf)-R1simint(tshare < tf))./R1interp(tshare < tf)));
avg_relerr_p2 = mean((abs(R2interp(tshare < tf)-R2simint(tshare < tf))./R2interp(tshare < tf)));

0.5*(avg_relerr_p1+avg_relerr_p2)

function sol = runFullNonsphericalModelForBasicSimulation(opts, Rmax, Req, ...
    ep0, epd0, epeq, modes, mu, G, alph, sig, p_a, f_a, rho, p8, ...
    tfNd, tsteps, ultra, ani)
sol = struct('success', false, 't', [], 'R', [], 'epnm', [], ...
    'modes', [], 'message', "");

if ~opts.Enabled
    sol.message = "disabled";
    return
end

tol = opts.IsotropicTolerance;
if opts.RunOnlyWhenIsotropic && any(abs(ani(:)) > tol)
    sol.message = "skipped because ani is not [0 0]";
    return
end

fullModelRoot = resolveFullModelRoot(opts.Root);
if strlength(string(fullModelRoot)) == 0
    sol.message = "IMR_nonspherical_dynamics root was not found";
    warning('s_basic_simulation:FullModelRootMissing', '%s', sol.message);
    return
end

commonDir = fullfile(fullModelRoot, 'common');
radialSolverDir = resolveFullModelRadialSolverDir(fullModelRoot);
if strlength(string(radialSolverDir)) == 0
    sol.message = "sibling Code/IMRv2/src/forward_solver was not found";
    warning('s_basic_simulation:FullModelRadialMissing', '%s', sol.message);
    return
end

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(radialSolverDir, '-begin')
addpath(commonDir, '-begin')
clear f_imr_fd f_call_params f_odesolve compute_rotational_perturbation_evolution

try
    tstepsFull = tsteps;
    if isfield(opts, 'MaxSteps') && isfinite(opts.MaxSteps) && opts.MaxSteps > 0
        tstepsFull = min(tsteps, opts.MaxSteps);
    end

    [tRadial, RRadial, RdRadial, RddRadial] = runFullModelRadialHistory( ...
        radialSolverDir, Rmax, Req, mu, G, alph, sig, p_a, f_a, ...
        rho, p8, tfNd, tstepsFull, ultra);

    Lmax = Rmax/Req;
    tReq = tRadial(:).' .* Lmax;
    RReq = RRadial(:).' .* Lmax;
    RdReq = RdRadial(:).';
    RddReq = RddRadial(:).' ./ Lmax;

    nMode = numel(modes);
    T0 = zeros(nMode, opts.xN);
    Td0 = T0;

    Lc = Req;
    rhoc = rho;
    tc = sqrt(rhoc/p8)*Lc;
    Uc = Lc/tc;
    pc = rhoc*Uc^2;
    Ca = pc/G;
    Re = Lc*sqrt(rhoc*pc)/mu;
    We = pc*Lc/(2*sig);

    [epFull, ~, ~, ~, RReqOut, ~, tReqOut] = ...
        compute_rotational_perturbation_evolution(opts.xN, opts.L, ...
        modes, ep0, epd0, epeq, T0, Td0, 1, RReq, RdReq, RddReq, ...
        Ca, alph, Re, We, tReq, opts.TimeSteppingMethod, ...
        opts.ForcedEp, opts.Model, "rot", 'Verbose', opts.Verbose);

    sol.success = true;
    sol.t = tReqOut(:)./Lmax;
    sol.R = RReqOut(:)./Lmax;
    sol.epnm = epFull.';
    sol.modes = modes(:).';
    sol.message = "ok";
    fprintf('Full nonspherical model completed with %d time samples.\n', ...
        numel(sol.t));
catch ME
    sol.success = false;
    sol.message = string(ME.message);
    warning('s_basic_simulation:FullModelFailed', ...
        'Full nonspherical model failed: %s', ME.message);
end

clear pathCleanup
end

function [t, R, Rd, Rdd] = runFullModelRadialHistory(radialSolverDir, ...
    Rmax, Req, mu, G, alph, sig, p_a, f_a, rho, p8, tfNd, tsteps, ultra)
startDir = pwd;
dirCleanup = onCleanup(@() cd(startDir));
cd(radialSolverDir)
clear f_imr_fd f_call_params f_odesolve

kappa = 1;
T8 = 298.15;
radial = 2;
vapor = 0;
collapse = 0;
bubtherm = 0;
medtherm = 0;
masstrans = 0;
stress = 2;

if ultra
    pa = p_a;
    omega = 2*pi*f_a;
    wavetype = 4;
else
    pa = 0;
    omega = 0;
    wavetype = 0;
end

tcLIC = Rmax*sqrt(rho/p8);
tvector = linspace(0, tfNd*tcLIC, tsteps);
varin = {'progdisplay', 0, 'radial', radial, 'bubtherm', bubtherm, ...
    'tvector', tvector, 'vapor', vapor, 'medtherm', medtherm, ...
    'masstrans', masstrans, 'method', 23, 'stress', stress, ...
    'collapse', collapse, 'mu', mu, 'g', G, 'lambda1', 0e-7, ...
    'lambda2', 0, 'alphax', alph, 'surft', sig, 'r0', Rmax, ...
    'req', Req, 'kappa', kappa, 't8', T8, 'rho8', rho, ...
    'p8', p8, 'pa', pa, 'omega', omega, 'wave_type', wavetype};

[t, R, Rd, ~, ~, ~, ~, Rdd] = f_imr_fd(varin{:}, 'Nt', 75);
clear dirCleanup
end

function fullModelRoot = resolveFullModelRoot(rootIn)
scriptDir = "";
stack = dbstack('-completenames');
for kk = 1:numel(stack)
    [candidateDir, candidateName] = fileparts(stack(kk).file);
    if strcmp(candidateName, 's_basic_simulation')
        scriptDir = string(candidateDir);
        break
    end
end
if strlength(scriptDir) == 0
    scriptPath = which('s_basic_simulation');
    if strlength(string(scriptPath)) > 0
        scriptDir = string(fileparts(scriptPath));
    else
        scriptDir = string(fileparts(mfilename('fullpath')));
    end
end
if strlength(scriptDir) == 0
    scriptDir = string(pwd);
end

currentDir = string(pwd);
repoDir = string(fileparts(char(scriptDir)));
candidates = strings(0, 1);
if strlength(string(rootIn)) > 0
    candidates(end + 1) = string(rootIn);
end
candidates(end + 1) = string(fullfile(char(scriptDir), '..', '..', ...
    'IMR_nonspherical_dynamics'));
candidates(end + 1) = string(fullfile(char(scriptDir), '..', ...
    'IMR_nonspherical_dynamics'));
candidates(end + 1) = string(fullfile(char(repoDir), '..', ...
    'IMR_nonspherical_dynamics'));
candidates(end + 1) = string(fullfile(char(currentDir), '..', '..', ...
    'IMR_nonspherical_dynamics'));
candidates(end + 1) = string(fullfile(char(currentDir), '..', ...
    'IMR_nonspherical_dynamics'));
candidates(end + 1) = string(fullfile(fileparts(char(currentDir)), ...
    'IMR_nonspherical_dynamics'));
candidates(end + 1) = string(fullfile(fileparts(fileparts(char(currentDir))), ...
    'IMR_nonspherical_dynamics'));

fullModelRoot = "";
for ii = 1:numel(candidates)
    candidate = char(candidates(ii));
    solverFile = fullfile(candidate, 'common', ...
        'compute_rotational_perturbation_evolution.m');
    if exist(candidate, 'dir') == 7 && exist(solverFile, 'file') == 2
        fullModelRoot = candidate;
        return
    end
end
end

function radialSolverDir = resolveFullModelRadialSolverDir(fullModelRoot)
codeDir = fileparts(fullModelRoot);
candidate = fullfile(codeDir, 'IMRv2', 'src', 'forward_solver');
if exist(fullfile(candidate, 'f_imr_fd.m'), 'file') == 2
    radialSolverDir = candidate;
else
    radialSolverDir = "";
end
end

function [values, hLine] = plotFullModelModeIfAvailable(ax, fullSol, ...
    modeNumber, ms)
values = [];
hLine = gobjects(0);
if ~fullSol.success
    return
end

modeIdx = find(fullSol.modes == modeNumber, 1, 'first');
if isempty(modeIdx)
    return
end

values = fullSol.epnm(:, modeIdx);
hLine = plot(ax, fullSol.t, values, ':', 'Color', [0.08 0.08 0.08], ...
    'LineWidth', ms.LineWidthAlt);
end



