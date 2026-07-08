%% plot_bubble_step5_time_mirrored.m
% Plot time-resolved Abaqus bubble shapes from one arbitrary ODB/job name.
%
%
% Expected input file:
%
%  Job1_Step5_bubble_time.dat
%
% This script assumes the extracted data contains ONLY the true quarter
% bubble surface. It mirrors that quarter surface across x=0 and y=0 to
% reconstruct the full 2D bubble cross-section.
% 
 
clear all; clc; close all;

% ===================== USER SETTINGS =====================

% Folder containing the extracted bubble collapse time-history .dat file.
DATA_DIR = fullfile(pwd);

% Job name WITHOUT .odb.
% Example:
%   Job1.odb -> JOB_NAME = 'Job1'
% for a non-spherical bubble in an NHKV medium, set JOB_NAME='Job1'
% for a non-spherical bubble in an anisotropic medium, set JOB_NAME='Job_QSR'
% QSR = Quadratic Standard fiber Reincorcement model

% Uncomment line below for non-spherical bubble in NHKV
% JOB_NAME = 'Job1';
JOB_NAME = 'Job';

% Uncomment line below for non-spherical bubble in fiber reinofced material
% JOB_NAME = 'Job_QSR';

% Number of Step-5 (bubble collapse step) frames to plot from beginning to end.
NUM_FRAMES_TO_PLOT = 39;

% Plot full mirrored bubble or only quarter surface.
MIRROR_FULL_BUBBLE = true;

% Normalize each curve by its own max radius.
% Usually keep false if you want actual size evolution.
normalizeByMaxRadius = true;

% Show node markers.
showMarkers = true;

% Save the plotted figure.
SAVE_FIGURE = true;
IMAGE_DIR = "BUBBLE_MIRRORED_IMAGES";
IMAGE_RESOLUTION = 450;

% Export full mirrored time history to DAT.
EXPORT_FULL_MIRRORED_DAT = true;
EXPORT_DIR = "BUBBLE_FULL_MIRRORED_DAT";

% Appearance.
lineWidth = 1.8;
markerSize = 4;
COLOR_MAP_NAME = 'turbo';

% ==========================================================


% ===================== LOAD DATA =====================

% fileName = sprintf('../data/Sims_Brown_Surya/%s_bubble_time.dat', JOB_NAME);

% fileName = sprintf('../data/Sims_Brown_Surya/%s_spherical_QSR_Step5_bubble_time.dat', JOB_NAME);
fileName = sprintf('../data/Sims_Brown_Surya/%s_anisotropic_new.dat', JOB_NAME);
filePath = fullfile(DATA_DIR, fileName);

if ~isfile(filePath)
    error('Could not find file:\n%s\n\nCheck DATA_DIR and JOB_NAME.', filePath);
end

T = readStep5BubbleDat(filePath);

requiredVars = {'frame_id','step_time','node_label','x','y'};
for k = 1:numel(requiredVars)
    if ~ismember(requiredVars{k}, T.Properties.VariableNames)
        error('Missing required column "%s" in %s', requiredVars{k}, filePath);
    end
end

allFrames = unique(T.frame_id);
nFrames = numel(allFrames);

fprintf('Loaded file: %s\n', filePath);
fprintf('Number of frames: %d\n', nFrames);
fprintf('Rows: %d\n', height(T));

if nFrames == 0
    error('No frames found in file.');
end

% Select evenly spaced frames.
nPlot = nFrames;
frameIdx = round(linspace(1, nFrames, nPlot));
%frameIdx = unique(frameIdx, 'stable');
framesToPlot = allFrames(frameIdx);

fprintf('Plotting %d frames.\n', numel(framesToPlot));


% ===================== PLOT STEP-5 SHAPE EVOLUTION =====================

addpath src\common\

colors = getColorMap(numel(framesToPlot), COLOR_MAP_NAME);

for i = 1:nFrames

    frameID = framesToPlot(i);
    Tf = T(T.frame_id == frameID, :);

    % Raw quarter surface from Abaqus.
    xq = Tf.x;
    yq = Tf.y;

    % Sort by polar angle so the quarter arc is ordered.
    [xq, yq] = sortQuarterCurveByAngle(xq, yq);

    % Remove duplicate consecutive points, if any.
    [xq, yq] = uniqueXY(xq, yq, 1e-10);

    % Mirror if requested.
    if MIRROR_FULL_BUBBLE
        [xPlot, yPlot] = mirrorQuarterCurveClosed(xq, yq);
    else
        xPlot = xq;
        yPlot = yq;
    end

    % Rotate clockwise so theta = 0 corresponds to the original +y axis.
    xBeforeRotation = xPlot;
    xPlot = yPlot;
    yPlot = -xBeforeRotation;
    
    radius = sqrt(xPlot.^2 + yPlot.^2);
    theta = acos(xPlot./radius);
    simdata = [theta, radius];
    if i == 1
        figure
        plot(theta, radius, 'o')
    end
    [mode_extractf(:,i), amp_extractf(:,i), phase_extractf(:,i)] = fft_extract_axissym(simdata, 26);
    tStep(i) = Tf.step_time(1);

end

%% ===================== COMPARE SIMULATED AND EXTRACTED SURFACES =====================

surfaceDiffL2 = zeros(1, nFrames);
surfaceDiffRelativeL2 = zeros(1, nFrames);
surfaceDiffRMSE = zeros(1, nFrames);
surfaceDiffMaxAbs = zeros(1, nFrames);

thetaComparisonGrid = linspace(0, pi, 361).';
surfaceRadiusResidualGrid = nan(numel(thetaComparisonGrid), nFrames);
simSurfaceX = cell(1, nFrames);
simSurfaceY = cell(1, nFrames);
extractedSurfaceX = cell(1, nFrames);
extractedSurfaceY = cell(1, nFrames);

for i = 1:nFrames

    frameID = framesToPlot(i);
    Tf = T(T.frame_id == frameID, :);

    xq = Tf.x;
    yq = Tf.y;

    [xq, yq] = sortQuarterCurveByAngle(xq, yq);
    [xq, yq] = uniqueXY(xq, yq, 1e-10);

    if MIRROR_FULL_BUBBLE
        [xPlot, yPlot] = mirrorQuarterCurveClosed(xq, yq);
    else
        xPlot = xq;
        yPlot = yq;
    end

    % Rotate clockwise so theta = 0 corresponds to the original +y axis.
    xBeforeRotation = xPlot;
    xPlot = yPlot;
    yPlot = -xBeforeRotation;

    radiusSim = hypot(xPlot, yPlot);
    cosTheta = xPlot ./ radiusSim;
    cosTheta = min(max(cosTheta, -1), 1);
    thetaSim = acos(cosTheta);

    radiusExtracted = reconstructAxisymSurfaceRadius(thetaSim, amp_extractf(:,i));
    radiusResidual = radiusExtracted - radiusSim;

    surfaceDiffL2(i) = norm(radiusResidual, 2);
    surfaceDiffRelativeL2(i) = surfaceDiffL2(i) ./ norm(radiusSim, 2);
    surfaceDiffRMSE(i) = sqrt(mean(radiusResidual.^2));
    surfaceDiffMaxAbs(i) = max(abs(radiusResidual));

    [thetaUnique, ~, thetaGroup] = unique(thetaSim);
    radiusResidualUnique = accumarray(thetaGroup, radiusResidual, [], @mean);
    surfaceRadiusResidualGrid(:,i) = interp1(thetaUnique, radiusResidualUnique, ...
        thetaComparisonGrid, 'linear', NaN);

    simSurfaceX{i} = xPlot;
    simSurfaceY{i} = yPlot;
    extractedSurfaceX{i} = radiusExtracted .* xPlot ./ radiusSim;
    extractedSurfaceY{i} = radiusExtracted .* yPlot ./ radiusSim;

end

surfaceComparisonTable = table(framesToPlot(:), tStep(:), ...
    surfaceDiffL2(:), surfaceDiffRelativeL2(:), surfaceDiffRMSE(:), ...
    surfaceDiffMaxAbs(:), ...
    'VariableNames', {'Frame', 'StepTime', 'L2RadiusError', ...
    'RelativeL2RadiusError', 'RMSERadiusError', 'MaxAbsRadiusError'});
disp(surfaceComparisonTable)

figure('Color','w');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile
plot(tStep, surfaceDiffL2, 'o-', 'LineWidth', 1.4)
hold on
plot(tStep, surfaceDiffRMSE, 's-', 'LineWidth', 1.4)
plot(tStep, surfaceDiffMaxAbs, '^-', 'LineWidth', 1.4)
grid on
box on
xlabel('Step time')
ylabel('Radius error')
title('Extracted surface error norms')
legend({'L2', 'RMSE', 'max abs'}, 'Location', 'best')

nexttile
plot(tStep, surfaceDiffRelativeL2, 'o-', 'LineWidth', 1.4)
grid on
box on
xlabel('Step time')
ylabel('Relative L2 radius error')
title('Relative extracted surface error')

figure('Color','w');
imagesc(tStep, thetaComparisonGrid, surfaceRadiusResidualGrid)
set(gca, 'YDir', 'normal')
colorbar
box on
xlabel('Step time')
ylabel('\theta')
title('Extracted - simulated radius over time')

figure('Color','w');
hold on
axis equal
box on
grid on
xlabel('x')
ylabel('y')
title('Simulated and extracted surfaces')
comparisonFrameIdx = unique(round(linspace(1, nFrames, min(8, nFrames))), 'stable');
comparisonColors = getColorMap(numel(comparisonFrameIdx), COLOR_MAP_NAME);
for k = 1:numel(comparisonFrameIdx)
    i = comparisonFrameIdx(k);
    plot(simSurfaceX{i}, simSurfaceY{i}, 'o', ...
        'Color', comparisonColors(k,:), 'MarkerSize', markerSize, ...
        'LineStyle', 'none', ...
        'DisplayName', sprintf('sim frame %d', framesToPlot(i)));
    plot(extractedSurfaceX{i}, extractedSurfaceY{i}, '-', ...
        'Color', comparisonColors(k,:), 'LineWidth', 1.6, ...
        'DisplayName', sprintf('extract frame %d', framesToPlot(i)));
end
legend('Location', 'bestoutside', 'Interpreter', 'none')

%%
% load("../data/Sims_Brown_Surya/aniso_sim_FEA_sphequil.mat")

figure
nmodes = size(amp_extractf,1);

for i = 1:nmodes
    plotl = ceil(sqrt(nmodes));
    subplot(plotl, plotl, i)
    if i > 1
        amp = amp_extractf(i,:)./amp_extractf(1,:);
    else
        amp = amp_extractf(i,:);
    end
    plot(tStep, amp, 'o-')
end

%%

figure
plot(tStep, amp_extractf(1,:), 'o-')

figure
plot(tStep, amp_extractf(2,:), 'o-')







%% ===================== LOCAL FUNCTIONS =====================

function T = readStep5BubbleDat(filePath)
    % Read comma-delimited .dat file written by:
    %   extract_bubble_step5_time_singlejob.py
    %
    % Expected header:
    % frame_id,step_time,node_label,X0,Y0,Z0,U1,U2,U3,x,y,z,angle_rad

    fid = fopen(filePath, 'r');
    if fid < 0
        error('Could not open %s', filePath);
    end

    headerLineNumber = 0;
    lineNumber = 0;

    while ~feof(fid)
        line = fgetl(fid);
        lineNumber = lineNumber + 1;

        if ischar(line)
            lineTrim = strtrim(line);
            if startsWith(lineTrim, 'frame_id')
                headerLineNumber = lineNumber;
                break;
            end
        end
    end

    fclose(fid);

    if headerLineNumber == 0
        error('Could not find header line beginning with frame_id in %s', filePath);
    end

    opts = detectImportOptions(filePath, 'FileType', 'text');
    opts.Delimiter = ',';
    opts.VariableNamesLine = headerLineNumber;
    opts.DataLines = [headerLineNumber + 1, Inf];

    T = readtable(filePath, opts);
end


function [xSorted, ySorted] = sortQuarterCurveByAngle(x, y)
    % Sort quarter-bubble points by polar angle.
    %
    % For a first-quadrant bubble surface, this orders the curve from
    % the +x-axis endpoint toward the +y-axis endpoint.

    x = x(:);
    y = y(:);

    ang = atan2(y, x);
    [~, idx] = sort(ang);

    xSorted = x(idx);
    ySorted = y(idx);
end


function [xFull, yFull] = mirrorQuarterCurveClosed(xq, yq)
    % Mirror first-quadrant curve across x=0 and y=0.
    %
    % Assumes xq,yq contain ONLY the real quarter bubble surface.
    %
    % The output is a continuous closed loop:
    %   Q1: +x,+y
    %   Q2: -x,+y
    %   Q3: -x,-y
    %   Q4: +x,-y

    xq = xq(:);
    yq = yq(:);

    [xq, yq] = sortQuarterCurveByAngle(xq, yq);
    [xq, yq] = uniqueXY(xq, yq, 1e-10);

    if numel(xq) < 3
        error('Quarter curve has fewer than 3 points. Check BubbleSet extraction.');
    end

    % Q1 goes from +x axis to +y axis.
    Q1x = xq;
    Q1y = yq;

    % Q2 should continue from +y axis to -x axis.
    Q2x = -flipud(xq(1:end-1));
    Q2y =  flipud(yq(1:end-1));

    % Q3 should continue from -x axis to -y axis.
    Q3x = -xq(2:end);
    Q3y = -yq(2:end);

    % Q4 should continue from -y axis back to +x axis.
    Q4x =  flipud(xq(2:end-1));
    Q4y = -flipud(yq(2:end-1));

    xFull = [Q1x; Q2x; Q3x; Q4x; Q1x(1)];
    yFull = [Q1y; Q2y; Q3y; Q4y; Q1y(1)];
end


function [xu, yu] = uniqueXY(x, y, tol)
    % Remove consecutive near-duplicate points.

    x = x(:);
    y = y(:);

    keep = true(size(x));

    for i = 2:numel(x)
        if hypot(x(i) - x(i-1), y(i) - y(i-1)) < tol
            keep(i) = false;
        end
    end

    xu = x(keep);
    yu = y(keep);
end


function colors = getColorMap(n, mapName)
    % Return n visually distinct colors.
    % Falls back to lines() if turbo() is unavailable.

    if nargin < 2
        mapName = 'turbo';
    end

    switch lower(mapName)
        case 'turbo'
            try
                colors = turbo(n);
            catch
                colors = lines(n);
            end

        case 'lines'
            colors = lines(n);

        case 'parula'
            colors = parula(n);

        otherwise
            warning('Unknown color map "%s". Using lines().', mapName);
            colors = lines(n);
    end
end


function radiusExtracted = reconstructAxisymSurfaceRadius(theta, ampExtract)
    % Reconstruct r(theta) from the extracted axisymmetric modal amplitudes.

    theta = theta(:);
    ampExtract = ampExtract(:);

    maxMode = numel(ampExtract) - 1;
    radiusExtracted = ampExtract(1) .* ones(size(theta));
    x = cos(theta);

    for n = 1:maxMode
        radiusExtracted = radiusExtracted + ampExtract(n + 1) .* ...
            extract_real_spherical_harmonic(x, 0, n, 0);
    end
end


function exportFullMirroredDat(T, outDat)
    % Export full mirrored bubble boundary for every frame.
    %
    % Output columns:
    %   frame_id, step_time, point_id, x, y
    %
    % This writes one closed mirrored curve per frame.

    allFrames = unique(T.frame_id);

    fid = fopen(outDat, 'w');
    if fid < 0
        error('Could not open output file: %s', outDat);
    end

    fprintf(fid, '# Full mirrored bubble boundary time history\n');
    fprintf(fid, '# Generated by plot_bubble_step5_time_mirrored.m\n');
    fprintf(fid, '# Assumes input BubbleSet contains only true quarter bubble surface.\n');
    fprintf(fid, '# columns:\n');
    fprintf(fid, 'frame_id,step_time,point_id,x,y\n');

    for iframe = 1:numel(allFrames)

        frameID = allFrames(iframe);
        Tf = T(T.frame_id == frameID, :);

        xq = Tf.x;
        yq = Tf.y;

        [xq, yq] = sortQuarterCurveByAngle(xq, yq);
        [xq, yq] = uniqueXY(xq, yq, 1e-10);

        [xFull, yFull] = mirrorQuarterCurveClosed(xq, yq);

        tStep = Tf.step_time(1);

        for ip = 1:numel(xFull)
            fprintf(fid, '%d,%.16e,%d,%.16e,%.16e\n', ...
                frameID, tStep, ip, xFull(ip), yFull(ip));
        end
    end

    fclose(fid);
end
