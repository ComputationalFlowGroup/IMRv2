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

load('../data/SicongJinChicken/chicken_Rt_data/Jin_15_33_11/ellipse_fitting_results.mat')
% load('../data/SicongJinChicken/chicken_Rt_data/Jin_15_43_39/ellipse_fitting_results.mat')
% load('../data/SicongJinChicken/PVA_Rt_data/Jin_17_25_03/ellipse_fitting_results.mat')

addpath src\common\

nFrames = size(CircleEdgePtSave,2);

% Frames used to estimate one global rotation.
% Use frames with visible deformation, not purely spherical/noisy frames.
frameIDsForRotation = 1:nFrames;
maxDeflectionAxis = 'x'; % align maximum deflection with extraction +x

[globalAxisAngleRaw, rotInfo] = estimateGlobalThetaRotation( ...
    CircleEdgePtSave, frameIDsForRotation, ...
    'NumAngles', 720, ...
    'CenterMode', 'circle');

[globalAxisAngle, poleInfo] = orientAxisAngleToMaxDeflection( ...
    CircleEdgePtSave, frameIDsForRotation, globalAxisAngleRaw, ...
    'CenterMode', 'circle', ...
    'NumTheta', 1001, ...
    'NumBins', 500, ...
    'SmoothWindow', 21, ...
    'TargetAxis', maxDeflectionAxis);

fprintf('Raw global axis angle      = %.3f deg\n', rad2deg(globalAxisAngleRaw));
fprintf('Corrected global axis angle = %.3f deg\n', rad2deg(globalAxisAngle));
fprintf('Applied correction          = %.3f deg\n', rad2deg(poleInfo.appliedCorrection));
fprintf('Maximum deflection aligned with extraction +%s axis.\n', ...
    maxDeflectionAxis);

%%
Nmax = 26;

for i = 1:nFrames

    data = CircleEdgePtSave{i};
    if isempty(data)
        continue
    end

    xq = data(:,1);
    yq = data(:,2);

    % Get cleaned theta-radius data after rotating max deflection to +x.

    [theta, radiusSmooth] = axisymThetaRadiusFromXY( ...
        xq, yq, globalAxisAngle, ...
        'CenterMode', 'circle', ...
        'NumTheta', 1001, ...
        'NumBins', 500, ...
        'SmoothWindow', 31, ...
        'LegendreCleanN', []);

    simdata = [theta, radiusSmooth];

    [modeFFT, ampFFT, phaseFFT] = fft_extract_axissym(simdata, Nmax);

    mode_extract_fft(:,i)  = modeFFT(:);
    amp_extract_fft(:,i)   = ampFFT(:);
    phase_extract_fft(:,i) = phaseFFT(:);
end





%% DEBUGGING
Nmax = 26;
debugFrame = nFrames;%min(15, nFrames);

for i = debugFrame % or 1:nFrames

    data = CircleEdgePtSave{i};
    if isempty(data)
        warning('Debug frame %d is empty; skipping orientation plot.', i);
        continue
    end

    xq = data(:,1);
    yq = data(:,2);

    % Get cleaned theta-radius data using the global rotation

    [theta, radiusSmooth, shapeInfo] = axisymThetaRadiusFromXY( ...
        xq, yq, globalAxisAngle, ...
        'CenterMode', 'circle', ...
        'NumTheta', 1001, ...
        'NumBins', 500, ...
        'SmoothWindow', 31, ...
        'LegendreCleanN', []);

    simdata = [theta, radiusSmooth];

    [modeFFT, ampFFT, phaseFFT] = fft_extract_axissym(simdata, Nmax);

    mode_extract_fft(:,i)  = modeFFT(:);
    amp_extract_fft(:,i)   = ampFFT(:);
    phase_extract_fft(:,i) = phaseFFT(:);

    % Reconstruction from the current script's axisymmetric SH extraction.
    modeCurrent = modeFFT;
    ampCurrent = ampFFT;
    radiusCurrent = reconstruct_axisym_spharm(theta, modeCurrent, ...
        ampCurrent);

    % Method 2: same SH basis, but sin(theta)-weighted

    [modeWeighted, ampWeighted, radiusWeighted] = ...
        fitAxisymSHModes(theta, radiusSmooth, Nmax, ...
        'WeightBySinTheta', true, ...
        'Ridge', 1e-10);

    % Errors

    errCurrent = sqrt(mean((radiusCurrent - radiusSmooth).^2)) / ...
        mean(radiusSmooth);
    errWeighted = sqrt(mean((radiusWeighted - radiusSmooth).^2)) / mean(radiusSmooth);

    fprintf('Frame %d\n', i);
    fprintf('  Current SH extraction RMSE    = %.4e\n', errCurrent);
    fprintf('  Weighted SH extraction RMSE   = %.4e\n', errWeighted);

    % Convert reconstructions to aligned x-y branches. In this debug view,
    % the maximum-deflection direction should be horizontal along +x.

    [xCurrentUpper, yCurrentUpper, xCurrentLower, yCurrentLower] = ...
        thetaRadiusToXYBranches(theta, radiusCurrent, 0);

    [xWeightedUpper, yWeightedUpper, xWeightedLower, yWeightedLower] = ...
        thetaRadiusToXYBranches(theta, radiusWeighted, 0);

    % Plot 1: aligned x-y reconstruction comparison

    figure
    plot(shapeInfo.xAlignedRaw, shapeInfo.yAlignedRaw, 'o', ...
        'DisplayName', 'raw aligned data')
    hold on

    plot(xCurrentUpper, yCurrentUpper, 'r-', 'LineWidth', 1.6, ...
        'DisplayName', 'current SH reconstruction')
    plot(xCurrentLower, yCurrentLower, 'r-', 'LineWidth', 1.6, ...
        'HandleVisibility', 'off')

    plot(xWeightedUpper, yWeightedUpper, 'k--', 'LineWidth', 1.6, ...
        'DisplayName', 'weighted SH reconstruction')
    plot(xWeightedLower, yWeightedLower, 'k--', 'LineWidth', 1.6, ...
        'HandleVisibility', 'off')

    Rplot = max(radiusSmooth);
    plot([-Rplot, Rplot], [0, 0], ...
         'b:', 'LineWidth', 1.2, ...
         'DisplayName', 'extraction x-axis')
    plot([0, 0], [-Rplot, Rplot], ...
         'c:', 'LineWidth', 1.0, ...
         'DisplayName', 'extraction y-axis')

    axis equal
    grid on
    xlabel('aligned x')
    ylabel('aligned y')
    title(sprintf('Frame %d: aligned x-y reconstruction comparison', i))
    legend('Location', 'best')

    % Plot 2: theta-radius reconstruction comparison

    figure
    plot(shapeInfo.thetaRaw, shapeInfo.radiusRaw, 'o', ...
        'DisplayName', 'raw folded data')
    hold on

    plot(shapeInfo.thetaBinned, shapeInfo.radiusBinned, '.', ...
        'MarkerSize', 12, ...
        'DisplayName', 'binned data')

    plot(theta, radiusSmooth, 'b-', 'LineWidth', 1.4, ...
        'DisplayName', 'smoothed/interpolated data')

    plot(theta, radiusCurrent, 'r-', 'LineWidth', 1.6, ...
        'DisplayName', sprintf('current SH reconstruction, RMSE %.3g', errCurrent))

    plot(theta, radiusWeighted, 'k--', 'LineWidth', 1.6, ...
        'DisplayName', sprintf('weighted SH reconstruction, RMSE %.3g', errWeighted))

    grid on
    xlabel('\theta')
    ylabel('radius')
    title(sprintf('Frame %d: theta-radius reconstruction comparison', i))
    legend('Location', 'best')

    % Plot 3: apples-to-apples amplitude comparison

    figure
    stem(modeCurrent, abs(ampCurrent), 'r', 'LineWidth', 1.3, ...
        'DisplayName', 'current unweighted SH |a_n|')
    hold on

    stem(modeWeighted, abs(ampWeighted), 'k--', 'LineWidth', 1.3, ...
        'DisplayName', 'weighted SH |a_n|')

    grid on
    xlabel('mode number')
    ylabel('amplitude magnitude')
    title(sprintf('Frame %d: extracted axisymmetric SH amplitudes', i))
    legend('Location', 'best')

end





%%
% load('../data/SicongJinChicken/chicken_Rt_data/processed_11.mat')

figure
nmodes = size(amp_extract_fft,1);
tStep = 1:size(amp_extract_fft,2);

for i = 1:nmodes
    plotl = ceil(sqrt(nmodes));
    subplot(plotl, plotl, i)
    if i > 1
        amp = amp_extract_fft(i,:)./amp_extract_fft(1,:);
    else
        amp = amp_extract_fft(i,:);
    end
    plot(tStep, amp, 'o-')
end

%%

figure
plot(tStep, amp_extract_fft(1,:), 'o-')

figure
plot(tStep, amp_extract_fft(2,:), 'o-')







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

function [axisAngle, info] = estimateGlobalThetaRotation(CircleEdgePtSave, frameIDs, varargin)
% estimateGlobalThetaRotation
%
% Finds one global in-plane rotation angle so that, after rotation, the
% upper and lower halves of the bubble collapse as closely as possible onto
% a single axisymmetric radius(theta) curve.
%
% The returned angle is the physical polar axis angle in the original x-y
% coordinates. After applying this angle, theta = 0 is placed on the side
% with the larger average radius.

    p = inputParser;
    addParameter(p, 'NumAngles', 720);
    addParameter(p, 'CenterMode', 'circle'); % 'circle', 'mean', or 'none'
    addParameter(p, 'ThetaGrid', linspace(0.03, pi-0.03, 401).');
    addParameter(p, 'MinValidFraction', 0.40);
    parse(p, varargin{:});

    numAngles = p.Results.NumAngles;
    centerMode = p.Results.CenterMode;
    thetaGrid = p.Results.ThetaGrid(:);
    minValidFraction = p.Results.MinValidFraction;

    alphas = linspace(0, pi, numAngles + 1);
    alphas(end) = [];

    scores = inf(size(alphas));
    directionMetric = zeros(size(alphas));

    for ia = 1:numel(alphas)

        alpha = alphas(ia);

        frameErrors = [];
        frameDirections = [];

        for kk = 1:numel(frameIDs)

            i = frameIDs(kk);

            if i < 1 || i > numel(CircleEdgePtSave)
                continue
            end

            data = CircleEdgePtSave{i};

            if isempty(data) || size(data,2) < 2
                continue
            end

            x = data(:,1);
            y = data(:,2);

            good = isfinite(x) & isfinite(y);
            x = x(good);
            y = y(good);

            if numel(x) < 20
                continue
            end

            [x, y, ~, ~] = centerXYLocal(x, y, centerMode);

            [rTop, rBot, valid] = mirroredBranchRadiiLocal(x, y, alpha, thetaGrid);

            if mean(valid) < minValidFraction
                continue
            end

            rt = rTop(valid);
            rb = rBot(valid);

            scale = var([rt; rb], 1);
            if scale < eps
                continue
            end

            err = mean((rt - rb).^2) / scale;
            frameErrors(end+1,1) = err; %#ok<AGROW>

            rMean = 0.5 * (rTop + rBot);
            validMean = isfinite(rMean);

            if nnz(validMean) > 10
                firstID = find(validMean, 1, 'first');
                lastID  = find(validMean, 1, 'last');

                % Positive means theta = 0 side is larger than theta = pi side.
                frameDirections(end+1,1) = rMean(firstID) - rMean(lastID); %#ok<AGROW>
            end
        end

        if ~isempty(frameErrors)
            scores(ia) = median(frameErrors);
        end

        if ~isempty(frameDirections)
            directionMetric(ia) = median(frameDirections);
        end
    end

    [bestScore, idxBest] = min(scores);
    alphaBest = alphas(idxBest);

    % alpha and alpha + pi have the same symmetry error.
    % Choose the direction such that theta = 0 corresponds to the larger
    % average radius side.
    if directionMetric(idxBest) < 0
        axisAngle = mod(alphaBest + pi, 2*pi);
    else
        axisAngle = mod(alphaBest, 2*pi);
    end

    info.axisAngle = axisAngle;
    info.axisAngleDegrees = rad2deg(axisAngle);
    info.alphas = alphas;
    info.scores = scores;
    info.bestScore = bestScore;
    info.directionMetric = directionMetric;
end


function [xCentered, yCentered, cx, cy] = centerXYLocal(x, y, centerMode)

    switch lower(centerMode)

        case 'circle'
            [cx, cy] = fitCircleKasaLocal(x, y);

        case 'mean'
            cx = mean(x);
            cy = mean(y);

        case 'none'
            cx = 0;
            cy = 0;

        otherwise
            error('Unknown CenterMode: %s', centerMode);
    end

    xCentered = x - cx;
    yCentered = y - cy;
end


function [cx, cy, R] = fitCircleKasaLocal(x, y)
% Algebraic least-squares circle fit.

    x = x(:);
    y = y(:);

    A = [2*x, 2*y, ones(size(x))];
    b = x.^2 + y.^2;

    sol = A \ b;

    cx = sol(1);
    cy = sol(2);
    c  = sol(3);

    R = sqrt(max(c + cx^2 + cy^2, 0));
end


function [rTopGrid, rBotGrid, valid] = mirroredBranchRadiiLocal(x, y, alpha, thetaGrid)

    ca = cos(alpha);
    sa = sin(alpha);

    % Rotate coordinates so the guessed symmetry axis is the new +x axis.
    xr =  ca*x + sa*y;
    yr = -sa*x + ca*y;

    rr = hypot(xr, yr);
    beta = atan2(yr, xr);

    top = beta >= 0;
    bot = beta <= 0;

    thetaTop = abs(beta(top));
    thetaBot = abs(beta(bot));

    rTop = rr(top);
    rBot = rr(bot);

    [thetaTop, rTop] = averageDuplicateThetaLocal(thetaTop, rTop);
    [thetaBot, rBot] = averageDuplicateThetaLocal(thetaBot, rBot);

    rTopGrid = nan(size(thetaGrid));
    rBotGrid = nan(size(thetaGrid));

    if numel(thetaTop) >= 5
        rTopGrid = interp1(thetaTop, rTop, thetaGrid, 'pchip', nan);
    end

    if numel(thetaBot) >= 5
        rBotGrid = interp1(thetaBot, rBot, thetaGrid, 'pchip', nan);
    end

    valid = isfinite(rTopGrid) & isfinite(rBotGrid);
end


function [thetaOut, rOut] = averageDuplicateThetaLocal(thetaIn, rIn)

    thetaIn = thetaIn(:);
    rIn = rIn(:);

    good = isfinite(thetaIn) & isfinite(rIn);
    thetaIn = thetaIn(good);
    rIn = rIn(good);

    [thetaIn, idx] = sort(thetaIn);
    rIn = rIn(idx);

    tol = 1e-8;
    thetaRounded = round(thetaIn / tol) * tol;

    [thetaOut, ~, ic] = unique(thetaRounded);
    rOut = accumarray(ic, rIn, [], @median);

    thetaOut = thetaOut(:);
    rOut = rOut(:);
end


function [thetaGrid, radiusSmooth, info] = axisymThetaRadiusFromXY(x, y, axisAngle, varargin)
% axisymThetaRadiusFromXY
%
% Takes raw x-y bubble edge points, applies one fixed global axis rotation,
% folds the data into theta in [0, pi], bins/averages duplicate top-bottom
% data, interpolates onto a clean theta grid, and optionally cleans the
% profile with a Legendre-series fit.

    p = inputParser;
    addParameter(p, 'CenterMode', 'circle'); % 'circle', 'mean', or 'none'
    addParameter(p, 'NumTheta', 1001);
    addParameter(p, 'NumBins', 500);
    addParameter(p, 'SmoothWindow', 31);
    addParameter(p, 'LegendreCleanN', []);
    addParameter(p, 'Ridge', 1e-8);
    parse(p, varargin{:});

    centerMode = p.Results.CenterMode;
    numTheta = p.Results.NumTheta;
    numBins = p.Results.NumBins;
    smoothWindow = p.Results.SmoothWindow;
    legendreCleanN = p.Results.LegendreCleanN;
    ridge = p.Results.Ridge;

    x = x(:);
    y = y(:);

    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    [xCentered, yCentered, cx, cy] = centerXYLocal(x, y, centerMode);

    ca = cos(axisAngle);
    sa = sin(axisAngle);

    % Rotate coordinates so the global polar axis is the new +x axis.
    xr =  ca*xCentered + sa*yCentered;
    yr = -sa*xCentered + ca*yCentered;

    rRaw = hypot(xr, yr);

    % Full signed polar angle first.
    beta = atan2(yr, xr);

    % Then fold top and bottom halves into theta in [0, pi].
    thetaRaw = abs(beta);

    % Bin-average noisy/duplicate folded data.
    [thetaBinned, radiusBinned] = binAverageThetaLocal(thetaRaw, rRaw, numBins);

    thetaGrid = linspace(0, pi, numTheta).';

    % Interpolate onto a clean uniform theta grid.
    radiusInterp = interp1(thetaBinned, radiusBinned, thetaGrid, 'pchip', 'extrap');

    radiusSmooth = radiusInterp;

    % Light smoothing. Use odd window length.
    if ~isempty(smoothWindow) && smoothWindow > 2
        smoothWindow = round(smoothWindow);
        if mod(smoothWindow,2) == 0
            smoothWindow = smoothWindow + 1;
        end

        radiusSmooth = smoothdata(radiusSmooth, 'movmedian', smoothWindow);
        radiusSmooth = smoothdata(radiusSmooth, 'sgolay', smoothWindow);
    end

    % Optional final cleanup by projecting onto P_n(cos theta).
    % This is usually better than FFT smoothing for axisymmetric harmonics.
    if ~isempty(legendreCleanN)
        [coeffClean, radiusClean] = fitAxisymLegendreModes(thetaGrid, radiusSmooth, legendreCleanN, ...
            'Ridge', ridge);
        radiusSmooth = radiusClean;
    else
        coeffClean = [];
    end

    % Useful plotting outputs in aligned coordinates.
    xUpper = radiusSmooth .* cos(thetaGrid);
    yUpper = radiusSmooth .* sin(thetaGrid);

    xLower = radiusSmooth .* cos(thetaGrid);
    yLower = -radiusSmooth .* sin(thetaGrid);

    info.cx = cx;
    info.cy = cy;
    info.axisAngle = axisAngle;
    info.axisAngleDegrees = rad2deg(axisAngle);

    info.xRawCentered = xCentered;
    info.yRawCentered = yCentered;

    info.xAlignedRaw = xr;
    info.yAlignedRaw = yr;

    info.thetaRaw = thetaRaw;
    info.radiusRaw = rRaw;

    info.thetaBinned = thetaBinned;
    info.radiusBinned = radiusBinned;

    info.xUpperAligned = xUpper;
    info.yUpperAligned = yUpper;
    info.xLowerAligned = xLower;
    info.yLowerAligned = yLower;

    info.legendreCleanCoeff = coeffClean;
end

function [coeff, radiusFit, ampNorm] = fitAxisymLegendreModes(theta, radius, Nmax, varargin)
% fitAxisymLegendreModes
%
% Fits
%
%   radius(theta) = sum_{n=0}^{Nmax} coeff(n+1) P_n(cos(theta))
%
% using sin(theta)-weighted least squares.
%
% coeff(1) is the spherical/base radius contribution.
% coeff(n+1) is the absolute radial amplitude of P_n(cos theta).
% ampNorm(n+1) is coeff(n+1)/coeff(1).

    p = inputParser;
    addParameter(p, 'Ridge', 0);
    addParameter(p, 'PenaltyPower', 2);
    parse(p, varargin{:});

    ridge = p.Results.Ridge;
    penaltyPower = p.Results.PenaltyPower;

    theta = theta(:);
    radius = radius(:);

    good = isfinite(theta) & isfinite(radius);
    theta = theta(good);
    radius = radius(good);

    x = cos(theta);

    B = zeros(numel(theta), Nmax + 1);

    for n = 0:Nmax
        L = legendre(n, x.');
        B(:, n+1) = L(1,:).';
    end

    % Spherical-axisymmetric projection weight.
    w = sin(theta);

    % Avoid zero weights exactly at the poles.
    w = max(w, 1e-6 * max(w));

    sw = sqrt(w);

    Bw = B .* sw;
    rw = radius .* sw;

    if ridge > 0
        penalties = (0:Nmax).^penaltyPower;
        penalties(1) = 0; % do not penalize mean radius
        D = diag(penalties);

        coeff = (Bw.'*Bw + ridge*(D.'*D)) \ (Bw.'*rw);
    else
        coeff = Bw \ rw;
    end

    radiusFit = B * coeff;

    if abs(coeff(1)) > eps
        ampNorm = coeff / coeff(1);
    else
        ampNorm = nan(size(coeff));
    end
end

function [thetaBinned, radiusBinned] = binAverageThetaLocal(thetaRaw, radiusRaw, numBins)

    thetaRaw = thetaRaw(:);
    radiusRaw = radiusRaw(:);

    good = isfinite(thetaRaw) & isfinite(radiusRaw) & ...
           thetaRaw >= 0 & thetaRaw <= pi;

    thetaRaw = thetaRaw(good);
    radiusRaw = radiusRaw(good);

    edges = linspace(0, pi, numBins + 1);
    binID = discretize(thetaRaw, edges);

    good = isfinite(binID);
    binID = binID(good);
    thetaRaw = thetaRaw(good);
    radiusRaw = radiusRaw(good);

    thetaBinned = accumarray(binID, thetaRaw, [numBins, 1], @mean, nan);
    radiusBinned = accumarray(binID, radiusRaw, [numBins, 1], @median, nan);

    good = isfinite(thetaBinned) & isfinite(radiusBinned);

    thetaBinned = thetaBinned(good);
    radiusBinned = radiusBinned(good);

    [thetaBinned, idx] = sort(thetaBinned);
    radiusBinned = radiusBinned(idx);

    % Remove repeated theta values, just in case.
    [thetaBinned, ia] = unique(thetaBinned, 'stable');
    radiusBinned = radiusBinned(ia);

    if numel(thetaBinned) < 5
        error('Too few valid theta bins. Check axis rotation, centering, or input points.');
    end
end


function radiusRecon = reconstructFftAxissym(theta, modeFFT, ampFFT, phaseFFT, radiusReference)
% reconstructFftAxissym
%
% Reconstructs a theta-radius profile from the output of fft_extract_axissym.
%
% Assumed convention:
%
%   r(theta) = a0 + sum_k amp_k cos(mode_k theta + phase_k)
%
% If your fft_extract_axissym uses the opposite phase convention, change
% the plus sign to a minus sign below.

    theta = theta(:);
    modeFFT = modeFFT(:);
    ampFFT = ampFFT(:);
    phaseFFT = phaseFFT(:);

    if nargin < 5 || isempty(radiusReference)
        radiusReference = [];
    else
        radiusReference = radiusReference(:);
    end

    radiusRecon = zeros(size(theta));

    hasZeroMode = any(modeFFT == 0);

    if hasZeroMode
        idx0 = find(modeFFT == 0, 1, 'first');
        radiusRecon = radiusRecon + ampFFT(idx0);
    else
        % If fft_extract_axissym does not return the mean radius,
        % use the mean of the input profile.
        if isempty(radiusReference)
            warning('No zero mode found and no radius reference supplied. Using zero mean.');
        else
            radiusRecon = radiusRecon + mean(radiusReference, 'omitnan');
        end
    end

    for k = 1:numel(modeFFT)

        n = modeFFT(k);

        if n == 0
            continue
        end

        radiusRecon = radiusRecon + ampFFT(k) .* cos(n .* theta + phaseFFT(k));

        % If reconstruction looks phase-flipped, try this instead:
        % radiusRecon = radiusRecon + ampFFT(k) .* cos(n .* theta - phaseFFT(k));
    end
end




function radiusRecon = reconstruct_axisym_spharm(theta, mode_extractf, amp_extractf)
% reconstruct_axisym_spharm
%
% Reconstructs the radius profile using the SAME convention as your
% fft_extract_axissym.m:
%
%   r(theta) = a0 + sum_{n=1}^{N} a_n Y_n^0(theta)
%
% where:
%
%   Y_n^0(theta) = sqrt((2n+1)/(4*pi)) P_n(cos(theta))
%
% Important:
%   This does NOT use cos(n*theta + phase).
%   The phase output from fft_extract_axissym should be ignored.

    theta = theta(:);
    mode_extractf = mode_extractf(:);
    amp_extractf = amp_extractf(:);

    radiusRecon = zeros(size(theta));

    x = cos(theta);

    for k = 1:numel(mode_extractf)

        n = mode_extractf(k);
        a = amp_extractf(k);

        if n == 0
            radiusRecon = radiusRecon + a;
        else
            Yn0 = extract_real_spherical_harmonic(x, 0, n, 0);
            radiusRecon = radiusRecon + a .* Yn0;
        end
    end
end



function [mode, amp, radiusFit, ampNorm] = fitAxisymSHModes(theta, radius, max_mode, varargin)
% fitAxisymSHModes
%
% Fits the same basis as fft_extract_axissym:
%
%   r(theta) = a0 + sum_{n=1}^{N} a_n Y_n^0(theta)
%
% but allows sin(theta)-weighted least squares.

    p = inputParser;
    addParameter(p, 'WeightBySinTheta', true);
    addParameter(p, 'Ridge', 0);
    addParameter(p, 'PenaltyPower', 2);
    parse(p, varargin{:});

    weightBySinTheta = p.Results.WeightBySinTheta;
    ridge = p.Results.Ridge;
    penaltyPower = p.Results.PenaltyPower;

    theta = theta(:);
    radius = radius(:);

    good = isfinite(theta) & isfinite(radius);
    theta = theta(good);
    radius = radius(good);

    x = cos(theta);

    basis = ones(numel(theta), max_mode + 1);

    for n = 1:max_mode
        basis(:, n+1) = extract_real_spherical_harmonic(x, 0, n, 0);
    end

    if weightBySinTheta
        w = sin(theta);
        w = max(w, 1e-8 * max(w));
    else
        w = ones(size(theta));
    end

    sw = sqrt(w);

    Bw = basis .* sw;
    rw = radius .* sw;

    if ridge > 0
        penalties = (0:max_mode).^penaltyPower;
        penalties(1) = 0;
        D = diag(penalties);

        amp = (Bw.'*Bw + ridge*(D.'*D)) \ (Bw.'*rw);
    else
        amp = Bw \ rw;
    end

    radiusFit = basis * amp;

    mode = (0:max_mode).';

    if abs(amp(1)) > eps
        ampNorm = amp ./ amp(1);
    else
        ampNorm = nan(size(amp));
    end
end

function [xUpper, yUpper, xLower, yLower] = thetaRadiusToXYBranches(theta, radius, axisAngle)
% thetaRadiusToXYBranches
%
% Converts an axisymmetric theta-radius profile into upper and lower x-y
% branches, then rotates them back into the original centered x-y frame.

    theta = theta(:);
    radius = radius(:);

    xUpperA = radius .* cos(theta);
    yUpperA = radius .* sin(theta);

    xLowerA = radius .* cos(theta);
    yLowerA = -radius .* sin(theta);

    ca = cos(axisAngle);
    sa = sin(axisAngle);

    xUpper = ca*xUpperA - sa*yUpperA;
    yUpper = sa*xUpperA + ca*yUpperA;

    xLower = ca*xLowerA - sa*yLowerA;
    yLower = sa*xLowerA + ca*yLowerA;
end

function [axisAngleOut, info] = orientAxisAngleToMaxDeflection(CircleEdgePtSave, frameIDs, axisAngleIn, varargin)
% orientAxisAngleToMaxDeflection
%
% Takes a candidate symmetry-axis angle and resolves the 90-degree ambiguity
% so that the maximum radius/deflection is placed on the requested aligned
% coordinate axis.
%
% It tests:
%
%   axisAngleIn
%   axisAngleIn + pi/2
%   axisAngleIn + pi
%   axisAngleIn + 3*pi/2
%
% and chooses the one with the largest radius near the target axis. For
% TargetAxis='x', the winner puts the maximum at theta = 0 (+x). For
% TargetAxis='y', the winner puts the maximum at theta = pi/2.

    p = inputParser;
    addParameter(p, 'CenterMode', 'circle');
    addParameter(p, 'NumTheta', 1001);
    addParameter(p, 'NumBins', 500);
    addParameter(p, 'SmoothWindow', 21);
    addParameter(p, 'PoleWidth', 0.10);      % radians near theta = 0
    addParameter(p, 'EquatorWidth', 0.10);   % radians near theta = pi/2
    addParameter(p, 'UseDeflectionFromMean', true);
    addParameter(p, 'TargetAxis', 'x');
    parse(p, varargin{:});

    centerMode = p.Results.CenterMode;
    numTheta = p.Results.NumTheta;
    numBins = p.Results.NumBins;
    smoothWindow = p.Results.SmoothWindow;
    poleWidth = p.Results.PoleWidth;
    equatorWidth = p.Results.EquatorWidth;
    useDeflectionFromMean = p.Results.UseDeflectionFromMean;
    targetAxis = validatestring(p.Results.TargetAxis, {'x', 'y'}, ...
        mfilename, 'TargetAxis');

    corrections = [0, pi/2, pi, 3*pi/2];
    candidateAngles = mod(axisAngleIn + corrections, 2*pi);

    scores = nan(size(candidateAngles));
    poleValues = nan(size(candidateAngles));
    equatorValues = nan(size(candidateAngles));
    oppositePoleValues = nan(size(candidateAngles));

    for j = 1:numel(candidateAngles)

        alpha = candidateAngles(j);

        frameScores = [];
        framePoleVals = [];
        frameEqVals = [];
        frameOppVals = [];

        for kk = 1:numel(frameIDs)

            i = frameIDs(kk);

            if i < 1 || i > numel(CircleEdgePtSave)
                continue
            end

            data = CircleEdgePtSave{i};

            if isempty(data) || size(data,2) < 2
                continue
            end

            xq = data(:,1);
            yq = data(:,2);

            try
                [theta, radius] = axisymThetaRadiusFromXY( ...
                    xq, yq, alpha, ...
                    'CenterMode', centerMode, ...
                    'NumTheta', numTheta, ...
                    'NumBins', numBins, ...
                    'SmoothWindow', smoothWindow, ...
                    'LegendreCleanN', []);
            catch
                continue
            end

            theta = theta(:);
            radius = radius(:);

            if useDeflectionFromMean
                signal = radius - mean(radius, 'omitnan');
            else
                signal = radius;
            end

            poleID = theta <= poleWidth;
            eqID = abs(theta - pi/2) <= equatorWidth;
            oppPoleID = theta >= pi - poleWidth;

            if nnz(poleID) < 3 || nnz(eqID) < 3 || nnz(oppPoleID) < 3
                continue
            end

            poleVal = mean(signal(poleID), 'omitnan');
            eqVal = mean(signal(eqID), 'omitnan');
            oppVal = mean(signal(oppPoleID), 'omitnan');

            switch targetAxis
                case 'x'
                    % theta = 0 is the extraction +x direction.
                    score = poleVal - max(eqVal, oppVal);
                case 'y'
                    % theta = pi/2 is the extraction y direction.
                    score = eqVal - max(poleVal, oppVal);
            end

            frameScores(end+1,1) = score; %#ok<AGROW>
            framePoleVals(end+1,1) = poleVal; %#ok<AGROW>
            frameEqVals(end+1,1) = eqVal; %#ok<AGROW>
            frameOppVals(end+1,1) = oppVal; %#ok<AGROW>
        end

        if ~isempty(frameScores)
            scores(j) = median(frameScores, 'omitnan');
            poleValues(j) = median(framePoleVals, 'omitnan');
            equatorValues(j) = median(frameEqVals, 'omitnan');
            oppositePoleValues(j) = median(frameOppVals, 'omitnan');
        end
    end

    [~, idxBest] = max(scores);

    axisAngleOut = candidateAngles(idxBest);

    info.axisAngleIn = axisAngleIn;
    info.axisAngleOut = axisAngleOut;
    info.axisAngleInDegrees = rad2deg(axisAngleIn);
    info.axisAngleOutDegrees = rad2deg(axisAngleOut);
    info.targetAxis = targetAxis;

    info.corrections = corrections;
    info.correctionsDegrees = rad2deg(corrections);
    info.appliedCorrection = corrections(idxBest);
    info.appliedCorrectionDegrees = rad2deg(corrections(idxBest));

    info.candidateAngles = candidateAngles;
    info.candidateAnglesDegrees = rad2deg(candidateAngles);
    info.scores = scores;

    info.poleValues = poleValues;
    info.equatorValues = equatorValues;
    info.oppositePoleValues = oppositePoleValues;
end
