%% graded geometry images
% Parameters
l1 = 1.5;
l2 = 3.0;
a = 2.0;
n = 0.3;
maxr0 = 3.5;
res = 1000;
%%
% Polar grid
r = linspace(0, maxr0, res);
theta = linspace(0, 2*pi, res);
[R, T] = meshgrid(r, theta);
% Convert to cartesian
X = R .* cos(T);
Y = R .* sin(T);

% Radial gradient function
G = zeros(size(R));
mask1 = R > l1 & R < l2;
mask2 = R >= l2;
f_lambda = (l2 - R(mask1)) ./ (R(mask1) - l1);
G(mask1) = (1 + f_lambda.^a).^((n - 1)/a);
G(mask2) = 1;

% Plot
figure;

% Step 1: Plot rectangle background with max color
colormap(sky(1000)); %turbo, jet, parula, sky
cmap = colormap;
maxRGB = cmap(end,:); % The maximum value of the gradient
hold on;
fill([-maxr0 maxr0 maxr0 -maxr0], [-maxr0 -maxr0 maxr0 maxr0], maxRGB, ...
    'EdgeColor', 'none'); 
% Rest of your overlays...
axis equal off;

% Step 2: Overlay the gradient
p = pcolor(X, Y, G);
shading interp;
%uistack(p,'top');
cb = colorbar;
cb.FontSize = 14;
%title('2D Radial Gradient Cross-Section','FontSize',10);

% opacity for colors
alpha(0.8);

% Draw polar grid manually
hold on;
% for rt = 0.5:0.3:maxr0
%     theta_grid = linspace(0, 2*pi, res/2);
%     plot(rt*cos(theta_grid), rt*sin(theta_grid), '-', 'Color', [0.6 0.6 0.6]);
% end
% dtheta = pi/8;
% for th = 0:dtheta:2*pi
%     plot([0 maxr0*cos(th)], [0 maxr0*sin(th)], '-', 'Color', [0.6 0.6 0.6]);
% end

% Radial dashed lines for l1 and l2
l_theta = linspace(0, 2*pi, res);
plot(l1*cos(l_theta), l1*sin(l_theta), '--k', 'LineWidth', 1.2);
plot(l2*cos(l_theta), l2*sin(l_theta), '--k', 'LineWidth', 1.2);
%plot(maxr0*cos(l_theta), maxr0 *sin(l_theta), '-k', 'LineWidth', 1.5);

% Overlay boundary rings
%hold on;
%viscircles([0, 0], l1, 'LineStyle', '--', 'Color', 'k');
%viscircles([0, 0], l2, 'LineStyle', '--', 'Color', 'k');

% Add arrows and labels
%quiver(0, 0, 0.8, 0.1, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
%text(1, 0.2, 'Near Field', 'FontSize', 12);

%quiver(0, 0, 1.5, 0.5, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
%text(1.8, 0.6, 'Graded Region', 'FontSize', 12);

%quiver(0, 0, 2.5, 0.8, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
%text(2.8, 0.9, 'Far Field', 'FontSize', 12);

% Arrows + Labels
axis equal
%quiver(0, 0, l1, 0, 0, 'k', 'LineWidth', 1.2,'ShowArrowHead','on','AutoScale','on');
%text(0.75, 0.1, 'Near', 'FontSize', 15);

%quiver(0, 0, l2, 1.23, 0.92, 'k', 'LineWidth', 1.2, 'ShowArrowHead','on','AutoScale','on');
%text(2, 1.15, 'Graded', 'FontSize', 15);

%quiver(0, 0, maxr0, 2.75, 0, 'k', 'LineWidth', 1.2, 'ShowArrowHead','on','AutoScale','on');
%text(3.15, 2.5, 'Far', 'FontSize', 12);
%quiver(0, 0, maxr0, 2.75, 1.0, 'k', 'LineWidth', 1.2, 'ShowArrowHead','on','AutoScale','on');
%text(3.15, 2.0, 'Far', 'FontSize', 15);

% rectangle boundary
%rectangle('Position',[-maxr0,-maxr0,maxr0*2,maxr0*2], 'EdgeColor','k','LineWidth',1.5);

print('-dpng','-r600','2D_geo')
print('-depsc','-r600','2D_geo')

%%
% annotation
xlim([-maxr0 maxr0]);
ylim([-maxr0 maxr0]);
ax = gca;
ax.Units = 'normalized';
ax_pos = ax.Position;

% convert data coordinates to normalized figure coordinates
data2norm = @(x, y) [(x - ax.XLim(1)) / diff(ax.XLim) * ax_pos(3) + ax_pos(1), ...
                     (y - ax.YLim(1)) / diff(ax.YLim) * ax_pos(4) + ax_pos(2)];

% Arrow 1: "Near" (point to l1)
[start1] = data2norm(0, 0);
[end1]   = data2norm(l1 * cos(dtheta), l1 * sin(dtheta));
annotation('arrow', [start1(1), end1(1)], [start1(2), end1(2)], ...
    'LineWidth', 1.2);
textPos1 = data2norm(l1 * cos(dtheta) - 0.1, l1 * sin(dtheta) - 0.1);
annotation('textbox', [textPos1(1), textPos1(2), 0.05, 0.05], ...
    'String', 'Near', 'EdgeColor', 'none', 'FontSize', 12);

% Arrow 2: "Gradient" (point to l2)
[start2] = data2norm(0, 0);
[end2]   = data2norm(l2 * cos(dtheta), l2 * sin(dtheta));
annotation('arrow', [start2(1), end2(1)], [start2(2), end2(2)], ...
    'LineWidth', 1.2);
textPos2 = data2norm(l2 * cos(dtheta) + 0.2, l2 * sin(dtheta) - 0.05);
annotation('textbox', [textPos2(1), textPos2(2), 0.05, 0.05], ...
    'String', 'Gradient', 'EdgeColor', 'none', 'FontSize', 9);

% Arrow 3: "Far" (point to edge of domain)
[start3] = data2norm(0, 0);
[end3]   = data2norm(maxr0 * cos(dtheta), maxr0 * sin(dtheta));
annotation('arrow', [start3(1), end3(1)], [start3(2), end3(2)], ...
    'LineWidth', 1.2);
textPos3 = data2norm(maxr0 * cos(dtheta) + 0.2, maxr0* sin(dtheta) + 0.05);
annotation('textbox', [textPos3(1), textPos3(2), 0.05, 0.05], ...
    'String', 'Far', 'EdgeColor', 'none', 'FontSize', 9);

%%
% Arrows + Labels
quiver(0, 0, l1, 0, 0, 'k', 'LineWidth', 1.2,'ShowArrowHead','on','AutoScale','on','AutoScaleFactor',l1,'Alignment','tail');
text(0.75, 0.1, 'Near', 'FontSize', 12);

%quiver(0, 0, l2, 0.75, 0, 'k', 'LineWidth', 1.2, 'ShowArrowHead','on','AutoScale','on','AutoScaleFactor',l2,'Alignment','tail');
%text(1.75, 0.68, 'Graded', 'FontSize', 12);
quiver(0, 0, l2, 1.5, 0, 'k', 'LineWidth', 1.2, 'ShowArrowHead','on','AutoScale','on','AutoScaleFactor',l2,'Alignment','tail');
text(1.75, 0.68, 'Graded', 'FontSize', 12);

quiver(0, 0, maxr0, 2.75, 0, 'k', 'LineWidth', 1.2, 'ShowArrowHead','on','AutoScale','on','AutoScaleFactor',maxr0,'Alignment','tail');
text(3.0, 1.0, 'Far', 'FontSize', 12);

% debuggin
quiver(0,0,0,0,'g','LineWidth',1.2);

% rectangle boundary
rectangle('Position',[-maxr0,-maxr0,maxr0*2,maxr0*2], 'EdgeColor','k','LineWidth',1.5);
%%
% Rectangular block dimensions
Lx = 7; Ly = 5; Lz = 3; % Block size (length, width, height)

% Grid on top face
x = linspace(-Lx/2, Lx/2, res);
y = linspace(-Ly/2, Ly/2, res);
[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);
% prefunction G (homogeneous)
G = zeros(size(R));
G0 = 100; 
G1 = 1000;

% function
mask1 =  R > l1 & R < l2;
mask2 = R >= l2;
f = (l2 - R(mask1)) ./ (R(mask1) - l1);
G(mask1) = (1+f.^a).^((n-1)/a);
G(mask2) = 1;
% G = G0 + (G1 - G0) *G;

% Plot
figure; hold on;

% % Front face (y = -Ly/2) with real gradient
% surf(XF, -Ly/2 * ones(size(XF)), ZF, GF, ...
%      'EdgeColor', 'none', 'FaceAlpha', 0.8); hold on;
% 
% % Side face (x = Lx/2) with real gradient
% surf(Lx/2 * ones(size(YF)), YF, ZS, GS, ...
%      'EdgeColor', 'none', 'FaceAlpha', 0.8);
% 
% % Top face (X-Y plane) with real gradient
% surf(X, Y, Lz * ones(size(G)), G, ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.7);

% NOW color the **outside edge of the front face** with the max gradient
% (Top edge of front face = y = -Ly/2, z = Lz)
% fill3(x, -Ly/2 * ones(size(x)), Lz * ones(size(x)), ...
%       maxColor * ones(size(x)), ...
%       'EdgeColor', 'none', 'FaceAlpha', 1);
% 
% % Same for side face vertical edge (x = Lx/2, z)
% fill3(Lx/2 * ones(size(y)), y, Lz * ones(size(y)), ...
%       maxColor * ones(size(y)), ...
%       'EdgeColor', 'none', 'FaceAlpha', 1);

% Base of block
%patch([-Lx/2 Lx/2 Lx/2 -Lx/2], [-Ly/2 -Ly/2 Ly/2 Ly/2], [0 0 0 0], [0.9 0.9 0.9]);
% xtop = [-Lx/2 Lx/2 Lx/2 -Lx/2];
% ytop = [-Ly/2 -Ly/2 Ly/2 Ly/2];
% ztop = Lz * ones(size(xtop));
% plot3(xtop, ytop, ztop, 'k-','LineWidth',1.5);
% hold on;

% Settings
colormap(sky(res));
% cb = colorbar;
% cb.FontSize = 14;
axis equal;
axis off;

% Front face (y = -Ly/2)
[XF, ZF] = meshgrid(x, linspace(0, Lz, res));
RF = sqrt(XF.^2 + (Ly/2).^2);
GF = zeros(size(RF));
mask1f = RF > l1 & RF < l2;
mask2f = RF >= l2;
f_lambda_f = (l2 - RF(mask1f)) ./ (RF(mask1f) - l1);
GF(mask1f) = (1 + f_lambda_f.^a).^((n - 1)/a);
GF(mask2f) = 1;
surf(XF, -Ly/2 * ones(size(XF)), ZF, GF, 'EdgeColor', 'none'); hold on;

% Front face border (y = -Ly/2)
front_y = -Ly/2;
x_edge = [-Lx/2, Lx/2, Lx/2, -Lx/2, -Lx/2];
z_edge = [0, 0, Lz, Lz, 0];
y_edge = front_y * ones(size(x_edge));
plot3(x_edge, y_edge, z_edge, 'k-', 'LineWidth', 1.5);

% Side face (x = Lx/2)
[YF, ZS] = meshgrid(y, linspace(0, Lz, res));
RS = sqrt((Lx/2).^2 + YF.^2);
GS = zeros(size(RS));
mask1s = RS > l1 & RS < l2;
mask2s = RS >= l2;
f_lambda_s = (l2 - RS(mask1s)) ./ (RS(mask1s) - l1);
GS(mask1s) = (1 + f_lambda_s.^a).^((n - 1)/a);
GS(mask2s) = 1;
surf(Lx/2 * ones(size(YF)), YF, ZS, GS, 'EdgeColor', 'none'); hold on;
% Side face border (x = Lx/2)
side_x = Lx/2;
y_edge = [-Ly/2, Ly/2, Ly/2, -Ly/2, -Ly/2];
z_edge = [0, 0, Lz, Lz, 0];
x_edge = side_x * ones(size(y_edge));
plot3(x_edge, y_edge, z_edge, 'k-', 'LineWidth', 1.5);

% Top face
% surf(X, Y, Lz * ones(size(G)), G, 'EdgeColor', 'none');
% [ZF, ZT] = meshgrid(z, linspace(0, Lz, res));
% RT = sqrt((Lx/2).^2 + ZF.^2);
% GS = zeros(size(RT));
% mask1t = RT > l1 & RT < l2;
% mask2t = RT >= l2;
% f_lambda_t = (l2 - RT(mask1st)) ./ (RT(mask1st) - l1);
% GS(mask1st) = (1 + f_lambda_s.^a).^((n - 1)/a);
% GS(mask2st) = 1;
% surf(Lx/2 * ones(size(ZF)), ZF, ZT, GS, 'EdgeColor', 'none');

% Draw sides of the block (light gray)
% patch([-Lx/2 Lx/2 Lx/2 -Lx/2], [-Ly/2 -Ly/2 Ly/2 Ly/2], [0 0 0 0], [0.9 0.9 0.9]); hold on;
% patch([-Lx/2 -Lx/2 -Lx/2 -Lx/2], [-Ly/2 -Ly/2 Ly/2 Ly/2], [0 Lz Lz 0], [0.9 0.9 0.9]);
% patch([Lx/2 Lx/2 Lx/2 Lx/2], [-Ly/2 -Ly/2 Ly/2 Ly/2], [0 Lz Lz 0], [0.9 0.9 0.9]);
% xfront = [-Lx/2 Lx/2 Lx/2 -Lx/2];
% yfront = [-Ly/2 -Ly/2 Ly/2 Ly/2];
% zfront = Lz * ones(size(xedge));
% plot3(xedge, yedge, zedge, 'k-','LineWidth',1.5);

% Top face with radial gradient
surf(X, Y, Lz*ones(size(G)), G, 'EdgeColor', 'none');
% Top face border
top_z = Lz;
x_edge = [-Lx/2, Lx/2, Lx/2, -Lx/2, -Lx/2];
y_edge = [-Ly/2, -Ly/2, Ly/2, Ly/2, -Ly/2];
z_edge = top_z * ones(size(x_edge));
plot3(x_edge, y_edge, z_edge, 'k-', 'LineWidth', 1.5);

view(5, 25);
axis off;
axis equal;
clim([G0/G1, G1/G1]);
print('-dpng','-r600','3D_geo')
print('-depsc','-r600','3D_geo')
%title('Radial Gradient on Top Face of Block', 'FontSize', 10);



%%
% isosurface_modulus_shells.m
% Concentric isosurfaces representing radial modulus field G(r0)

% Parameters
G0 = 100;
G1 = 1000;
l1 = 5;
l2 = 10;
a = 2;
n = 0.3;

% Grid setup
N = 100;    % resolution
L = 15;     % domain size
[x, y, z] = meshgrid(linspace(-L, L, N));
r0 = sqrt(x.^2 + y.^2 + z.^2);

% Initialize modulus field
G = zeros(size(r0));

% Region 1: r0 <= l1
G(r0 <= l1) = G0;

% Region 2: l1 < r0 < l2
idx2 = (r0 > l1) & (r0 < l2);
f = (l2 - r0(idx2)) ./ (r0(idx2) - l1);
G(idx2) = G0 + (G1 - G0) .* (1 + f.^a).^((n - 1)/a);

% Region 3: r0 >= l2
G(r0 >= l2) = G1;

% Define isosurface levels
num_levels = 4;  % number of shells
iso_values = linspace(G0, G1, num_levels + 2);  % skip extreme edges

% Plot
figure;
hold on;
for i = 2:(num_levels + 1)  % skip first and last to avoid flat regions
    fv = isosurface(x, y, z, G, iso_values(i));
    p = patch(fv);
    isonormals(x, y, z, G, p);
    set(p, 'FaceColor', 'interp', ...
           'EdgeColor', 'none', ...
           'FaceVertexCData', iso_values(i)*ones(size(fv.vertices,1),1), ...
           'FaceAlpha', 0.5);  % semi-transparent
end

colormap(sky);
colorbar;
clim([G0, G1]);
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;
title('Concentric Isosurfaces of G(r_0)');
view(3);
camlight; lighting flat; %gouraud;
%%

% animate_modulus_transition.m
% Animated radial buildup of modulus field G(r0)

% Parameters
G0 = 100;
G1 = 1000;
l1 = 5;
l2 = 10;
a = 2;
n = 3;

% 2D grid (cross-section)
N = 300;
L = 15;
[x, y] = meshgrid(linspace(-L, L, N));
r0 = sqrt(x.^2 + y.^2);

% Initialize figure
figure;
set(gcf, 'Color', 'w');
filename = 'modulus_transition.gif';  % Output file

% Time steps simulate radial expansion
num_frames = 60;
r_max_list = linspace(0, L, num_frames);

for k = 1:num_frames
    r_cutoff = r_max_list(k);
    
    % Compute G only up to current radius
    G = nan(size(r0));  % NaNs for "not yet formed" regions
    
    idx = r0 <= r_cutoff;
    r_vals = r0(idx);
    
    G_local = zeros(size(r_vals));
    % Region 1
    G_local(r_vals <= l1) = G0;
    
    % Region 2
    idx2 = (r_vals > l1) & (r_vals < l2);
    f = (l2 - r_vals(idx2)) ./ (r_vals(idx2) - l1);
    G_local(idx2) = G0 + (G1 - G0) .* (1 + f.^a).^((n - 1)/a);
    
    % Region 3
    G_local(r_vals >= l2) = G1;
    
    G(idx) = G_local;

    % Plot
    imagesc(linspace(-L, L, N), linspace(-L, L, N), G);
    axis equal;
    set(gca, 'YDir', 'normal');
    xlabel('x'); ylabel('y');
    title(sprintf('Radial Modulus Field – Frame %d / %d', k, num_frames));
    colorbar;
    colormap(jet);
    clim([G0, G1]);
    drawnow;
    
    % Capture frame for GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end

%%
% surface_plot_modulus.m
% 3D Surface plot of G(r0) as height over radial domain

% Parameters
G0 = 100;
G1 = 1000;
l1 = 5;
l2 = 10;
a = 2;
n = 3;

% Grid setup (2D radial domain)
N = 200;
R = 15;  % Max radius
[x, y] = meshgrid(linspace(-R, R, N));
r0 = sqrt(x.^2 + y.^2);

% Compute G(r0)
G = zeros(size(r0));

% Region 1
G(r0 <= l1) = G0;

% Region 2
idx2 = (r0 > l1) & (r0 < l2);
f = (l2 - r0(idx2)) ./ (r0(idx2) - l1);
G(idx2) = G0 + (G1 - G0) .* (1 + f.^a).^((n - 1)/a);

% Region 3
G(r0 >= l2) = G1;

% Plot
figure;
surf(x, y, G, 'EdgeColor', 'none');
axis equal;
xlabel('x'); ylabel('y'); zlabel('G(r_0)');
title('Surface Plot of Radial Modulus G(r_0)');
colormap(jet);
colorbar;
clim([G0, G1]);
view(30, 40);  % angled view
lighting gouraud;
camlight headlight;
