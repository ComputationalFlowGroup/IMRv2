%% bicomposite stress field
addpath('./examples/')
clc; close all;
%tickrange = ceil(min(t)) : floor(max(t));
tickrange = linspace(ceil(min(t)),floor(max(t)),8);
lR = length(R);
lr_N = lR;
lr_length = 5;
nt = length(t);
r_coord = ones(lR,lr_N).*linspace(0,lr_length,lr_N);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
Ca = Pref/G0;
Ca1 = Pref/G1;
[taurr_near,taurr_far,mask_near,mask_far,r1] = f_bicomp_stress(r_coord,R,Req/R0,Ca,Ca1,l1);
% 
taurr = NaN(size(r_coord));
taurr(mask_near) = taurr_near;
taurr(mask_far) = taurr_far;
ntaurr = taurr;

maxtaurr = max(max(abs(taurr))); 
ntaurr = taurr/maxtaurr;

% diverging color map
rgb = [ ...0
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;

% Interpolate to 200 colors
n_col = 200;
old_x = linspace(0,1,size(rgb,1));
new_x = linspace(0,1,n_col);
rgb200 = interp1(old_x, rgb, new_x);
colormap(rgb200);
cbar = colorbar;
clevels = n_col;

% contour figure
figure(1);
hold on;
xlabel('$t / t_{\mathrm{c}}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$R/R_{\mathrm{max}}$','Interpreter','latex','FontSize',24);
cbar.Label.String = '$\tau_{rr}/\max(\tau_{rr})$';
%cbar.Label.String = '$\frac{\tau_{rr}}{\mathrm{max}(\tau_{rr})}$';
set(cbar,'TickLabelInterpreter','latex','FontSize',18);
cbar.Label.Rotation = 0;
%cbar.Label.Interpreter = 'latex';
pos = get(cbar, 'Position');
cbar.Label.Position = [pos(1) + pos(3), pos(2) - 0.15,0];
pos = get(cbar,'Position');
%clim([-1 1]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
ya = gca;
ya.TickLength = [.015 .015];
ya.LineWidth = 1.5;
box on;
plot(t,R,'LineWidth',3,'Color','k');
contourf(xcon,ycon,ntaurr',clevels,'edgecolor','none')
hold on;
plot(t,r1','c--','LineWidth',2,'DisplayName','l1')
ylim([0 lr_length])
hold off;


function [taurr_near,taurr_far,mask_near,mask_far,r1] = f_bicomp_stress(rcoord,R,Req,Ca,Ca1,l1)

r1 = (l1^3 + R.^3 - Req^3).^(1/3);
%x1 = @(x) (el1.^3+lam(x).^3-1).^(1/3);
% reference coordinate calculation
r0coord = (rcoord.^3 - R.^3 + Req.^3).^(1/3);
r0shift = (rcoord.^3 + Req^3 - (R+0.001).^3).^(1/3);
% removing the data within the bubble and slightly away from the bubble wall
r0coord(r0coord < r0shift) = NaN;
% 
% % masks for each region
mask_near = r0coord < l1 + 1e-12;
mask_far = r0coord >= l1;
% 
% function for graded region
tau = @(r,r0) (2/3)*(r0.^4 ./ r.^4 - r.^2 ./ r0.^2);

% near field
r0near = r0coord(mask_near);
rnear = rcoord(mask_near); %????
% diff_rnear = abs(rnear - r0near);
% fprintf('Min |r - r0|: %g, Max |r - r0|: %g\n', min(diff_rnear(:)), max(diff_rnear(:)))
% taurr(mask_near) = 1/Ca *tau(rnear,r0near);
% taurr_near = taurr(mask_near);
taurr_near = 1/Ca *tau(rnear,r0near);

% far field
r0far = r0coord(mask_far);
rfar = rcoord(mask_far);
% diff_rfar = abs(rfar - r0far);
% fprintf('Min |r - r0|: %g, Max |r - r0|: %g\n', min(diff_rfar(:)), max(diff_rfar(:)))
% taurr(mask_far) = 1/Ca1 * tau(rfar,r0far);
% taurr_far = taurr(mask_far);
taurr_far = 1/Ca1 * tau(rfar,r0far);

taurr_near(isinf(taurr_near)) = 0;

end
