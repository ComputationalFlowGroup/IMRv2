%% bicomposite stress fields
% function f_stress_field(R,t)
addpath('./examples/')
clc;
%tickrange = ceil(min(t)) : floor(max(t));
tickrange = linspace(ceil(min(t)),floor(max(t)),8);
lR = length(R);
lr_N = lR;
lr_length =3;
nt = length(t);
r_coord = ones(lR,lr_N).*linspace(0,lr_length,lr_N);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
%addpath('./common/')
%[taurr,r1,r2] = f_graded_stress(r_coord,R,Req/R0,Ca,Ca1,l1/R0,l2/R0,v_nc,v_a);
[taurr_near,taurr_far,mask_near,mask_far,r1] = f_bi_stress(r_coord,R,Req/R0,Ca,Ca1,l1);

taurr = NaN(size(r_coord));
taurr(mask_near) = taurr_near;
taurr(mask_far) = taurr_far;

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
% inverted diverging colormap
% rgb200 = rgb200(end:-1:1,:);
colormap(rgb200);
cbar = colorbar;
clevels = n_col;

% contour figure
figure(1);
hold on;
xlabel('$t / t_{\mathrm{c}}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$r/R_{\mathrm{max}}$','Interpreter','Latex','FontSize',24);
%cbar.Label.String = '$\tau_{rr}/\max(\tau_{rr})$';
cbar.Label.String = '$\frac{\tau_{rr}}{\mathrm{max}(\tau_{rr})}$';
set(cbar,'TickLabelInterpreter','Latex','FontSize',18);
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
pos = get(cbar, 'Position');
cbar.Label.Position = [pos(1) + pos(3), pos(2) - 0.15,0];
pos = get(cbar,'Position');
clim([0 1]);
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
% contourf(xcon,ycon,taurr_near,clevels,'edgecolor','none')
% contourf(xcon,ycon,taurr_mid,clevels,'edgecolor','none')
% contourf(xcon,ycon,taurr_far,clevels,'edgecolor','none')
hold on;
plot(t,r1','c--','LineWidth',2,'DisplayName','l1')
ylim([0 lr_length])
hold off;

function [taurr_near,taurr_far,mask_near,mask_far,r1,r2] = f_bi_stress(rcoord,R,Req,Ca,Ca1,l1)

r1 = (l1^3 + R.^3 - Req^3).^(1/3);

% reference coordinate calculation
r0coord = (rcoord.^3 - R.^3 + Req.^3).^(1/3);
r0shift = (rcoord.^3 + Req^3 - (R+0.001).^3).^(1/3);
% removing the data within the bubble and slightly away from the bubble wall
r0coord(r0coord < r0shift) = NaN;

% masks for each region
mask_near = r0coord <= l1 - 1e-12;
mask_far = r0coord >= l1 + 1e-12;

% stress
tau = @(r,r0) (2/3)*(r0.^4 ./ r.^4 - r.^2 ./ r0.^2);

% near field
r0near = r0coord(mask_near);
rnear = rcoord(mask_near); %????
diff_rnear = abs(rnear - r0near);
fprintf('Min |r - r0|: %g, Max |r - r0|: %g\n', min(diff_rnear(:)), max(diff_rnear(:)))
% taurr(mask_near) = 1/Ca *tau(rnear,r0near);
% taurr_near = taurr(mask_near);
taurr_near = 1/Ca *tau(rnear,r0near);

% far field
r0far = r0coord(mask_far);
rfar = rcoord(mask_far);
diff_rfar = abs(rfar - r0far);
fprintf('Min |r - r0|: %g, Max |r - r0|: %g\n', min(diff_rfar(:)), max(diff_rfar(:)))
% taurr(mask_far) = 1/Ca1 * tau(rfar,r0far);
% taurr_far = taurr(mask_far);
taurr_far = 1/Ca1 * tau(rfar,r0far);


% % zero out regions larger than the near field
% r0near(r0near > l1) = 0;
% taurr1 = (2/(3*Ca))*((r0near.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0near.^2));
% % zero out the stress above and below the near field
% taurr1(isinf(taurr1)) = 0;
% % far field
% r0far = r0coord;
% % zero out regions less than the far field
% r0far(r0far < l2) = 0;
% taurr3 = (2/(3*Ca1))*((r0far.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0far.^2));
% % zero out the stress below the far field
% taurr3(isinf(taurr3)) = 0;
% % compute the graded region coordinate
% r0mid = r0coord - r0far - r0near;
% % graded stress
% taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).* ...
%     ((r0mid.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0mid.^2));
% % removing the negative infinities
% taurr2(isinf(taurr2)) = 0;
% % sum of the near, far, and graded stress fields
% %taurr = taurr1 + taurr2 + taurr3;

% figure
% strain = abs(rcoord - r0coord) ./ r0coord;
% imagesc(strain); colorbar; title("Relative radial strain");

end
