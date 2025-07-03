% function f_stress_field(R,t)
clc;
%tickrange = ceil(min(t)) : floor(max(t));
tickrange = linspace(ceil(min(t)),floor(max(t)),8);
lR = length(R);
lr_N = lR;
lr_length = 3;
nt = length(t);
r_coord = ones(lR,lr_N).*linspace(0,lr_length,lr_N);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
%addpath('./common/')
[taurr,r1,r2] = f_graded_stress(r_coord,R,Req/R0,Ca,Ca1,l1/R0,l2/R0,v_nc,v_a);
maxtaurr = -max(max(abs(taurr))); 
ntaurr = taurr/maxtaurr;

% diverging color map
rgb = [ ...
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
hold on;
plot(t,r1','c--','LineWidth',2,'DisplayName','l1')
plot(t,r2','m--','LineWidth',2,'DisplayName','l2')
ylim([0 lr_length])
hold off;