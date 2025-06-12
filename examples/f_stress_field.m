% function f_stress_field(R,t)
fig_tend = 4;
%tickrange= 0:2:fig_tend;
tickrange = ceil(min(t)) : floor(max(t));
lR = length(R);
%lr_max = 100;
lr_N = 500;
lr_length = 5;
nt = length(t);
r_coord = ones(lR,lr_N).*logspace(-1,lr_length,lr_N);
%r_coord = ones(lR,lr_N).*linspace(0.1,3, lr_N);
%R_coord = (ones(lR,lr_N).*R);

[xcon,ycon] = meshgrid(t,r_coord(1,:));
ycon = log10(ycon);
clevels = 1000;
addpath('./')
[taurr1,taurr,r1,r2] = f_gradedstress(r_coord',R',Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);
maxtaurr = 1; %max(max(taurr)); %max(abs(taurr(:))); %1; %max(max(abs(taurr)));
ntaurr = taurr/maxtaurr;

% f contour figure
figure
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}/p_{\infty}$';
set(cbar,'TickLabelInterpreter','latex','FontSize',18);
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -21];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
clim([min(ntaurr(:)) max(ntaurr(:))]);
%clim([-0.75 .25])
%xlim([0 fig_tend]);
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
plot(t, log10(R),'LineWidth',3,'Color','k');
plot(t, R,'LineWidth',3,'Color','k');
plot(t, R,'LineWidth',3,'Color','k');
contourf(xcon,ycon,log(abs(ntaurr)),clevels,'edgecolor','none')
%contourf(xcon,ycon,ntaurr,clevels,'edgecolor','none')
hold on;
plot(t,log10(r1'),'b--','LineWidth',2,'DisplayName','l1')
plot(t,log10(r2'),'r--','LineWidth',2,'DisplayName','l2')
%ylim([log10(min(r_coord(:))), log10(max(r_coord(:)))])
% end