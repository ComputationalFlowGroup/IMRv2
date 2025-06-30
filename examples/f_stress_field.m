% function f_stress_field(R,t)
fig_tend = 12;
tickrange = ceil(min(t)) : floor(max(t));
lR = length(R);
lr_N = 500;
lr_length = 12;
nt = length(t);
%%
r_coord = ones(lR,lr_N).*linspace(0.1,lr_length,lr_N);
%r_coord = linspace(0.1,1,lr_N);

[xcon,ycon] = meshgrid(t,r_coord(1,:));
clevels = 1000;
%addpath('./common/')
[taurr1,taurr,r1,r2] = f_gradedstress(r_coord,R,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);
%[~,~,taurr1,taurr,r1,r2] = f_graded(stress,r_coord',R',Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);
maxtaurr = -max(max(abs(taurr))); %max(abs(taurr(:))); %1; %max(max(abs(taurr)));
ntaurr = taurr/maxtaurr;

% f contour figure
figure
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$\it{r}/R_o$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}/\mathrm{max}(\tau_{rr})$';
set(cbar,'TickLabelInterpreter','Latex','FontSize',18);
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -1];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
%clim([min(real(ntaurr(:))) max(real(ntaurr(:)))]);
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
plot(t,R,'LineWidth',3,'Color','k');
contourf(xcon,ycon,ntaurr',clevels,'edgecolor','none')
hold on;
plot(t,r1','b--','LineWidth',2,'DisplayName','l1')
plot(t,r2','r--','LineWidth',2,'DisplayName','l2')
ylim([0 fig_tend])
hold off;


%% now in logspace
lr_length = 4;
logr_coord = ones(lR,lr_N).*logspace(-1,lr_length,lr_N);

[xcon,ycon] = meshgrid(t,logr_coord(1,:));
ycon = log10(ycon);
clevels = 1000;
addpath('./')
[logtaurr1,logtaurr,logr1,logr2] = f_gradedstress(logr_coord,log10(R),Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);
maxlogtaurr = -max(max(logtaurr));
logntaurr = logtaurr/maxlogtaurr;

% f contour figure
figure
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$log(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}/\mathrm{max}(\tau_{rr})$';
set(cbar,'TickLabelInterpreter','latex','FontSize',18);
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -1];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
%clim([min(logntaurr(:)) max(logntaurr(:))]);
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
contourf(xcon,ycon,logntaurr',clevels,'edgecolor','none')
hold on;
plot(t,log10(logr1)','b--','LineWidth',2,'DisplayName','l1')
plot(t,log10(logr2)','r--','LineWidth',2,'DisplayName','l2')
%ylim([log10(min(r_coord(:))), log10(max(r_coord(:)))])
%ylim([0 lr_length])
hold off;