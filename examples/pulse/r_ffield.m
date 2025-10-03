% fig_tend = 30;
% tickrange= [0:5:fig_tend];
% lR = length(R);
% lr_max = 100;
% lr_N = 200;
% lr_length = 4;
% r_coord = ones(lR,lr_N).*logspace(-0.5,lr_length,lr_N);
% calculating the shear
% varsigmadot_r = f_gammadot_r(r_coord,R,U,lR,lr_N);
% f_r = sf_carreau(v_nc,v_lambda,varsigmadot_r);
% [xcon,ycon] = meshgrid(t,r_coord(1,:));
% ycon = log10(ycon);
% clevels = 50;
% feps_r = f_f_filter(f_r,lR,lr_N);
% tau_r = 2*(feps_r./DRe+1./Re8).*varsigmadot_r;
% min_tau_r = min(min(tau_r))
% max_tau_r = max(max(tau_r))
% ntau_r = tau_r/max_tau_r;

%s_graded_stress in examples
addpath('./examples')
tickrange = linspace(ceil(min(t)),floor(max(t)),8);
lR = length(R);
lr_N = lR;
lr_length = 3;
nt = length(t);
r_coord = ones(lR,lr_N).*linspace(0,lr_length,lr_N);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
%addpath('./common/')
%[taurr,r1,r2] = f_graded_stress(r_coord,R,Req/R0,Ca,Ca1,l1/R0,l2/R0,v_nc,v_a);
[taurr_near,taurr_mid,taurr_far,mask_near,mask_grad,mask_far,r1,r2] = f_graded_stress(r_coord,R,Req/R0,Ca,Ca1,l1,l2,v_nc,v_a);

taurr = NaN(size(r_coord));
taurr(mask_near) = taurr_near;
taurr(mask_grad) = taurr_mid;
taurr(mask_far) = taurr_far;

maxtaurr = max(max(abs(taurr))); 
ntaurr = taurr/maxtaurr;

% f contour figure
figure(5)
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$m(\dot{\varsigma})$';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -0.04]; 
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([0 1]);
xlim([0 fig_tend]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
box on;
%plot(t, log10(R),lm,'LineWidth',3); 
plot(t, R,lm,'LineWidth',3); 
contourf(xcon',ycon',f_r,clevels,'edgecolor','none')
%saveas(gcf,'./figs/baseline/fcon_T','png')

% shear stress contour figure
figure(6)
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}/\mathrm{max}(\tau_{rr})$';
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [1.5*pos(1) -1.2];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
caxis([-1 1])
xlim([0 fig_tend]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
box on;
%plot(t, log10(R),lm,'LineWidth',3); 
plot(t,R,lm,'LineWidth',3); 
contourf(xcon',ycon',ntau_r,clevels,'edgecolor','none');
%saveas(gcf,'./figs/baseline/taucon_T','png')