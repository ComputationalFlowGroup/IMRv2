fig_tend = 55;
tickrange= 0:5:fig_tend;
lma = 'r';
lm ='k-';

figure(1)  
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{R} / R_o$', 'Interpreter', 'Latex', 'FontSize', 20); 
%set(leg,'Interpreter','latex','Location','northeast');
%set(leg,'FontSize',18);
set(gcf,'color','w'); %Changes background to white 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
box on;
plot(t,R, lm,'LineWidth',2); 
%saveas(gcf,'./figs/baseline/R_T','png')

figure(2)  
hold on
xlabel('$\it{t}$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$\it{R}$', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gcf,'color','w');
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
box on;
plot(tdim,Rdim, lm,'LineWidth',2); 

figure(3)
hold on
xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$m$', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gcf,'color','w');
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
% ylim([-1.5 1.5]*1E-5)
box on;
lambda = (R.^3 - Req^3).^(1/3);
Lambda = R./Req;
l = ((lambda.^3 - 1) ./ (Lambda.^3 - 1)).^(1/3);
f = (l2.*l - 1) ./ (1 - l1.*l);
m = (1+f.^v_a).^((v_nc-1)/v_a);
plot(t,f,lma,'LineWidth',2);
hold on;
plot(t,m,lm,'LineWidth',2); 
hold off;
%%
% Q: HOW TO EXTRACT r0 and r from R0 and R?
lR = length(R);
lr_N = lR;
lr_length = 3;
nt = length(t);
r_coord = ones(lR,lr_N).*linspace(0,lr_length,lr_N);
r0coord = (r_coord.^3 - R.^3 + Req.^3).^(1/3);
figure(4)
hold on
xlabel('\it{r_0} / $R_0$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$m$', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gcf,'color','w');
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(tickrange)
xlim([0 fig_tend])
% ylim([-1.5 1.5]*1E-5)
box on;
lambda = (R.^3 - Req^3).^(1/3);
Lambda = R./Req;
l = ((lambda.^3 - 1) ./ (Lambda.^3 - 1)).^(1/3);
f = (l2.*l - 1) ./ (1 - l1.*l);
m = (1+f.^v_a).^((v_nc-1)/v_a);
plot(r0coord,m,lm,'LineWidth',2); 


% figure(4)
% hold on
% xlabel('\it{t} / $t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
% ylabel('$\tau_{rr}|_R / p_{\infty}$', 'Interpreter', 'Latex', 'FontSize', 20); 
% leg = legend('$\mu_{water}$','$\mu_{blood,\infty}$', '$\mu_{blood}$', '$\mu_{blood,0}$' ); 
% set(leg,'Interpreter','latex','Location','northeast');
% set(leg,'FontSize',18);
% set(gcf,'color','w'); %Changes background to white 
% set(gca, 'FontName', 'Times', 'FontSize',20); 
% set(gca,'TickLabelInterpreter','latex')
% xa = gca;
% xa.TickLength = [.03 .03];
% xa.LineWidth = 1.5;
% xticks(tickrange)
% xlim([0 fig_tend])
% box on;
% tau = 2*(m_slope./DRe+1./Re8).*gammadot_R;
% plot(t,tau,lm,'LineWidth',2);
% saveas(gcf,'./figs/baseline/tau_T','png')