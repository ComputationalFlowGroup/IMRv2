% graded parameters
graded = 1;
v_nc = 0.3; 
v_a = 2;
R0 = 100e-6;
l1 = 1.2*R0; %trying dim
l2 = 2*R0;
G0 = 1E3;
G1 = 10E3;

eps = 1E-8;
r0_l1 = l1+eps;
r0 = linspace(r0_l1, l2-eps, 500);
x = r0./R0;
f = (l2 - r0) ./ (r0 - l1);
m = (1+ f.^v_a).^((v_nc - 1)/v_a);

f_l1 = (l2-r0_l1) / (r0_l1 - l1);
m_l1 = (1+ f_l1^v_a)^((v_nc-1)/v_a);

f_l2 = (l2 - l2) /(l2-l1);
m_l2 = (1+f_l2^v_a)^((v_nc-1)/v_a);

figure;
hold on;
xlabel('$r_0$ / $R_0$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$m$', 'Interpreter','Latex', 'FontSize', 20);
set(gcf,'color','w');
%set(gca,'YScale','log');
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
plot(x,m,'k','LineWidth',2);
% xval = [min(x), max(x)];
% plot(xval,[m_l1, m_l1],'r--')
% plot(xval,[m_l2, m_l2],'b--')
xline(l1/R0,'r--','LineWidth',2);
xline(l2/R0,'b--','LineWidth',2);
%plot(l1/R0,m_l1,'r--')
%plot(l2/R0,m_l2,'b--')
hold off;


%r0_tanh = linspace(l1-0.2*(l2-l2),l2+0.2*(l2-l1),500);
r0_tanh = linspace(l1,l2,500);
b = 0.5 *(l1+l2) / R0;
f_tanh = (r0 - (l1+l2)/2) / (l2-l1);
m_tanh = 0.5*(1+tanh(2*b*f_tanh));
x_tanh = r0_tanh/R0;

figure;
hold on;
xlabel('$r_0$ / $R_0$', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('$m$', 'Interpreter','Latex', 'FontSize', 20);
set(gcf,'color','w');
%set(gca,'YScale','log');
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
plot(x_tanh,m_tanh,'k','LineWidth',2);
xline(l1/R0,'r--');
xline(l2/R0,'b--');
%xlim([l1+eps,l2+eps])

% Define b values to test (from 0.1*(l1+l2) to 1*(l1+l2))
%b_vals = linspace(0.1*(l1+l2)/R0, 1*(l1+l2)/R0, 10);  % 5 steps, adjust as needed
b_fac = linspace(0.1,1,10);
b_vals = b_fac*(l1+l2)/R0;
figure; 
hold on;
color = turbo(length(b_vals));
for i = 1:length(b_vals)
    b = b_vals(i);
    m_tanh = 0.5 * (1 + tanh(2 * b * f_tanh));
    plot(x_tanh, m_tanh, 'Color', color(i,:), 'LineWidth', 2, 'DisplayName', ['b = ' num2str(b_fac(i),'%.1f')]);
end
xline(l1 / R0, 'k--','HandleVisibility','off');
xline(l2 / R0, 'k--','HandleVisibility','off');
xlim([l1 / R0, l2 / R0]);
xlabel('$r_0 / R_0$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$m$', 'Interpreter', 'latex', 'FontSize', 20);
lgd = legend('show', 'Location', 'best');
%set(gca, 'FontName', 'Times', 'FontSize', 20);
lgd.FontSize = 10;
box on;
hold off;

% xticks(tickrange)
% xlim([0 fig_tend])
% % ylim([-1.5 1.5]*1E-5)=
% lambda = (R.^3 - Req^3).^(1/3);
% Lambda = R./Req;
% l = ((lambda.^3 - 1) ./ (Lambda.^3 - 1)).^(1/3);
% f = (l2.*l - 1) ./ (1 - l1.*l);
% m = (1+f.^v_a).^((v_nc-1)/v_a);
% plot(t,f,lma,'LineWidth',2);
% hold on;
% plot(t,m,lm,'LineWidth',2); 
% hold off;
