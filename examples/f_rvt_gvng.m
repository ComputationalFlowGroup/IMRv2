% generate R v t curves
addpath('~/IMRv2/src/forward_solver/');
%addpath('./common/')

% options
kappa = 1.4;
T8 = 298.15;
rho8 = 1064;
mu = 0;
lambda1 = 0;
lambda2 = 0;
alphax = 0;
Pref = 101325;

R0 = 100e-6;
Req = R0/8; 
tfin = 75E-6;
%tfin = 1.25*R0*sqrt(rho8/Pref);
tvector = linspace(0,tfin,10000);

collapse = 0;
radial = 2;
vapor = 0;
bubtherm = 0;
medtherm = 0;
masstrans = 0;
stress = 1;

% graded parameters
graded = 1;
v_nc = 0.3; 
v_a = 2;
l1 = 1.1; 
l2 = 2.2;
gfun = 1;

% for loop
isgraded = [0 0 1];
% positive gradient: soft-to-stiff
G0_all = [1000 10000 1000]; 
G1_all = [0 0 10000];
Ca_G0 = Pref./G0_all;
Ca_G1 = Pref./G1_all;

figure;
hold on;
%colors = flipud(turbo(length(G0_all)));
colors = [1 0 0; %red
          0 0 0; %black
          0 1 0]; %green or 0.5 for forest green

for i = 1:length(isgraded)
    G0 = G0_all(i);
    G1 = G1_all(i);
    graded = isgraded(i);

varin = {'progdisplay',0,...
    'radial',radial,...
    'bubtherm',bubtherm,...
    'tvector',tvector,...
    'vapor',vapor,...
    'medtherm',medtherm,...
    'masstrans',masstrans,...
    'method',23,...
    'stress',stress,...
    'collapse',collapse,...
    'mu',mu,...
    'g',G0,...
    'graded',graded,...
    'g1',G1,...
    'l1',l1,...
    'l2',l2,...
    'v_a',v_a,...
    'v_nc',v_nc,...
    'gfun',gfun,...
    'lambda1',lambda1,...
    'lambda2',lambda2,...
    'alphax',alphax,...
    'r0',R0,...
    'req',Req,...
    'kappa',kappa,...
    't8',T8,...
    'rho8',rho8};
% generate R v t data
[t,R,~] = f_imr_fd(varin{:},'Nt',50,'Mt',150);

if i == length(G0_all)
    linestyle = '--';
else
    linestyle = '-';
end

plot(t,R,'Color',colors(i,:),'LineStyle',linestyle,'DisplayName',sprintf('G0 = %g, G1 = %g',G0,G1));
end
ylim([0 1.2]);
xlim([0 4]);
ylabel('$R$ / $R_0$', 'Interpreter', 'Latex', 'FontSize', 20); 
xlabel('$t$ / $t_c$', 'Interpreter','Latex', 'FontSize', 20);
set(gcf,'color','w');
%set(gca,'YScale','log');
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;
legend show;
hold off;

%%
