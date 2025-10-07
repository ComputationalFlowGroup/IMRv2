
isgraded = [0 1 0];
% positive gradient: soft-to-stiff
G0_mat = [1000 1000 10000];
G1_mat = [0 10000 0];

% pulse
addpath('~/IMRv2/src/forward_solver/');

% Hersey parameters
wave_type = 1; %1 gaussian, 0 impulse, 2 histotripsy (FU)
Pref = 101325;
rho8 = 1064;
% Initital Radii
R0 = 0.5E-6;
Req = R0;

% computing nondim params for waveform
P8 = Pref;
Uc = sqrt(P8/rho8);
t0 = R0/Uc;
tw = 1;
dt = 10;
TW = tw*t0; % Gaussian width (s)
DT = dt*t0; % delay (s)

% nondimensional pressure amplitude (Pa) Hersey
pA = 1.15E6;
% f frequency (rad/s) Hersey
omega = 2E6;
% power shift /exponent for waveform
mn = 0;
tend = 30 * TW;
tvector = linspace(0,tend,1000);

% options
kappa = 1.4;
T8 = 298.15;
mu = 0; %5E-3;
lambda1 = 0;
lambda2 = 0;
alphax = 0;

collapse = 0;
radial = 2;
vapor = 0;
bubtherm = 0;
medtherm = 0;
masstrans = 0;
stress = 1;

% graded parameters
v_nc = 0.3;
v_a = 2;
l1 = 1.2;
l2 = 2.2;
for i = 1:length(isgraded)
    graded = isgraded(i);
    G0 = G0_mat(i);
    G1 = G1_mat(i);
    %Ca = G0/Pref;
    %Ca1 = G1/Pref;

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
        'lambda1',lambda1,...
        'lambda2',lambda2,...
        'alphax',alphax,...
        'r0',R0,...
        'req',Req,...
        'kappa',kappa,...
        't8',T8,...
        'rho8',rho8,...
        'wave_type',wave_type,...
        'pA',pA,...,
        'omega',omega,...
        'TW',TW,...
        'DT',DT,...
        'mn',mn
        };

    % % generate R v t data
    [t,R,U,P,~] = f_imr_fd(varin{:},'Nt',16,'Mt',64);
    hold on
    fig_tend = 30;
    tickrange= 0:5:fig_tend;

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
plot(t,R,'LineWidth',2);
end