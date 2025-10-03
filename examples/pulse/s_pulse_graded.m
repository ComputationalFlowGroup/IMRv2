%close all; clear; clc;
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
TW = tw*t0;
DT = dt*t0;


% nondimensional pressure amplitude (Pa) Hersey
pA = 5E5; %500E3*Pref; %.5*Pref; %500E3*Pref;
% f frequency (rad/s) Hersey
omega = 2E6; 
% Gaussian width (s)
%TW = 3E-6; %1.01*R0*sqrt(rho8/Pref);
% delay (s)
%DT = 10*TW; %10*TW
% power shift /exponent for waveform
mn = 0;
tstart = 8E-6;
tend = 30 * TW; %1.5E-6;
%tvector = linspace(tstart,tend,1000);
tvector = linspace(0,tend,1000);

% options
kappa = 1.4;
T8 = 298.15;
mu = 5E-3;
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
graded = 0;
v_nc = 0.3; 
v_a = 2;
l1 = 1.2; 
l2 = 2.2;
G0 = 1E3;
G1 = 10E3;
Ca = G0/Pref;
Ca1 = G1/Pref;

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
    'g',0,...
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
% dimensional R
Rdim = R.*R0;
tchar = sqrt(rho8/Pref)*R0;
tdim = t.*tchar;


figure
hold on;
%plot(tdim,Rdim,'bx')
semilogy(t,R)
xlabel('t/tc')
ylabel('R/R0')
hold off

figure 
plot(tdim,P)
xlabel('t [s]')
ylabel('pressure')

figure
plot(t,P)
xlabel('t/t0')
ylabel('pressure')

figure
tplot = linspace(0,max(tdim),1000);
TW = 3E-6;
Pplot = pA * exp(-((tplot - DT)./TW).^2);
plot(tplot,Pplot)
xlabel('Time [\mu s]')
ylabel('Pressure [Pa]')

lambda = (Rdim.^3 - Req^3).^(1/3);
Lambda = Rdim./Req;
l = ((lambda.^3 - 1) ./ (Lambda.^3 - 1)).^(1/3);
f = (l2.*l - 1) ./ (1 - l1.*l);
m = (1+f.^v_a).^((v_nc-1)/v_a);

% post_processing
%lm = 'k';
r_fig;
%r_ffield;


%%
% Baseline values
A = 1E6; 
f = 2E6;
R0= 0.5E-6;

% RUNNING FREQUENCY VARIATION
p_range = A*2.^(-2.5:0.5:2.5);
p_l_mu_max = zeros(length(p_range),1);
p_l_mu_avg = zeros(length(p_range),1);
R_max = zeros(length(p_range),1);
R_ave = zeros(length(p_range),1);
tfactor = p_range./p_range(1);

t0 = tstart;
Pinf = Pf8;
xrange = (Pinf./p_range*f*t0)';
% figure(2); hold on;
for i=1:length(p_range)
    Pext_Amp_Freq = [p_range(i), f]; 
    tend = 0.5E-6;%*tfactor(i);
    [ t , R ,U ,P, T,C, Tm,tdel,Tdel,Cdel] = m_cavitation...
    (tend,R0,NT,NTM,Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,...
    Dim,comp,vmaterial,vmodel);
    vec = f_char_length(R,U,v_lambda);
    p_l_mu_max(i) = max(vec);
    p_l_mu_avg(i) = mean(vec); 
    R_max(i) = max(R);
    R_ave(i) = mean(R);
%     figure(2)
%     plot(t,R/R_max(i))
end

figure(1)
hold on;
plot(xrange,p_l_mu_max,'or','MarkerFaceColor','r','MarkerSize',8);
plot(xrange,R_max,'ob','MarkerFaceColor','b','MarkerSize',8);
% plot(p_range/Pinf,p_l_mu_avg,'^r','MarkerFaceColor','r','MarkerSize',8);
% plot(p_range/Pinf,R_ave,'^b','MarkerFaceColor','b','MarkerSize',8);
xlabel('$(R_{o}/R_{o,b})(f R_o\sqrt{\rho_{\ell}/p_{\infty}})(p_{\infty}/A)$', 'Interpreter', 'Latex', 'FontSize', 20); 
% ylabel('$ \ell $', 'Interpreter', 'Latex', 'FontSize', 20); 
set(gca, 'YScale', 'log','XScale', 'log')
set(gcf,'color','w'); 
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex');
ylim([1E0 1E6])
xlim([1E-3 1E-1])
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
box on;


%% sandbox
t = linspace(0, 30e-6, 1000);  % longer time
TW = 3e-6;
DT = 10e-6;
pA = 500e3 / 101325 / 10; % normalized pressure amplitude

p_t = pA * exp(-((t - DT)/TW).^2);

figure;
plot(t*1e6, p_t);
xlabel('Time [\mus]');
ylabel('Pressure');
title('Gaussian Pressure Pulse');


