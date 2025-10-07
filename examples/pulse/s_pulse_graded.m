close all; clear; clc;
% pulse
addpath('~/IMRv2/src/forward_solver/');

% Hersey parameters
wave_type = 1; %1 gaussian, 0 impulse, 2 histotripsy (FU)
Pref = 101325;
rho8 = 1064;
% Initital Radii 
R0 = 0.5E-6;
Req = R0 + 1E-10;

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
mu = 0;%5E-3;
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
graded = 1;
v_nc = 0.3; 
v_a = 2;
l1 = 2*R0; %trying dim
l2 = 4*R0;
G0 = 1E3;
G1 = 10E3;

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
r_fig;

%%
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

p_t = pA * exp(-((tvector - DT)/TW).^2);
figure;
plot(tvector, p_t);
xlabel('Time [s]');
ylabel('Pressure [Pa]');
title('Gaussian Pressure Pulse');

lambda = (Rdim.^3 - Req^3).^(1/3);
Lambda = Rdim./Req;
l = ((lambda.^3 - 1) ./ (Lambda.^3 - 1)).^(1/3);
f = (l2.*l - 1) ./ (1 - l1.*l);
m = (1+f.^v_a).^((v_nc-1)/v_a);

% post_processing
r_fig;
%r_ffield;
