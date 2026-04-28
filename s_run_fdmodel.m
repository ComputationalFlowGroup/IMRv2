clc;
clear;
%close;

addpath(genpath('src'));

% equation options
% ------- Initial condition and equilibrium ----------%
R0 = 50e-6;
Req = R0/1;

% ------- Material Properties ------------------------%
kappa = 1.4;
T8 = 298.15;
rho8 = 1048;
mu = 1e-2;
Gelastic = 1e4;

% ------- Simulation settings ------------------------%
radial = 2;
vapor = 1;
collapse = 0;
bubtherm = 1;
medtherm = 0;
masstrans = 1;
stress = 2;

% --------- Ultrasound settins -----------------------%
pa = 100e3;
omega = 2*pi*150e3;
wavetype = 4;

% ------ Simulation time ---------------------------- %
tc = R0*sqrt(rho8/101325);
tfin = 25*tc;
tvector = linspace(0,tfin,1000);
varin = {'progdisplay',0,'radial',radial,'bubtherm',bubtherm,'tvector',tvector,...
         'vapor',vapor,'medtherm',medtherm,'masstrans',masstrans,'method',23,...
         'stress',stress,'collapse',collapse,'mu',mu,'g',Gelastic,'lambda1',0e-7,...
         'lambda2',0,'alphax',10,'r0',R0,'req',Req,'kappa',kappa,'t8',T8,...
         'rho8',rho8, 'pa',pa 'omega', omega, 'wave_type', wavetype};

[t,R,Rd,~,~,~,~,Rddot] = f_imr_fd(varin{:},'Nt',75);

%figure
hold on
plot(t, R)

