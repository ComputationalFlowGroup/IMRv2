clc;
clear;
close;

addpath('src');

% equation options
R0 = 50e-6;
Req = R0/10;
tfin = 15E-6;
kappa = 1.4;
Lheat = 2.378193575129533e+04;
T8 = 298.15;
rho8 = 1064;
tvector = linspace(0,tfin,500);
radial = 1;
vapor = 1;
bubtherm = 0;
medtherm = 0;
masstrans = 0;
stress = [8,2];
G = 2500;
dGdhs = 1E+4;
for i = 1:length(stress)
varin = {'progdisplay',0,...
    'radial',radial,...
         'bubtherm',bubtherm,...
         'tvector',tvector,...
         'vapor',vapor,...
         'medtherm',medtherm,...
         'method',23,...
         'stress',stress(i),...
         'collapse',1,...
         'mu',1E-3,...
         'g',G,...
         'lambda1',1e-7,...
         'lambda2',0,...
         'alphax',1,...
         'dgdr0',dGdr0,...
         'r0',R0,...
         'req',Req,...
         'kappa',kappa,...
         't8',T8,...
         'rho8',rho8};
Pref = 101325;
Ca = Pref/G;
%[t,R,~] = m_imr_fd(varin{:},'Nt',70,'Mt',70);
[t,R,~] = m_imr_fd(varin{:},'Nt',50,'Mt',300); %mu = 1E-4, alphax = 1e-3
G = 2500; dGdhs = 1E+4; 

figure(1)
hold on;
plot(t,R,'x')
ylim([0 1.2])
legend('nHKV-gmlin-fd','nHKV-fd')
Lambda = R0/Req;
r_ffield_baseline


% figure(1)
% hold on;
% [t2,R2,~] = m_imr_spectral(varin{:},'Nt',12,'Mt',12);
% plot(t2,R2,'^')
% legend('nHKV-gmlin-sp','nHKV-sp')
% %legend('nHqKV-gmlin-sp')
end