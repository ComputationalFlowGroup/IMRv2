clc;
clear;
%close;

addpath(genpath('src'));

% equation options
R0 = 230e-6;
Req = 34e-6;
tc = R0*sqrt(1048/101325);
tfin = 2.5*tc; %160E-6;
kappa = 1.4;
Lheat = 2.378193575129533e+04;
T8 = 298.15;
rho8 = 1048;
ST = 0.072;
alphax = 0;
mu = 1e-1;
Gqs = 8.66e3;
rho = 1048;
modes = 7;
epnm0 = 0.01.*modes;
epnmd0 = 0.*epnm0;
  % equation options
    kappa = 1.4;
    T8 = 298.15;


    % simulation equation options
    radial = 2;
    vapor = 1;
    collapse = 0;
    bubtherm = 1;
    medtherm = 0;
    masstrans = 1;
    perturbed = 1;
    stress = 2;

    % combine all inputs into varin
    varin = {'progdisplay',0,...
        'radial',radial,...
        'bubtherm',bubtherm,...
        'perturbed', perturbed, ...
        'tvector',linspace(0, tfin, round(1e7*tfin)),...
        'vapor',vapor,...
        'medtherm',medtherm,...
        'masstrans',masstrans,...
        'method',45,...
        'stress',stress,...
        'collapse',collapse,...
        'mu',mu,...
        'alphax', alphax, ...
        'g',Gqs,...
        'lambda1',0,...
        'lambda2',0,...
        'surft', ST, ...
        'r0',R0,...
        'req',Req,...
        'kappa',kappa,...
        't8',T8,...
        'rho8', rho, 'modes', modes,'epnm0', epnm0, ...
        'epnmd0', epnmd0, 'reltol', 1e-5, 'abstol', 1e-7, 'Nt', 100};
tic
if perturbed
    [tfd,Rfd,Rfddot,Pfd,Tfd,Tmfd,kvfd, epnm, epnmd] = f_imr_fd(varin{:},'Nt',100,'Mt',100);
else
    [tfd,Rfd,Rfddot,Pfd,Tfd,Tmfd,kvfd] = f_imr_fd(varin{:},'Nt',100,'Mt',100);
end
toc

%%
figure(1)
hold on;
plot(tfd,Rfd,'-');
% plot(tsp,Rsp,'r^');
yline(Req/R0)
ylim([0 1]);
hold on
%%
addpath ../../../Common_functions/

for i = 1:length(modes)
    x_all{i} = struct('time', tfd.*tc, 'n', modes(i), 'Ro', R0, 'sig', ST);
end

eppred = pert_relax_multi_dataset_plot([Gqs, mu], x_all);

%%

figure(2)
plot(tfd, epnm,'-.')
hold on
for i = 1:length(eppred)
    plot(tfd, epnm0(i).*eppred{i}, '-')
end
ylim([-.1 .1])
box on
grid on
