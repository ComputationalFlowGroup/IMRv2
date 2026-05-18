clc;
clear;
close;

addpath(genpath('src'));

% equation options
R0 = 20*1.25e-6;
Req = 20e-6;
tc = R0*sqrt(1048/101325);
tfin = 10*tc; %160E-6;
kappa = 1.4;
Lheat = 2.378193575129533e+04;
T8 = 298.15;
rho8 = 1048;
ST = 0.04;
alphax = 0;
mu = 1e-3;
Gqs = 7.81e3;
rho = 1048;
modes = 8;
epnm0 = 0.1;
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
        'pertmod', 0, ...
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
plot(tfd, epnm,'-.', LineWidth=6)
% plot(tsp,Rsp,'r^');
yline(Req/R0)
ylim([0 1]);
hold on
%%
addpath ../../Common_functions/
addpath cmap\


Cmap = viridis(5);

for i = 1:length(modes)
    x_all{i} = struct('time', tfd.*tc, 'n', modes(i), 'Ro', R0, 'sig', ST, 'dim', 0);
end



eppred = pert_relax_multi_dataset([Gqs, mu], x_all);



figure(2)
plot(tfd, epnm,'-', 'Color', Cmap(1, :) , LineWidth=3)
hold on
plot(tfd, epnm0(i).*eppred, '--', 'Color', Cmap(2,:), LineWidth=3)
% ylim([-.1 .1])


% equation options
R0 = 20*1.25e-6;
Req = 20e-6;
tc = R0*sqrt(1048/101325);
tfin = 10*tc; %160E-6;
kappa = 1.4;
Lheat = 2.378193575129533e+04;
T8 = 298.15;
rho8 = 1048;
ST = 0.04;
alphax = 0;
mu = 1e-2;
Gqs = 7.81e3;
rho = 1048;
modes = 8;
epnm0 = 0.1;
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
        'pertmod', 0, ...
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
plot(tfd, epnm,'-','Color', Cmap(3,:), LineWidth=3)

for i = 1:length(modes)
    x_all{i} = struct('time', tfd.*tc, 'n', modes(i), 'Ro', R0, 'sig', ST, 'dim', 0);
end

eppred = pert_relax_multi_dataset([Gqs, 1e-2], x_all);
plot(tfd, epnm0(i).*eppred, '-.','Color', Cmap(4,:), LineWidth=3)



xlabel('$t/t_c^{\rm LIC}$', 'Interpreter','latex', FontSize=16)
ylabel('$\epsilon_7(t)$', 'Interpreter','latex', FontSize=16)
legend('IMR forward solve $\mu = 1$ mPa s', 'Exact $\mu = 1$ mPa s','IMR forward solve $\mu = 10$ mPa s','Exact $\mu = 10$ mPa s', 'Interpreter', 'latex')
box on
grid on
