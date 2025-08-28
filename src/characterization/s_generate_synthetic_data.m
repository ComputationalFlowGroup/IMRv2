% file s_generate_synthetic_data.m
% brief contains a loop to run the froward solver a specified number of
% times to create datasets that resemble experimental data

clear all
close all
clc
addpath ../common/
addpath ../forward_solver/

ndata = 10;
space0 = [275e-6 0.25];
lb = [200e-6 1/8];
ub = [375e-6 1/2.5];
dim = numel(space0);
array = repmat(lb, ndata, 1) + lhsdesign(ndata, dim) .* repmat(ub-lb, ndata, 1);

    
%%
tic
for i = 1:ndata
    fps = 2e6;
    % equation options
    R0 = array(i,1);
    Req = R0*array(i,2);
    tfin = 192/fps;
    tc = R0*sqrt(1048/101325);
    kappa = 1.4;
    Lheat = 2.378193575129533e+04;
    T8 = 298.15;
    rho8 = 1048;
    ST = 0.032;
    alphaxs = 1;
    alphax = alphaxs*(1-0.1+0.2*rand(1,1));
    mus = 5e-2;
    mu = mus*(1-0.1+0.2*rand(1,1));
    G = 1e3;
    tvector = linspace(0,tfin,192);
    radial = 2;
    vapor = 1;
    collapse = 0;
    bubtherm = 0;
    medtherm = 0;
    masstrans = 0;
    perturbed = 1;
    stress = 2;
    modes = 2:9;
    orders = 2:9;
    epnm0 = -1e-4+2e-4*rand(length(modes), 1);
    epnmd0 = -0.05+0.1*rand(length(modes), 1);
    varin = {'progdisplay',0,...
        'radial',radial,...
        'bubtherm',bubtherm,...
        'perturbed', perturbed, ...
        'tvector',tvector,...
        'vapor',vapor,...
        'medtherm',medtherm,...
        'masstrans',masstrans,...
        'method',23,...
        'stress',stress,...
        'collapse',collapse,...
        'mu',mu,...
        'alphax', alphax, ...
        'g',G,...
        'lambda1',0,...
        'lambda2',0,...
        'surft', ST, ...
        'r0',R0,...
        'req',Req,...
        'kappa',kappa,...
        't8',T8,...
        'rho8',rho8, 'modes', modes, 'orders', orders, 'epnm0', epnm0, ...
        'epnmd0', epnmd0, 'abstol', 1e-7, 'reltol', 1e-5};
    [tfd,Rfd,Rfddot,Pfd,Tfd,Tmfd,kvfd, epnm, epnmd] = f_imr_fd(varin{:},'Nt',100);
    Rout = zeros(length(tfd), 1);
    epnmout = 0.*epnm;
    for s = 1:length(tfd)
        Rout(s) = Rfd(s).*R0.*(1-0.05+0.1*rand(1,1));
        for j = 1:length(modes)
            epnmout(s,j) = epnm(s,j).*(1-0.05+0.1*rand(1,1));
        end
    end
    kindata{i} = struct('time', tfd.*tc, 'R', Rout, 'Req', Req, 'epnm', epnmout,...
        'epnmd0', epnmd0, 'n', modes, 'm', orders, 'G', G, 'alpha', alphax, ...
        'mu', mu, 'ST', ST, 'fps', fps,'rho', rho8);
    figure
    plot(tfd, Rout)
    hold on
    yline(Req/R0)
end
toc

%%

figure(1)
plot(tfd, Rout, '--')
hold on

figure(2)
plot(tfd, epnmout, '--')
hold on
