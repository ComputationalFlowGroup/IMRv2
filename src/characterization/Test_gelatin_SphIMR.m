clear all
clc
load("../Experimental_data/UM_data/dataset2_IMR.mat")

gel = [1:10 21:34]; % 5% gel
gel = 11:20; % 10% gel

gelidxgood = [15 21 22 23 24];
for i = gel

    expt1 = expts(i);
    R = expt1.Roft;
    t = expt1.t;
    Req = expt1.Req;

    [Rmax, idxmax] = max(R);

    tc = Rmax*sqrt(1048/101325);
    R = R(idxmax:end)./Rmax;
    t = t(idxmax:end)./tc;
    t = t - t(1);
    %
    % file s_generate_synthetic_data.m
    % brief contains a loop to run the froward solver a specified number of
    % times to create datasets that resemble experimental data

    addpath src/common/
    addpath src/forward_solver/

    fps = 2e6;
    % equation options
    R0 = Rmax;
    Req = Req;
    tfin = 192/fps;
    kappa = 1.4;
    Lheat = 2.378193575129533e+04;
    T8 = 298.15;
    rho8 = 1048;
    ST = 0.04;
    alphax = 1.25;
    mu = 1.5e-2;
    G = 7.81e3;
    tvector = linspace(0,tfin,192);
    radial = 2;
    vapor = 1;
    collapse = 0;
    bubtherm = 1;
    medtherm = 0;
    masstrans = 1;
    perturbed = 0;
    pertmod = 0;
    stress = 2;
    modes = 2;
    orders = 2;
    epnm0 = 0;
    epnmd0 =0;
    varin = {'progdisplay',0,...
        'radial',radial,...
        'bubtherm',bubtherm,...
        'tvector',tvector,...
        'vapor',vapor,...
        'perturbed', perturbed, ...
        'pertmod', pertmod, ...
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
        'rho8',rho8, 'abstol', 1e-7, 'reltol', 1e-5};
    [tfd,Rfd,Rfddot,Pfd,Tfd,Tmfd,kvfd] = f_imr_fd(varin{:},'Nt',100);
    Rout = zeros(length(tfd), 1);
    for s = 1:length(tfd)
        Rout(s) = Rfd(s);
    end
    figure
    hold on
    plot(tfd, Rout)
    hold on
    plot(t,R, 'o')
    xlim([0 3])
end