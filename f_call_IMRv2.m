function [t, R, epnm] = f_call_IMRv2(Rmax, Req, epnm0, epnmd0, epeq, modes, orders, mu, G, alph, ani, sig, p_a, f_a, tf_nd, tsteps, ultra)

addpath ../IMRv2/src/forward_solver/

% equation options
% ------- Material Properties ------------------------%
kappa = 1.4;
T8 = 298.15;
rho8 = 1048;

% ------- Simulation settings ------------------------%
radial = 2;
vapor = 1;
collapse = 0;
bubtherm = 1;
medtherm = 0;
masstrans = 1;
stress = 2;
perturbed = 1;
pertmod = 0;

if ultra
    % --------- Ultrasound settins -----------------------%
    % sinusoidal
    pa = p_a;
    omega = 2*pi*f_a;
    wavetype = 4;


%     % histo
%     p = ee*(0.5 + 0.5*cos(om*(t - dt))).^mn;
%     pa = p_a;
%     omega = 2*pi*f_a;
%     mn = 3.7;
%     dt = pi/omega;
% 
else
    pa = 0;
    omega = 0;
    wavetype = 0;
end

% ------ Simulation time ---------------------------- %
tc = Rmax*sqrt(rho8/101325);
tfin = tf_nd*tc;
tvector = linspace(0,tfin,tsteps);
varin = {'progdisplay',0,'radial',radial,'bubtherm',bubtherm,'tvector',tvector,...
    'vapor',vapor,'medtherm',medtherm,'masstrans',masstrans,'method',23,...
    'stress',stress,'collapse',collapse,'mu',mu,'g',G,'lambda1',0e-7,...
    'lambda2',0,'alphax', alph, 'ani', ani, 'surft', sig,'r0',Rmax,'req',Req,'kappa',kappa,'t8',T8,...
    'rho8',rho8, 'pa',pa 'omega', omega, 'wave_type', wavetype, 'perturbed', perturbed, ...
    'modes', modes, 'orders', orders, 'epnm0', epnm0, 'pertmod', pertmod, ...
    'epnmd0', epnmd0, 'epnmeq',epeq, 'reltol', 1e-4, 'abstol', 1e-6, 'Nt', 75};
% run the forward solver
[t,R,~,~,~,~,~,epnm, ~] = f_imr_fd(varin{:});


end