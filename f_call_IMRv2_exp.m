function [t, R, epnm] = f_call_IMRv2_exp(Rmax, Req, epnm0, epnmd0, epeq, modes, orders, mu, G, alph, ani, sig, p_a, f_a, tf_nd, tsteps, ultra, varargin)

thisDir = fileparts(mfilename('fullpath'));
persistent imrPathsAdded
if isempty(imrPathsAdded)
    addpath(fullfile(thisDir, 'src', 'forward_solver'));
    addpath(fullfile(thisDir, 'src', 'common'));
    imrPathsAdded = true;
end
opts = parseSolverOptions(varargin{:});

% equation options
% ------- Material Properties ------------------------%
kappa = 1;
T8 = 298.15;
rho8 = 1000;

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
    % pa = p_a;
    % omega = 2*pi*f_a;
    % wavetype = 4;


    % histo
    pa = p_a;
    omega = 2*pi*f_a;
    mn = 3.7;
    dt = pi/omega;
    wavetype = 2;

else
    pa = 0;
    omega = 0;
    wavetype = 0;
    dt = 0; mn = 1;
end


% ------ Simulation time ---------------------------- %
tc = Rmax*sqrt(rho8/101325);
if isempty(opts.OutputTimeNd)
    tfin = tf_nd*tc;
    tvector = linspace(0,tfin,tsteps);
else
    % For optimization this is the dimensional experimental time vector.
    % MATLAB ODE solvers return output at every entry in this tspan vector.
    tvector = opts.OutputTimeNd(:).' .* tc;
end
varin = {'progdisplay',0,'radial',radial,'bubtherm',bubtherm,'tvector',tvector,...
    'vapor',vapor,'medtherm',medtherm,'masstrans',masstrans,'method',opts.Method,...
    'stress',stress,'collapse',collapse,'mu',mu,'g',G,'lambda1',0e-7,...
    'lambda2',0,'alphax', alph, 'ani', ani, 'surft', sig,'r0',Rmax,'req',Req,'kappa',kappa,'t8',T8,...
    'rho8',rho8, 'pa',pa, 'omega', omega, 'wave_type', wavetype, 'perturbed', perturbed, ...
    'modes', modes, 'orders', orders, 'epnm0', epnm0, 'pertmod', pertmod, ...
    'epnmd0', epnmd0, 'epnmeq',epeq, 'reltol', opts.RelTol, 'abstol', opts.AbsTol, 'Nt', opts.Nt, ...
    'dt', dt, 'mn', mn};
if opts.MaxWallTime > 0
    varin = [varin, {'maxwalltime', opts.MaxWallTime}];
end
[t,R,~,~,~,~,~,epnm, ~] = f_imr_fd(varin{:});
end

function opts = parseSolverOptions(varargin)
opts = struct('RelTol', 1e-6, 'AbsTol', 1e-7, 'Nt', 75, ...
    'Method', 45, 'MaxWallTime', 0, 'OutputTimeNd', []);

if mod(numel(varargin), 2) ~= 0
    error('Optional solver inputs must be name-value pairs.');
end

for ii = 1:2:numel(varargin)
    name = lower(varargin{ii});
    value = varargin{ii+1};
    switch name
        case 'reltol'
            opts.RelTol = value;
        case 'abstol'
            opts.AbsTol = value;
        case 'nt'
            opts.Nt = value;
        case 'method'
            opts.Method = value;
        case {'maxwalltime', 'maxtime', 'walltime'}
            opts.MaxWallTime = value;
        case {'outputtimend', 'timevectornd', 'toutnd'}
            opts.OutputTimeNd = value;
        otherwise
            error('Unknown f_call_IMRv2_exp option "%s".', varargin{ii});
    end
end
end
