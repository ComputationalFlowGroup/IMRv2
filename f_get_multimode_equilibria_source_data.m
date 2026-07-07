function analysis = f_get_multimode_equilibria_source_data()
%F_GET_MULTIMODE_EQUILIBRIA_SOURCE_DATA Build source-term diagnostics.

imr_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(imr_dir);
data_dir = fullfile(repo_dir, 'data', 'Sims_Brown_Surya');
addpath(fullfile(imr_dir,'src','forward_solver'));
addpath(fullfile(imr_dir,'src','common'));

equil_file = fullfile(data_dir,'FEM_Equil_analysis.mat');
sim_file = fullfile(data_dir,'aniso_sim_FEA.mat');
require_file(equil_file);
require_file(sim_file);

equil_data = load(equil_file,'amp_extractf');
Req = equil_data.amp_extractf(1,1).*1e-6;
epnmeq_all = equil_data.amp_extractf(3:size(equil_data.amp_extractf,1),1).*1e-6./Req;

sim_data = load(sim_file,'amp_extractf','mode_extractf');
Rmax = sim_data.amp_extractf(1,1).*1e-6;

mode_rows = 3:13;
n = sim_data.mode_extractf(mode_rows,1).';
m = zeros(size(n));
ep0 = (sim_data.amp_extractf(mode_rows,1)./sim_data.amp_extractf(1,1)).';
epd0 = zeros(size(ep0));
epeq_fem = epnmeq_all(mode_rows-2).';
epeq_zero = zeros(size(epeq_fem));

mu = 0.2625;
G = 105e3;
alph = 0.0;
ani = [5 0];
sig = 0.0;
rho = 1000;
p8 = 101325;
tc = Rmax*sqrt(rho/p8);
tvector = linspace(0,3.5*tc,1500);

cases = struct([]);
cases(1).name = 'fem-epeq far root';
cases(1).short_name = 'fem';
cases(1).epeq = epeq_fem;
cases(1).root = [0.981729377481; ...
    0.403240673; 0; 0.133016399; 0; 0.0572179797; 0; ...
    0.0241291673; 0; 0.00905540539; 0; 0.00272857706];
cases(2).name = 'zero-epeq far root';
cases(2).short_name = 'zero';
cases(2).epeq = epeq_zero;
cases(2).root = [0.789112464273; ...
    0.166775097; 0; 0.0679752349; 0; 0.0303516056; 0; ...
    0.0124513588; 0; 0.00421969508; 0; 0.000935443187];

for c = 1:numel(cases)
    args = {'progdisplay',0,'radial',2,'bubtherm',0,'tvector',tvector, ...
        'vapor',0,'medtherm',0,'masstrans',0,'method',23,'stress',2, ...
        'collapse',0,'mu',mu,'g',G,'lambda1',0,'lambda2',0, ...
        'alphax',alph,'ani',ani,'surft',sig,'r0',Rmax,'req',Req, ...
        'kappa',1,'t8',298.15,'rho8',rho,'pa',0,'omega',0, ...
        'wave_type',0,'perturbed',1,'modes',n,'orders',m, ...
        'epnm0',ep0,'pertmod',0,'epnmd0',epd0,'epnmeq',cases(c).epeq, ...
        'reltol',1e-4,'abstol',1e-5,'Nt',75,'dt',0,'mn',1};

    R = cases(c).root(1);
    ep = cases(c).root(2:end);

    D = f_imr_rhs_diagnostics([R;0;NaN;ep;zeros(size(ep))], args{:});
    B = f_imr_anisotropic_projection_breakdown(R, ep, args{:});

    xeq = [Req/Rmax; cases(c).epeq(:)];
    xfar = cases(c).root;
    [Aeq, JA] = linearize_anisotropy(xeq, args);
    Afin = anisotropy_vector(xfar, args);
    Alin = Aeq + JA*(xfar - xeq);
    Dlinear = f_imr_rhs_diagnostics([xfar(1);0;NaN;xfar(2:end); ...
        zeros(numel(n),1)], args{:}, 'rad_mod_override', Alin(1), ...
        'ep_mod_override', Alin(2:end));

    cases(c).args = args;
    cases(c).R = R;
    cases(c).ep = ep(:);
    cases(c).diag = D;
    cases(c).breakdown = B;
    cases(c).radial_without_aniso = D.P - 1 + D.S;
    cases(c).radial_with_aniso = cases(c).radial_without_aniso + D.rad_mod;
    cases(c).Afin = Afin;
    cases(c).Alin = Alin;
    cases(c).linearized_diag = Dlinear;
    cases(c).linearized_residual_norm = norm([Dlinear.Rddot; Dlinear.epddot]);
    cases(c).finite_residual_norm = norm([D.Rddot; D.epddot]);
end

analysis = struct();
analysis.imr_dir = imr_dir;
analysis.Req = Req;
analysis.Rmax = Rmax;
analysis.Req_nd = Req/Rmax;
analysis.n = n(:);
analysis.m = m(:);
analysis.ep0 = ep0(:);
analysis.epeq_fem = epeq_fem(:);
analysis.epeq_zero = epeq_zero(:);
analysis.cases = cases;

end

function [A0, JA] = linearize_anisotropy(x0, args)
A0 = anisotropy_vector(x0, args);
JA = zeros(length(A0), length(x0));
for j = 1:length(x0)
    h = 1e-6*(1 + abs(x0(j)));
    xp = x0;
    xm = x0;
    xp(j) = xp(j) + h;
    xm(j) = xm(j) - h;
    JA(:,j) = (anisotropy_vector(xp, args) - ...
        anisotropy_vector(xm, args))./(2*h);
end
end

function A = anisotropy_vector(x, args)
B = f_imr_anisotropic_projection_breakdown(x(1), x(2:end), args{:});
A = [B.rad_mod; B.ep_mod(:)];
end

function require_file(pathname)
if ~isfile(pathname)
    error('Required file not found: %s', pathname);
end
end
