% Parsimonious IMR (pIMR) toward rapid first-estimate of viscoelastic
% properties: this main script reads a datafile containing micro-cavitation
% experiments and fits for a suitable set of viscoelastic properties

% acquiring required data input
infile = 'data.mat';   % Name of file to read
load(infile);

infile2 = 'Rdata.mat';
load(infile2);

% reading required data input
nX = size(data,1);   % # of experiments
RX = data(:,1);      % All Rmax
LX = data(:,2);      % All amplification Lmax
T1X = data(:,3);     % All collapse time t1

% physical parameters
% Atmospheric Pressure (Pa)
p_inf = 101325;
% Density of characterized material (kg/m^3)
rho = 998.2;
% Surface tension (N/m)
gam = 0.070;
% Wave speed in characterized material (m/s)
cwave = 1484;
% Saturated vapor pressure at far-field temperature (Pa)
pvsat = 3116.7757;
% pIMR non-viscoelastic parameters
C_kap = 1.4942;

% setting up model parameters
pbar = p_inf;
We = pbar*RX/(2*gam);
uc = sqrt(pbar/rho);
Ma = cwave/uc;
% ARC ~= 1/0.9147 = 1/trc
ARC = 1/( sqrt(pi/6)*gamma(5/6)/gamma(4/3) );
f_Ma = (2/sqrt( 1 + (1 + 4*(Ma/ARC)^2 ))) * ones(nX,1);
f_We = - pi*ARC./(sqrt(6) * We); % This varies with Rmax
Req = (RX./LX);
PG0 = (p_inf).*(LX.^(-3));
Alpha = 1 + pvsat./PG0;
f_gas = C_kap*(LX.^(-3)).*Alpha;
Ca_scale = pbar*ones(nX,1); % Ca*G
Re_scale = rho*uc*RX; % = Re*mu
De_scale = RX/uc; % characteristic time scale

% Construct new matrix to pass to solver:
data_fit = [RX,LX,f_Ma,f_We,f_gas,Ca_scale,Re_scale,De_scale];

% Convert T1X to dimensionless, relative to t_{RC} at each scale.
tRC = (RX/uc)/ARC;
T1_ND = T1X./tRC;

% graded NH fit
t1approx_vals = @(X) f_fit_gNH(X(1),X(2),X(3),X(4),X(5),X(6),data_fit);
err_gNH = @(X) (T1_ND./t1approx_vals(X)).^2 - 1;
err_fn_gNH = @(X) log10((err_gNH(X))'*(err_gNH(X))/nX);

% define search range and initial start for each graded parameter (dim)
G0_min = 1E-3;
G0_max = 1E5;
G0_start = 1E3;

G1_min = 1E-3;
G1_max = 1E5;
G1_start = 10E3;

% l1_min is at least 1*max(Req)
Req_max = max(RX)/min(LX);
l1_min = Req_max;
l1_start = 1.2E-4;
l1_max = l1_start + 1e-6;

% l2 must be greater than l1
l2_min = l1_start + 1e-6;
l2_max = 5*Req_max;
l2_start = 1.8e-4;

% min should be greater than 1, upperbound?
va_min = 1;
va_max = 4;
va_start = 2;

% should be less than 0.5, lowerbound?
vnc_min = 0;
vnc_max = 0.5;
vnc_start = 0.3;

%%
% contains all X values the optimizer tries
X_log = [];
% add options to set optimization tolerances
opts = optimset('Display','iter','TolX',1e-1,'TolFun',1e-1,'MaxIter',50,'MaxFunEvals',120,'OutputFcn', @output_log);

% optimizing based on search ranges
opt_fit_gNH = f_minsearchbnd(err_fn_gNH,...
                [G0_start,G1_start,l1_start,l2_start,va_start,vnc_start],...
                [G0_min,G1_min,l1_min,l2_min,va_min,vnc_min],...
                [G0_max,G1_max,l1_max,l2_max,va_max,vnc_max],...
                opts);

% extract best fit of each parameter
G0_opt_NH = opt_fit_gNH(1);
G1_opt_NH = opt_fit_gNH(2);
l1_opt_NH = opt_fit_gNH(3);
l2_opt_NH = opt_fit_gNH(4);
va_opt_NH = opt_fit_gNH(5);
vnc_opt_NH = opt_fit_gNH(6);
disp("Overall Best Fit: G0 = " + G0_opt_NH + " Pa, G1 = " + G1_opt_NH + " Pa, l1 = " + l1_opt_NH + "m, l2 = " + l2_opt_NH + "m, v_a = " + va_opt_NH + ", v_nc = " + vnc_opt_NH)
%%
% reconstruct all t1approx_vals
t1_approx_all = zeros(size(X_log,1),length(T1_ND));
for i = 1:size(X_log)
    X = X_log(i,:);
    err_vec = (T1_ND ./ fit_gNH(X(1),X(2),X(3),X(4),X(5),X(6),data_fit)).^2 - 1;
    format long;
    t1_approx_all(i, :) = T1_ND ./ sqrt(err_vec + 1);  % exact inverse
end
%%
% Also find single-parameter fits:
opt_fit_NHG0 = fminsearchbnd(err_fn_gNH,[G0_start,0,0,0,0,0],[G0_min,0,0,0,0,0],[G0_max,0,0,0,0,0]);
opt_fit_NHG1 = fminsearchbnd(err_fn_gNH,[0,G1_start,0,0,0,0],[0,G1_min,0,0,0,0],[0,G1_max,0,0,0,0]);
opt_fit_NHl1 = fminsearchbnd(err_fn_gNH,[0,0,l1_start,0,0,0],[0,0,l1_min,0,0,0],[0,0,l1_max,0,0,0]);
opt_fit_NHl2 = fminsearchbnd(err_fn_gNH,[0,0,0,l2_start,0,0],[0,0,0,l2_min,0,0],[0,0,0,l2_max,0,0]);
opt_fit_NHva = fminsearchbnd(err_fn_gNH,[0,0,0,0,va_start,0],[0,0,0,0,va_min,0],[0,0,0,0,va_max,0]);
opt_fit_NHvnc = fminsearchbnd(err_fn_gNH,[0,0,0,0,0,vnc_start],[0,0,0,0,0,vnc_min],[0,0,0,0,0,vnc_max]);

disp("NH Best Fit for near field: G = " + opt_fit_NHG0(1) + " Pa.")
disp("NH Best Fit for far field: G1 = " + opt_fit_NHG1(1) + " Pa.")
disp("NH Best Fit for near field region: l1 = " + opt_fit_NHl1(1) + "m")
disp("NH Best Fit for far field region: l2 = " + opt_fit_NHl2(1) + "m")
disp("Best Fit for C-Y parameter: v_a = " + opt_fit_NHva(1))
disp("Best Fit for C-Y parameter: v_nc = " + opt_fit_NHvnc(1))
disp("Overall Best Fit: G0 = " + G0_opt_NH + " Pa, G1 = " + G1_opt_NH + " Pa, l1 = " + l1_opt_NH + "m, l2 = " + l2_opt_NH + "m, v_a = " + va_opt_NH + ", v_nc = " + vnc_opt_NH)
