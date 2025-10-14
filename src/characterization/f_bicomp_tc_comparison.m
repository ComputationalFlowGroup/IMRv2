%% bicomposite collapse time comparison

% --- Load simulation data ---
load('bicomp_sim_data.mat');
%load('R_bicomp_sim_data.mat');

% reading required data input
nX = size(data,1);   % # of experiments
RX = data(:,1);      % All Rmax
LX = data(:,2);      % All amplification Lmax
T1X = data(:,3);     % All collapse time t1 

% physical parameters
% Atmospheric Pressure (Pa)
p_inf = 101325; pbar = p_inf; Pref = p_inf;
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
%
uc = sqrt(pbar/rho);
Ma = cwave/uc;
% ARC ~= 1/0.9147 = 1/trc
ARC = 1/( sqrt(pi/6)*gamma(5/6)/gamma(4/3) );
f_Ma = (2/sqrt( 1 + (1 + 4*(Ma/ARC)^2 )));

tRC = (RX ./uc)/ARC;
T1_ND = T1X .* ARC; % / tRC

rel_err_eb = zeros(size(RX));
rel_err_pimr = zeros(size(RX));
% tg_all = zeros(size(RX));
% tc_pimr_all = zeros(size(RX));
n = size(data,1);
tg_all = zeros(n,1);
tc_pimr_all = zeros(n,1);

for k = 1:nX
    % nX = size(data,1);   % # of experiments
    Rmax_k = data(k,1);    % All Rmax
    Lambda_k = data(k,2);  % All amplification Lmax
    T1_k = data(k,3);      % All collapse time t1

    % setting up model parameters
    We_k = pbar*Rmax_k/(2*gam);

    tRC_k = (Rmax_k/uc)/ARC;
    T1_ND_k = T1_k * ARC; % / tRC
    Req_k = Rmax_k / Lambda_k;

    f_We_k = - pi*ARC./(sqrt(6) * We_k); % This varies with Rmax

    %Req = (RX./LX);
    PG0_k = (p_inf)*(Lambda_k^(-3));
    Alpha_k = 1 + pvsat/PG0_k;
    f_gas_k = C_kap*(Lambda_k^(-3))*Alpha_k;

    %Ca_scale = pbar*ones(nX,1); % Ca*G
    Re_scale_k = rho*uc*Rmax_k; % = Re*mu
    De_scale_k = Rmax_k/uc; % characteristic time scale

    data_fit_k = [Rmax_k,Lambda_k,f_Ma,f_We_k,f_gas_k,Re_scale_k,De_scale_k];

    % call pimr and energy balance collapse time solvers
    %tg = f_bicomp_tc(Req_k,R,G0,G1,l1,Pref,rho8);
    tg_dim = f_bicomp_predict_tg(Req_k,Rmax_k,G0,G1,l1,Pref,rho8);
    % nondimensionalizing
    tg = tg_dim / Rmax_k; %(Rmax_k/sqrt(Pref/rho8));
    tg_all(k) = tg;
  
    tc_pimr = f_pimr_bicomp(G0,G1,l1,Pref,data_fit_k);
    tc_pimr_all(k) = tc_pimr;

    rel_err_eb(k) = abs(tg-T1_ND_k) / T1_ND_k; %* 100;
    rel_err_pimr(k) = abs(tc_pimr-T1_ND_k) / T1_ND_k;% * 100;
end
figure;
scatter(T1_ND,tg_all,60,'filled');
% scatter(T1_ND,rel_err_eb,60,'filled'); hold on;
% scatter(T1_ND,rel_err_pimr,60,'d','filled');
hold on;
scatter(T1_ND(:),tc_pimr_all(:),60,'d','filled');
plot([min(T1_ND),max(T1_ND)],[min(T1_ND),max(T1_ND)],'k--');
xlabel('simulation t_c');
ylabel('theory t_c');
xlim([0 1])
ylim([0 1])
grid on;
axis equal;
