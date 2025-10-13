%% bicomposite collapse time comparison

% --- Load simulation data ---
load('bicomp_sim_data.mat');  % Should give R_sim, t_sim

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

%Ca_scale = pbar*ones(nX,1); % Ca*G
Re_scale = rho*uc*RX; % = Re*mu
De_scale = RX/uc; % characteristic time scale

data_fit = [RX,LX,f_Ma,f_We,f_gas,Re_scale,De_scale];

% Convert T1X to dimensionless, relative to t_{RC} at each scale.
tRC = (RX/uc)/ARC;
T1_ND = T1X./tRC; %or multiply by tc-sim because that's what it's saved by


rel_err_eb = zeros(size(RX));
rel_err_pimr = zeros(size(RX));
for k = 1:nX
    nX = size(data,1);   % # of experiments
    RX = data(k,1);      % All Rmax
    LX = data(k,2);      % All amplification Lmax
    T1X = data(k,3);     % All collapse time t1

    tRC = (RX/uc)/ARC;
    T1_ND = T1X./tRC;

    % call pimr and energy balance collapse time solvers
    [tg] = f_bicomp_tc(Req,R,G0,G1,l1,Pref,rho8);
    [tc_pimr] = f_pimr_bicomp(G0,G1,l1,data_fit);

    rel_err_eb(k) = abs(tg-T1_ND) / T1_ND * 100;
    rel_err_pimr(k) = abs(tc_pimr-T1_ND) / T1_ND * 100;
end
figure;
scatter(T1_ND(:),tg(:),60,'filled');
hold on;
scatter(T1_ND(:),tc_pimr(:),60,'d','filled');
plot([min(T1_ND),max(T1_ND)],[min(T1_ND),max(T1_ND)],'k--');
xlabel('simulation t_c');
ylabel('theory t_c');
grid on;
axis equal;
%legend('simulations','perfect agreement line');

% reconstruct all t1approx_vals
% t1_approx_all = zeros(size(X_log,1),length(T1_ND));
% for i = 1:size(X_log)
%     X = X_log(i,:);
%     err_vec = (T1_ND ./ fit_gNH(X(1),X(2),X(3),X(4),X(5),X(6),data_fit)).^2 - 1;
%     format long;
%     t1_approx_all(i, :) = T1_ND ./ sqrt(err_vec + 1);  % exact inverse
% end
