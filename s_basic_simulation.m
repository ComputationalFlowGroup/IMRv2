%% basic simulation input
clear all
clc
close all

addpath src\common\
% % addpath ../Anisotropic_material_IMR/IMRv2/
% %%
% % load data and process
% load("../data/high_stretch/R1R2_data.mat")
% R1 = table2array(R1_Lamb(:,2));
% t1 = table2array(R1_Lamb(:,1));
% R2 = table2array(R2_Lamb(:,2));
% t2 = table2array(R2_Lamb(:,1));
% 
% 
% 
% % load("../data/stiff/R1_and_R_2.mat")
% % R1 = table2array(R1_anastas(:,2));
% % t1 = table2array(R1_anastas(:,1));
% % R2 = table2array(R2_anastas(:,2));
% % t2 = table2array(R2_anastas(:,1));
% 
% % load("../data/R1_and_R2_soft.mat")
% % R1 = table2array(R1_soft(:,2));
% % t1 = table2array(R1_soft(:,1));
% % R2 = table2array(R2_soft(:,2));
% % t2 = table2array(R2_soft(:,1));
% 
% % Remove duplicate time points, keeping the first occurrence
% [t1_unique, idx1] = unique(t1, 'stable');
% R1_unique = R1(idx1);
% 
% [t2_unique, idx2] = unique(t2, 'stable');
% R2_unique = R2(idx2);
% 
% % Optional: overwrite originals
% t1 = t1_unique;
% R1 = R1_unique;
% 
% t2 = t2_unique;
% R2 = R2_unique;
% 
% % tshare = (t1+t2)./2;
% 
% tshare = t2;
% 
% R1interp = interp1(t1, R1, tshare);
% R2interp = interp1(t2, R2, tshare);
% 
% theta = [0, pi/2]; Y20 = sqrt(5/(16*pi))*(3*cos(theta).^2 - 1);
% 
% M = [1 Y20(1); 1 Y20(2)];
% for i = 1:length(R1interp)
%     b = [R1interp(i); R2interp(i)];
%     x = M \ b;
%     Rbar(i) = x(1); ep2(i) = x(2)./Rbar(i);
% end
% 
% Rmax = max(Rbar).*1e-6;
% 
% R1interp = R1interp./max(Rbar);
% R2interp = R2interp./max(Rbar);
% Rbar = Rbar./max(Rbar);
% 
% tc =  Rmax*sqrt(1048/101325);
% tshare = tshare.*1e-6./tc;
% 
% figure
% % plot(tshare, R1interp, '^--')
% hold on
% % plot(tshare, R2interp, '^--')
% plot(tshare, Rbar, 'o')
% plot(tshare, ep2, 'o')


load("../data/Sims_Brown_Surya/FEM_Equil_analysis.mat")

% Req = amp_extractf(1,:).*1e-6;
% epnmeq =  0.*amp_extractf(3:end,:).*1e-6./Req;
epnmeq =  0.*amp_extractf(3:end,:).*1e-6;



% load("../data/Sims_Brown_Surya/aniso_sim_FEA_sphequil.mat")
load("../data/Sims_Brown_Surya/iso_sim_FEA_slight.mat")
% load("../data/Sims_Brown_Surya/aniso_sim_FEA.mat")
Req =  amp_extractf(1,end).*1e-6;
Rexp = amp_extractf(1,:).*1e-6;
epnmeq =  amp_extractf(3:end,end).*1e-6./Req;

Rmax = Rexp(1);
k = 0;
idxs = 3:2:9;%size(amp_extractf,1);
for i = idxs
    k = k+1;
    amp(k,:) = amp_extractf(i,:)./amp_extractf(1,:);
end
texp = tStep;

%%
addpath src/common/
tic
% -------- Radial Solver ----------------------------------------%
% Rmax = 150e-6;
% Rmax = 50e-6;
% Req = Rmax;
% Rmax = 100e-6;
% Req = Rmax;
mu =  0.2625;
G = 105e3;
alph = 0.0;
ani = [0 0];
sig = 0.0;
p_a = -1.15*101325; f_a = 50e3;
rho = 1000;
p8 = 101325;
tcLIC = Rmax*sqrt(rho/p8);
tf_nd = 3;
tsteps = 30000; ultra = false;


% -------- perturbation solver initial conditions ---------------%
% Mode numbers
n = mode_extractf(idxs,1)';
m = zeros(size(n));
N = n;
ep0 = amp(:,1);
epd0 = zeros(size(ep0));
epeq = epnmeq(idxs-2);
% ep0 = epeq.*0;
% Rmax = Req.*1.5;



t = linspace(0, tf_nd, tsteps);
[t, R, epnm] = f_call_IMRv2(Rmax, Req, ep0, epd0, epeq, n, m, mu, G, alph, ani, sig, p_a, f_a, tf_nd, tsteps, ultra);

[t, Riso, epnmiso] = f_call_IMRv2(Rmax, Req, ep0, epd0, epeq, n, m, mu, G, alph, [0 0], sig, p_a, f_a, tf_nd, tsteps, ultra);



% Rsiminterp = interp1(t, R, tshare(tshare < max(t)), 'linear');
% epnmsiminterp = interp1(t, epnm, tshare(tshare < max(t)), 'linear');
% 
% Rbar(isnan(Rbar)) = 0;
% ep2(isnan(ep2)) = 0;
% 
% rmseR = rmse(Rsiminterp,Rbar(tshare < max(t))')
% rmseep = rmse(epnmsiminterp,ep2(tshare < max(t))')

%%
% load("../data/Sims_Brown_Surya/aniso_model_FEA_rad_aniso.mat")
figure
plotl = ceil(sqrt(size(epnm,2)+1));

subplot(plotl, plotl, 1)
hold on
plot(t, R, '-')
plot(t, Riso, '--')
plot(texp./tcLIC, Rexp./Rexp(1), 'o')
% yline(Req/Rmax, 'k--', 'LineWidth',2)
xlabel("t^*")
ylabel("R")


for i = 1:size(epnm,2)
    subplot(plotl, plotl, i + 1)
    hold on
    plot(t, epnm(:,i), '-')
    plot(t, epnmiso(:,i), '--')
    plot(texp./tcLIC, amp(i,:), 'r^')
    xlabel("t^*")
    str = sprintf('$\\epsilon_{%.0f}$', n(i));
    ylabel(str, 'Interpreter', 'latex')
    % yline(epeq(i), 'k--', 'LineWidth',2)
end

% xlim([0 1])

%%
clear all
close all
addpath ../cmap/

cmap = viridis(5);

clc
load("../data/soft/data_and_sim.mat")
figure
plot(tshare, Rbar, 'o', 'MarkerFaceColor',cmap(1,:), 'Color',cmap(1,:))
hold on
plot(tshare, ep2, '^', 'MarkerFaceColor',cmap(2,:), 'Color',cmap(3,:))
plot(t,R, '-','LineWidth',2, 'Color',cmap(3,:))
plot(t, epnm(:,1), '-', 'LineWidth',2, 'Color',cmap(4,:))


% load("../data/soft/fullmodel_noaniso.mat")
% plot(t./max(R), R./max(R), '--')
% hold on
% plot(t./max(R), ep, 'r--')
% hold on
% plot(t./max(R), epirr(1, 1:length(t)), '-.')
ylim([-0.25 1])
xlim([0 1])
ax = gca;

ylabel('Interface Quantity', 'Interpreter','latex', 'FontSize',16)
xlabel('Time', 'Interpreter','latex', 'FontSize',16)
ax.TickLabelInterpreter = 'latex';
legend('$\overline{R}$ Tzoumaka 2022', '$\epsilon_2$ Tzoumaka 2022', ...
    '$\overline{R}$ Current model', '$\epsilon_2$ Current model', 'Interpreter', ...
    'latex', 'FontSize', 14)
ax = gca;
ax.FontSize = 14;
ax.LabelFontSizeMultiplier = 1.5;

%%
clear all
clc
close all

load("../data/high_stretch/data_and_sim.mat")

R1sim = R;
R2sim = R;
tf = max(t);

x = cos(theta);

for i = 1:length(n)
    Pnt = legendre(n(i),x);
    Rmod = Pnt(1,:);
    R1sim = R1sim + Rmod(1).*R.*epnm(:,i).*sqrt((2*n(i)+1)/(4*pi));
    R2sim = R2sim + Rmod(2).*R.*epnm(:,i).*sqrt((2*n(i)+1)/(4*pi));
end

figure(1)
plot(tshare, R1interp, 'o')
hold on
plot(tshare, R2interp, 'o')
plot(t, R1sim)
plot(t, R2sim)

R1simint = interp1(t, R1sim, tshare);
R2simint = interp1(t, R2sim, tshare);


figure(2)
semilogy(tshare,abs(R1interp-R1simint)./R1interp, 'LineWidth',2)
hold on
figure(3)
semilogy(tshare,abs(R2interp-R2simint)./R2interp, 'LineWidth',2)
hold on

avg_relerr_pa1 = mean((abs(R1interp(tshare < tf)-R1simint(tshare < tf))./R1interp(tshare < tf)));
avg_relerr_pa2 = mean((abs(R2interp(tshare < tf)-R2simint(tshare < tf))./R2interp(tshare < tf)));

0.5*(avg_relerr_pa1+avg_relerr_pa2)


load("../data/high_stretch/fullmodel_noaniso.mat")
t = t./max(R);
R = R./max(R);


R1sim = R;
R2sim = R;

x = cos(theta);

for i = 1:length(n)
    Pnt = legendre(n(i),x);
    Rmod = Pnt(1,:);
    R1sim = R1sim + Rmod(1).*R.*ep(:).*sqrt((2*n(i)+1)/(4*pi));
    R2sim = R2sim + Rmod(2).*R.*ep(:).*sqrt((2*n(i)+1)/(4*pi));
end
figure(1)
plot(t, R1sim, '--')
hold on
plot(t, R2sim, '--')
ylim([-1 1.5])


R1simint = interp1(t, R1sim, tshare);
R2simint = interp1(t, R2sim, tshare);



figure(2)
semilogy(tshare,abs(R1interp-R1simint)./R1interp, '--')
hold on
figure(3)
semilogy(tshare,abs(R2interp-R2simint)./R2interp,  '--')
hold on

avg_relerr_f1 = mean((abs(R1interp(tshare < tf)-R1simint(tshare < tf))./R1interp(tshare < tf)));
avg_relerr_f2 = mean((abs(R2interp(tshare < tf)-R2simint(tshare < tf))./R2interp(tshare < tf)));

0.5*(avg_relerr_f1+avg_relerr_f2)

load("../data/high_stretch/fullmodel_noaniso.mat")
t = t./max(R);
R = R./max(R);

R1sim = R;
R2sim = R;

x = cos(theta);

for i = 1:length(n)
    Pnt = legendre(n(i),x);
    Rmod = Pnt(1,:);
    R1sim = R1sim + Rmod(1).*R.*epirr(:,i).*sqrt((2*n(i)+1)/(4*pi));
    R2sim = R2sim + Rmod(2).*R.*epirr(:,i).*sqrt((2*n(i)+1)/(4*pi));
end
figure(1)
plot(t, R1sim, '-.')
plot(t, R2sim, '-.')
ylim([-1 1.5])


R1simint = interp1(t, R1sim, tshare);
R2simint = interp1(t, R2sim, tshare);



figure(2)
semilogy(tshare,abs(R1interp-R1simint)./R1interp, 'k-.')
hold on
figure(3)
semilogy(tshare,abs(R2interp-R2simint)./R2interp, 'k-.')
hold on


avg_relerr_p1 = mean((abs(R1interp(tshare < tf)-R1simint(tshare < tf))./R1interp(tshare < tf)));
avg_relerr_p2 = mean((abs(R2interp(tshare < tf)-R2simint(tshare < tf))./R2interp(tshare < tf)));

0.5*(avg_relerr_p1+avg_relerr_p2)




