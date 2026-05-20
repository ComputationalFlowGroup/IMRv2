%% basic simulation input
clear all
clc
% close all
%%

addpath src\common\


% load data and process
% load("../data/R1_and_R2_soft.mat")
load("../data/R1_and_R_2.mat")
R1 = table2array(R1_anastas(:,2));
t1 = table2array(R1_anastas(:,1));
R2 = table2array(R2_anastas(:,2));
t2 = table2array(R2_anastas(:,1));
% R1 = table2array(R1_soft(:,2));
% t1 = table2array(R1_soft(:,1));
% R2 = table2array(R2_soft(:,2));
% t2 = table2array(R2_soft(:,1));
% Remove duplicate time points, keeping the first occurrence
[t1_unique, idx1] = unique(t1, 'stable');
R1_unique = R1(idx1);

[t2_unique, idx2] = unique(t2, 'stable');
R2_unique = R2(idx2);

% Optional: overwrite originals
t1 = t1_unique;
R1 = R1_unique;

t2 = t2_unique;
R2 = R2_unique;

% tshare = (t1+t2)./2;

tshare = t1;

R1interp = interp1(t1, R1, tshare);
R2interp = interp1(t2, R2, tshare);

theta = [0, pi/2]; Y20 = sqrt(5/(16*pi))*(3*cos(theta).^2 - 1);

M = [1 Y20(1); 1 Y20(2)];
for i = 1:length(R1interp)
    b = [R1interp(i); R2interp(i)];
    x = M \ b;
    Rbar(i) = x(1); ep2(i) = x(2)./Rbar(i);
end

Rmax = max(Rbar).*1e-6;

R1interp = R1interp./max(Rbar);
R2interp = R2interp./max(Rbar);
Rbar = Rbar./max(Rbar);

tc =  Rmax*sqrt(1048/101325);
tshare = tshare.*1e-6./tc;

figure
% plot(tshare, R1interp, '^--')
hold on
% plot(tshare, R2interp, '^--')
plot(tshare, Rbar, 'o')
plot(tshare, ep2, 'o')


%%


tic
% -------- Radial Solver ----------------------------------------%
% Rmax = 150e-6;
% Rmax = 50e-6;
% Req = Rmax;
Req = Rmax/3;
mu =  0.1;
G = 15.09e3;
alph = 0;
ani = [3 0];
sig = 0.0;
p_a = -10*101325; f_a = 1e6/pi;
rho = 1000;
p8 = 101325;
tcLIC = Rmax*sqrt(rho/p8);
tf_nd = 4/3;
tsteps = 5000; ultra = false;


% -------- perturbation solver initial conditions ---------------%
% Mode numbers
n = [2];
m = [0];
N = n;
ep0 = [ep2(2)];
epd0 = [0];
epeq = [0];


t = linspace(0, tf_nd, tsteps);
[t, R, epnm] = f_call_IMRv2(Rmax, Req, ep0, epd0, epeq, n, m, mu, G, alph, ani, sig, p_a, f_a, tf_nd, tsteps, ultra);

% Rsiminterp = interp1(t, R, tshare(tshare < max(t)), 'linear');
% epnmsiminterp = interp1(t, epnm, tshare(tshare < max(t)), 'linear');
% 
% Rbar(isnan(Rbar)) = 0;
% ep2(isnan(ep2)) = 0;
% 
% rmseR = rmse(Rsiminterp,Rbar(tshare < max(t))')
% rmseep = rmse(epnmsiminterp,ep2(tshare < max(t))')


%%
figure
plot(t.*3, epnm, '-')
hold on
plot(t.*3,R.*3, '-')
plot(tshare.*3, Rbar.*3, 'o')
plot(tshare.*3, ep2, 'o')
 ylim([-1 3])




