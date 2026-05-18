%% basic simulation input
clear all
clc
%close all

addpath src\common\
% load('data/Expt_1_EML_IC.mat')


tic
% -------- Radial Solver ----------------------------------------%
Rmax = 172e-6;
Req = Rmax/3;
mu =  0.15;
G = 15.07e3;
alph = 0;
sig = 0.056;
p_a = -50e3; f_a = 28e3;
rho = 1048;
p8 = 101325;
tcLIC = Rmax*sqrt(rho/p8);
tf_nd = 2;
tsteps = 5000; ultra = false;


% -------- perturbation solver initial conditions ---------------%
% Mode numbers
n = 2;
N = n;
ep0 = 0.25;
epd0 = 0;


t = linspace(0, tf_nd, tsteps);
[t, R, epnm] = f_call_IMRv2(Rmax, Req, ep0, epd0, n, mu, G, alph, sig, p_a, f_a, tf_nd, tsteps, ultra);

%%
figure
plot(t, epnm)
hold on
plot(t,R)




