% call f_tcol_calc_graded given each Rmax, Lambdamax val
close all; clear all; clc;

% input synthetic experiments
infile = 'data.mat';   % Name of file to read
load(infile);
infile2 = 'Rdata.mat';
load(infile2);

% read synthetic experiment data
nX = size(data,1);   % # of experiments
RX = data(:,1);      % All Rmax 
LX = data(:,2);      % All amplification Lmax
T1X = data(:,3);     % All collapse time tc [approx tg]

% EXAMPLE DATA - TO MATCH SYNTHETIC DATA
stress = 1; G0 = 1E3; G1 = 1E4; l1 = 1.2e-4; l2 = 1.8e-4; v_a = 2; v_nc = 0.3; 
rho8 = 1064; Pref = 101325;
format long;

% nondimensionalize approx collapse time
tc_all = f_tc_nondim(nX,RX,T1X,rho8,Pref);
%%

% generate tg from energy balance based on {Rmax, Lambdamax} [exact tg]
tg_all = f_tg_generate(nX,RX,LX,Rdata,stress,G0,G1,l1,l2,v_a,v_nc,rho8,Pref);

%%
% compare approx and exact tg
% 1. scatter plot of collapse time vs Rmax for 
figure;
plot(RX, tc_all, 'm^-', 'DisplayName', 'Experimental');
hold on;
plot(RX, tg_all, 'bo-', 'DisplayName', 'Analytical');
xlabel('R_{max}');
ylabel('Collapse time');
legend show;
title('Collapse Time vs R_{max}');
grid on;
residuals = tc_all - tg_all;
figure;
plot(RX, residuals, 'ko-');
xlabel('R_{max}');
ylabel('t_{c} - t_{g}');
title('Residuals between Experimental and Analytical Collapse Times');
grid on;
figure;
plot(tc_all, tg_all, 'ms');
hold on;
plot([min(tc_all), max(tc_all)], [min(tg_all), max(tg_all)], 'k--', 'DisplayName', 'y = x');
xlabel('Experimental t_{g} (T1X)');
ylabel('Analytical t_{g} (tg\_all)');
title('Parity Plot');
legend show;
axis equal;
grid on;
% Errors
errors = tc_all - tg_all;
% MAE
MAE = mean(abs(errors));
% RMSE
RMSE = sqrt(mean(errors.^2));
% R^2
SS_res = sum(errors.^2);
SS_tot = sum((tc_all - mean(tc_all)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('MAE: %.4f\n', MAE);
fprintf('RMSE: %.4f\n', RMSE);
fprintf('R^2: %.4f\n', R2);


function tg_all = f_tg_generate(nX,RX,LX,Rdata,stress,G0,G1,l1,l2,v_a,v_nc,rho8,Pref)
% for each {Rmax, Lambdamax} call f_tcol_calc_graded
tg_all = zeros(1,nX);
for i = 1:nX
    % WAIT IS THIS DIMENSIONAL? because f_tcol takes in and outputs NONDIM
    % Req_dim_i = RX(i) / LX(i); %dim
    % Req_i = Req_dim_i / RX(i);
    Req_i = 1/LX(i); %Req/R0
    R_i = Rdata{i}; %nondim already
    % nondim RX
    nd_R0 = RX(i) / RX(i);
    % dim Req
    d_Req = Req_i * RX(i);
    Ca = Pref/G0; Ca1 = Pref/G1;
    % el1 = l1/Req_i; el2 = l2/Req_i;
    el1 = l1/d_Req; el2 = l2/d_Req;
    %tg_all(i) = f_tcol_calc_graded(stress,Req_i,R_i,nd_R0,Ca,Ca1,Pref,el1,el2,v_a,v_nc,rho8);
    tg_all(i) = f_tg_calc(stress,Req_i,nd_R0,Ca,Ca1,Pref,el1,el2,v_a,v_nc,rho8);
end
end

function tc_all = f_tc_nondim(nX,RX,T1X,rho8,Pref)
tc_all = zeros(1,nX);
for i = 1:nX
    tchar = sqrt(rho8/Pref)*RX(i);
    tc_all(i) = T1X(i)/tchar;
end
end