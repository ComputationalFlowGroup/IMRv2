clc; clear; close;

% generate R v t curves
addpath('src/forward_solver/');
%addpath('./common/')

% options
R0 = 100e-6;
Req = R0/8; 
tfin = 75E-6;
tvector = linspace(0,tfin,1000);

collapse = 0;
radial = 1;
vapor = 1;
bubtherm = 1;
medtherm = 0;
masstrans = 0;
stress = 1;

% graded parameters
graded = 1;
v_nc = 0.3; 
v_a = 2;
l1 = 1.2*R0; 
l2 = 2.2*R0;
G0 = 1E3;
G1 = 10E3;

kappa = 1.4;
T8 = 298.15;
rho8 = 1064;
mu = 5E-2;
%%
varin = {'progdisplay',0,...
    'radial',radial,...
    'bubtherm',bubtherm,...
    'tvector',tvector,...
    'vapor',vapor,...
    'medtherm',medtherm,...
    'masstrans',masstrans,...
    'method',23,...
    'stress',stress,...
    'collapse',collapse,...
    'mu',mu,...
    'g',G0,...
    'graded',graded,...
    'g1',G1,...
    'l1',l1,...
    'l2',l2,...
    'v_a',v_a,...
    'v_nc',v_nc,...
    'lambda1',1e-7,...
    'lambda2',0,...
    'alphax',1e-3,...
    'r0',R0,...
    'req',Req,...
    'kappa',kappa,...
    't8',T8,...
    'rho8',rho8};
% generate R v t data
[t,R,~] = f_imr_fd(varin{:},'Nt',16,'Mt',64);

figure
hold on;
plot(t,R,'bx')
ylim([0 1.2])

% addpath('./src/characterization/')
% tg = f_tcol_calc_graded(stress,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
% plot(tg,R,'r--')
% hold off;

%% graded and ungraded stress vs time

addpath('examples/');
% G0 = [1000 1000 0];
% G1 = [0 10000 10000];
legendName = {'G0 = 1000', 'G0 = 1000, G1 = 10000', 'G1 = 10000'};

isgraded = [0 1 0];
% positive gradient: soft-to-stiff
G0 = [1000 1000 5000]; 
G1 = [0 5000 0];
% nondimensionalize
Pref=101325;
Ca_G0 = Pref./G0;
Ca_G1 = Pref./G1;

% % parameters
nt = length(tvector);
lr_N = nt;
lr_length = 6;
r_coord = linspace(0,lr_length,lr_N); %fixed Eulerian grid
% last location to evaluate stress
r_far = max(r_coord)*0.8;

% locations to evaluate
nloc=4;
location_labels = {sprintf('At R_0 = %.3f',R0/R0),...
                   sprintf('Beginning of graded region (l_1 = %.3f)',l1/R0),...
                   sprintf('End of graded region (l_2 = %.3f)',l2/R0),...
                   sprintf('Far field (r=%.3f)',r_far)};

% looping parameters
nmat = length(isgraded); %number of materials
stress_all = zeros(nt, nloc, nmat); %[time, location, material]
R_all = zeros(nt,nmat);

% for each material
for i = 1:nmat
    Ca = Ca_G0(i);
    Ca1 = Ca_G1(i);
    % run solver to get R
    varin = {'progdisplay',0,...
    'radial',radial,...
    'bubtherm',bubtherm,...
    'tvector',tvector,...
    'vapor',vapor,...
    'medtherm',medtherm,...
    'masstrans',masstrans,...
    'method',23,...
    'stress',stress,...
    'collapse',collapse,...
    'mu',mu,...
    'g',G0(i),...
    'graded',isgraded(i),...
    'g1',G1(i),...
    'l1',l1,...
    'l2',l2,...
    'v_a',v_a,...
    'v_nc',v_nc,...
    'lambda1',1e-7,...
    'lambda2',0,...
    'alphax',1e-3,...
    'r0',R0,...
    'req',Req,...
    'kappa',kappa,...
    't8',T8,...
    'rho8',rho8};
    [t,R,~] = f_imr_fd(varin{:},'Nt',150,'Mt',150);

    % store all R data for each material
    R_all(:,i) = R;
    % compute stresses
    stress_all(:,:,i) = f_stress_v_time(isgraded(i),nt,nloc,r_coord,R,R0/R0,Req,r_far,Ca,Ca1,l1,l2,v_nc,v_a);
    stress_all(:,:,i) = stress_all(:,:,i) / -max(max(abs(stress_all(:,:,i))),[],'all');
    % normalized R0?
end

% preparing plots
labels = {'Soft','Graded','Stiff'};
colors = {'r','c--','k--'};

figure
hold on;
plot(t,R_all(:,1),'LineWidth',2,'DisplayName',labels{1})
plot(t,R_all(:,2),'LineWidth',2,'DisplayName',labels{2})
plot(t,R_all(:,3),'LineWidth',2,'DisplayName',labels{3})
legend show;
ylim([0 1.2])
hold off;

% plotting
for loc = 1:nloc
    figure
    hold on;
    for mat = 1:nmat
        plot(t,stress_all(:,loc,mat),colors{mat},'LineWidth',2,'DisplayName',labels{mat})
    end
    xlim([0 lr_length]);
    xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
    ylabel('$\tau_{rr} / \max(\tau_{rr})$','Interpreter','Latex','FontSize',20);
    title(location_labels{loc},'FontSize',18)
    legend show;
    hold off;
end
% filename = sprintf('./stress_%s.png', legendName{i});
% saveas(gcf,filename)
