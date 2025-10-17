function f_graded_checks(R,R0,Req,t,l1,l2,v_nc,v_a,rho8)
%% check 1: recover classical Rayleigh collapse time
G = 0.1; G1 = 0.1;
Pref = 101325;
Ca = Pref/G; Ca1 = Pref/G1;
%%[Ca,Ca1] = shearmod(G,G1);
% NEED TO RERUN IMR TO GET NEW R SINCE CHANGED  [t,R,~] = f_imr_fd(varin{:},'Nt',16,'Mt',64);
addpath('../src/characterization/')
%graded stress should approach classical
tg = f_tcol_calc_graded(1,Req/R0,R,R0/R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
%classical Rayleigh collapse time
trc = 0.91468*R0/R0*sqrt(rho8/Pref); %had to nondim Rmax
if tg == trc
    disp("check 1 correct: Rayleigh collapse time recovered")
else
    disp("check 1 incorrect: Rayleigh collapse time not recovered")
end

%% check 2: collapse not energetically permissible
G = 1E9; G1 = 1E9;
Pref = 101325;
Ca = Pref/G; Ca1 = Pref/G1;
%%[Ca,Ca1] = shearmod(G,G1);
% NEED TO RERUN IMR TO GET NEW R SINCE CHANGED G
tg = f_tcol_calc_graded(1,Req/R0,R,R0/R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
if abs(imag(tg)) > 0 %% WHAT IS THIS TRYING TO SAY? DO YOU MEAN IF IMAGINARY WHY NEED ABS HERE?
    disp("check 2 expected: collapse not energetically permissible")
else
    disp("check 1 incorrect: collapse occurred")
end

%% check 3: uniform shear modulus of homogeneous elastic media
G1 = 1E4; G3 = G;
Ca = Pref/G1; Ca1 = Pref/G3;
%[Ca,Ca1] = shearmod(G,G1);

% time check
tg = f_tcol_calc_graded(1,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8)
%classical Rayleigh collapse time
trc = 0.915*Req*sqrt(rho8/Pref)
if tg == trc
    disp("check 1 correct: Rayleigh collapse time recovered with only G")
else
    disp("check 1 incorrect: Rayleigh collapse time not recovered with only G")
end

% no graded region, should restore nH stress integral and stress
%%
% stress integral check
G1 = 1E4; G3 = G;
Ca = Pref/G1; Ca1 = Pref/G3;
Rst = linspace(0.2,10,200); %R_max/R_0
l1 = 1.5; l2 = 3; a = 2; n = 0.3;
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2
ycy = @(x,Rst) (G3+(G1-G3)*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
      (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
reltol = 1e-8; abstol = 1e-8;
S2 = zeros(1,length(Rst));
for i = 1:length(Rst)
    rst = Rst(i);
    S2(i) = 2*integral(@(x) ycy(x,rst),x1(rst),x2(rst),...
            'RelTol',reltol,'AbsTol',abstol);
end
S1 = (G1/2)*(1./Rst.^4 + 4./Rst - (1./x1(Rst).^4 + 4./x1(Rst)));
S3 = -(G3/2)*(5 - 4./x2(Rst) - 1./x2(Rst).^4);
SG1 = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);
SG3 = -(G3/2)*(5 - 4./Rst - 1./Rst.^4);
figure(1)
hold on;
% plot(Rst,S1,'m')
% plot(Rst,S2,'k')
% plot(Rst,S3,'b')
plot(Rst,SG1/G1,'r','LineWidth',3)
plot(Rst,SG3/G1,'k--','LineWidth',3)
plot(Rst,(S1+S2+S3)/G1,'-.g','LineWidth',3)
ylim([-5 5])
xlabel('$R_{\mathrm{max}}/R_{0}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$S/G_1$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
tickrange= 0:2:10;
xticks(tickrange)
tickrange= -5:2:5;
yticks(tickrange)
box on;


%%
% looping parameters
nmat = length(isgraded); %number of materials
stress_all = zeros(nt,nmat); %[time, location, material]
R_all = zeros(nt,nmat);
max_stress_each = zeros([1,nmat]);

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
    % DO l1 and l2 NEED TO BE DIVIDED BY REQ???


    % store all R data for each material
    R_all(:,i) = R;
    % compute stresses
    stress_all(:,:,i) = f_stress_v_time(isgraded(i),nt,nloc,R,R0/R0,Req/R0,r_far,Ca,Ca1,l1/R0,l2/R0,v_nc,v_a);
    % per material normalization of stress (temporal features)
    % max_stress_each(i) = max(max(abs(stress_all(:,:,i)))); % across each location and each time
    % stress_all(:,:,i) = stress_all(:,:,i) / max_stress_each(i);
    % alternatively: global normalization of stress (relative strength)
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


%% check 4: uniform shear modulus with G1
G1 = 1E4; G = G1;
Ca = Pref/G; Ca1 = Pref/G1;
%[Ca,Ca1] = shearmod(G,G1);
tg = f_tcol_calc_graded(1,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
%classical Rayleigh collapse time
trc = 0.915*Req*sqrt(rho8/Pref);
if tg == trc
    disp("check 1 correct: Rayleigh collapse time recovered with only G1")
else
    disp("check 1 incorrect: Rayleigh collapse time not recovered with only G1")
end

%% is taurr increasing with r?
nt = 1;
lr_N = 500;
r_coord = linspace(0.1,3,lr_N);

nloc=4; %locations to evaluate
for time_idx = 1:nt
    Rnow = R(time_idx);
    r0_coord = r_coord.^3 - Rnow^3 + Req^3;
    taurr_check = (r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2);
    figure;
    plot(r_coord,taurr_check)
    xlabel('r');
    ylabel('\tau_{rr}')
end

%% at R_max
Rnow = max(R);
r_coord = linspace(0.1,3,500);
r0_coord = (r_coord.^3 - Rnow^3 + Req^3).^3;
valid = r0_coord > 0 & isreal(r0_coord); %where r0_coord values are negative (inside bubble)
r0_coord_correct = r0_coord.*valid;

f_cy = (l2 - r0_coord_correct) ./ (r0_coord_correct - l1);
taurr_base = (r0_coord_correct.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord_correct.^2);

G0 = 1000; G1 = 5000;
v_a = 2; v_nc = 0.3;
taurr1 = (2*G0/3) * taurr_base;
taurr2 = ( G0 + (G1-G0)*(1+f_cy.^v_a).^((v_nc-1)/v_a) ).*taurr_base;
taurr3 = (2*G1/3) * taurr_base;
taurrtot = taurr1 + taurr2 + taurr3;
figure;
plot(r_coord,taurrtot)
xlabel('r'); ylabel('\tau_{rr}');
%%
    function [Ca,Ca1] = shearmod(G,G1)
        Pref = 101325;
        Ca = Pref/G; Ca1 = Pref/G1;
    end
end