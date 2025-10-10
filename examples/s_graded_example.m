clc; clear; close all;
% figure 1 JMPS

Rst = linspace(0.2,10,200); %R_max/R_0

G1 = 0.5;
G3 = 1;

l1 = 1.5;
l2 = 3;
a = 2.5;
n = 0.3;
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2

% ycy = @(x,Rst) (1+(G1-1)*(1+((x-x1(Rst))./(x2(Rst)-x)).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
ycy = @(x,Rst) (G1+(G3-G1)*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
      (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
% ype = @(x,Rst) (G3+(G1-G3)*asinh((x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ympe = @(x,Rst) (G3+(G1-G3)*log((x-x1(Rst))./(x2(Rst)-x)+1)./((x-x1(Rst))./(x2(Rst)-x)).^a).*(1./x.^5+1./x.^2);
% ycr = @(x,Rst) (G3+(G1-G3)*1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).*(1./x.^5+1./x.^2);
% yscr = @(x,Rst) (G3+(G1-G3)*1./(1+(x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ymcr = @(x,Rst) (G3+(G1-G3)*(1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).^a).*(1./x.^5+1./x.^2);

f = @(x,Rst) (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
      (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3));
ytanh = @(x,Rst) (G1+(G3-G1)*(1/2)*(1+tanh(a*(2*f(x,Rst) -1)))).*(1./x.^5+1./x.^2);

%
reltol = 1e-8;
abstol = 1e-8;
S2 = zeros(1,length(Rst));
for i = 1:length(Rst)
    rst = Rst(i);
    S2(i) = 2*integral(@(x) ytanh(x,rst),x1(rst),x2(rst),...
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
ylabel('$S/G_0$','Interpreter','Latex','FontSize',24);
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
%saveas(gcf,'./fig_graded_stress_integral','png')
%%
% figure 2 JMPS
% dimensional, in terms of lambda
tau = @(x) (2/3)*(1./x.^4 - x.^2);
%gtau = @(x,Rst) G3+((G1-G3)*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
%        (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a));
gtau = @(x,Rst) (G1+(G3-G1)*(1/2)*(1+tanh(a*(2*f(x,Rst) -1))));

gtau1 = zeros(size(Rst));
gtau2 = zeros(size(Rst));
gtau3 = zeros(size(Rst));
for i = 1:length(Rst)
    rst = Rst(i);
    gtau1(i) = G1*tau(rst);
    gtau2(i) = gtau(rst,rst)*tau(rst);
    gtau3(i) = G3*tau(rst);
end

figure(2)
hold on;
plot(Rst,gtau1,'r','LineWidth',3)
plot(Rst,gtau3,'k--','LineWidth',3)
plot(Rst,gtau2,'c--','LineWidth',3)
%plot(Rst,gtau1+gtau2+gtau3,'-.g','LineWidth',3)
ylim([-5 5])
xlabel('$R_{\mathrm{max}}/R_{0}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$\tau_{rr}/p_{\infty}$','Interpreter','Latex','FontSize',24);
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
%saveas(gcf,'./fig_graded_stress_Rst','png')

%% stress integral as a function of Lambda
% Lambda = linspace(0.1,1,200); %R(t)/R0
% Rt = Lambda * Req;
% Rst = R0/Req; %Rmax/Req
R0=1; %initial radius Rmax
Req=0.125; % equilibrium radius R0
Rst = R0/Req; % Rmax/Req (constant)
Lambda = linspace(1,Rst,200); % R(t)/Req, decreasing from Rmax to Req

G1 = 0.5;
G3 = 1;
l1 = 1.5;
l2 = 3;
a = 2.5;
n = 0.3;
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2


% ycy = @(x,Rst) (G1+(G3-G1)*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
%       (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
% 
% f = @(x,Rst) (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
%       (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3));
% ytanh = @(x,Rst) (G1+(G3-G1)*(1/2)*(1+tanh(a*(2*f(x,Rst) -1)))).*(1./x.^5+1./x.^2);

%
reltol = 1e-8;
abstol = 1e-8;
S2 = zeros(1,length(Lambda));
for i = 1:length(Lambda)
    lambda= Rst(i);
    S2(i) = 2*integral(@(x) ytanh(x,rst),x1(rst),x2(rst),...
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
ylabel('$S/G_0$','Interpreter','Latex','FontSize',24);
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
%saveas(gcf,'./fig_graded_stress_integral','png')
%% stress as a function of ratio l1/l2
% at an instant, how does graded material's stress in graded region change with ratio
G0 = 1000; G1 = 5000;
v_nc = 0.3; v_a = 2;
R0 = 300e-5; %Rmax
Req = R0/3; %R_0

l1_l2 = linspace(0.1,1,100);
r_coord = 1; %linspace(0.1,3,500); %fixed Eulerian grid
Rnow = R0;
r0coord = (r_coord.^3 - Rnow.^3 + Req^3);
% valid = r0coord > 0; %where r0_coord values are negative (inside bubble)
% r0_coord = (r0coord.*valid).^(1/3);
comp = zeros(size(l1_l2));
for i = 1:length(l1_l2)
    % r0_l2 = r0_coord / l2; %this still relies on a FIXED value of l2
    % fcy = (1 - r0_l2) ./ (r0_l2 - l1_l2);
    % mcy = @(fcy) (1 + fcy.^v_a).^((v_nc-1)/v_a);

    % at a fixed location - midpoint of graded region
    gamma = l1_l2(i);
    beta = 0.5*(gamma + 1);
    fcy_mid = ((1 - beta) ./ (beta - gamma));
    %mcy_mid = mcy(fcy_mid);
    mcy_mid = (1 + fcy_mid.^v_a).^((v_nc-1)/v_a);

   tau = 2/3 * ((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr1 = G0*tau;
   taurr2 = (G0 + (G1 - G0)*mcy_mid).*tau;
  % taurr3 = G1*tau;

  % taurr = taurr1 + taurr2 + taurr3;

   %still have to fix l_2 to use l_1 or vice versa
   l2 = 1;
   l1 = gamma*l2;
   comp(i) = 1./ ((r0_coord./l1).* (((1/(G1-G0)).*((3/2).*(taurr2./taurr1) -1)).^(v_a/(v_nc-1)) -1).^(1/v_a) + (r0_coord/l1) - 1);
end
figure;
plot(l1_l2,comp,'k','LineWidth',3)
xlabel('$l_1 / l_2$', 'Interpreter', 'Latex', 'FontSize', 20);
%ylabel('$\tau_{rr} / \mathrm{max}(\tau_{rr})$','Interpreter','Latex','FontSize',20);
ylabel('$\tau_{rr}$','Interpreter','Latex','FontSize',20);

%% compare other non-Newtonian graded functions

Rst = linspace(0.2,10,200); %R_max/R_0

G0 = 0.5;
G1 = 1;
%should G1 and G0 be switched since G0 < G1 for positive gradient
l1 = 1.5;
l2 = 3;
a = 2;
n = 0.3;
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2

G = G1 + (G1 - G0);

% ycy = @(x,Rst) (1+(G1-1)*(1+((x-x1(Rst))./(x2(Rst)-x)).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
ycy = @(x,Rst) (G*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
      (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
% NEED TO CONSIDER HOW THESE CHANGE FOR EACH GRADED FUNCTION WITH l1 AND l2
ype = @(x,Rst) (G*asinh((x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
ympe = @(x,Rst) (G*log((x-x1(Rst))./(x2(Rst)-x)+1)./((x-x1(Rst))./(x2(Rst)-x)).^a).*(1./x.^5+1./x.^2);
ycr = @(x,Rst) (G*1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).*(1./x.^5+1./x.^2);
yscr = @(x,Rst) (G*1./(1+(x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
ymcr = @(x,Rst) (G*(1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).^a).*(1./x.^5+1./x.^2);
functions = {@ycy, @ympe, @ycr, @yscr, @ymcr};


reltol = 1e-8;
abstol = 1e-8;
%S2 = zeros(1,length(Rst));
S2 = zeros(length(Rst),length(functions));
for i = 1:length(Rst)
    rst = Rst(i);
    for f = 1:length(functions)
        func = functions{f};
        S2(i,f) = 2*integral(@(x) func(x,rst),x1(rst),x2(rst),...
            'RelTol',reltol,'AbsTol',abstol);
    end
end

%%
S1 = (G1/2)*(1./Rst.^4 + 4./Rst - (1./x1(Rst).^4 + 4./x1(Rst)));
S3 = -(G3/2)*(5 - 4./x2(Rst) - 1./x2(Rst).^4);
SG1 = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);
SG3 = -(G3/2)*(5 - 4./Rst - 1./Rst.^4);
figure(1)
hold on;
% plot(Rst,S1,'m')
% plot(Rst,S2,'k')
% plot(Rst,S3,'b')
plot(f, ype, 'r', 'LineWidth', 2); hold on;
plot(f, ympe, 'g', 'LineWidth', 2);
plot(f, ycr, 'b', 'LineWidth', 2);
plot(f, yscr, 'm', 'LineWidth', 2);
plot(f, ymcr, 'c', 'LineWidth', 2);
plot(Rst,SG1/G1,'r','LineWidth',3)
plot(Rst,SG3/G1,'k--','LineWidth',3)
plot(Rst,(S1+S2+S3)/G1,'-.g','LineWidth',3)

% Plotting the results of different functions (S2)
figure;
hold on;
for f = 1:length(functions)
    plot(Rst, S2(:, f), 'DisplayName', func2str(functions{f}));
end
xlabel('Rst');
ylabel('Integral Result');
legend show;
title('Integral Results for Different Functions');
hold off;

% Plotting the functions themselves for a specific rst value
rst_example = Rst(1); % Example rst to plot for
x_vals = linspace(x1(rst_example), x2(rst_example), 100); % x-values for plotting

figure;
hold on;
for f = 1:length(functions)
    y_vals = arrayfun(@(x) functions{f}(x, rst_example), x_vals);
    plot(x_vals, y_vals, 'DisplayName', func2str(functions{f}));
end
xlabel('x');
ylabel('Function Value');
legend show;
title('Functions for a Specific Rst');
hold off;