Rst = linspace(0.2,10,200); %R_max/R_0

G0 = 500;
G1 = 1000;
l1 = 1.5; %nondim

%% bicomposite material stress integral as a function of max stretch ratio
Rst = linspace(0.2,10,200); 
Rmax = 100E-6; %Rst.*Req;
Req = Rst./Rmax;

x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
S0 = (G0/2)*(1./Rst.^4 + 4./Rst - (1./x1(Rst).^4 + 4./x1(Rst)));
S1 = -(G1/2)*(5 - 4./x1(Rst) - 1./x1(Rst).^4);
SG0 = -(G0/2)*(5 - 4./Rst - 1./Rst.^4);
SG1 = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);
Pref = 101325;
%figure(1)
figure
hold on;
plot(Rst,abs(SG0)/G1,'r','LineWidth',3) %no abs
plot(Rst,abs(SG1)/G1,'k--','LineWidth',3) %no abs
plot(Rst,abs(S0+S1)/G1,'-.g','LineWidth',3) %no nabs
ylim([-1 10])
xlabel('$R_{\mathrm{max}}/R_{0}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$|S|/G1$','Interpreter','Latex','FontSize',24);
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

Rdot = @(Rst) cos(Rst);
Rstdot = Rdot(Rst) ./ Req;
x1dot = @(Rst) Rst.^2 .* Rstdot(Rst) ./ x1(Rst).^2;
S0dot = (G0/2)*(4.*x1dot(Rst)./x1(Rst).^5 + 4.*x1dot(Rst)./x1(Rs).^2 -4.*Rstdot./Rst.^5 - 4.*Rstdot./Rst.^2);
S1dot = -(G1/2)*(-4.*x1dot(Rst)./x1(Rst).^5 - 4.*x1dot(Rst)./x1(Rst).^2);
SG0dot = -(G0/2)*(- 4.*Rstdot(Rst)./Rst.^5 - 4.*Rstdot(Rst)./Rst.^2);
SG1dot = -(G1/2)*(- 4.*Rstdot(Rst)./Rst.^5 - 4.*Rstdot(Rst)./Rst.^2);
figure
hold on;
plot(Rst,abs(SG0dot)/G1,'r','LineWidth',3) %no abs
plot(Rst,abs(SG1dot)/G1,'k--','LineWidth',3) %no abs
plot(Rst,abs(S0dot+S1dot)/G1,'-.g','LineWidth',3) %no nabs
ylim([-1 10])
xlabel('$R_{\mathrm{max}}/R_{0}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$|\dot{S}|/G1$','Interpreter','Latex','FontSize',24);
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
