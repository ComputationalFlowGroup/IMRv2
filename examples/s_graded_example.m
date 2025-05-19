clc; clear; close all;


Rst = linspace(0.5,10,200);

G1 = 1;
G3 = 2;

l1 = 1.2;
l2 = 2;
a = 1.2;
n = 0.3;
x1 = @(Rst) (1+(Rst.^3-1)./(1+l1).^3).^(1/3);
x2 = @(Rst) (1+(Rst.^3-1)./(1+l2).^3).^(1/3);

ycy = @(x,Rst) (G3+(G1-G3)*(1+((x-x1(Rst))./(x2(Rst)-x)).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
% ype = @(x,Rst) (G3+(G1-G3)*asinh((x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ympe = @(x,Rst) (G3+(G1-G3)*log((x-x1(Rst))./(x2(Rst)-x)+1)./((x-x1(Rst))./(x2(Rst)-x)).^a).*(1./x.^5+1./x.^2);
% ycr = @(x,Rst) (G3+(G1-G3)*1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).*(1./x.^5+1./x.^2);
% yscr = @(x,Rst) (G3+(G1-G3)*1./(1+(x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ymcr = @(x,Rst) (G3+(G1-G3)*(1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).^a).*(1./x.^5+1./x.^2);

reltol = 1e-8;
abstol = 1e-8;
S2 = zeros(length(Rst),1);
for i = 1:length(Rst)
    rst = Rst(i);
    S2(i) = 2*integral(@(x) ycy(x,rst),x1(rst),x2(rst),...
            'RelTol',reltol,'AbsTol',abstol);
end


S1 = (G1/2)*(1./Rst.^4 + 1./Rst - (1./x1(Rst).^4 + 1./x1(Rst)));
S3 = -(G3/2)*(5 - 4./x2(Rst) - 1./x2(Rst).^4);
SG1 = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);
SG3 = -(G3/2)*(5 - 4./Rst - 1./Rst.^4);
figure(1)
hold on;
% plot(Rst,S1,'m')
% plot(Rst,S2,'k')
% plot(Rst,S3,'b')
plot(Rst,S1+S2+S3,'-.g','LineWidth',3)
plot(Rst,SG1,'r','LineWidth',3)
plot(Rst,SG3,'k--','LineWidth',3)
ylim([-6 6])
xlabel('$R_{\mathrm{max}}/R_{0}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$S$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
tickrange= 0:2:10;
xticks(tickrange)
box on;
saveas(gcf,'./fig_graded_stress_integral','png')