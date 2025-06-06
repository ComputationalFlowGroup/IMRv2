fig_tend = 4;
%tickrange= 0:2:fig_tend;
tickrange = ceil(min(t)) : floor(max(t));
lR = length(R);
%lr_max = 100;
lr_N = 500;
lr_length = 1.5;
nt = length(t);
r_coord = ones(lR,lr_N).*logspace(-0.5,lr_length,lr_N);
%r_coord = ones(lR,lr_N).*linspace(0.1,3, lr_N);
%R_coord = (ones(lR,lr_N).*R);

[xcon,ycon] = meshgrid(t,r_coord(1,:));
ycon = log10(ycon);
clevels = 1000;
[taurr,r1,r2] = gradedstress(r_coord',R',Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);
maxtaurr = 1; %max(abs(taurr(:))); %1; %max(max(abs(taurr)));
ntaurr = taurr/maxtaurr;

% f contour figure
figure(2)
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
cbar.Label.String = '$\tau_{rr}/p_{\infty}$';
set(cbar,'TickLabelInterpreter','latex','FontSize',18);
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -0.9];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
%clim([min(ntaurr(:)) max(ntaurr(:))]);
%clim([-0.75 .25])
xlim([0 fig_tend]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
ya = gca;
ya.TickLength = [.015 .015];
ya.LineWidth = 1.5;
box on;
plot(t, log10(R),'LineWidth',3,'Color','k');
%plot(t, R,'LineWidth',3,'Color','k');
%plot(t, R,'LineWidth',3,'Color','k');
contourf(xcon,ycon,log(abs(ntaurr)),clevels,'edgecolor','none')
hold on;
plot(t,log10(r1'),'b--','LineWidth',2,'DisplayName','l1')
plot(t,log10(r2'),'r--','LineWidth',2,'DisplayName','l2')
%ylim([log10(min(r_coord(:))), log10(max(r_coord(:)))])

function [taurr,r1,r2] = gradedstress(r_coord,R,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a)
   aa = r_coord.^3 - R.^3 + Req^3;
   aa = (1./(1-(aa<0))).*aa;
   r0_coord = real((aa).^(1/3));
    
   f_cy = (l2 - r0_coord) ./ (r0_coord - l1);
        
   taurr1 = (2/(3*Ca))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr3 = (2/(3*Ca1))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));

   taurr = taurr1 + taurr2 + taurr3;

   r1 = ((l1/R0)^3 + R.^3 - Req^3).^(1/3);
   r2 = ((l2/R0)^3 + R.^3 - Req^3).^(1/3);
end