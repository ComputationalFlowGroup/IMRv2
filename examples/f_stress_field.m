clc;
fig_tend = 4;
%tickrange= 0:2:fig_tend;
tickrange = ceil(min(t)) : floor(max(t));
lR = length(R);
%lr_max = 100;
lr_N = 500;
lr_length = 1;
nt = length(t);
r_coord = ones(lR,lr_N).*logspace(-0.5,lr_length,lr_N);
R_coord = (ones(lR,lr_N).*R);

% calculating the shear
%varsigmadot_r = f_gammadot_r(r_coord,R,U,lR,lr_N);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
ycon = log10(ycon);
clevels = 1000;
% feps_r = f_f_filter(f_r,lR,lr_N);
% tau_r = 2*(feps_r./DRe+1./Re8).*varsigmadot_r;
%ntau_r = tau_r/max_tau_r;
[taurr,r1,r2] = gradedstress(r_coord',R_coord',Req,Ca,Ca1,l1,l2,v_nc,v_a);
maxtaurr = 1; %max(abs(taurr(:))); %1; %max(max(abs(taurr)));
ntaurr = taurr/maxtaurr;
% N = length(taurr);
% M = lr_N;
% %feps_r = f_f_filter(taurr,lR,lr_N);
% %feps_r = f_f_filters(taurr,r_coord,R,N,M);
% tau_rr = f_f_filter(taurr,N,M);
% maxtau_rr = max(tau_rr);
% ntau_rr = tau_rr/maxtau_rr;

% f contour figure
figure(2)
hold on;
xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
colormap jet;
cbar = colorbar;
%cbar.Label.String = '$m(\dot{\varsigma})$';
cbar.Label.String = '$\tau_{rr}/p_{\infty}$';
set(cbar,'TickLabelInterpreter','latex','FontSize',18);
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -0.9];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
%clim([min(ntaurr(:)) max(ntaurr(:))]);
clim([-0.75 .25])
%clim([0.1 0]);
xlim([0 fig_tend]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
ytickrange = round(linspace(log10(min(r_coord(:))),log10(max(r_coord(:))),3),2);
yticks(ytickrange)
ya = gca;
ya.TickLength = [.015 .015];
ya.LineWidth = 1.5;
box on;
plot(t, log10(R),'LineWidth',3,'Color','k');
%plot(t, R,'LineWidth',3,'Color','k');
contourf(xcon,ycon,ntaurr,clevels,'edgecolor','none')
hold on;
plot(t,log10(r1'),'b--','LineWidth',2,'DisplayName','l1')
plot(t,log10(r2'),'r--','LineWidth',2,'DisplayName','l2')
%saveas(gcf,'./figs/baseline/fcon_T','png')
%ylim([0 1])
ylim([log10(min(r_coord(:))), log10(max(r_coord(:)))])

function [taurr,r1,r2] = gradedstress(r_coord,R_coord,Req,Ca,Ca1,l1,l2,v_nc,v_a)
   aa = r_coord.^3 - R_coord.^3 + Req^3;
   aa = (1./(1-(aa<0))).*aa;
   r0_coord = real((aa).^(1/3));
    
   f_cy = (l2 - r0_coord) ./ (r0_coord - l1);
        
   taurr1 = (2/(3*Ca))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr3 = (2/(3*Ca1))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));

   taurr = taurr1 + taurr2 + taurr3;

   r1 = (l1^3 + R_coord.^3 - Req^3).^(1/3);
   r2 = (l2^3 + R_coord.^3 - Req^3).^(1/3);
end

% shear as a function of r (radial coordinate) calculation
function sofr = f_gammadot_r(r,R,Rdot,N,M)
    
    sofr = zeros(N,M);
    for i = 1:N
        for j = 1:M
            if r(i,j) >=  R(i,1)
                sofr(i,j) = -2*Rdot(i,1)*(R(i,1)^2)/((r(i,j))^3);
            else
                sofr(i,j) = NaN;
            end
        end
    end
end
%is this supposed to be the velocity? shouldn't it be Rdot*R^2/r^2 then?

function [feps_r] = f_f_filter(f_r,N,M)
    eps = 0.01;
    feps_r = zeros(size(f_r));
    for i = 1:N
        for j = 1:M
            if f_r(i,j) > 1-eps
                feps_r(i,j) = NaN;
            elseif f_r(i,j) < eps
                feps_r(i,j) = NaN;
            else
                feps_r(i,j) = f_r(i,j);
            end
        end
    end   
end

function [taurr] = f_f_filters(taurr,r,R,N,M)
    for i = 1:N
        for j = 1:M
            if r(i,j) < R(i,j)
                taurr(i,j) = NaN;
            end
        end
    end
end