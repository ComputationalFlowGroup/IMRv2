fig_tend = 5;
tickrange= 0:5:fig_tend;
lR = length(R);
lr_max = 100;
lr_N = 200;
lr_length = 1;
r_coord = ones(lR,lr_N).*logspace(-0.5,lr_length,lr_N);
R_coord = (ones(lR,lr_N).*R);
%r0_coord = (r_coord.^3 - R_coord.^3 + (Req/R_coord).^3).^(1/3);

% calculating the shear
%varsigmadot_r = f_gammadot_r(r_coord,R,U,lR,lr_N);
[xcon,ycon] = meshgrid(t,r_coord(1,:));
ycon = log10(ycon);
clevels = 50;
% feps_r = f_f_filter(f_r,lR,lr_N);
% tau_r = 2*(feps_r./DRe+1./Re8).*varsigmadot_r;
%ntau_r = tau_r/max_tau_r;
taurr = gradedstress(lR,R_coord,Req,Ca,Ca1,l1,l2,v_nc,v_a);
maxtaurr = max(taurr);
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
set(cbar,'TickLabelInterpreter','latex');
pos = get(cbar,'Position');
cbar.Label.Position = [pos(1) -0.1];
cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';
clim([0 1]);
%clim([0.1 0]);
xlim([0 fig_tend]);
xticks(tickrange)
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
box on;
plot(t, log10(R),'LineWidth',3,'Color','k');
contourf(xcon',ycon',ntaurr,clevels,'edgecolor','none')
%saveas(gcf,'./figs/baseline/fcon_T','png')

function [taurr] = gradedstress(lR,R_coord,Req,Ca,Ca1,l1,l2,v_nc,v_a)
    reltol = 1e-8;
    abstol = 1e-8;
    Se1 = zeros(1,lR); Se2 = zeros(1,lR); Se3 = zeros(1,lR); taurr = zeros(1,lR);
    for i = 1:lR
        Rst = R_coord(i)/Req;
        x1 = (1 + (Rst.^3 - 1)./(l1^3)).^(1/3);
        x2 = (1 + (Rst.^3 - 1)./(l2^3)).^(1/3);
    
        f_cy = @(x) (l2*((x.^3 - 1)/(Rst.^3 - 1)).^(1/3) - 1)/(1-l1*((x.^3 - 1)/(Rst.^3 - 1)).^(1/3));
        g = @(x) (1/Ca + (1/Ca1 - 1/Ca)*(1+f_cy(x(i)).^v_a).^((v_nc-1)/v_a)).*((1./x.^5) + (1./x.^2));
    
        %Sv = - 4/Re8*Rdot/R - 6*intfnu*iDRe;
        
        Se1(i) = (1/(2*Ca))*(1/(Rst.^4) + 4/Rst - (1./x1.^4 + 4./x1));
        Se2(i) = 2*integral(@(x) g(x),x1,x2,'RelTol',reltol,'AbsTol',abstol);
        Se3(i) = -(1/(2*Ca1))*(5 - 4./x2 - 1./x2.^4);
        
        taurr(i) = Se1(i) + Se2(i) + Se3(i);
    end
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