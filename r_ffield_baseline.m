
    %%
    fig_tend = t(end);
    tickrange= 0:5:fig_tend;
    lR = length(R);
    lr_max = 500;
    lr_N = 500;
    lr_length = 1;
    r_coord = ones(lR,lr_N).*logspace(-1,lr_length,lr_N);
    R_coord = (ones(lR,lr_N).*R);
    r0_coord = (r_coord.^3 - R_coord.^3 + (1/Lambda).^3).^(1/3);
    % calculating the shear
   
    % varsigmadot_r = f_gammadot_r(r_coord,R,U,lR,lr_N);
    % f_r = sf_carreau(v_nc,v_lambda,varsigmadot_r);
    tau_rr = 2/(3*Ca)*((r0_coord./r_coord).^4 - (r_coord./r0_coord).^2);% + 2*dGdhsnd/(Ca*3)*((r0_coord./r_coord).^3 - (r0_coord./r_coord).^3);
    [xcon,ycon] = meshgrid(t,r_coord(1,:));
    ycon = log10(ycon);
    clevels = 500;
    
    % feps_r = f_f_filter(f_r,lR,lr_N);
    % tau_r = 2*(feps_r./DRe+1./Re8).*varsigmadot_r;
    % min_tau_r = min(min(tau_r));
    % max_tau_r = max(max(tau_r));
    % ntau_r = tau_r/max_tau_r;
    
    
    % call filter
    N = length(tau_rr); M = lr_N;

    [tau_rr] = f_f_filter(tau_rr,r_coord,R_coord,N,M);
    
    % f contour figure
    figure(5)
    hold on;
    xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
    ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
    colormap jet;
    cbar = colorbar;
    cbar.Label.String = ['$\tau_{rr}/p_{\infty}$'];
    set(cbar,'TickLabelInterpreter','latex');
    pos = get(cbar,'Position');
    cbar.Label.Position = [pos(1) -0.11]; 
    cbar.Label.Rotation = 0;
    cbar.Label.Interpreter = 'latex';
    clim([-0.1 0]);
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
    contourf(xcon',ycon',tau_rr,clevels,'edgecolor','none')
    %saveas(gcf,'./figs/baseline/fcon_T','png')
    % %%
    % % shear stress contour figure
    % figure(6)
    % hold on;
    % xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20); 
    % ylabel('log$_{10}(\it{r}/R_o)$','Interpreter','Latex','FontSize',24);
    % colormap jet;
    % cbar = colorbar;
    % cbar.Label.String = '$\tau_{rr}/\mathrm{max}(\tau_{rr})$';
    % set(cbar,'TickLabelInterpreter','latex');
    % pos = get(cbar,'Position');
    % cbar.Label.Position = [1.5*pos(1) -1.2];
    % cbar.Label.Rotation = 0;
    % cbar.Label.Interpreter = 'latex';
    % clim([-1 1])
    % xlim([0 fig_tend]);
    % xticks(tickrange)
    % set(gcf,'color','w');
    % set(gca,'FontName','Times','FontSize',20);
    % set(gca,'TickLabelInterpreter','latex')
    % xa = gca;
    % xa.TickLength = [.015 .015];
    % xa.LineWidth = 1.5;
    % box on;
    % plot(t, log10(R),lm,'LineWidth',3); 
    % contourf(xcon',ycon',ntau_r,clevels,'edgecolor','none');
    % saveas(gcf,'./figs/baseline/taucon_T','png')