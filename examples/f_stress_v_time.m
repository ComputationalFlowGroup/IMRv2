function [gstress_all, ungstress_all] = f_stress_v_time(R,R0,Req,t,Ca,Ca1,l1,l2,v_nc,v_a)
nt = length(t);
lr_N = 500;
%lr_length = 5;
%r_coord = ones(lR,lr_N).*logspace(-1,lr_length,lr_N);
r_coord = linspace(0.1,3, lr_N); %fixed Eulerian grid
%[taurr1,taurr,~,~] = f_gradedstress(r_coord,R,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);

nloc=4; %locations to evaluate
gstress_all = zeros(nt,nloc);
ungstress_all = zeros(nt,nloc);
for time_idx = 1:nt
    %r_spat = r_coord(time_idx,:);
    %stress_slice = taurr(time_idx,:); unstress_slice = taurr1(time_idx,:);
    Rnow = R(time_idx);

    [taurr1_now,taurr_now,~,~] = f_gradedstress(r_coord,Rnow,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);

    r1 = ((l1/R0)^3 + Rnow.^3 - Req^3).^(1/3);
    r2 = ((l2/R0)^3 + Rnow.^3 - Req^3).^(1/3);
    r_far = max(r_coord)*0.8;
    r_eval = [R0,r1,r2,r_far];

    for loc = 1:nloc
        r_current = r_eval(loc);
        gstress_all(time_idx,loc) = interp1(r_coord,taurr_now,r_current,'pchip');
        ungstress_all(time_idx,loc) = interp1(r_coord,taurr1_now,r_current,'pchip');
    end
end
location_labels = {sprintf('At R_0',R0),...
                   sprintf('Beginning of graded region (l_1 = %.3f)',l1/R0),...
                   sprintf('End of graded region (l_2 = %.3f)',l2/R0),...
                   sprintf('Far field (r=%.3f)',max(r_coord(:))*0.8)};
for loc=1:nloc
    figure
    hold on;
    plot(t,gstress_all(:,loc),'r--','LineWidth',3,'DisplayName','graded')
    plot(t,ungstress_all(:,loc),'b--','LineWidth',3,'DisplayName','ungraded')
    %formatting plot
    xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
    ylabel('$\tau_{rr} / \mathrm{max}(\tau_{rr})$','Interpreter','Latex','FontSize',20);
    %title(sprintf('Stress Comparison at %s',location_labels{loc},'Interpreter','Latex','FontSize',18))
    title(location_labels{loc},'FontSize',18)
    legend show;
    hold off;
end
end