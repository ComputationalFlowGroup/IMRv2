function f_stress_v_time
r_fixed = [min(R), l1/R0, l2/R0, max(r_coord(:))*0.8];
location_labels = {sprintf('Bubble Surface (R_{min} = %.3f',min(R)),...
                   sprintf('Beginning of graded region (l_1 = %.3f)',l1/R0),...
                   sprintf('End of graded region (l_2 = %.3f)',l2/R0),...
                   sprintf('Far field (r=%.3f)',max(r_coord(:))*0.8)};
gstress_all = zeros(nt,length(r_fixed));
ungstress_all = zeros(nt,length(r_fixed));

% for loc = 1:length(r_fixed)
%     r_current = r_fixed(loc);
%     for time_idx = 1:nt
%         r_spat = r_coord(time_idx,:);
%         gstress_all(time_idx,loc) = interp1(r_spat,taurr(:,time_idx),r_current,'spline');
%         ungstress_all(time_idx,loc) = interp1(r_spat,taurr1(:,time_idx),r_current,'spline');
% 
%     end
%         figure
%         hold on;
%         plot(t,gstress_all(:,loc),'r--','LineWidth',3,'DisplayName','Graded')
%         plot(t,ungstress_all(:,loc),'b--','LineWidth',3,'DisplayName','Ungraded')
%         %formatting plot
%         xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
%         ylabel('$\tau_{rr}/p_{\infty}$','Interpreter','Latex','FontSize',20);
%         title(sprintf('Stress Comparison at %s',location_labels{loc},'Interpreter','Latex','FontSize',18))
%         legend show;
% end
% hold off;
%%
nloc=2; %locations to evaluate
gstress_all = zeros(nt,nloc);
ungstress_all = zeros(nt,nloc);
r_current = r_fixed(loc);
for time_idx = 1:nt
    r_spat = r_coord(time_idx,:);
    stress_slice = taurr(time_idx,:); unstress_slice = taurr1(time_idx,:);
    Rnow = R(time_idx);

    r1 = ((l1/R0)^3 + Rnow.^3 - Req^3).^(1/3);
    r2 = ((l2/R0)^3 + Rnow.^3 - Req^3).^(1/3);
    r_far = max(r_spat)*0.8;
    r_eval = [r1,r2,r_far];

    for loc = 1:nloc
        r_current = r_eval(loc);
        gstress_all(time_idx,loc) = interp1(r_spat,taurr(:,time_idx),r_current,'spline');
        ungstress_all(time_idx,loc) = interp1(r_spat,taurr1(:,time_idx),r_current,'spline');

        figure
        hold on;
        plot(t,gstress_all(:,loc),'r--','LineWidth',3,'DisplayName','Graded')
        plot(t,ungstress_all(:,loc),'b--','LineWidth',3,'DisplayName','Ungraded')
        %formatting plot
        xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
        ylabel('$\tau_{rr}/p_{\infty}$','Interpreter','Latex','FontSize',20);
        title(sprintf('Stress Comparison at %s',location_labels{loc},'Interpreter','Latex','FontSize',18))
        legend show;
    end
end
hold off;
end