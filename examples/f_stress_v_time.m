function f_stress_v_time(R,R0,Req,t,Ca,Ca1,l1,l2,v_nc,v_a)
nt = length(t);
lr_N = 500;
%lr_length = 5;
%r_coord = ones(lR,lr_N).*logspace(-1,lr_length,lr_N);
r_coord = ones(nt,lr_N).*linspace(0.1,3, lr_N);
[taurr1,taurr,~,~] = f_gradedstress(r_coord,R,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);

nloc=3; %locations to evaluate
gstress_all = zeros(nt,nloc);
ungstress_all = zeros(nt,nloc);
for time_idx = 1:nt
    r_spat = r_coord(time_idx,:);
    %stress_slice = taurr(time_idx,:); unstress_slice = taurr1(time_idx,:);
    Rnow = R(time_idx);

    r1 = ((l1/R0)^3 + Rnow.^3 - Req^3).^(1/3);
    r2 = ((l2/R0)^3 + Rnow.^3 - Req^3).^(1/3);
    r_far = max(r_spat)*0.8;
    r_eval = [r1,r2,r_far];

    for loc = 1:nloc
        r_current = r_eval(loc);
        gstress_all(time_idx,loc) = interp1(r_spat,taurr(time_idx,:),r_current,'spline');
        ungstress_all(time_idx,loc) = interp1(r_spat,taurr1(time_idx,:),r_current,'spline');
    end
end
location_labels = {sprintf('Beginning of graded region (l_1 = %.3f)',l1/R0),...
                   sprintf('End of graded region (l_2 = %.3f)',l2/R0),...
                   sprintf('Far field (r=%.3f)',max(r_coord(:))*0.8)};
for loc=1:nloc
    figure
    hold on;
    plot(t,gstress_all(:,loc),'r--','LineWidth',3,'DisplayName','Graded')
    plot(t,ungstress_all(:,loc),'b--','LineWidth',3,'DisplayName','Ungraded')
    %formatting plot
    xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
    ylabel('$\tau_{rr}/p_{\infty}$','Interpreter','Latex','FontSize',20);
    title(sprintf('Stress Comparison at %s',location_labels{loc},'Interpreter','Latex','FontSize',18))
    legend show;
    hold off;
end
end
%%
function f_stress_time(R, R0, Req, t, Ca, Ca1, l1, l2, v_nc, v_a)
    nt = length(t);
    lr_N = 500;

    % Spatial domain: fixed in space, uniform for all time steps
    r_vec = linspace(0.1, 3, lr_N);         % [1 x lr_N]
    r_coord = repmat(r_vec, nt, 1);         % [nt x lr_N]

    % Compute stress data at all time steps and all radii
    [taurr1, taurr, ~, ~] = f_gradedstress(r_coord, R, Req, R0, Ca, Ca1, l1, l2, v_nc, v_a);

    % Locations: beginning, end of graded region, far field
    nloc = 3;
    gstress_all = zeros(nt, nloc);
    ungstress_all = zeros(nt, nloc);

    for time_idx = 1:nt
        r_spat = r_coord(time_idx, :);       % [1 x lr_N]
        Rnow = R(time_idx);

        % Compute spatial locations for this time
        r1 = ((l1/R0)^3 + Rnow^3 - Req^3)^(1/3);
        r2 = ((l2/R0)^3 + Rnow^3 - Req^3)^(1/3);
        r_far = max(r_spat) * 0.8;
        r_eval = [r1, r2, r_far];

        % Interpolate at fixed spatial locations
        for loc = 1:nloc
            r_current = r_eval(loc);
            gstress_all(time_idx, loc) = interp1(r_spat, taurr(time_idx, :), r_current, 'spline', 'extrap');
            ungstress_all(time_idx, loc) = interp1(r_spat, taurr1(time_idx, :), r_current, 'spline', 'extrap');
        end
    end

    % Labels for legend
    location_labels = {
        sprintf('Beginning of graded region (l_1 = %.3f)', l1/R0), ...
        sprintf('End of graded region (l_2 = %.3f)', l2/R0), ...
        sprintf('Far field (r = %.3f)', max(r_coord(:)) * 0.8)
    };

    % Plot for each location
    for loc = 1:nloc
        figure;
        hold on;
        plot(t, gstress_all(:, loc), 'r--', 'LineWidth', 3, 'DisplayName', 'Graded');
        plot(t, ungstress_all(:, loc), 'b--', 'LineWidth', 3, 'DisplayName', 'Ungraded');
        xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
        ylabel('$\tau_{rr}/p_{\infty}$', 'Interpreter', 'Latex', 'FontSize', 20);
        title(sprintf('Stress Comparison at %s', location_labels{loc}), 'Interpreter', 'Latex', 'FontSize', 18);
        legend show;
        hold off;
    end
end

