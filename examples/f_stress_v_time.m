% clc; clear; close all;

% isgraded = [0 1 0];
% % positive gradient: soft-to-stiff
% G0 = [1000 1000 0];
% G1 = [0 10000 10000];
% % nondimensionalize
% Pref=101325;
% Ca_G0 = Pref./G0;
% Ca_G1 = Pref./G1;
% 
% % additional graded parameters
% l1 = 1.2;
% l2 = 2.2;
% v_a = 2;
% v_nc = 0.3;
% 
% % IMR parameters
% R0 = 100e-6;
% Req = R0/8; 
% tfin = 75E-6;
% t = linspace(0,tfin,1000);
% 
% % parameters
% nt = length(t);
% lr_N = nt;
% lr_length = 3;
% r_coord = linspace(0,lr_length,lr_N); %fixed Eulerian grid
% % last location to evaluate stress
% r_far = max(r_coord)*0.8;

% locations to evaluate
% nloc=4;
% location_labels = {sprintf('At R_0',R0),...
%                    sprintf('Beginning of graded region (l_1 = %.3f)',l1/R0),...
%                    sprintf('End of graded region (l_2 = %.3f)',l2/R0),...
%                    sprintf('Far field (r=%.3f)',r_far)};
% 
% % looping parameters
% nmat = length(isgraded); %number of materials
% stress_all = zeros(length(t), nloc, nmat); %[time, location, material]
% nt = length(t);
% 
% % for each material
% for i = 1:nmat
%     Ca = Ca_G0(i);
%     Ca1 = Ca_G1(i);
%     stress_all(:,:,i) = f_stress_v_time(nt,nloc,r_coord,R,R0,Req,t,Ca,Ca1,l1,l2,v_nc,v_a);
%     %stress_all(:,:,i) = stress_all(:,:,i) / -max(max(stress_all(:,:,i)),[],'all');
% end
% 
% % preparing plots
% labels = {'Soft','Graded','Stiff'};
% colors = {'r','c--','k--'};
% 
% % plotting
% for loc = 1:nloc
%     figure
%     hold on;
%     for mat = 1:nmat
%         plot(t,stress_all(:,loc,mat),colors{mat},'LineWidth',3,'DisplayName',labels{mat})
%     end
%     xlabel('$t / t_c$', 'Interpreter', 'Latex', 'FontSize', 20);
%     ylabel('$\tau_{rr} / \max(\tau_{rr})$','Interpreter','Latex','FontSize',20);
%     title(location_labels{loc},'FontSize',18)
%     legend show;
%     hold off;
% end

function [gstress_all] = f_stress_v_time(isgraded,nt,nloc,r_coord,R,R0,Req,r_far,Ca,Ca1,l1,l2,v_nc,v_a)

% for first two columns of stress_all
gstress_all = zeros(nt,nloc);
for time_idx = 1:nt
    % current time
    Rnow = R(time_idx);
    
    % if not graded, compute homogeneous stress
    if isgraded
        % get stress profile at current time
        [taurr_now,r1_now,r2_now] = f_graded_stress(r_coord,Rnow,Req,Ca,Ca1,l1,l2,v_nc,v_a);
    else
        [taurr_now,r1_now,r2_now] = f_ng_stress(r_coord,Rnow,Req,Ca,l1,l2);
    end

    % location of each tracer (at each time)
    r_eval = [R0,r1_now,r2_now,r_far];

    % interpolate stress at each location of each tracer (at each time)
    for loc=1:nloc
        r_current = r_eval(loc);
        warning('off','MATLAB:interp1:NaNstrip')
        gstress_all(time_idx,loc) = interp1(r_coord,taurr_now,r_current,'pchip');
    end
end
end