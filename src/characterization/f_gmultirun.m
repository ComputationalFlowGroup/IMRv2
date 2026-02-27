% file f_gmultirun.m
% brief contains function f_gmultirun

% brief This function generates latin hypercube samples for R_max and
% Lambda_max for pIMR. Note all inputs received in f_gmultirun must be
% dimensional (except v_a and v_nc, which are scalars)
function [data,Rdata] = f_gmultirun(n,R_min,R_max,Lambda_min,Lambda_max,...
    G0,G1,l1,l2,v_a,v_nc)

format long;
% preallocate output
data = zeros(n,3);
Rdata = cell(1,n);
% Generate Latin Hypercube Samples
lhs = lhsdesign(n, 2);  % n samples in 2 dimensions

% Scale samples so instead of [0,1], they're from [R_min,R_max] and
% [Lambda_min,Lambda_max] respectively
Rmax_vals = R_min + (R_max - R_min) * lhs(:, 1);
Lambdamax_vals = Lambda_min + (Lambda_max - Lambda_min) * lhs(:, 2);

% Plot
figure;
scatter(Rmax_vals.*1e6, Lambdamax_vals, 50, 'k', 'filled');
xlabel('$R_{\mathrm{max}}$ [$\mu$m]', 'FontName','TimesNewRoman','FontSize',20,'Interpreter', 'latex');
ylabel('$\Lambda_{\mathrm{max}}$', 'FontName','TimesNewRoman','FontSize',20,'Interpreter', 'latex');
%titlestr = sprintf('Latin Hypercube Sampling of $R_{\mathrm{max}}$, $\Lambda_{\mathrm{max}}$ $(n = %d)$',n);
%title(titlestr, 'Interpreter', 'latex');
%title(sprintf('Latin Hypercube Sampling of \\R_{\\mathrm{max}}, \\Lambda_{\\mathrm{max}} (n = %d)',n),'Interpreter', 'latex');
xlim([R_min.*1e6 - 10, R_max.*1e6 + 10]);
ylim([Lambda_min - 0.5, Lambda_max + 0.5]);
grid on;
axis square;
fname = sprintf('syn_dat_dist_n%d.png',n);
saveas(gcf, fname);
fprintf('Figure saved as %s\n',fname);

counter = 0;
% call IMR, with KM, keep vapor on but turn off medtherm and masstherm
% turn off mu, lambda1, lambda2, alphax
addpath('../forward_solver/');
for i = 1:length(Rmax_vals)
    % options
    R0 = Rmax_vals(i);
    Req = R0/Lambdamax_vals(i);
    
    kappa = 1.4;
    T8 = 298.15;
    rho8 = 1064;
    mu = 0;
    lambda1 = 0;
    lambda2 = 0;
    alphax = 0;
    Pref = 101325;
    
    % compute tfin
    tfin = 1.25*R0*sqrt(rho8/Pref);
    tvector = linspace(0,tfin,1000);
    
    collapse = 0;
    radial = 1;
    vapor = 1;
    bubtherm = 0;
    %1; can't directly compare bubble pressure/temp in energy balance
    medtherm = 0;
    masstrans = 0;
    stress = 1;
    
    % graded parameters (nondim)
    el1 = l1/Req;
    el2 = l2/Req;
    
    varin = {'progdisplay',0,...
        'radial',radial,...
        'bubtherm',bubtherm,...
        'tvector',tvector,...
        'vapor',vapor,...
        'medtherm',medtherm,...
        'masstrans',masstrans,...
        'method',23,...
        'stress',stress,...
        'collapse',collapse,...
        'mu',mu,...
        'g',G0,...
        'graded',1,...
        'g1',G1,...
        'l1',el1,...
        'l2',el2,...
        'v_a',v_a,...
        'v_nc',v_nc,...
        'lambda1',lambda1,...
        'lambda2',lambda2,...
        'alphax',alphax,...
        'r0',R0,...
        'req',Req,...
        'kappa',kappa,...
        't8',T8,...
        'rho8',rho8};
    
    % generate R v t data
    [t,R,~] = f_imr_fd(varin{:},'Nt',16,'Mt',64);
    [~,idx_minR] = min(R);
    Rdata{i} = R;
    % need to dimensionalize time
    tchar = sqrt(rho8/Pref)*R0;
    tc = t(idx_minR)*tchar;
    data(i,:) = [Rmax_vals(i), Lambdamax_vals(i), tc];
    counter = counter + 1;
    
end
% save the output
save('data.mat','data');
save('Rdata.mat','Rdata');
% Outputs
% nX = size(data,1);   % # of experiments
% RX = data(:,1);      % All Rmax
% LX = data(:,2);      % All amplification Lmax
% T1X = data(:,3);     % All collapse time t1
% Ri = Rdata{i};       % spatial info of bubble wall
end
