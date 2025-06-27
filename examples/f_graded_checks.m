function f_graded_checks(R,R0,Req,t,l1,l2,v_nc,v_a,rho8)
%% check 1: recover classical Rayleigh collapse time
G = 0.1; G1 = 0.1;
Pref = 101325;
Ca = Pref/G; Ca1 = Pref/G1;
%%[Ca,Ca1] = shearmod(G,G1);
%graded stress should approach classical
tg = f_tcol_calc_graded(1,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
%classical Rayleigh collapse time
trc = 0.915*Req*sqrt(rho8/Pref);
if tg == trc
    disp("check 1 correct: Rayleigh collapse time recovered")
else
    disp("check 1 incorrect: Rayleigh collapse time not recovered")
end

%% check 2: collapse not energetically permissible
G = 1E9; G1 = 1E9;
Pref = 101325;
Ca = Pref/G; Ca1 = Pref/G1;
%%[Ca,Ca1] = shearmod(G,G1);
tg = f_tcol_calc_graded(1,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
if abs(imag(tg)) > 0
    disp("check 2 expected: collapse not energetically permissible")
else
    disp("check 1 incorrect: collapse occurred")
end

%% check 3: uniform shear modulus of homogeneous elastic media
G = 1E4; G1 = G;
Ca = Pref/G; Ca1 = Pref/G1;
%[Ca,Ca1] = shearmod(G,G1);
tg = f_tcol_calc_graded(1,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
%classical Rayleigh collapse time
trc = 0.915*Req*sqrt(rho8/Pref);
if tg == trc
    disp("check 1 correct: Rayleigh collapse time recovered with only G")
else
    disp("check 1 incorrect: Rayleigh collapse time not recovered with only G")
end

%% check 4: uniform shear modulus with G1
G1 = 1E4; G = G1;
Ca = Pref/G; Ca1 = Pref/G1;
%[Ca,Ca1] = shearmod(G,G1);
tg = f_tcol_calc_graded(1,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8);
%classical Rayleigh collapse time
trc = 0.915*Req*sqrt(rho8/Pref);
if tg == trc
    disp("check 1 correct: Rayleigh collapse time recovered with only G1")
else
    disp("check 1 incorrect: Rayleigh collapse time not recovered with only G1")
end

%% is taurr increasing with r?
nt = 1;
lr_N = 500;
r_coord = linspace(0.1,3,lr_N);

nloc=4; %locations to evaluate
for time_idx = 1:nt
    Rnow = R(time_idx);
    r0_coord = r_coord.^3 - Rnow^3 + Req^3;
    taurr_check = (r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2);
    figure;
    plot(r_coord,taurr_check)
    xlabel('r');
    ylabel('\tau_{rr}')
end

%% at R_max
Rnow = max(R);
r_coord = linspace(0.1,3,500);
r0_coord = (r_coord.^3 - Rnow^3 + Req^3).^3;
valid = r0_coord > 0 & isreal(r0_coord); %where r0_coord values are negative (inside bubble)
r0_coord_correct = r0_coord.*valid;

f_cy = (l2 - r0_coord_correct) ./ (r0_coord_correct - l1);
taurr_base = (r0_coord_correct.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord_correct.^2);

G0 = 1000; G1 = 5000;
v_a = 2; v_nc = 0.3;
taurr1 = (2*G0/3) * taurr_base;
taurr2 = ( G0 + (G1-G0)*(1+f_cy.^v_a).^((v_nc-1)/v_a) ).*taurr_base;
taurr3 = (2*G1/3) * taurr_base;
taurrtot = taurr1 + taurr2 + taurr3;
figure;
plot(r_coord,taurrtot)
xlabel('r'); ylabel('\tau_{rr}');
%%
    function [Ca,Ca1] = shearmod(G,G1)
        Pref = 101325;
        Ca = Pref/G; Ca1 = Pref/G1;
    end
end