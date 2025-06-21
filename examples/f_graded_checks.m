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


%%
    function [Ca,Ca1] = shearmod(G,G1)
        Pref = 101325;
        Ca = Pref/G; Ca1 = Pref/G1;
    end

end