% file f_pvsat.m
% brief contains function f_pvsat

% brief This function computes the saturated water vapor pressure per
% (A. Preston 2004, Thesis), Pv = c1*exp(c2/T). The constants c1 (Pa) and
% c2 (K) are returned as optional outputs so that this file is their single
% source.
function [ Pv, c1, c2 ] = f_pvsat( T )
    c1 = 1.17e11;
    c2 = -5200;
    Pv = c1*exp(c2./(T));
end
