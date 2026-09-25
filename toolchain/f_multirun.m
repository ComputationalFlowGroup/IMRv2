% file f_multirun.m
% brief contains function f_multirun

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function sweeps the material parameters of the Newtonian,
% neo-Hookean Kelvin-Voigt, UCM and Oldroyd-B models with the nodal solver
% and returns every radius history as [t, R] (dimensionless), in cell
% arrays shaped like the parameter grids. The parameters are given in SI
% units and converted with f_dim2nondim; the relaxation times are
% lambda1 = lambda1_factor*mu/G and lambda2 = lambda2_factor*mu/G. The
% sweeps run with f_imr_sweep on a parallel pool.
function [newtonian, neohook, voigt, ucm, oldb] = f_multirun()
    
    N = 8^4;
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'src'));
    addpath(fullfile(root,'toolchain'));
    sweep = @(base, varargin) f_imr_sweep(@f_imr_nodal, base, varargin{:}, ...
        'transform', @f_si_inputs, 'parallel', true);
    
    % one parameter: Newtonian and neo-Hookean models
    newtonian = f_tr(sweep({'stress',1,'g',0}, 'mu', logspace(-6,-1,N)))';
    neohook = f_tr(sweep({'stress',1}, 'g', logspace(-6,-1,N)))';
    
    % two parameters: Kelvin-Voigt model
    n2 = round(N^(1/2));
    voigt = f_tr(sweep({'stress',1}, 'mu', logspace(-6,-1,n2), ...
        'g', logspace(-6,-1,n2)));
    
    % three parameters: upper-convected Maxwell model
    n3 = ceil(N^(1/3));
    ucm = f_tr(sweep({'stress',5,'lambda2',0}, 'mu', logspace(-6,-1,n3), ...
        'g', logspace(-6,-1,n3), 'lambda1_factor', logspace(-6,-1,n3)));
    
    % four parameters: Oldroyd-B model
    n4 = round(N^(1/4));
    oldb = f_tr(sweep({'stress',5}, 'mu', logspace(-6,-1,n4), ...
        'g', logspace(-6,-1,n4), 'lambda1_factor', logspace(-6,-1,n4), ...
        'lambda2_factor', logspace(-6,-1,n4)));
    
end

% SI inputs of one run: the relaxation times from their factors, then the
% conversion to the solver's dimensionless inputs
function varin = f_si_inputs(varargin)
    names = lower(varargin(1:2:end));
    vals = varargin(2:2:end);
    keep = true(size(names));
    extra = {};
    for k = 1:numel(names)
        if any(strcmp(names{k}, {'lambda1_factor','lambda2_factor'}))
            % mu/G from the last mu and g given
            mu = vals{find(strcmp(names,'mu'),1,'last')};
            G = vals{find(strcmp(names,'g'),1,'last')};
            extra = [extra, {strrep(names{k},'_factor',''), vals{k}*mu/G}]; %#ok<AGROW>
            keep(k) = false;
        end
    end
    pairs = reshape(varargin, 2, []);
    pairs = pairs(:, keep);
    varin = f_dim2nondim(pairs{:}, extra{:});
end

% sweep results as [t, R] per run ([] for a failed run)
function c = f_tr(res)
    c = reshape(arrayfun(@(r) [r.out{1}, r.out{2}], res, ...
        'UniformOutput', false), size(res));
end
