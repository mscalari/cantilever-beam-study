function p = viga_params()
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  Geometric and material properties   %
%  of the cantilever beam model.       %
%--------------------------------------%
%
% OUTPUT
%   p : struct with the beam properties
%
% Rectangular steel cantilever, clamped at x = 0, free at x = L.
% Slenderness L/h = 100 places the model well inside the
% Euler-Bernoulli (shear-free) regime.

    % --- Material ---
    p.E   = 210e9;      % Young's modulus            [Pa]
    p.rho = 7850;       % Mass density               [kg/m^3]

    % --- Geometry ---
    p.L = 1.00;         % Length                     [m]
    p.b = 0.050;        % Cross-section width        [m]
    p.h = 0.010;        % Cross-section height       [m]

    % --- Derived cross-section properties ---
    p.A = p.b*p.h;              % Area                        [m^2]
    p.I = p.b*p.h^3/12;         % Second moment of area       [m^4]
    p.slenderness = p.L/p.h;    % Slenderness ratio           [-]

    % --- Discretization ---
    p.N = 10;           % Number of 2-node beam elements

    % --- Damping ---
    p.xi_target = 0.2;  % Target modal damping ratio (same as the SDOF study)
    p.mode_i    = 1;    % First anchor mode for Rayleigh damping
    p.mode_j    = 3;    % Second anchor mode for Rayleigh damping

    % --- Impulse initial condition ---
    p.V0 = 1.0;         % Initial tip transverse velocity   [m/s]

    % ======================================================================
    %  ONDE FICAM OS PARAMETROS NUMERICOS DOS ESTUDOS B-F
    % ======================================================================
    % Este arquivo descreve apenas O QUE A VIGA E: material, geometria,
    % malha, o amortecimento de referencia do Estudo A e a CI de impulso.
    % Mexer aqui muda o PROBLEMA FISICO sendo resolvido.
    %
    % Tudo que define COMO a simulacao e rodada - passo de tempo, numero
    % de periodos, quais casos cada estudo varre (xi, modos excitados,
    % razoes h/T_modo, coeficientes de Rayleigh alternativos, niveis de
    % reducao modal) e o visual das figuras - esta em UM UNICO arquivo:
    %
    %       viga_config_estudos.m
    %
    % Esses campos moravam aqui como "Extended study parameters" ate
    % 20/09/2026; foram movidos para la porque definem o EXPERIMENTO, nao
    % a viga. Equivalencias, caso algum script antigo ainda os use:
    %
    %       p.xi_list         ->  cfg.B.xi_list
    %       p.modal_modes     ->  cfg.C.modos
    %       p.modal_modes_D   ->  cfg.D.modos
    %       p.hD_ratios       ->  cfg.D.hRatios
    %       p.rayleigh_alt    ->  cfg.E.rayleigh
    %       p.reduction_modes ->  cfg.F.reduction_modes
    %
    % com cfg = viga_config_estudos().

    % --- Sanity check ---
    if p.slenderness < 10
        warning('viga_params:slenderness', ...
            'L/h = %.1f < 10: Euler-Bernoulli theory may be inaccurate.', ...
            p.slenderness);
    end
end
