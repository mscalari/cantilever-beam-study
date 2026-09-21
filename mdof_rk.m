function [t, U, V, A] = mdof_rk(K, M, C, u0, v0, h, nSteps, order)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Explicit Runge-Kutta methods RK1 to %
%  RK4 - matrix form for MDOF systems. %
%--------------------------------------%
%
% Matrix counterpart of euler_last.m, euler_mod_last.m, RK3_last.m and
% RK4_last.m.  The second-order system is written in state-space form
%
%       yd = Ast*y ,        y = {u ; ud}
%
%       Ast = [   0        I   ]
%             [ -M\K     -M\C  ]
%
% and the approximation operator of an order-j explicit RK method is the
% truncated Taylor expansion of exp(h*Ast), exactly as in Eq. (14) of the
% paper:
%
%       Ahat = I + sum_{i=1}^{j} (h^i/i!) * Ast^i
%
% INPUT
%   K, M, C  : [n x n] system matrices
%   u0, v0   : [n x 1] initial displacement and velocity
%   h        : time step                                          [s]
%   nSteps   : number of steps to advance
%   order    : 1, 2, 3 or 4  (RK1 = Euler, RK2 = Modified Euler,
%                             RK3 = Heun, RK4 = classical RK4)
%
% OUTPUT
%   t        : [1 x nSteps+1] time vector
%   U, V, A  : [n x nSteps+1] displacement, velocity, acceleration
%
% Ahat is assembled once and applied repeatedly, so the cost per step is
% a single dense matrix-vector product of size 2n.

    if nargin < 8 || isempty(order), order = 4; end
    if ~ismember(order, 1:4)
        error('mdof_rk:order','order must be 1, 2, 3 or 4.');
    end

    n = size(K,1);

    % --- State matrix ---
    MiK = M \ K;
    MiC = M \ C;
    Ast = [ zeros(n), eye(n) ;
            -MiK    , -MiC   ];

    % --- Approximation operator: truncated Taylor series ---
    Ahat = eye(2*n);
    term = eye(2*n);
    for i = 1:order
        term = (term*(h*Ast))/i;
        Ahat = Ahat + term;
    end

    % RK step counts are set by explicit-scheme STABILITY, not accuracy,
    % so nSteps here can reach into the millions (see viga_h_critico.m) -
    % at n = ndof this is where U, V and A stop being "small" arrays.
    % Estudos B-F only ever ask for [t,U] (they only need the tip row of
    % the displacement history), so V and A are computed and stored ONLY
    % when the caller actually requests them via nargout, exactly
    % preserving the original 4-output behaviour (e.g. viga_main_simulacao.m's
    % [t,U,V,A] = mdof_rk(...) calls, Estudo A's reference case) while
    % cutting peak memory by up to 3x for the new callers.
    computeV = nargout >= 3;
    computeA = nargout >= 4;

    U = zeros(n, nSteps+1);
    if computeV, V = zeros(n, nSteps+1); end
    if computeA, A = zeros(n, nSteps+1); end

    y = [u0(:); v0(:)];
    U(:,1) = y(1:n);
    if computeV, V(:,1) = y(n+1:end); end
    if computeA, A(:,1) = -MiK*y(1:n) - MiC*y(n+1:end); end

    % --- Time loop ---
    for k = 1:nSteps
        y = Ahat*y;
        U(:,k+1) = y(1:n);
        if computeV, V(:,k+1) = y(n+1:end); end
        if computeA, A(:,k+1) = -MiK*y(1:n) - MiC*y(n+1:end); end
    end

    t = (0:nSteps)*h;
end
