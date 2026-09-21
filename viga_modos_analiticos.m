function [wn, betaL, phi] = viga_modos_analiticos(p, nModes)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Analytical modes of the continuous  %
%  Euler-Bernoulli cantilever beam.    %
%--------------------------------------%
%
% Used ONLY to validate the finite element mesh (see
% viga_main_convergencia.m).  The reference solution employed in the
% time-domain comparisons is the exact modal solution of the
% semi-discrete FE system (viga_solucao_exata.m).
%
% INPUT
%   p       : struct from viga_params()
%   nModes  : number of modes requested
%
% OUTPUT
%   wn      : [nModes x 1] natural circular frequencies       [rad/s]
%   betaL   : [nModes x 1] roots of cos(bL)*cosh(bL) + 1 = 0  [-]
%   phi     : function handle phi(x, n) -> mode shape n at x
%
% The characteristic equation of the clamped-free beam is
%       cos(beta*L)*cosh(beta*L) + 1 = 0
% and the natural frequencies follow from
%       wn = (beta*L)^2 * sqrt(E*I/(rho*A*L^4)).

    if nargin < 2, nModes = 6; end

    f = @(x) cos(x).*cosh(x) + 1;

    betaL = zeros(nModes,1);
    for k = 1:nModes
        % Asymptotic estimate (2k-1)*pi/2 is excellent for k >= 2 and
        % acceptable as a bracketing seed for k = 1.
        guess = (2*k-1)*pi/2;
        lo = max(guess - 0.9, 1e-3);
        hi = guess + 0.9;
        % Guarantee a sign change before calling fzero
        while f(lo)*f(hi) > 0
            lo = max(lo - 0.2, 1e-3);
            hi = hi + 0.2;
        end
        betaL(k) = fzero(f, [lo hi]);
    end

    wn = betaL.^2 * sqrt(p.E*p.I/(p.rho*p.A*p.L^4));

    % Mass-normalised-in-shape mode functions (unit scaling at the tip
    % is NOT enforced; only the shape matters for validation).
    phi = @(x,k) modeShape(x, betaL(k), p.L);
end

function y = modeShape(x, bL, L)
    beta  = bL/L;
    sigma = (sinh(bL) + sin(bL)) / (cosh(bL) + cos(bL));
    y = cosh(beta*x) - cos(beta*x) - sigma*(sinh(beta*x) - sin(beta*x));
end
