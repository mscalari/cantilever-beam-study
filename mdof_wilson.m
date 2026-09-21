function [t, U, V, A] = mdof_wilson(K, M, C, u0, v0, h, nSteps, theta)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Wilson-theta method - matrix form   %
%  for MDOF free vibration.            %
%--------------------------------------%
%
% Matrix counterpart of wilson_last.m, with the default theta = 1.4 used
% throughout the paper.  In free vibration the extrapolated load at
% t + theta*h vanishes, so only the inertia and damping contributions
% remain on the right-hand side.
%
% Solves   M*udd + C*ud + K*u = 0 .
%
% INPUT
%   K, M, C  : [n x n] system matrices
%   u0, v0   : [n x 1] initial displacement and velocity
%   h        : time step                                          [s]
%   nSteps   : number of steps to advance
%   theta    : Wilson parameter (optional, default 1.4)
%
% OUTPUT
%   t        : [1 x nSteps+1] time vector
%   U, V, A  : [n x nSteps+1] displacement, velocity, acceleration
%
% NOTE - the Wilson-theta scheme is known to OVERSHOOT in the first
% steps.  On a finite element mesh the overshoot is driven by the modes
% with the largest w*h, so it becomes dramatic for an impulsive initial
% condition, which excites the whole spectrum.

    if nargin < 8 || isempty(theta), theta = 1.4; end

    n = size(K,1);

    % --- Integration constants ---
    a0 = 6/(theta*h)^2;
    a1 = 3/(theta*h);
    a2 = 2*a1;
    a3 = (theta*h)/2;
    a4 = a0/theta;
    a5 = -a2/theta;
    a6 = 1 - 3/theta;
    a7 = h/2;
    a8 = h^2/6;

    % --- Effective stiffness, factorised once ---
    Keff = K + a0*M + a1*C;
    dKe  = decomposition(Keff);

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);
    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    % --- Time loop ---
    for k = 1:nSteps
        Reff  = M*(a0*u + a2*v + 2*a) + C*(a1*u + 2*v + a3*a);
        uTh   = dKe \ Reff;                 % displacement at t + theta*h

        aNew = a4*(uTh - u) + a5*v + a6*a;
        vNew = v + a7*(aNew + a);
        uNew = u + h*v + a8*(aNew + 2*a);

        u = uNew;  v = vNew;  a = aNew;
        U(:,k+1) = u;  V(:,k+1) = v;  A(:,k+1) = a;
    end

    t = (0:nSteps)*h;
end
