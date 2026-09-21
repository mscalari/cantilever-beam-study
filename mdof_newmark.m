function [t, U, V, A] = mdof_newmark(K, M, C, u0, v0, h, nSteps, alpha, delta)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Newmark method - matrix form for    %
%  MDOF free vibration.                %
%--------------------------------------%
%
% Matrix counterpart of newmark_last.m.  Default parameters
% alpha = 1/4, delta = 1/2 give the unconditionally stable constant
% average acceleration (trapezoidal) rule used in the paper.
%
% Solves   M*udd + C*ud + K*u = 0 .
%
% INPUT
%   K, M, C       : [n x n] system matrices
%   u0, v0        : [n x 1] initial displacement and velocity
%   h             : time step                                     [s]
%   nSteps        : number of steps to advance
%   alpha, delta  : Newmark parameters (optional)
%
% OUTPUT
%   t             : [1 x nSteps+1] time vector
%   U, V, A       : [n x nSteps+1] displacement, velocity, acceleration
%
% The effective stiffness is constant in a linear analysis, therefore it
% is factorised ONCE before the time loop.

    if nargin < 8 || isempty(alpha), alpha = 1/4; end
    if nargin < 9 || isempty(delta), delta = 1/2; end

    n = size(K,1);

    % --- Integration constants (Bathe, Finite Element Procedures) ---
    a0 = 1/(alpha*h^2);
    a1 = delta/(alpha*h);
    a2 = 1/(alpha*h);
    a3 = 1/(2*alpha) - 1;
    a4 = delta/alpha - 1;
    a5 = (h/2)*(delta/alpha - 2);
    a6 = h*(1 - delta);
    a7 = delta*h;

    % --- Effective stiffness, factorised once ---
    Keff = K + a0*M + a1*C;
    dKe  = decomposition(Keff);

    % --- Allocation and initial state ---
    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);

    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    % --- Time loop ---
    for k = 1:nSteps
        Reff = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a);
        uNew = dKe \ Reff;
        aNew = a0*(uNew - u) - a2*v - a3*a;
        vNew = v + a6*a + a7*aNew;

        u = uNew;  v = vNew;  a = aNew;
        U(:,k+1) = u;  V(:,k+1) = v;  A(:,k+1) = a;
    end

    t = (0:nSteps)*h;
end
