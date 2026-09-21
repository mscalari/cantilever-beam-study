function [t, U, V, A] = mdof_dc(K, M, C, u0, v0, h, nSteps)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Central Differences method -        %
%  matrix form for MDOF free vibration.%
%--------------------------------------%
%
% Matrix counterpart of dc_last.m.  Conditionally stable: the critical
% step of the DISCRETE system is h_cr = 2/w_max, where w_max is the
% HIGHEST frequency of the finite element mesh - not the fundamental
% frequency.  See viga_h_critico.m.
%
% Solves   M*udd + C*ud + K*u = 0 .
%
% INPUT
%   K, M, C  : [n x n] system matrices
%   u0, v0   : [n x 1] initial displacement and velocity
%   h        : time step                                          [s]
%   nSteps   : number of steps to advance
%
% OUTPUT
%   t        : [1 x nSteps+1] time vector
%   U, V, A  : [n x nSteps+1] displacement, velocity, acceleration
%
% Velocities and accelerations are central estimates and therefore
% require the neighbouring displacements; the first and last columns are
% filled with the exact initial state and a one-sided estimate.

    n = size(K,1);

    % --- Integration constants ---
    a0 = 1/h^2;
    a1 = 1/(2*h);
    a2 = 2*a0;
    a3 = 1/a2;

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);

    % --- Fictitious displacement at t = -h (Taylor expansion) ---
    uPrev = u - h*v + a3*a;

    % --- Effective mass, factorised once ---
    Meff = a0*M + a1*C;
    dMe  = decomposition(Meff);

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);
    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    uCurr = u;
    for k = 1:nSteps
        Reff  = -(K - a2*M)*uCurr - (a0*M - a1*C)*uPrev;
        uNext = dMe \ Reff;

        % Central estimates at the CURRENT step
        if k > 1
            V(:,k) = a1*(uNext - uPrev);
            A(:,k) = a0*(uNext - 2*uCurr + uPrev);
        end

        uPrev = uCurr;
        uCurr = uNext;
        U(:,k+1) = uCurr;
    end

    % One-sided closure of the last column
    if nSteps >= 1
        V(:,end) = (U(:,end) - U(:,max(end-1,1)))/h;
        A(:,end) = M \ (-C*V(:,end) - K*U(:,end));
    end

    t = (0:nSteps)*h;
end
