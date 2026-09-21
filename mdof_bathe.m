function [t, U, V, A] = mdof_bathe(K, M, C, u0, v0, h, nSteps)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Bathe method (two sub-steps) -      %
%  matrix form for MDOF free vibration.%
%--------------------------------------%
%
% Matrix counterpart of bathe_last.m, with the standard splitting ratio
% gamma = 1/2 (two equal sub-steps).
%
%   Sub-step 1 : trapezoidal rule over h/2, equilibrium at t + h/2
%   Sub-step 2 : three-point Euler backward, equilibrium at t + h
%
% Solves   M*udd + C*ud + K*u = 0 .
%
% INPUT
%   K, M, C  : [n x n] system matrices
%   u0, v0   : [n x 1] initial displacement and velocity
%   h        : FULL time step (both sub-steps included)            [s]
%   nSteps   : number of full steps to advance
%
% OUTPUT
%   t        : [1 x nSteps+1] time vector (full steps only)
%   U, V, A  : [n x nSteps+1] displacement, velocity, acceleration
%
% Both effective stiffness matrices are constant in a linear analysis
% and are factorised once.  With gamma = 1/2 they are DIFFERENT, so two
% factorisations are required - this is the well known extra cost of the
% Bathe scheme relative to a single-step method.

    n = size(K,1);

    % --- Sub-step 1 constants (trapezoidal rule over h/2) ---
    b0 = 16/h^2;
    b1 = 4/h;
    b2 = 8/h;

    % --- Sub-step 2 constants (three-point Euler backward) ---
    c1 =  1/h;
    c2 = -4/h;
    c3 =  3/h;

    Keff1 = K + b0*M + b1*C;
    Keff2 = K + c3^2*M + c3*C;
    dK1   = decomposition(Keff1);
    dK2   = decomposition(Keff2);

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);
    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    % --- Time loop ---
    for k = 1:nSteps
        % ----- First sub-step: trapezoidal rule, t -> t + h/2 -----
        R1 = M*(b0*u + b2*v + a) + C*(b1*u + v);
        uH = dK1 \ R1;
        vH = b1*(uH - u) - v;
        aH = b1*(vH - v) - a;

        % ----- Second sub-step: Euler 3-point backward, t + h/2 -> t + h --
        R2 = -M*(c1*v + c2*vH + c3*c1*u + c3*c2*uH) - C*(c1*u + c2*uH);
        uNew = dK2 \ R2;
        vNew = c1*u + c2*uH + c3*uNew;
        aNew = c1*v + c2*vH + c3*vNew;

        u = uNew;  v = vNew;  a = aNew;
        U(:,k+1) = u;  V(:,k+1) = v;  A(:,k+1) = a;
    end

    t = (0:nSteps)*h;
end
