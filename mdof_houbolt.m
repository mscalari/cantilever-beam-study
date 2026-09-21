function [t, U, V, A] = mdof_houbolt(K, M, C, u0, v0, h, nSteps)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Houbolt method - matrix form for    %
%  MDOF free vibration.                %
%--------------------------------------%
%
% Matrix counterpart of houbolt_last.m.  The Houbolt scheme is a
% three-step method and is therefore NOT self-starting: the
% displacements at t = h and t = 2h must be supplied by another
% integrator.  Here they are obtained with the Newmark trapezoidal rule
% using the SAME step size, which is the standard starting procedure.
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

    n = size(K,1);

    % --- Integration constants ---
    a0 = 2/h^2;
    a1 = 11/(6*h);
    a2 = 5/h^2;
    a3 = 3/h;
    a4 = -2*a0;
    a5 = -a3/2;
    a6 = a0/2;
    a7 = a3/9;

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    % --- Starting procedure: two Newmark steps ---
    nStart = min(2, nSteps);
    [~, Us, Vs, As] = mdof_newmark(K, M, C, u0, v0, h, nStart);
    U(:,1:nStart+1) = Us;
    V(:,1:nStart+1) = Vs;
    A(:,1:nStart+1) = As;

    if nSteps < 3
        t = (0:nSteps)*h;
        return
    end

    % --- Effective stiffness, factorised once ---
    Keff = K + a0*M + a1*C;
    dKe  = decomposition(Keff);

    % --- Time loop (index k corresponds to time t = (k-1)*h) ---
    for k = 3:nSteps
        uc = U(:,k);      % t
        u1 = U(:,k-1);    % t - h
        u2 = U(:,k-2);    % t - 2h

        Reff = M*(a2*uc + a4*u1 + a6*u2) + C*(a3*uc + a5*u1 + a7*u2);
        uNew = dKe \ Reff;

        U(:,k+1) = uNew;
        A(:,k+1) = a0*uNew - a2*uc - a4*u1 - a6*u2;
        V(:,k+1) = a1*uNew - a3*uc - a5*u1 - a7*u2;
    end

    t = (0:nSteps)*h;
end
