function [K, M, info] = viga_matrizes_EF(p, N)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Finite element matrices of the      %
%  Euler-Bernoulli cantilever beam.    %
%--------------------------------------%
%
% Two-node beam element, four degrees of freedom (w1, theta1, w2, theta2),
% cubic Hermite interpolation, CONSISTENT mass matrix.
%
% INPUT
%   p : struct from viga_params()
%   N : number of elements (optional, defaults to p.N)
%
% OUTPUT
%   K    : [n x n] global stiffness matrix, clamped DOFs removed
%   M    : [n x n] global consistent mass matrix, clamped DOFs removed
%   info : struct with mesh bookkeeping
%            .N        number of elements
%            .Le       element length
%            .ndof     number of free DOFs (n = 2*N)
%            .tipW     index of the tip transverse DOF
%            .tipTh    index of the tip rotation DOF
%            .xNodes   nodal coordinates
%
% The clamped boundary condition w(0) = w'(0) = 0 is imposed by deleting
% the two first rows/columns of the assembled system.

    if nargin < 2 || isempty(N), N = p.N; end

    Le   = p.L/N;
    ndof = 2*(N+1);

    K = zeros(ndof);
    M = zeros(ndof);

    % --- Element matrices (constant along the uniform mesh) ---
    Ke = (p.E*p.I/Le^3) * [ 12,      6*Le,   -12,      6*Le;
                             6*Le,   4*Le^2, -6*Le,    2*Le^2;
                           -12,     -6*Le,    12,     -6*Le;
                             6*Le,   2*Le^2, -6*Le,    4*Le^2 ];

    Me = (p.rho*p.A*Le/420) * [ 156,     22*Le,   54,    -13*Le;
                                 22*Le,   4*Le^2, 13*Le,  -3*Le^2;
                                 54,     13*Le,  156,    -22*Le;
                                -13*Le,  -3*Le^2,-22*Le,   4*Le^2 ];

    % --- Assembly ---
    for e = 1:N
        d = [2*e-1, 2*e, 2*e+1, 2*e+2];
        K(d,d) = K(d,d) + Ke;
        M(d,d) = M(d,d) + Me;
    end

    % --- Clamp at x = 0 : remove w1 and theta1 ---
    K = K(3:end, 3:end);
    M = M(3:end, 3:end);

    % --- Enforce exact symmetry (kills round-off asymmetry) ---
    K = 0.5*(K + K.');
    M = 0.5*(M + M.');

    info.N      = N;
    info.Le     = Le;
    info.ndof   = size(K,1);          % = 2*N
    info.tipW   = size(K,1) - 1;      % transverse DOF of the free end
    info.tipTh  = size(K,1);          % rotation   DOF of the free end
    info.xNodes = (0:N).'*Le;
end
