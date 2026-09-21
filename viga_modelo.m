function mdl = viga_modelo(p, N)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  One-call assembly of the complete   %
%  cantilever beam model.              %
%--------------------------------------%
%
% Convenience wrapper that builds, in the correct order, everything the
% time-domain simulations need.
%
% INPUT
%   p : struct from viga_params() (optional)
%   N : number of elements (optional, defaults to p.N)
%
% OUTPUT
%   mdl : struct with fields
%           .p      parameters
%           .info   mesh bookkeeping
%           .K,.M   FE matrices
%           .C      Rayleigh damping matrix (xi = p.xi_target case)
%           .Cz     zero damping matrix     (xi = 0 case)
%           .w      natural frequencies of the discrete system
%           .Phi    mass-normalised mode shapes
%           .damp   Rayleigh data (a0, a1, modal xi)
%           .u0,.v0,.a0v  impulse initial conditions
%           .T1     fundamental period
%           .wmax   highest frequency of the mesh

    if nargin < 1 || isempty(p), p = viga_params(); end
    if nargin < 2 || isempty(N), N = p.N;           end

    [K, M, info]  = viga_matrizes_EF(p, N);
    [w, Phi]      = viga_autovalores_EF(K, M);
    [C, damp]     = viga_amortecimento(K, M, w, p);
    [u0, v0, a0v] = viga_ci_impulso(K, M, C, info, p);

    mdl.p    = p;
    mdl.info = info;
    mdl.K    = K;
    mdl.M    = M;
    mdl.C    = C;
    mdl.Cz   = zeros(size(K));
    mdl.w    = w;
    mdl.Phi  = Phi;
    mdl.damp = damp;
    mdl.u0   = u0;
    mdl.v0   = v0;
    mdl.a0v  = a0v;
    mdl.T1   = 2*pi/w(1);
    mdl.wmax = w(end);
end
