function [C, damp] = viga_amortecimento(K, M, w, p, a0in, a1in)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Rayleigh (proportional) damping     %
%  matrix of the FE beam system.       %
%--------------------------------------%
%
% TWO MODES OF USE
%
%   1) Anchored (original behaviour, 4 input arguments): builds
%      C = a0*M + a1*K such that the modal damping ratio equals
%      p.xi_target exactly at the two anchor modes p.mode_i and p.mode_j.
%
%   2) Direct coefficients (Estudo E, advisor's revision): pass a0in and
%      a1in explicitly - e.g. a0=0 with a1=1e-4 or 1e-5, as suggested by
%      the advisor to see the RESULTING xi per mode instead of imposing
%      one. p.xi_target/.mode_i/.mode_j are then ignored and damp.anchors
%      is returned empty.
%
% In both cases,
%
%   xi_k = a0/(2*w_k) + a1*w_k/2
%
% INPUT
%   K, M       : FE matrices
%   w          : [n x 1] natural frequencies from viga_autovalores_EF()
%   p          : struct from viga_params()
%   a0in, a1in : (optional) explicit Rayleigh coefficients. When given,
%                the anchor-based solve is skipped entirely.
%
% OUTPUT
%   C    : [n x n] damping matrix
%   damp : struct with
%            .a0, .a1    Rayleigh coefficients (given or solved for)
%            .xi         [n x 1] resulting modal damping ratios (EVERY
%                        mode of the mesh, not just the anchors)
%            .anchors    indices of the anchor modes ([] in mode 2)
%            .xi_target  target xi at the anchors ([] in mode 2)
%
% NOTE - the stiffness-proportional term makes xi grow linearly with
% frequency, so the highest FE modes become strongly overdamped
% (xi >> 1).  This is the physically correct behaviour of Rayleigh
% damping and is precisely the regime in which physical damping is
% reported to interfere with the stability of time integration schemes.

    if nargin < 4 || isempty(p), p = viga_params(); end

    if nargin >= 6 && ~isempty(a0in) && ~isempty(a1in)
        % ---- Mode 2: direct coefficients ----
        a0 = a0in;
        a1 = a1in;
        damp.anchors    = [];
        damp.xi_target  = [];
    else
        % ---- Mode 1: anchored at two modes (original behaviour) ----
        i = p.mode_i;  j = p.mode_j;
        wi = w(i);     wj = w(j);
        xi = p.xi_target;

        a1 = 2*xi/(wi + wj);
        a0 = a1*wi*wj;

        damp.anchors    = [i j];
        damp.xi_target  = xi;
    end

    C = a0*M + a1*K;
    C = 0.5*(C + C.');

    damp.a0 = a0;
    damp.a1 = a1;
    damp.xi = a0./(2*w) + a1*w/2;      % every mode of the mesh
end
