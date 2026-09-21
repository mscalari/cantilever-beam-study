function [u0, v0, a0] = viga_ci_impulso(K, M, C, info, p)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Impulse initial conditions applied  %
%  at the free end of the beam.        %
%--------------------------------------%
%
% The impulse is imposed as an initial TRANSVERSE VELOCITY at the tip
% translation DOF, so that the problem reduces to free vibration
%
%       M*udd + C*ud + K*u = 0 ,   u(0) = 0 ,  ud(0) = v0
%
% which is the matrix counterpart of the SDOF free-vibration equation
% analysed in Sections 2 to 4 of the paper.
%
% INPUT
%   K, M, C : FE system matrices
%   info    : struct from viga_matrizes_EF() (supplies the tip DOF index)
%   p       : struct from viga_params()
%
% OUTPUT
%   u0 : [n x 1] initial displacement (zero)
%   v0 : [n x 1] initial velocity     (V0 at the tip translation DOF)
%   a0 : [n x 1] initial acceleration, from equilibrium at t = 0:
%                  a0 = M \ (-C*v0 - K*u0)
%                Note a0 is NOT zero when damping is present.
%
% Because the system is linear, the magnitude V0 only scales the
% response and has no effect on any stability, PE or AD conclusion.

    if nargin < 5, p = viga_params(); end

    n  = size(K,1);
    u0 = zeros(n,1);
    v0 = zeros(n,1);
    v0(info.tipW) = p.V0;

    a0 = M \ (-C*v0 - K*u0);
end
