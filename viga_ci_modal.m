function [u0, v0, a0] = viga_ci_modal(K, M, C, Phi, modeIdx, u0scale)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Modal (displacement) initial        %
%  conditions - advisor's revision.    %
%--------------------------------------%
%
% Instead of an impulsive velocity at the tip (viga_ci_impulso.m, which
% excites many modes at once), this displaces EVERY node according to a
% single mode shape and releases the beam from rest:
%
%       u(0) = u0scale * phi_j ,   ud(0) = 0
%
% Because Phi is mass-normalised (Phi'*M*Phi = I) and the damping is
% proportional, this initial condition projects onto ONLY mode j of the
% modal decomposition (viga_solucao_exata.m): q0(j) = u0scale, every
% other q0(k) = 0. The response is therefore, in principle, a pure
% single-frequency oscillation at w(j) - exactly what lets Estudos C
% and D isolate how each integrator behaves at a chosen natural
% frequency, instead of the mixture of modes excited by the impulse.
%
% INPUT
%   K, M, C   : FE system matrices (or a reduced model's Kr, Mr, Cr)
%   Phi       : [n x n] mass-normalised mode shapes (viga_autovalores_EF)
%               - or the reduced Phir from viga_reducao_modal.m
%   modeIdx   : index of the mode to excite (column of Phi)
%   u0scale   : (optional, default 1) modal amplitude q0(modeIdx).
%               The system is linear, so this only scales the response
%               and has no effect on any stability/PE/AD conclusion -
%               same remark as in viga_ci_impulso.m.
%
% OUTPUT
%   u0 : [n x 1] initial displacement = u0scale*Phi(:,modeIdx)
%   v0 : [n x 1] initial velocity     = 0
%   a0 : [n x 1] initial acceleration, from equilibrium at t = 0:
%                  a0 = M \ (-C*v0 - K*u0)
%                NOT zero when damping is present, exactly as in
%                viga_ci_impulso.m.

    if nargin < 6 || isempty(u0scale), u0scale = 1; end

    n = size(K,1);
    if modeIdx < 1 || modeIdx > size(Phi,2)
        error('viga_ci_modal:modeIdx', ...
            'modeIdx = %d is out of range for a model with %d retained mode(s).', ...
            modeIdx, size(Phi,2));
    end

    u0 = u0scale * Phi(:,modeIdx);
    v0 = zeros(n,1);
    a0 = M \ (-C*v0 - K*u0);
end
