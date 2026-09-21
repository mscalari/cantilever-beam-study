function red = viga_reducao_modal(mdl, m)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Modal reduction (truncation) of the %
%  cantilever beam model - Estudo F.   %
%--------------------------------------%
%
% Retains only the first m modes of the full model and expresses the
% reduced equations of motion in MODAL coordinates q:
%
%       I*qdd + diag(2*xi_k*w_k)*qd + diag(w_k^2)*q = 0 ,   k = 1..m
%
% This is EXACT, not an approximation on top of the FE discretisation:
% because Phi is mass-normalised (Phi'*M*Phi = I) and the damping is
% Rayleigh (proportional, hence diagonal in the modal basis), dropping
% modes m+1..n from the modal sum is the definition of modal truncation.
% The only approximation is physical - assuming modes m+1..n do not
% matter for the response of interest - and that assumption is exactly
% what Estudo F is meant to check (compare against the full model).
%
% Because the reduced system is itself a linear MDOF system (Kr, Mr, Cr
% below), every routine already written for the full beam - mdof_*.m,
% viga_h_critico.m, viga_pe_ad.m via viga_operador.m, viga_solucao_exata.m
% - applies to it UNCHANGED, simply by passing the reduced matrices in
% place of the full ones.
%
% INPUT
%   mdl : struct from viga_modelo() (the FULL model)
%   m   : number of modes retained (1 <= m <= n)
%
% OUTPUT
%   red : struct with
%           .m       number of modes retained
%           .K,.M,.C [m x m] reduced matrices, all diagonal:
%                      M = eye(m), K = diag(w(1:m).^2),
%                      C = diag(2*xi(1:m).*w(1:m))
%           .Phi     [n x m] retained mode shapes (columns of mdl.Phi)
%           .w,.xi   [m x 1] frequencies/damping ratios retained
%           .wmax    w(m) - the reduced model's highest frequency,
%                    used in place of mdl.wmax for stability/h_cr
%           .T1      fundamental period, 2*pi/w(1) (same as the full
%                    model's, since mode 1 is always retained)
%           .u0,.v0  impulse initial condition PROJECTED onto the
%                    retained modes: q0 = Phi'*mdl.M*mdl.u0 (analogous
%                    for v0). Feed these to mdof_*.m together with
%                    red.K, red.M, red.C exactly as u0/v0 are used with
%                    K/M/C in the full model.
%           .tipRow  1 x m row such that u_tip(t) = tipRow * q(t)
%                    recovers the physical tip displacement from the
%                    modal-coordinate response of any integrator.
%           .Mref    eye(m) - pass as "M" to viga_solucao_exata together
%                    with Phi = eye(m), since q is already the physical
%                    state of the reduced system.
%
% Example (reduced-model impulse response, Bathe method):
%
%   red = viga_reducao_modal(mdl, 4);
%   [t,Q] = mdof_bathe(red.K, red.M, red.C, red.u0, red.v0, h, nSteps);
%   u_tip = red.tipRow * Q;                      % physical tip response
%   hcr4  = viga_h_critico(red.K, red.M, red.C, 'RK4');

    n = size(mdl.K,1);
    if m < 1 || m > n
        error('viga_reducao_modal:m', ...
            'm = %d is out of range for a model with %d DOF.', m, n);
    end

    w  = mdl.w(1:m);
    xi = mdl.damp.xi(1:m);
    Phi_r = mdl.Phi(:,1:m);

    red.m    = m;
    red.M    = eye(m);
    red.K    = diag(w.^2);
    red.C    = diag(2*xi(:).*w(:));
    red.Phi  = Phi_r;
    red.w    = w;
    red.xi   = xi;
    red.wmax = w(end);
    red.T1   = 2*pi/w(1);

    red.u0 = Phi_r.' * (mdl.M * mdl.u0);
    red.v0 = Phi_r.' * (mdl.M * mdl.v0);

    red.tipRow = Phi_r(mdl.info.tipW, :);
    red.Mref   = eye(m);
end
