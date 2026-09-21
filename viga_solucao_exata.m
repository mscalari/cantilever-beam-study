function [U, V, A] = viga_solucao_exata(t, M, w, Phi, xi, u0, v0)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Exact (reference) response of the   %
%  semi-discrete FE beam system.       %
%--------------------------------------%
%
% Because the damping is proportional, the FE system decouples in the
% modal basis and each modal coordinate obeys exactly the SDOF equation
% analysed in Section 2 of the paper:
%
%       eta_k'' + 2*xi_k*w_k*eta_k' + w_k^2*eta_k = 0
%
% The response below is therefore EXACT for the semi-discrete system:
% it carries no time integration error whatsoever and serves as the
% reference against which every integrator is compared.  Spatial
% discretization error is assessed separately in
% viga_main_convergencia.m.
%
% INPUT
%   t      : [1 x nt] time vector                                 [s]
%   M      : mass matrix (used to project the initial conditions)
%   w      : [n x 1] natural frequencies                      [rad/s]
%   Phi    : [n x n] mass-normalised mode shapes (Phi'*M*Phi = I)
%   xi     : [n x 1] modal damping ratios
%   u0, v0 : [n x 1] initial displacement and velocity
%
% OUTPUT
%   U, V, A : [n x nt] displacement, velocity and acceleration histories
%
% All three damping regimes are handled: underdamped (xi < 1),
% critically damped (xi = 1) and overdamped (xi > 1).  The overdamped
% branch is written with decaying exponentials instead of sinh to avoid
% floating point overflow at the strongly damped high-frequency modes.

    t  = t(:).';
    n  = numel(w);
    nt = numel(t);

    U = zeros(n, nt);
    V = zeros(n, nt);
    A = zeros(n, nt);

    % --- Modal projection of the initial conditions ---
    % With Phi'*M*Phi = I the generalised mass is unity.
    q0  = Phi.' * (M * u0);
    dq0 = Phi.' * (M * v0);

    tolc = 1e-10;

    for k = 1:n
        wk = w(k);
        zk = xi(k);
        e0 = q0(k);
        d0 = dq0(k);

        if wk < eps
            % Rigid-body mode (not present in a cantilever, kept for safety)
            eta  = e0 + d0*t;
            deta = d0*ones(1,nt);
            dde  = zeros(1,nt);

        elseif zk < 1 - tolc
            % ---- Underdamped ----
            wd  = wk*sqrt(1 - zk^2);
            s   = zk*wk;
            ex  = exp(-s*t);
            cwd = cos(wd*t);
            swd = sin(wd*t);
            Ac  = e0;
            Bc  = (d0 + s*e0)/wd;
            eta  = ex.*(Ac*cwd + Bc*swd);
            deta = ex.*((-s*Ac + wd*Bc)*cwd + (-s*Bc - wd*Ac)*swd);
            dde  = -2*zk*wk*deta - wk^2*eta;

        elseif zk < 1 + tolc
            % ---- Critically damped ----
            ex   = exp(-wk*t);
            eta  = ex.*(e0 + (d0 + wk*e0)*t);
            deta = ex.*((d0 + wk*e0) - wk*(e0 + (d0 + wk*e0)*t));
            dde  = -2*wk*deta - wk^2*eta;

        else
            % ---- Overdamped (exponential form, overflow safe) ----
            sg = sqrt(zk^2 - 1);
            r1 = -(zk - sg)*wk;     % slow decaying root  (r1 < 0)
            r2 = -(zk + sg)*wk;     % fast decaying root  (r2 < r1 < 0)
            c1 = (d0 - r2*e0)/(r1 - r2);
            c2 = e0 - c1;
            E1 = exp(r1*t);
            E2 = exp(r2*t);
            eta  = c1*E1 + c2*E2;
            deta = c1*r1*E1 + c2*r2*E2;
            dde  = c1*r1^2*E1 + c2*r2^2*E2;
        end

        U = U + Phi(:,k)*eta;
        V = V + Phi(:,k)*deta;
        A = A + Phi(:,k)*dde;
    end
end
