function [hcr, rhoFun] = viga_h_critico(K, M, C, method, hHi, tolRho)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Critical time step of a conditionally%
%  stable scheme applied to an MDOF     %
%  finite element system.               %
%--------------------------------------%
%
% The stability of an integrator applied to the FE system is governed by
% the spectral radius of its MDOF approximation operator, exactly as in
% Section 3.2 of the paper but now for a block operator of size 2n.  The
% critical step is found by bisection on
%
%       rho(Ahat(h)) <= 1
%
% INPUT
%   K, M, C : [n x n] system matrices
%   method  : 'RK1','RK2','RK3','RK4' or 'DC'
%   hHi     : upper bracket for the bisection (optional, default 1e-1)
%   tolRho  : tolerance on the unit spectral radius (optional, 1e-8).
%             A strictly zero tolerance is unusable here because the
%             operator becomes extremely ill conditioned as h -> 0.
%
% OUTPUT
%   hcr     : critical time step [s].  Returns 0 when no stable step is
%             found within 200 bisection halvings, i.e. the method is
%             unstable for every practical h.
%   rhoFun  : function handle h -> rho(Ahat(h)), useful for plotting
%
% IMPORTANT - for a finite element mesh the critical step is dictated by
% the HIGHEST frequency of the discrete system, not by the fundamental
% frequency.  This is why the limits found here are orders of magnitude
% smaller than the h/Tn limits reported for the SDOF oscillator in
% Sections 4.1 to 4.3.
%
% Reference checks (undamped, consistent mass):
%   DC  : hcr = 2/w_max
%   RK3 : hcr = sqrt(3)/w_max
%   RK4 : hcr = sqrt(8)/w_max

    if nargin < 5 || isempty(hHi),    hHi    = 1e-1; end
    if nargin < 6 || isempty(tolRho), tolRho = 1e-8; end

    n = size(K,1);

    switch upper(method)
        case {'RK1','RK2','RK3','RK4'}
            order = str2double(method(3));
            Ast = [zeros(n), eye(n); -(M\K), -(M\C)];
            rhoFun = @(h) rhoRK(Ast, h, order);

        case 'DC'
            rhoFun = @(h) rhoCD(K, M, C, h);

        otherwise
            error('viga_h_critico:method', ...
                'Unsupported method "%s". Use RK1..RK4 or DC.', method);
    end

    % --- Already stable at the upper bracket? ---
    if rhoFun(hHi) <= 1 + tolRho
        hcr = hHi;
        warning('viga_h_critico:bracket', ...
            'Method %s still stable at h = %g; enlarge hHi.', method, hHi);
        return
    end

    % --- Search downwards for a stable lower bracket ---
    hLo   = hHi;
    found = false;
    for it = 1:200
        hLo = hLo/2;
        if rhoFun(hLo) <= 1 + tolRho
            found = true;
            break
        end
    end
    if ~found
        hcr = 0;                     % unstable for every practical step
        return
    end

    % --- Bisection ---
    for it = 1:200
        hMid = 0.5*(hLo + hHi);
        if rhoFun(hMid) <= 1 + tolRho
            hLo = hMid;
        else
            hHi = hMid;
        end
    end
    hcr = hLo;
end

% ----------------------------------------------------------------------
function r = rhoRK(Ast, h, order)
% Approximation operator of an explicit RK method: truncated Taylor
% expansion of exp(h*Ast) up to the method order.
    m    = size(Ast,1);
    Ahat = eye(m);
    term = eye(m);
    for i = 1:order
        term = (term*(h*Ast))/i;
        Ahat = Ahat + term;
    end
    r = max(abs(eig(Ahat)));
end

% ----------------------------------------------------------------------
function r = rhoCD(K, M, C, h)
% The two-step Central Differences recurrence
%   (a0*M + a1*C) u_{n+1} = (a2*M - K) u_n + (a1*C - a0*M) u_{n-1}
% written as a first-order block operator acting on {u_{n+1}; u_n}.
    n  = size(K,1);
    a0 = 1/h^2;  a1 = 1/(2*h);  a2 = 2*a0;
    Meff = a0*M + a1*C;
    B1 = Meff \ (a2*M - K);
    B2 = Meff \ (a1*C - a0*M);
    Ahat = [B1, B2; eye(n), zeros(n)];
    r = max(abs(eig(Ahat)));
end
