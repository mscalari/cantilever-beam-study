function Ahat = viga_operador(method, wn, xi, h)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Approximation operator Ahat of each  %
%  scheme for a single modal equation.  %
%--------------------------------------%
%
% Returns the operator of Eq. (4) of the paper for the decoupled modal
% equation
%
%       xdd + 2*xi*wn*xd + wn^2*x = 0
%
% Expressions taken from Section 2.2 and Appendix A of the paper.
%
% INPUT
%   method : 'DC','Houbolt','Newmark','Wilson','Bathe','RK1'..'RK4'
%   wn     : modal natural frequency        [rad/s]
%   xi     : modal damping ratio            [-]
%   h      : time step                      [s]
%
% OUTPUT
%   Ahat   : approximation operator (2x2 or 3x3)
%
% WHY THIS IS ENOUGH FOR THE BEAM
% Rayleigh damping is proportional, so M, C and K are simultaneously
% diagonalised by the mode shapes.  Every integrator considered here is
% built from M, C and K through rational expressions, hence its MDOF
% operator is block-diagonal in the modal basis.  The spectral radius of
% the MDOF operator is therefore the maximum of the scalar spectral radii
% over the modes of the mesh - verified numerically in
% viga_main_validacao.m.

    switch upper(method)

        case 'DC'       % ---- Central Differences ----
            Ahat = [ (2 - wn^2*h^2)/(1 + xi*wn*h), -(1 - xi*wn*h)/(1 + xi*wn*h);
                      1,                            0 ];

        case 'HOUBOLT'  % ---- Houbolt, Eqs. (A.1)-(A.2) ----
            bt = (2/(wn^2*h^2) + 11*xi/(3*wn*h) + 1)^(-1);
            ki = (xi*bt)/(wn*h);
            Ahat = [ 5*bt/(wn^2*h^2) + 6*ki, -(4*bt/(wn^2*h^2) + 3*ki), bt/(wn^2*h^2) + 2*ki/3;
                     1, 0, 0;
                     0, 1, 0 ];

        case 'NEWMARK'  % ---- Newmark, alpha = 1/4, delta = 1/2 ----
            a  = 1 + xi*wn*h + (wn^2*h^2)/4;
            c1 = -wn^2/a;
            c2 = -(2*xi*wn + wn^2*h)/a;
            c3 = -(xi*wn*h + (wn^2*h^2)/4)/a;
            Ahat = [ 1 + (h^2/4)*c1, h + (h^2/4)*c2, (h^2/4)*(1 + c3);
                         (h/2)*c1,   1 + (h/2)*c2,   (h/2)*(1 + c3);
                              c1,          c2,             c3        ];

        case 'WILSON'   % ---- Wilson-theta, theta = 1.4, Eqs. (A.3)-(A.4) ----
            th = 1.4;
            b3 = (th/(wn^2*h^2) + xi*th^2/(wn*h) + th^3/6)^(-1);
            k3 = (xi*b3)/(wn*h);
            Ahat = [ 1 - b3*th^2/3 - 1/th - k3*th,        (-b3*th - 2*k3)/h,      -b3/h^2;
                     h*(1 - 1/(2*th) - b3*th^2/6 - k3*th/2), 1 - b3*th/2 - k3,   -b3/(2*h);
                     h^2*(1/2 - 1/(6*th) - b3*th^2/18 - k3*th/6), h*(1 - b3*th/6 - k3/3), 1 - b3/6 ];

        case 'BATHE'    % ---- Bathe, Eqs. (A.5)-(A.6) ----
            b1 = 16 + 8*xi*wn*h + wn^2*h^2;
            b2 =  9 + 6*xi*wn*h + wn^2*h^2;
            Ahat = (1/(b1*b2)) * [ ...
                -4*wn*h*(24*xi + 7*wn*h), ...
                 wn*(-288*xi + 14*xi*wn^2*h^2 - 144*wn*h + 5*wn^3*h^3 + 48*xi^2*wn*h), ...
                 wn^2*(24*xi*wn*h + 19*wn^2*h^2 - 144); ...
                -4*h*(-12 + wn^2*h^2), ...
                 144 - 47*wn^2*h^2 - 8*xi*wn^3*h^3 - 24*xi*wn*h, ...
                 wn^2*h*(-96 - 24*xi*wn*h + wn^2*h^2); ...
                 4*h^2*(7 + 2*xi*wn*h), ...
                 h*(144 - 5*wn^2*h^2 + 80*xi*wn*h + 16*xi^2*wn^2*h^2), ...
                -19*wn^2*h^2 + 144 + 168*xi*wn*h + 48*xi^2*wn^2*h^2 - 2*xi*wn^3*h^3 ];

        case {'RK1','RK2','RK3','RK4'}   % ---- Explicit RK, Eq. (14) ----
            j = str2double(method(3));
            A = [0 1; -wn^2 -2*xi*wn];
            Ahat = eye(2);
            term = eye(2);
            for i = 1:j
                term = (term*(h*A))/i;
                Ahat = Ahat + term;
            end

        otherwise
            error('viga_operador:method','Unknown method "%s".', method);
    end
end
