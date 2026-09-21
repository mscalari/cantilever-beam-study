function [PE, AD, rho, lambda] = viga_pe_ad(method, wn, xi, h)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Period elongation and amplitude     %
%  decay, Section 3.3 of the paper.    %
%--------------------------------------%
%
% Applies the SAME definitions used in Sections 4.2 and 4.3, now to the
% decoupled modal equations of the finite element beam.
%
% INPUT
%   method : 'DC','Houbolt','Newmark','Wilson','Bathe','RK1'..'RK4'
%   wn     : modal natural frequency        [rad/s]
%   xi     : PHYSICAL modal damping ratio   [-]
%   h      : time step                      [s]
%
% OUTPUT
%   PE     : period elongation   (NaN when not defined, see below)
%   AD     : amplitude decay     (NaN when not defined, see below)
%   rho    : spectral radius of Ahat        (always returned)
%   lambda : dominant eigenvalue of Ahat    (always returned)
%
% UNDAMPED CASE (xi = 0), Eqs. (19)-(21).  PE and AD are returned in
% PERCENT, as in Figures 5 to 8:
%
%       xibar      = -ln(rho)/Omega_n
%       Omegabar_d = Omega_n*sqrt(1 - xibar^2)
%       PE = (Omega_n/Omegabar_d - 1)*100
%       AD = (1 - exp(-2*pi*xibar*Omega_n/Omegabar_d))*100
%
% DAMPED CASE (0 < xi < 1), Eqs. (23)-(26).  PE and AD are returned as
% RATIOS: PE = 1 and AD = 1 mean the integrator reproduces the physical
% period and the physical decay exactly.
%
%       PE = Omega_d/angle(lambda)
%       AD = exp(-2*pi*xibar*Omega_n/angle(lambda)) / exp(-xi*wn*T_d)
%
% NOTE - the paper gives two apparently different damped formulations,
% one for the Runge-Kutta family (Eqs. 23-24, built on the closed-form
% solution of the difference equation) and one for the direct methods
% (Eqs. 25-26, built on the dominant eigenvalue).  They are algebraically
% identical:  Tbar_d/h = 2*pi/angle(lambda), so
%
%       rho^(Tbar_d/h) = exp(ln(rho)*2*pi/angle(lambda))
%                      = exp(-2*pi*xibar*Omega_n/angle(lambda))
%
% since xibar = -ln(rho)/Omega_n.  A single branch therefore serves both
% families, which is worth stating explicitly in the paper.
%
% ---------------------------------------------------------------------
% WHEN PE AND AD ARE NOT DEFINED
%
% Both quantities compare a numerical PERIOD and a numerical per-period
% DECAY against their physical counterparts.  Three situations leave
% them undefined, and each returns NaN rather than a misleading number:
%
%   (i)  PHYSICALLY OVERDAMPED MODE, xi >= 1.
%        The reference motion does not oscillate at all: the damped
%        frequency wd = wn*sqrt(1-xi^2) is imaginary and the physical
%        period Td = 2*pi/wd does not exist.  There is nothing for the
%        numerical period to be compared with.
%        This case is NOT hypothetical for the beam: Rayleigh damping
%        anchored at xi = 0.2 on modes 1 and 3 drives xi past unity from
%        mode 7 onwards, reaching xi = 18.4 at the highest mesh mode.
%        Without this guard, wd, Td and Omega_d all turn complex and PE
%        and AD are returned as complex numbers.
%
%   (ii) NON-OSCILLATORY NUMERICAL RESPONSE, angle(lambda) ~ 0.
%        The dominant eigenvalue is real and positive, so the scheme
%        produces no rotation and no numerical period is defined.
%
%   (iii) NUMERICAL DAMPING RATIO OUTSIDE THE UNDERDAMPED INTERVAL,
%        |xibar| >= 1.  The scheme is so far outside its stable range
%        that Omegabar_d = Omega_n*sqrt(1-xibar^2) turns complex.
%
% In all three cases the spectral radius rho is still returned and is
% the quantity to report.

    % --- Approximation operator and its dominant eigenvalue -----------
    Ahat = viga_operador(method, wn, xi, h);

    e = eig(Ahat);
    [rho, idx] = max(abs(e));
    lambda = e(idx);

    Omega_n = wn*h;
    th      = angle(lambda);

    % Numerical damping ratio introduced by the scheme, Eq. (19)
    xibar = -log(rho)/Omega_n;

    PE = NaN;  AD = NaN;

    % --- (i) physically overdamped or critically damped mode ----------
    % Checked FIRST, before anything that uses sqrt(1-xi^2).
    if xi >= 1
        return
    end

    % --- (ii) numerical response does not oscillate -------------------
    if th <= eps
        return
    end

    % --- (iii) numerical damping ratio outside the underdamped range --
    if abs(xibar) >= 1
        return
    end

    % --- Definitions --------------------------------------------------
    if xi == 0
        % Undamped, Eqs. (20)-(21), in percent
        Omegabar_d = Omega_n*sqrt(1 - xibar^2);
        PE = (Omega_n/Omegabar_d - 1)*100;
        AD = (1 - exp(-2*pi*xibar*Omega_n/Omegabar_d))*100;
    else
        % Damped, Eqs. (23)-(26), as ratios.  Reaching this point
        % guarantees 0 < xi < 1, so wd and Td are real.
        wd      = wn*sqrt(1 - xi^2);
        Td      = 2*pi/wd;
        Omega_d = wd*h;
        PE = Omega_d/th;
        AD = exp(-2*pi*xibar*Omega_n/th) / exp(-xi*wn*Td);
    end
end