function r = viga_rho_beam(method, w, xi, h)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Spectral radius of the MDOF beam    %
%  operator, exact for proportional    %
%  damping (modal maximum).            %
%--------------------------------------%
%
% Extracted from the local function of the same purpose in
% viga_main_simulacao.m so it can be reused by the new Estudo B-F
% scripts and by viga_espectro.m, without touching that script.
%
% Because the damping is Rayleigh (proportional), M, C and K are
% simultaneously diagonalised by the mode shapes and every integrator's
% MDOF operator is block-diagonal in the modal basis (Eq. (23) of the
% paper, Section 4.4.1). The spectral radius of the full operator is
% therefore exactly the maximum of the SCALAR spectral radii over the
% modes - verified numerically in viga_main_validacao.m.
%
% INPUT
%   method : 'DC','Houbolt','Newmark','Wilson','Bathe','RK1'..'RK4'
%   w, xi  : [m x 1] modal frequencies and damping ratios of the case
%            (the full mesh, or a reduced model's retained modes)
%   h      : time step                                          [s]
%
% OUTPUT
%   r      : rho(Ahat_MDOF(h)) = max_k rho(Ahat_scalar(w_k, xi_k, h))

    r = 0;
    for k = 1:numel(w)
        rk = max(abs(eig(viga_operador(method, w(k), xi(k), h))));
        if rk > r, r = rk; end
    end
end
