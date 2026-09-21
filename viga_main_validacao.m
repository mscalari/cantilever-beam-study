%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  STEP 1 - Validation of the model    %
%  and of the matrix integrators.      %
%--------------------------------------%
%
% Three checks, all of which must pass before any result is trusted:
%
%   (a) SPATIAL - the FE frequencies converge to the analytical
%       Euler-Bernoulli frequencies, so the spatial discretization error
%       is negligible and everything seen later is time integration
%       error only.
%
%   (b) TEMPORAL - the matrix routines mdof_*.m reproduce the scalar
%       results of Sections 4.1 to 4.4 in the single-degree-of-freedom
%       limit, each to its own order of accuracy.
%
%   (c) MODAL - the spectral radius of the MDOF operator equals the
%       maximum of the scalar spectral radii over the modes of the mesh.
%       This is what allows viga_operador.m to be used for the stability,
%       PE and AD analysis of the beam.
%
% Run this script first.

clc; clear; close all;

p = viga_params();

fprintf('=====================================================\n');
fprintf(' (a) SPATIAL DISCRETIZATION\n');
fprintf('=====================================================\n');
fprintf('  Steel, E = %.3e Pa, rho = %.0f kg/m^3\n', p.E, p.rho);
fprintf('  L = %.3f m, b = %.3f m, h = %.4f m, L/h = %.0f\n\n', ...
    p.L, p.b, p.h, p.slenderness);

nModes = 3;
wAna   = viga_modos_analiticos(p, nModes);

fprintf('  %4s %6s %12s %12s %12s %14s\n', ...
    'N','ndof','err w1','err w2','err w3','w_max');
for N = [2 5 10 20 40]
    [Kf, Mf, infoF] = viga_matrizes_EF(p, N);
    wf = viga_autovalores_EF(Kf, Mf);
    nUse = min(nModes, numel(wf));
    fprintf('  %4d %6d', N, infoF.ndof);
    for k = 1:nModes
        if k <= nUse
            fprintf('%12.2e', abs(wf(k)-wAna(k))/wAna(k));
        else
            fprintf('%12s','-');
        end
    end
    fprintf('%14.4e\n', wf(end));
end

mdl = viga_modelo(p);
fprintf('\n  Mesh adopted: N = %d (%d DOF)\n', mdl.info.N, mdl.info.ndof);
fprintf('    w_1   = %12.4f rad/s   T_1 = %.6f s\n', mdl.w(1), mdl.T1);
fprintf('    w_2   = %12.4f rad/s\n', mdl.w(2));
fprintf('    w_3   = %12.4f rad/s\n', mdl.w(3));
fprintf('    w_max = %12.4e rad/s   (w_max/w_1 = %.1f)\n', ...
    mdl.wmax, mdl.wmax/mdl.w(1));
fprintf('\n  Rayleigh damping  C = a0*M + a1*K\n');
fprintf('    anchors: modes %d and %d at xi = %.2f\n', ...
    mdl.damp.anchors(1), mdl.damp.anchors(2), mdl.damp.xi_target);
fprintf('    a0 = %.6f     a1 = %.6e\n', mdl.damp.a0, mdl.damp.a1);
fprintf('    xi_1 = %.4f   xi_2 = %.4f   xi_3 = %.4f   xi_max = %.2f\n', ...
    mdl.damp.xi(1), mdl.damp.xi(2), mdl.damp.xi(3), mdl.damp.xi(end));

%% ==================================================================
fprintf('\n=====================================================\n');
fprintf(' (b) MATRIX INTEGRATORS IN THE SDOF LIMIT\n');
fprintf('=====================================================\n');

wnS = 0.314;  TnS = 2*pi/wnS;
hS  = 0.005*TnS;  nsS = round(10*TnS/hS);

names = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe', ...
         'RK1','RK2','RK3','RK4'};

for xiS = [0 0.2]
    Ks = wnS^2;  Ms = 1;  Cs = 2*xiS*wnS;
    wdS = wnS*sqrt(1-xiS^2);
    tS  = (0:nsS)*hS;
    uEx = exp(-xiS*wnS*tS).*(cos(wdS*tS) + (xiS*wnS/wdS)*sin(wdS*tS));

    fprintf('\n  xi = %.1f   (h/Tn = 0.005, 10 periods)\n', xiS);
    for k = 1:numel(names)
        [~, U] = run_scheme(names{k}, Ks, Ms, Cs, 1, 0, hS, nsS);
        fprintf('    %-20s max abs err = %.4e\n', names{k}, ...
            max(abs(U(1,:) - uEx)));
    end
end
fprintf('\n  Expected: RK4 < RK3 < RK2 < RK1, direct methods around 1e-4.\n');
fprintf('  Any value of order unity indicates an implementation error.\n');

%% ==================================================================
fprintf('\n=====================================================\n');
fprintf(' (c) MDOF SPECTRAL RADIUS = MAX OVER MODES\n');
fprintf('=====================================================\n');
fprintf('  Cross-check of the full 2n block operator against the modal\n');
fprintf('  maximum built from viga_operador.m.\n\n');

fprintf('  %-6s %12s %18s %18s %10s\n','method','h [s]','full operator','modal maximum','diff');
for meth = {'DC','RK2','RK4'}
    [~, rhoFull] = viga_h_critico(mdl.K, mdl.M, mdl.C, meth{1});
    for hTest = [1e-6 1e-5 5e-5]
        rFull  = rhoFull(hTest);
        rModal = max(arrayfun(@(k) max(abs(eig( ...
            viga_operador(meth{1}, mdl.w(k), mdl.damp.xi(k), hTest)))), ...
            1:numel(mdl.w)));
        fprintf('  %-6s %12.1e %18.10f %18.10f %10.1e\n', ...
            meth{1}, hTest, rFull, rModal, abs(rFull-rModal));
    end
end
fprintf('\n  Differences at round-off level confirm the modal decomposition.\n');
fprintf('=====================================================\n');

%% ==================================================================
%  LOCAL FUNCTIONS
%  ==================================================================
function [t, U] = run_scheme(name, K, M, C, u0, v0, h, nSteps)
    switch name
        case 'Central Differences', [t,U] = mdof_dc     (K,M,C,u0,v0,h,nSteps);
        case 'Houbolt',             [t,U] = mdof_houbolt(K,M,C,u0,v0,h,nSteps);
        case 'Newmark',             [t,U] = mdof_newmark(K,M,C,u0,v0,h,nSteps);
        case 'Wilson-Theta',        [t,U] = mdof_wilson (K,M,C,u0,v0,h,nSteps);
        case 'Bathe',               [t,U] = mdof_bathe  (K,M,C,u0,v0,h,nSteps);
        case 'RK1',                 [t,U] = mdof_rk     (K,M,C,u0,v0,h,nSteps,1);
        case 'RK2',                 [t,U] = mdof_rk     (K,M,C,u0,v0,h,nSteps,2);
        case 'RK3',                 [t,U] = mdof_rk     (K,M,C,u0,v0,h,nSteps,3);
        case 'RK4',                 [t,U] = mdof_rk     (K,M,C,u0,v0,h,nSteps,4);
        otherwise, error('run_scheme:name','Unknown method "%s".', name);
    end
end
