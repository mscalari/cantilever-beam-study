%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  ESTUDO E - Rayleigh damping NOT     %
%  anchored to a target xi.            %
%--------------------------------------%
%
% The advisor suggested a0 = 0 with a1 = 1e-4 or 1e-5 (item 18) instead
% of always solving for (a0,a1) that hit a chosen xi at two anchor
% modes. Both cases are run here (cfg.E.rayleigh). Because a0 = 0, only
% the stiffness-proportional term is present, so xi_k = a1*w_k/2 grows
% LINEARLY with frequency from xi_1 itself - there is no "anchor" and no
% two modes forced to share the same xi. The resulting xi per mode is
% COMPUTED (viga_amortecimento.m, direct-coefficient mode), not imposed,
% and reported for the first several modes (item 18's checklist).
%
% Same impulse initial condition as Estudo A, so the two Rayleigh
% coefficients can be compared directly against xi = 0.2 (Estudo A) and
% against each other.
%
% Run after viga_main_validacao.m. Results and figures go to
% ImgViga/EstudoE/.

clc; clear; close all;
tInicio = now;                       %#ok<TNOW1>

% Saida SEMPRE em <pasta do codigo>/ImgViga/EstudoE (nao em <pwd>); cada
% figura sai em .eps E .png (viga_salva_figura.m).
outDir = viga_pasta_saida('EstudoE');
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

% Toda a configuracao numerica vem de viga_config_estudos.m (secoes 1-4 e
% cfg.E) - nao ha numero de simulacao chumbado neste arquivo.
cfg = viga_config_estudos();
FS = cfg.FS; FIGPOS = cfg.FIGPOS;
hRat = cfg.hRat;
rayleigh = cfg.E.rayleigh;

p = viga_params();
[K, M, info] = viga_matrizes_EF(p);
[w, Phi]     = viga_autovalores_EF(K, M);
T1   = 2*pi/w(1);
iTip = info.tipW;

dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'};
dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'};
rkNames  = {'RK1','RK2','RK3','RK4'};
condM    = [rkNames, {'DC'}];

fprintf('=====================================================\n');
fprintf(' ESTUDO E - RAYLEIGH DAMPING a0 = 0, a1 GIVEN\n');
fprintf('=====================================================\n');

res = struct();
res.model = struct('K',K,'M',M,'w',w,'Phi',Phi,'T1',T1,'wmax',w(end),'info',info);

for ic = 1:numel(rayleigh)
    a0 = rayleigh(ic).a0;
    a1 = rayleigh(ic).a1;
    tag = sprintf('a1_%s', matlab.lang.makeValidName(sprintf('%.0e',a1)));
    fprintf('--- a0 = %.4g, a1 = %.4g ---\n', a0, a1);

    [C, damp] = viga_amortecimento(K, M, w, p, a0, a1);
    ximo = damp.xi;

    fprintf('   xi_1..xi_5 = %s\n', mat2str(round(ximo(1:5),4)));
    res.damp.(tag) = damp;

    [u0, v0] = viga_ci_impulso(K, M, C, info, p);

    % --- Spectral radius ---------------------------------------------
    xiLabel = sprintf('a_0=0, a_1=%.0e', a1);
    stab = viga_espectro(w, ximo, xiLabel, ['E_' tag], outDir, ...
        'hSweep', cfg.hSweep, 'FS', FS, 'FIGPOS', FIGPOS, ...
        'xScaleBy', 1/T1, 'xLabelStr', '$h/T_1$');
    res.stability.(tag) = stab;

    % --- Critical steps -------------------------------------------------
    hcr = zeros(numel(condM),1);
    for k = 1:numel(condM)
        hcr(k) = viga_h_critico(K, M, C, condM{k});
    end
    res.hcr.(tag) = cell2table([condM(:), num2cell(hcr), num2cell(hcr/T1)], ...
        'VariableNames', {'method','hcr','hcr_T1'});

    hImpl = hRat*T1;
    hRK   = cfg.fatorSegurancaRK*min(hcr(1:4));

    % --- PE/AD ----------------------------------------------------------
    allNames  = [dirNames, rkNames];
    allLabels = [dirLabel, rkNames];
    allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
    pead = viga_pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
        'nModesReport', cfg.nModesReport, ...
        'printHeader', sprintf('Estudo E - PE/AD, a0=0, a1=%.0e', a1));
    res.pead.(tag) = pead;

    % --- Time-domain response --------------------------------------------
    cfgT = struct();
    cfgT.K = K; cfgT.M = M; cfgT.C = C;
    cfgT.u0 = u0; cfgT.v0 = v0;
    cfgT.Mref = M; cfgT.wAll = w; cfgT.PhiAll = Phi; cfgT.xiAll = ximo;
    cfgT.tipExtract = @(U) U(iTip,:);
    cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = max(hRK, eps);
    cfgT.label = ['E_' tag]; cfgT.xiLabel = xiLabel;
    cfgT.outDir = outDir;
    cfgT = viga_aplica_config(cfgT, cfg);
    cfgT.runRK = hRK > 0;

    tresp = viga_resposta_tempo(cfgT);
    res.time.(tag) = tresp;
end

res.config = struct('hRat',hRat,'rayleigh_alt',rayleigh,'runDate',datestr(now));
res.cfg = cfg;   % configuracao completa usada nesta corrida (rastreabilidade)

save(fullfile(outDir,'estudoE_resultados.mat'), 'res', '-v7.3');
fprintf('\nEstudo E concluido. Resultados e figuras em %s\n', outDir);
viga_lista_obsoletos({outDir}, tInicio);
