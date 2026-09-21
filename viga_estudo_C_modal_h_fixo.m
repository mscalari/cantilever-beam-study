%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  ESTUDO C - Modal initial conditions,%
%  fixed time step (advisor's item 5). %
%--------------------------------------%
%
% "mantendo o passo de tempo constante, dá para visualizar como os
%  métodos se comportam para cada condição inicial diferente" - the
% advisor's own words (item 5). Displaces the beam according to a single
% mode shape at a time (cfg.C.modos = [1 2 3 5 10], viga_ci_modal.m)
% and holds h FIXED across every mode, so that any change observed comes
% from the mode's own frequency, not from a re-tuned time step.
%
% h is fixed in ABSOLUTE terms (h = 0.01*T1, the same value already used
% throughout Estudo A), not relative to each mode's own period - the
% point of the study is precisely to see how h/T_modo changes as the
% excited mode gets higher (item 24), so h/T_modo must be free to vary.
%
% Damping: xi = 0.01, the low-damping case the advisor asked to look at
% (item 20), so that any numerical dissipation stands out against a
% small physical damping instead of being masked by it.
%
% h_max of the mesh - and hence the Runge-Kutta stability limit - is set
% by w_max of the FULL mesh, not by which mode is used as the initial
% condition (K, M, C are unchanged; only u0 changes). The RK step is
% therefore computed ONCE, exactly as for Estudo B's xi = 0.01 case, and
% shared by every modal case here.
%
% Run after viga_main_validacao.m. Results and figures go to
% ImgViga/EstudoC/.

clc; clear; close all;

outDir = fullfile(pwd,'ImgViga','EstudoC');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Toda a configuracao numerica vem de viga_config_estudos.m (secoes 1-4 e
% cfg.C) - nao ha numero de simulacao chumbado neste arquivo.
cfg = viga_config_estudos();
FS = cfg.FS; FIGPOS = cfg.FIGPOS;
hRat = cfg.hRat;
xiL  = cfg.C.xi;
modos = cfg.C.modos;

p = viga_params();
[K, M, info] = viga_matrizes_EF(p);
[w, Phi]     = viga_autovalores_EF(K, M);
T1   = 2*pi/w(1);
iTip = info.tipW;

pp = p; pp.xi_target = xiL;
[C, damp] = viga_amortecimento(K, M, w, pp);
ximo = damp.xi;

dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'};
dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'};
rkNames  = {'RK1','RK2','RK3','RK4'};
condM    = [rkNames, {'DC'}];

hImpl = hRat*T1;
hcr   = zeros(numel(condM),1);
for k = 1:numel(condM)
    hcr(k) = viga_h_critico(K, M, C, condM{k});
end
hRK = cfg.fatorSegurancaRK*min(hcr(1:4));

fprintf('=====================================================\n');
fprintf(' ESTUDO C - MODAL ICs, FIXED h  (xi = %.3f)\n', xiL);
fprintf('=====================================================\n');
fprintf(' hImpl = %.4e s (h/T1 = %.4f), hRK = %.4e s\n', hImpl, hImpl/T1, hRK);
fprintf(' modes: %s\n\n', mat2str(modos));

res = struct();
res.model = struct('K',K,'M',M,'w',w,'Phi',Phi,'T1',T1,'wmax',w(end),'info',info);
res.damp  = damp;
res.hShared.hImpl = hImpl;
res.hShared.hRK   = hRK;

for im = 1:numel(modos)
    modeIdx = modos(im);
    tag = sprintf('modo%d', modeIdx);
    Tmodo = 2*pi/w(modeIdx);
    fprintf('--- Modo %d: w = %.4f rad/s, T = %.4e s, h/T_modo = %.4e ---\n', ...
        modeIdx, w(modeIdx), Tmodo, hImpl/Tmodo);

    [u0, v0] = viga_ci_modal(K, M, C, Phi, modeIdx);

    % --- PE/AD at the mode actually excited (and its neighbours) ------
    allNames  = [dirNames, rkNames];
    allLabels = [dirLabel, rkNames];
    allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
    % Relata pelo menos cfg.nModesReport modos, mas sempre o suficiente
    % para que o modo EXCITADO apareca na tabela.
    pead = viga_pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
        'nModesReport', max(modeIdx, cfg.nModesReport), ...
        'printHeader', sprintf('Estudo C - PE/AD, modo %d excitado', modeIdx));
    res.pead.(tag) = pead;
    res.Tmodo.(tag) = Tmodo;
    res.hT_modo.(tag) = hImpl/Tmodo;

    % --- Time-domain response ----------------------------------------
    cfgT = struct();
    cfgT.K = K; cfgT.M = M; cfgT.C = C;
    cfgT.u0 = u0; cfgT.v0 = v0;
    cfgT.Mref = M; cfgT.wAll = w; cfgT.PhiAll = Phi; cfgT.xiAll = ximo;
    cfgT.tipExtract = @(U) U(iTip,:);
    cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = hRK;
    cfgT.label = ['C_' tag];
    cfgT.xiLabel = sprintf('mode %d, \\xi = %.3f', modeIdx, xiL);
    cfgT.outDir = outDir;
    cfgT = viga_aplica_config(cfgT, cfg);

    tresp = viga_resposta_tempo(cfgT);
    res.time.(tag) = tresp;
end

res.config = struct('hRat',hRat,'hImpl',hImpl,'hRK',hRK,'xi',xiL, ...
    'modes',modos,'runDate',datestr(now));
res.cfg = cfg;   % configuracao completa usada nesta corrida (rastreabilidade)

save(fullfile(outDir,'estudoC_resultados.mat'), 'res', '-v7.3');
fprintf('\nEstudo C concluido. Resultados e figuras em %s\n', outDir);
