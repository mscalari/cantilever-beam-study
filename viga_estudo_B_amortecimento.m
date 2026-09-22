%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  ESTUDO B - Effect of physical       %
%  damping (advisor's revision).       %
%--------------------------------------%
%
% Repeats the impulse-response study of Estudo A (viga_main_simulacao.m)
% for xi = 0, 0.01 and 0.2 (cfg.B.xi_list), instead of only 0 and 0.2. The
% advisor identifies xi ~ 0.005-0.02 as typical of metallic structures
% and xi = 0.01 in particular as the case where numerical damping should
% become comparatively visible (items 16-20 and 33/Estudo B).
%
% METHODOLOGY NOTE - same shared Runge-Kutta step as Estudo A. Exactly as
% in viga_main_simulacao.m (Sections 3 and 5), a SINGLE Runge-Kutta step
% hRK, sized at 90% of the critical step of the MOST damped case
% (xi = max(cfg.B.xi_list)), is used across every xi value. This is what
% lets RK1/RK2 be shown diverging at low damping and stable at high
% damping WITH THE SAME TIME STEP - the central point of Section 4.1 of
% the paper, now extended to three damping levels instead of two. Only
% the implicit-scheme step (hImpl = 0.01*T1) and the per-case critical
% steps in the "hcr" table are computed independently for every xi.
%
% Run after viga_main_validacao.m. Does NOT modify viga_main_simulacao.m
% or its output (Estudo A stays the untouched reference case). Results
% and figures go to ImgViga/EstudoB/.

clc; clear; close all;
tInicio = now;                       %#ok<TNOW1>

% Saida SEMPRE em <pasta do codigo>/ImgViga/EstudoB (nao em <pwd>); cada
% figura sai em .eps E .png (viga_salva_figura.m).
outDir = viga_pasta_saida('EstudoB');
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

% Toda a configuracao numerica vem de viga_config_estudos.m - nao ha
% numero de simulacao chumbado neste arquivo. Para mudar passo, numero de
% periodos, quais xi varrer ou o visual das figuras, edite LA.
cfg = viga_config_estudos();
FS = cfg.FS; FIGPOS = cfg.FIGPOS;
hRat = cfg.hRat;

p = viga_params();
[K, M, info] = viga_matrizes_EF(p);
[w, Phi]     = viga_autovalores_EF(K, M);
T1   = 2*pi/w(1);
iTip = info.tipW;

dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'};
dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'};
rkNames  = {'RK1','RK2','RK3','RK4'};
condM    = [rkNames, {'DC'}];

xiList = cfg.B.xi_list;
[xiRef, iRef] = max(xiList);   %#ok<ASGLU>  % most-damped case sets the shared RK step

fprintf('=====================================================\n');
fprintf(' ESTUDO B - PHYSICAL DAMPING SWEEP  (xi = %s)\n', mat2str(xiList));
fprintf('=====================================================\n');

% --- Damping matrices and modal xi for every case, plus per-case h_cr --
Cs    = cell(numel(xiList),1);
ximoS = cell(numel(xiList),1);
dampS = cell(numel(xiList),1);
hcrS  = zeros(numel(condM), numel(xiList));

for ix = 1:numel(xiList)
    xiL = xiList(ix);
    if xiL == 0
        Cs{ix}    = zeros(size(K));
        ximoS{ix} = zeros(size(w));
        dampS{ix} = struct('a0',0,'a1',0,'xi',ximoS{ix});
    else
        pp = p; pp.xi_target = xiL;
        [Cs{ix}, dtmp] = viga_amortecimento(K, M, w, pp);
        ximoS{ix} = dtmp.xi;
        dampS{ix} = dtmp;
    end
    for k = 1:numel(condM)
        hcrS(k,ix) = viga_h_critico(K, M, Cs{ix}, condM{k});
    end
end

hcrTable = array2table(hcrS, 'VariableNames', ...
    matlab.lang.makeValidName(compose('xi_%.3g',xiList)), 'RowNames', condM);
disp(' Critical time steps per damping case [s]:'); disp(hcrTable);

hImpl = hRat*T1;
hRK   = cfg.fatorSegurancaRK*min(hcrS(1:4, iRef));
fprintf(' Shared Runge-Kutta step: hRK = %.4e s (from xi = %.3g, the most damped case)\n\n', hRK, xiRef);

res = struct();
res.model  = struct('K',K,'M',M,'w',w,'Phi',Phi,'T1',T1,'wmax',w(end),'info',info);
res.hcr    = hcrTable;
res.hShared.hImpl = hImpl;
res.hShared.hRK   = hRK;

for ix = 1:numel(xiList)
    xiL = xiList(ix);
    tag = matlab.lang.makeValidName(sprintf('xi_%.3g', xiL));
    fprintf('--- xi = %.3f (%s) ---\n', xiL, tag);

    C     = Cs{ix};
    ximo  = ximoS{ix};
    res.damp.(tag) = dampS{ix};
    [u0, v0] = viga_ci_impulso(K, M, C, info, p);

    % --- Spectral radius --------------------------------------------
    stab = viga_espectro(w, ximo, sprintf('\\xi = %.3g', xiL), ['B_' tag], outDir, ...
        'hSweep', cfg.hSweep, 'FS', FS, 'FIGPOS', FIGPOS, ...
        'xScaleBy', 1/T1, 'xLabelStr', '$h/T_1$');
    res.stability.(tag) = stab;

    % --- PE/AD, evaluated at the SAME steps used in Section 5 --------
    allNames  = [dirNames, rkNames];
    allLabels = [dirLabel, rkNames];
    allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
    pead = viga_pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
        'nModesReport', cfg.nModesReport, ...
        'printHeader', sprintf('Estudo B - PE/AD, xi = %.3g', xiL));
    res.pead.(tag) = pead;

    % --- Time-domain response, with a zoom for the low-damping case --
    cfgT = struct();
    cfgT.K = K; cfgT.M = M; cfgT.C = C;
    cfgT.u0 = u0; cfgT.v0 = v0;
    cfgT.Mref = M; cfgT.wAll = w; cfgT.PhiAll = Phi; cfgT.xiAll = ximo;
    cfgT.tipExtract = @(U) U(iTip,:);
    cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = hRK;
    cfgT.label = ['B_' tag]; cfgT.xiLabel = sprintf('\\xi = %.3g', xiL);
    cfgT.outDir = outDir;
    cfgT = viga_aplica_config(cfgT, cfg);
    % items 20-21: zoom on the damping case picked in cfg.B.xiZoom
    cfgT.doZoom = ~isempty(cfg.B.xiZoom) && (xiL == cfg.B.xiZoom);

    tresp = viga_resposta_tempo(cfgT);
    res.time.(tag) = tresp;
end

res.config = struct('hRat',hRat,'hImpl',hImpl,'hRK',hRK,'xiList',xiList, ...
    'runDate',datestr(now));
res.cfg = cfg;   % configuracao completa usada nesta corrida (rastreabilidade)

save(fullfile(outDir,'estudoB_resultados.mat'), 'res', '-v7.3');
fprintf('\nEstudo B concluido. Resultados e figuras em %s\n', outDir);
viga_lista_obsoletos({outDir}, tInicio);
