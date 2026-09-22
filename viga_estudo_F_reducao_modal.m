%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  ESTUDO F - Modal reduction          %
%  (truncation), items 11-15 and 32.   %
%--------------------------------------%
%
% PART 1 - full comparison. Runs the same three outputs as every other
% study (spectral radius, PE/AD, time response) for THREE models: the
% full mesh (20 DOF) and two modal truncations (viga_reducao_modal.m)
% keeping only the first 4 or 3 modes (cfg.F.reduction_modes). Each model
% uses its OWN critical step, exactly as Estudo A does for the full
% model - this is what answers items 11-14 ("a resposta relevante pode
% ser representada pelo sistema reduzido?").
%
% PART 2 - the advisor's specific item 32 experiment, isolating the
% effect of a LARGER h from the effect of the reduction itself, using
% RK4 (a representative conditionally-stable method) on the 4-mode
% reduction:
%
%   Caso 1 - full model,     h = 0.9*h_cr,full(RK4)
%   Caso 2 - reduced model,  h = SAME as Caso 1  (still far below its
%                            own, larger, critical step)
%   Caso 3 - reduced model,  h = 0.9*h_cr,red(RK4)  (exploits the
%                            reduced model's higher critical step)
%
% All three are compared against the FULL model's EXACT (time-
% integration-error-free) tip response, so the reported error mixes
% modal truncation error (present in Casos 2-3) and time-integration
% error (present in all three) exactly as the physical response would.
%
% Damping: xi = p.xi_target = 0.2 (Estudo A's value), so the comparison
% sits on the same physical case already in the paper.
%
% Run after viga_main_validacao.m. Results and figures go to
% ImgViga/EstudoF/.

clc; clear; close all;
tInicio = now;                       %#ok<TNOW1>

% Saida SEMPRE em <pasta do codigo>/ImgViga/EstudoF (nao em <pwd>); cada
% figura sai em .eps E .png (viga_salva_figura.m).
outDir = viga_pasta_saida('EstudoF');
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

% Toda a configuracao numerica vem de viga_config_estudos.m (secoes 1-4 e
% cfg.F) - nao ha numero de simulacao chumbado neste arquivo.
cfg = viga_config_estudos();
FS = cfg.FS; FIGPOS = cfg.FIGPOS; YCLIP = cfg.YCLIP;
nPerImpl = cfg.nPerImpl; nPerRK = cfg.nPerRK;
hRat = cfg.hRat;
reductionModes = cfg.F.reduction_modes;

p   = viga_params();                 % xi_target = 0.2 by default
mdl = viga_modelo(p);
K = mdl.K; M = mdl.M; C = mdl.C; T1 = mdl.T1; iTip = mdl.info.tipW;

dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'};
dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'};
rkNames  = {'RK1','RK2','RK3','RK4'};
condM    = [rkNames, {'DC'}];

fprintf('=====================================================\n');
fprintf(' ESTUDO F - MODAL REDUCTION  (xi = %.2f)\n', p.xi_target);
fprintf('=====================================================\n');
fprintf(' Full model: w_max = %.4e rad/s (mode %d)\n', mdl.wmax, mdl.info.ndof);

res = struct();
res.model_full = struct('K',K,'M',M,'w',mdl.w,'Phi',mdl.Phi,'xi',mdl.damp.xi, ...
    'T1',T1,'wmax',mdl.wmax);

% ======================================================================
%  PART 1 - full comparison: completo, 4 modos, 3 modos
% ======================================================================
% 'tag' e identificador (nome de campo no .mat e sufixo de arquivo);
% 'label' e o texto que vai para titulo/legenda das figuras, em ingles.
models = struct('tag',{},'label',{},'K',{},'M',{},'C',{},'w',{},'xi',{},'Phi',{},'Mref',{}, ...
    'u0',{},'v0',{},'tipExtract',{},'T1',{},'wmax',{});

models(1) = struct('tag','full', 'label','full (20 DOF)', ...
    'K',K,'M',M,'C',C, 'w',mdl.w,'xi',mdl.damp.xi, ...
    'Phi',mdl.Phi,'Mref',M, 'u0',mdl.u0,'v0',mdl.v0, ...
    'tipExtract',@(U) U(iTip,:), 'T1',T1,'wmax',mdl.wmax);

reds = cell(numel(reductionModes),1);
for i = 1:numel(reductionModes)
    m = reductionModes(i);
    red = viga_reducao_modal(mdl, m);
    reds{i} = red;
    models(end+1) = struct('tag',sprintf('red%d',m), 'label',sprintf('reduced, %d modes',m), ...
        'K',red.K,'M',red.M,'C',red.C, ...
        'w',red.w,'xi',red.xi, 'Phi',eye(m),'Mref',red.Mref, ...
        'u0',red.u0,'v0',red.v0, 'tipExtract',@(Q) red.tipRow*Q, ...
        'T1',red.T1,'wmax',red.wmax); %#ok<SAGROW>
    fprintf(' Modelo reduzido (%d modos): w_max = %.4e rad/s (%.4f x menor que o completo)\n', ...
        m, red.wmax, mdl.wmax/red.wmax);
end

for i = 1:numel(models)
    mo = models(i);
    fprintf('\n--- Modelo: %s ---\n', mo.tag);

    stab = viga_espectro(mo.w, mo.xi, sprintf('model: %s', mo.label), ['F_' mo.tag], outDir, ...
        'hSweep', cfg.hSweep, 'FS', FS, 'FIGPOS', FIGPOS, ...
        'xScaleBy', 1/mo.T1, 'xLabelStr', '$h/T_1$');
    res.stability.(mo.tag) = stab;

    hcr = zeros(numel(condM),1);
    for k = 1:numel(condM)
        hcr(k) = viga_h_critico(mo.K, mo.M, mo.C, condM{k});
    end
    res.hcr.(mo.tag) = cell2table([condM(:), num2cell(hcr), num2cell(hcr/mo.T1)], ...
        'VariableNames', {'method','hcr','hcr_T1'});

    hImpl = hRat*mo.T1;
    hRK   = cfg.fatorSegurancaRK*min(hcr(1:4));

    allNames  = [dirNames, rkNames];
    allLabels = [dirLabel, rkNames];
    allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
    % min(...) porque um modelo reduzido pode ter menos modos que
    % cfg.nModesReport - nao da para reportar o modo 3 de um modelo de 2.
    pead = viga_pead_tabela(mo.w, mo.xi, allNames, allLabels, allSteps, mo.T1, ...
        'nModesReport', min(cfg.nModesReport, numel(mo.w)), ...
        'printHeader', sprintf('Estudo F - PE/AD, modelo %s', mo.tag));
    res.pead.(mo.tag) = pead;

    cfgT = struct();
    cfgT.K = mo.K; cfgT.M = mo.M; cfgT.C = mo.C;
    cfgT.u0 = mo.u0; cfgT.v0 = mo.v0;
    cfgT.Mref = mo.Mref; cfgT.wAll = mo.w; cfgT.PhiAll = mo.Phi; cfgT.xiAll = mo.xi;
    cfgT.tipExtract = mo.tipExtract;
    cfgT.T1 = mo.T1; cfgT.hImpl = hImpl; cfgT.hRK = max(hRK,eps);
    cfgT.label = ['F_' mo.tag];
    cfgT.xiLabel = sprintf('model: %s', mo.label);
    cfgT.outDir = outDir;
    cfgT = viga_aplica_config(cfgT, cfg);
    cfgT.runRK = hRK > 0;

    tresp = viga_resposta_tempo(cfgT);
    res.time.(mo.tag) = tresp;
end

% --- Truncation-only error: reduced EXACT solution vs full EXACT solution
tCommon = (0:round(nPerImpl*T1/(hRat*T1)))*(hRat*T1);
UexFull = viga_solucao_exata(tCommon, M, mdl.w, mdl.Phi, mdl.damp.xi, mdl.u0, mdl.v0);
uTipFull = UexFull(iTip,:);
for i = 1:numel(reds)
    red = reds{i};
    UexRed = viga_solucao_exata(tCommon, red.Mref, red.w, eye(red.m), red.xi, red.u0, red.v0);
    uTipRed = red.tipRow*UexRed;
    errTrunc = max(abs(uTipRed - uTipFull))/max(abs(uTipFull));
    fprintf('\n Erro de truncamento modal (so a CI/resposta exata, sem integracao no tempo):\n');
    fprintf('   %d modos: erro relativo maximo na ponta = %.4e\n', red.m, errTrunc);
    res.truncation.(sprintf('red%d',red.m)).err = errTrunc;
    res.truncation.(sprintf('red%d',red.m)).t   = tCommon;
    res.truncation.(sprintf('red%d',red.m)).u_full = uTipFull;
    res.truncation.(sprintf('red%d',red.m)).u_red  = uTipRed;
end

% ======================================================================
%  PART 2 - item 32: mesmo h x h maior, RK4, reducao de 4 modos
% ======================================================================
fprintf('\n=====================================================\n');
fprintf(' PARTE 2 - Item 32: mesmo h x h maior (RK4, 4 modos)\n');
fprintf('=====================================================\n');

red4 = reds{1};                        % cfg.F.reduction_modes(1), hoje = 4
rkP2 = sprintf('RK%d', cfg.F.item32_ordemRK);
hcrFullRK4 = viga_h_critico(K, M, C, rkP2);
hcrRedRK4  = viga_h_critico(red4.K, red4.M, red4.C, rkP2);
h1 = cfg.fatorSegurancaRK*hcrFullRK4;  % Caso 1 e 2: mesmo h
h3 = cfg.fatorSegurancaRK*hcrRedRK4;   % Caso 3: h maior, baseado no wmax reduzido

fprintf(' h_cr,full(%s) = %.4e s -> h1=h2 = %.4e s\n', rkP2, hcrFullRK4, h1);
fprintf(' h_cr,red%d(%s) = %.4e s -> h3     = %.4e s  (%.2fx maior)\n', ...
    red4.m, rkP2, hcrRedRK4, h3, h3/h1);

nPerP2 = nPerRK;
% 'nome' e identificador (campo no .mat); 'label' vai para a legenda.
casos = struct('nome',{'Case1_full_h1','Case2_red4_h1','Case3_red4_h3'}, ...
    'label',{'Case 1: full, h_1','Case 2: reduced, same h_1','Case 3: reduced, larger h_3'}, ...
    'K',{K,red4.K,red4.K}, 'M',{M,red4.M,red4.M}, 'C',{C,red4.C,red4.C}, ...
    'u0',{mdl.u0,red4.u0,red4.u0}, 'v0',{mdl.v0,red4.v0,red4.v0}, ...
    'tipExtract',{@(U) U(iTip,:), @(Q) red4.tipRow*Q, @(Q) red4.tipRow*Q}, ...
    'h',{h1,h1,h3});

tRefFine = linspace(0, nPerP2*T1, cfg.nPontosExato);
UexFine  = viga_solucao_exata(tRefFine, M, mdl.w, mdl.Phi, mdl.damp.xi, mdl.u0, mdl.v0);
uexTipFine = UexFine(iTip,:);
uPeakP2  = max(abs(uexTipFine));    % pico global: normaliza o ERRO

% Mesma regra de janela dos demais estudos (viga_janela_plot.m): com xi =
% 0.2 a resposta morre em poucos ciclos, entao mostra-se a resposta
% inteira em vez dos nPerP2 ciclos, quase todos achatados. So o que e
% DESENHADO muda; a simulacao e o erro seguem cobrindo os nPerP2 ciclos.
[tIniP2, tFimP2, infoP2] = viga_janela_plot(tRefFine, uexTipFine, T1, ...
    cfg.X_ciclos_final_plot, cfg.tolDecaimento, true, cfg.ciclosDecaimentoRapido);
maskP2 = (tRefFine >= tIniP2) & (tRefFine <= tFimP2);
% Figura normalizada pelo pico DENTRO da janela (como nos outros estudos).
uScaleP2 = max(abs(uexTipFine(maskP2)));
txtJanP2 = '';
if infoP2.completa, txtJanP2 = ' (resposta inteira)'; end
fprintf(' Item 32: janela plotada t/T1 em [%.2f, %.2f]%s\n', ...
    infoP2.cicloIni, infoP2.cicloFim, txtJanP2);

fig = figure('Name','Estudo F - item 32','NumberTitle','off'); hold on;
plot(tRefFine(maskP2)/T1, uexTipFine(maskP2)/uScaleP2, 'k--', 'LineWidth', 2);
legP2 = {'Exact (full)'};
coresP2 = [0.85 0.1 0.1; 0.1 0.1 0.85; 0.1 0.6 0.1];

for ic = 1:numel(casos)
    cs = casos(ic);
    nSteps = round(nPerP2*T1/cs.h);
    try
        % Caso1 (full model, h1) can need millions of RK steps (stability
        % at the full mesh's w_max), so march with the memory-light
        % streaming helper (see viga_resposta_tempo.m) instead of
        % mdof_rk.m directly: a decimated series is all interp1 needs
        % below anyway, since it is only re-sampled onto the
        % cfg.nPontosExato-point tRefFine grid.
        saveIdx = unique(round(linspace(0, nSteps, min(nSteps+1, cfg.nMaxSalvar))));
        [t, uTip] = viga_rk_stream_tip(cs.K, cs.M, cs.C, cs.u0, cs.v0, cs.h, ...
            nSteps, cfg.F.item32_ordemRK, cs.tipExtract, saveIdx);
        % Interpolate onto the common fine grid for a fair error metric.
        uInterp = interp1(t, uTip, tRefFine, 'linear', 'extrap');
        errRel  = max(abs(uInterp - uexTipFine))/uPeakP2;
        fprintf('   %-16s h=%.3e s  N=%d passos  erro relativo max = %.4e\n', ...
            cs.nome, cs.h, nSteps, errRel);
        res.item32.(cs.nome) = struct('h',cs.h,'N',nSteps,'err',errRel);

        plot(tRefFine(maskP2)/T1, uInterp(maskP2)/uScaleP2, '-', 'Color', coresP2(ic,:), 'LineWidth', 1.3);
        legP2{end+1} = cs.label; %#ok<AGROW>
    catch ME
        warning('Item 32, %s: falhou (%s).', cs.nome, ME.message);
        res.item32.(cs.nome) = struct('h',cs.h,'failed',true,'message',ME.message);
    end
end

grid on;
xlabel('$t/T_1$', 'Interpreter','latex');
ylabel('$w_{tip}/w_{ref}$', 'Interpreter','latex');
title(sprintf('%s: full vs. reduced model (same h vs. larger h)', rkP2), ...
    'Interpreter','tex');
legend(legP2, 'Location','northeast', 'FontSize', 12);
xlim([tIniP2 tFimP2]/T1);  ylim([-YCLIP YCLIP]);
set(fig,'Units','centimeters','Position',FIGPOS);
set(gca,'FontName','Times New Roman','FontSize',FS);
viga_salva_figura(fig, outDir, 'viga_item32_completo_vs_reduzido');

res.config = struct('hRat',hRat,'reduction_modes',reductionModes, ...
    'xi',p.xi_target,'runDate',datestr(now));
res.cfg = cfg;   % configuracao completa usada nesta corrida (rastreabilidade)

save(fullfile(outDir,'estudoF_resultados.mat'), 'res', '-v7.3');
fprintf('\nEstudo F concluido. Resultados e figuras em %s\n', outDir);
viga_lista_obsoletos({outDir}, tInicio);
