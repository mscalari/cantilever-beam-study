%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  ESTUDO D - Modal initial conditions,%
%  time step swept per mode.           %
%--------------------------------------%
%
% Second half of the modal-IC experiment (item 25 of the advisor's
% document): for a few representative modes (cfg.D.modos = [1 3 10],
% low/medium/high), sweep h/T_modo through cfg.D.hRatios = [0.01 0.10 0.90]
% and see how PE, AD and stability move as h grows relative to the
% excited mode's OWN period - as opposed to Estudo C, which holds h
% fixed. Damping xi = 0.01, as in Estudo C.
%
% COST CONTROL (see the plan's risk section). PE/AD and the spectral
% radius are closed-form (viga_operador.m / viga_pe_ad.m): computed for
% EVERY (mode, ratio) combination, no time integration involved, so this
% part is cheap regardless of how small h gets.
%
% Full time-domain simulation is different: h/T_modo = 0.90 for a high
% mode can land far above the GLOBAL stability limit (set by w_max of
% the mesh, not by the excited mode - see viga_h_critico.m), especially
% for RK1/RK2/DC, and the point of this study is precisely to let that
% happen. To keep the run bounded, the time-domain figures are produced
% ONLY for the smallest and largest ratio of each mode (skipping the
% middle one) and with few periods (nPerD = 3, not 12), and every run is
% wrapped in try/catch so a genuine numerical blow-up in one case is
% reported and skipped instead of stopping the whole script.
%
% Run after viga_main_validacao.m. Results and figures go to
% ImgViga/EstudoD/.

clc; clear; close all;
tInicio = now;                       %#ok<TNOW1>

% Saida SEMPRE em <pasta do codigo>/ImgViga/EstudoD (nao em <pwd>); cada
% figura sai em .eps E .png (viga_salva_figura.m).
outDir = viga_pasta_saida('EstudoD');
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

% Toda a configuracao numerica vem de viga_config_estudos.m (secoes 1-4 e
% cfg.D) - nao ha numero de simulacao chumbado neste arquivo. Note que os
% periodos aqui sao os de cfg.D (menores que os globais), pelo motivo
% explicado no cabecalho e documentado no arquivo de configuracao.
cfg = viga_config_estudos();
FS = cfg.FS; FIGPOS = cfg.FIGPOS;
nPerD = cfg.D.nPer; ciclosPlotD = cfg.D.ciclosPlot;
xiL  = cfg.D.xi;

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
allNames  = [dirNames, rkNames];
allLabels = [dirLabel, rkNames];

modes  = cfg.D.modos;
ratios = cfg.D.hRatios;

fprintf('=====================================================\n');
fprintf(' ESTUDO D - MODAL ICs, h/T_modo VARYING  (xi = %.3f)\n', xiL);
fprintf('=====================================================\n');
fprintf(' modos: %s   razoes h/T_modo: %s\n\n', mat2str(modes), mat2str(ratios));

res = struct();
res.model = struct('K',K,'M',M,'w',w,'Phi',Phi,'T1',T1,'wmax',w(end),'info',info);
res.damp  = damp;

% --- Closed-form PE/AD/rho for EVERY (mode, ratio) combination --------
peadGrid = struct();
for im = 1:numel(modes)
    modeIdx = modes(im);
    Tmodo = 2*pi/w(modeIdx);
    tagM = sprintf('modo%d', modeIdx);
    for ir = 1:numel(ratios)
        h = ratios(ir)*Tmodo;
        tagR = matlab.lang.makeValidName(sprintf('r_%.3g', ratios(ir)));
        allSteps = repmat(h, 1, numel(allNames));
        pead = viga_pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
            'nModesReport', max(modeIdx, cfg.nModesReport), ...
            'printHeader', sprintf('Estudo D - modo %d, h/T_modo = %.3g (h=%.3e s)', ...
                modeIdx, ratios(ir), h));
        peadGrid.(tagM).(tagR) = pead;
    end
end
res.pead_grid = peadGrid;

% --- Time-domain figures: which ratios actually get simulated ----------
% cfg.D.ratiosSimuladas = 'extremos' (padrao, contem o custo) ou 'todas'.
% Em qualquer caso, PE/AD em forma fechada ja foram calculados acima para
% TODAS as razoes - isto controla so a simulacao no tempo.
switch lower(cfg.D.ratiosSimuladas)
    case 'todas'
        ratiosPlot = 1:numel(ratios);
    case 'extremos'
        [~, iMin] = min(ratios);
        [~, iMax] = max(ratios);
        ratiosPlot = unique([iMin iMax]);
    otherwise
        error('viga_estudo_D:ratiosSimuladas', ...
            'cfg.D.ratiosSimuladas deve ser ''extremos'' ou ''todas'' (recebido: ''%s'').', ...
            cfg.D.ratiosSimuladas);
end

for im = 1:numel(modes)
    modeIdx = modes(im);
    Tmodo = 2*pi/w(modeIdx);
    tagM = sprintf('modo%d', modeIdx);
    [u0, v0] = viga_ci_modal(K, M, C, Phi, modeIdx);

    for ir = ratiosPlot
        h = ratios(ir)*Tmodo;
        tagR = matlab.lang.makeValidName(sprintf('r_%.3g', ratios(ir)));
        fprintf('--- modo %d, h/T_modo = %.3g (h = %.3e s) ---\n', modeIdx, ratios(ir), h);

        cfgT = struct();
        cfgT.K = K; cfgT.M = M; cfgT.C = C;
        cfgT.u0 = u0; cfgT.v0 = v0;
        cfgT.Mref = M; cfgT.wAll = w; cfgT.PhiAll = Phi; cfgT.xiAll = ximo;
        cfgT.tipExtract = @(U) U(iTip,:);
        cfgT.T1 = T1; cfgT.hImpl = h; cfgT.hRK = h;
        cfgT.label = sprintf('D_%s_%s', tagM, tagR);
        cfgT.xiLabel = sprintf('mode %d, h/T_{mode} = %.3g', modeIdx, ratios(ir));
        cfgT.outDir = outDir;
        % Estudo D usa contagens proprias de ciclos (cfg.D), bem menores
        % que as globais - ver o cabecalho e a secao 5 da configuracao.
        cfgT = viga_aplica_config(cfgT, cfg, ...
            'nPerImpl', nPerD, 'nPerRK', nPerD, ...
            'nPerImpl_sem_xi', nPerD, 'nPerRK_sem_xi', nPerD, ...
            'X_ciclos_final_plot', ciclosPlotD, ...
            'X_ciclos_final_plot_sem_xi', ciclosPlotD);
        % Only ONE family is actually meaningful per figure here (both
        % use the same h, so run both but keep them on separate plots).

        try
            tresp = viga_resposta_tempo(cfgT);
            res.time.(tagM).(tagR) = tresp;
        catch ME
            warning('Estudo D: modo %d, razao %.3g divergiu ou falhou (%s). Pulando.', ...
                modeIdx, ratios(ir), ME.message);
            res.time.(tagM).(tagR) = struct('failed', true, 'message', ME.message);
        end
    end
end

res.config = struct('modes',modes,'ratios',ratios,'xi',xiL, ...
    'nPerD',nPerD,'ciclosPlotD',ciclosPlotD,'runDate',datestr(now));
res.cfg = cfg;   % configuracao completa usada nesta corrida (rastreabilidade)

save(fullfile(outDir,'estudoD_resultados.mat'), 'res', '-v7.3');
fprintf('\nEstudo D concluido. Resultados e figuras em %s\n', outDir);
viga_lista_obsoletos({outDir}, tInicio);
