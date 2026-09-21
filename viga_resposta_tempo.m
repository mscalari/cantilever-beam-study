function tresp = viga_resposta_tempo(cfg)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Tip response in time, generalised   %
%  for the Estudos B-F (with zoom and  %
%  axis-cut support - advisor's        %
%  revision, items 21-22).             %
%--------------------------------------%
%
% Generalises Sections 4 and 5 of viga_main_simulacao.m to an arbitrary
% case (any K,M,C / initial condition / modal reference), and adds:
%
%   - a DISPLAY WINDOW anchored on the decay of the response. With
%     damping the response dies out, and showing the last cycles of the
%     SIMULATION would show a flat line; so the window is placed on the
%     last cycles BEFORE the exact response drops below cfg.tolDecaimento
%     of its own peak. Without damping there is nothing to detect and the
%     window is simply the end of the simulation. See janela_plot below
%     and section 2 of viga_config_estudos.m;
%   - a "zoom" inset on a late oscillation window (item 21: "dar um zoom
%     lá na oscilação 20"), so that curves which look coincident at the
%     start can be told apart after many periods;
%   - an optional hard cap on the main figure's time axis, "tMaxShow"
%     (item 22: the response often visually vanishes well before the
%     window's right edge);
%   - automatic Central-Differences inclusion/exclusion in the figure,
%     decided from its OWN critical step at the case's hImpl (instead of
%     the hard-coded skip used in Estudo A, which was correct only
%     because DC always exceeded its limit at h/T1 = 0.01);
%   - N = t_final/h recorded per family, for item 31 (error vs. number
%     of steps).
%
% viga_main_simulacao.m is left untouched - Estudo A's figures are
% unaffected by anything in this file.
%
% INPUT - cfg (struct), required fields
%   .K, .M, .C     system matrices used by the integrators (mdof_*.m):
%                  the full model's or a reduced model's (see
%                  viga_reducao_modal.m)
%   .u0, .v0       initial condition in the SAME basis as K,M,C
%   .Mref          mass matrix used by viga_solucao_exata for the modal
%                  reference: mdl.M for a full model, eye(m) for a
%                  reduced one
%   .wAll, .PhiAll, .xiAll   full modal data for viga_solucao_exata, in
%                  the SAME basis as Mref (mdl.w/mdl.Phi/mdl.damp.xi for
%                  the full model; red.w/eye(m)/red.xi for a reduced one)
%   .tipExtract    function handle, state history [rows x nt] -> tip
%                  response [1 x nt]. For the full model this is
%                  @(U) U(tipIdx,:); for a reduced model,
%                  @(Q) red.tipRow*Q (see viga_reducao_modal.m)
%   .T1            reference period (mdl.T1, or red.T1 for a reduced
%                  model) used to scale h and the time axis
%   .hImpl, .hRK   time steps for the direct / Runge-Kutta families
%   .nPerImpl, .nPerRK   cycles of T1 COMPUTED, damped case
%   .X_ciclos_final_plot cycles SHOWN, damped case (counted backwards
%                  from the detected decay instant)
%   .YCLIP         y-axis clip of the normalised response
%   .label, .xiLabel, .outDir   filename tag / title string / figure folder
%
% INPUT - cfg (struct), optional fields
%   .nPerImpl_sem_xi, .nPerRK_sem_xi, .X_ciclos_final_plot_sem_xi
%               counterparts used when the case has NO damping; each
%               falls back to its damped equivalent when absent. Which
%               set applies is decided here, from cfg.xiAll - no caller
%               needs to branch on it.
%   .tolDecaimento  default 1e-3 - fraction of its own peak below which
%               the exact response counts as "zero" for the window
%   .verboseJanela  default true - prints one diagnostic line per family
%               saying where the window landed and whether cycles were
%               computed past the decay and thrown away
%   .dirNames   default {'DC','Houbolt','Newmark','Wilson','Bathe'}
%   .dirLabel   default {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'}
%   .rkOrders   default [1 2 3 4]
%   .rkNames    default {'RK1','RK2','RK3','RK4'}
%   .doZoom     default false. When true, adds a zoomed inset to BOTH
%               figures on the window cfg.zoomCycles.
%   .zoomCycles [c1 c2] in periods of T1 (default: the LAST 1.5 periods
%               of the displayed window)
%   .tMaxShow   caps the main axis at this many periods of T1 (default:
%               the full displayed window, i.e. no cut - same as before)
%   .doPlots    default true
%   .runDirect,.runRK   default true,true - set false to skip a family
%   .nPontosExato  default 4000 - points on the exact reference curve
%   .nMaxSalvar    default 2e5  - max samples STORED per RK method
%   .nMaxPlot      default 5000 - max points actually drawn per curve
%                  (these three are memory/plot knobs only; every step is
%                  still integrated - see viga_config_estudos.m, sec. 4)
%
% OUTPUT
%   tresp : struct with .exact, .(method) time histories (direct and
%           RK), .config (steps, N = t_final/h per family)

    if ~isfield(cfg,'dirNames'), cfg.dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'}; end
    if ~isfield(cfg,'dirLabel'), cfg.dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'}; end
    if ~isfield(cfg,'rkOrders'), cfg.rkOrders = [1 2 3 4]; end
    if ~isfield(cfg,'rkNames'),  cfg.rkNames  = {'RK1','RK2','RK3','RK4'}; end
    if ~isfield(cfg,'doZoom'),   cfg.doZoom   = false; end
    if ~isfield(cfg,'tMaxShow'), cfg.tMaxShow = []; end
    if ~isfield(cfg,'doPlots'),  cfg.doPlots  = true; end
    if ~isfield(cfg,'runDirect'),cfg.runDirect= true; end
    if ~isfield(cfg,'runRK'),    cfg.runRK    = true; end
    if ~isfield(cfg,'FS'),       cfg.FS       = 18; end
    if ~isfield(cfg,'FIGPOS'),   cfg.FIGPOS   = [0 0 20 14]; end
    % Amostragem/armazenamento - ver secao 4 de viga_config_estudos.m. Os
    % defaults abaixo reproduzem os valores que estavam chumbados aqui, de
    % modo que quem chamar sem passar esses campos nao muda de resultado.
    if ~isfield(cfg,'nPontosExato'), cfg.nPontosExato = 4000;  end
    if ~isfield(cfg,'nMaxSalvar'),   cfg.nMaxSalvar   = 2e5;   end
    if ~isfield(cfg,'nMaxPlot'),     cfg.nMaxPlot     = 5000;  end
    if ~isfield(cfg,'tolDecaimento'),cfg.tolDecaimento = 1e-3; end
    if ~isfield(cfg,'verboseJanela'),cfg.verboseJanela = true; end

    if cfg.doPlots && ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

    coresDir = [1.00 0.00 0.00; 0.00 0.20 1.00; 0.10 0.70 0.10; 0.60 0.00 0.60; 1.00 0.50 0.00];
    coresRK  = [1.00 0.00 0.00; 0.00 0.20 1.00; 0.10 0.70 0.10; 0.60 0.00 0.60];

    T1 = cfg.T1;
    tresp = struct();

    % ------------------------------------------------------------------
    %  Amortecido ou nao? A escolha entre os dois conjuntos de contagens
    %  de ciclos (secao 2 de viga_config_estudos.m) e feita AQUI, a partir
    %  do xi do caso recebido, para que nenhum viga_estudo_*.m precise
    %  repetir essa decisao.
    % ------------------------------------------------------------------
    amortecido = any(cfg.xiAll > 0);

    if amortecido
        nPerImpl   = cfg.nPerImpl;
        nPerRK     = cfg.nPerRK;
        ciclosPlot = cfg.X_ciclos_final_plot;
    else
        nPerImpl   = pick(cfg,'nPerImpl_sem_xi', cfg.nPerImpl);
        nPerRK     = pick(cfg,'nPerRK_sem_xi',   cfg.nPerRK);
        ciclosPlot = pick(cfg,'X_ciclos_final_plot_sem_xi', cfg.X_ciclos_final_plot);
    end

    tresp.config.amortecido = amortecido;
    tresp.config.ciclosPlot = ciclosPlot;

    % ======================================================================
    %  DIRECT METHODS
    % ======================================================================
    if cfg.runDirect
        nS = round(nPerImpl*T1/cfg.hImpl);

        tRef = (0:nS)*cfg.hImpl;
        Uex  = viga_solucao_exata(tRef, cfg.Mref, cfg.wAll, cfg.PhiAll, cfg.xiAll, cfg.u0, cfg.v0);
        uexTip = cfg.tipExtract(Uex);

        tresp.exact.t     = tRef;
        tresp.exact.u_tip = uexTip;

        % Janela mostrada: ancorada no decaimento quando ha amortecimento,
        % no fim da simulacao quando nao ha (ver janela_plot abaixo).
        [tIni, tFim, infoJan] = janela_plot(tRef, uexTip, T1, ciclosPlot, ...
            cfg.tolDecaimento, amortecido);
        tresp.config.janelaDireto = infoJan;
        if cfg.verboseJanela
            print_janela('direct schemes', infoJan, T1, nPerImpl);
        end

        mask   = (tRef >= tIni - 0.5*cfg.hImpl) & (tRef <= tFim + 0.5*cfg.hImpl);
        uScale = max(abs(uexTip(mask)));
        if uScale == 0, uScale = 1; end

        % Decide whether DC belongs in the figure: skip only if it is
        % actually unstable at this hImpl (generalises the hard-coded
        % skip of Estudo A, which held only because DC always exceeded
        % its own limit at h/T1 = 0.01 in that study).
        hcrDC   = viga_h_critico(cfg.K, cfg.M, cfg.C, 'DC');
        skipDC  = (hcrDC > 0) && (cfg.hImpl > hcrDC);

        if cfg.doPlots
            fig = figure('Name',sprintf('Tip disp direct %s',cfg.label),'NumberTitle','off');
            plot(tRef(mask)/T1, uexTip(mask)/uScale, 'k--', 'LineWidth', 2); hold on;
            leg = {'Exact'};
        end

        for k = 1:numel(cfg.dirNames)
            [t,U] = run_scheme_direct(cfg.dirNames{k}, cfg.K, cfg.M, cfg.C, cfg.u0, cfg.v0, cfg.hImpl, nS);
            uTip = cfg.tipExtract(U);

            tresp.(cfg.dirNames{k}).t     = t;
            tresp.(cfg.dirNames{k}).h     = cfg.hImpl;
            tresp.(cfg.dirNames{k}).N     = nS;
            tresp.(cfg.dirNames{k}).u_tip = uTip;
            tresp.(cfg.dirNames{k}).err   = max(abs(uTip - uexTip))/max(uScale,eps);

            if strcmpi(cfg.dirNames{k},'DC') && skipDC, continue; end
            if cfg.doPlots
                plot(t(mask)/T1, uTip(mask)/uScale, '-', 'Color', coresDir(k,:), 'LineWidth', 1.3);
                leg{end+1} = cfg.dirLabel{k}; %#ok<AGROW>
            end
        end

        if cfg.doPlots
            grid on;
            xlabel('$t/T_1$', 'Interpreter','latex');
            ylabel('$w_{tip}/w_{ref}$', 'Interpreter','latex');
            % 'tex', not 'latex': cfg.xiLabel shape varies by caller (plain
            % text, '\xi = ...', 'modo N, \xi = ...', ...) and 'latex'
            % requires the whole string to be valid LaTeX, which a bare
            % '_' or '\xi' outside '$...$' is not (see viga_espectro.m).
            title(sprintf('Tip displacement, direct schemes (%s)', cfg.xiLabel), 'Interpreter','tex');
            legend(leg, 'Location','northeast', 'FontSize', 12);
            xHi = tFim/T1;
            if ~isempty(cfg.tMaxShow), xHi = min(xHi, cfg.tMaxShow); end
            xlim([tIni/T1 xHi]);  ylim([-cfg.YCLIP cfg.YCLIP]);
            set(fig,'Units','centimeters','Position',cfg.FIGPOS);
            set(gca,'FontName','Times New Roman','FontSize',cfg.FS);

            if cfg.doZoom
                add_zoom_inset(gca, tRef, uexTip, uScale, cfg.dirNames, tresp, T1, coresDir, cfg, skipDC);
            end

            print(fig, fullfile(cfg.outDir, sprintf('viga_desloc_direto_%s.eps',cfg.label)), '-depsc2','-painters');
        end

        tresp.config.hImpl = cfg.hImpl;
        tresp.config.nStepsImpl = nS;
        tresp.config.N_impl = nS;   % N = t_final/h
    end

    % ======================================================================
    %  RUNGE-KUTTA METHODS
    % ======================================================================
    if cfg.runRK
        nRK = round(nPerRK*T1/cfg.hRK);

        % A janela precisa ser decidida ANTES de montar a grade fina da
        % solucao exata (que so cobre a janela). Para isso, avalia-se a
        % solucao exata numa grade grosseira sobre a simulacao INTEIRA,
        % barata porque e forma fechada, so para localizar o decaimento.
        tVarre  = linspace(0, nPerRK*T1, max(cfg.nPontosExato, 2000));
        UexV    = viga_solucao_exata(tVarre, cfg.Mref, cfg.wAll, cfg.PhiAll, cfg.xiAll, cfg.u0, cfg.v0);
        uexTipV = cfg.tipExtract(UexV);

        [tIniRK, tFimRK, infoJanRK] = janela_plot(tVarre, uexTipV, T1, ciclosPlot, ...
            cfg.tolDecaimento, amortecido);
        tresp.config.janelaRK = infoJanRK;
        if cfg.verboseJanela
            print_janela('Runge-Kutta', infoJanRK, T1, nPerRK);
        end

        tFine = linspace(tIniRK, tFimRK, cfg.nPontosExato);
        UexF  = viga_solucao_exata(tFine, cfg.Mref, cfg.wAll, cfg.PhiAll, cfg.xiAll, cfg.u0, cfg.v0);
        uexTipF = cfg.tipExtract(UexF);
        uScaleRK = max(abs(uexTipF));
        if uScaleRK == 0, uScaleRK = 1; end

        tresp.exactRK.t     = tFine;
        tresp.exactRK.u_tip = uexTipF;

        if cfg.doPlots
            fig = figure('Name',sprintf('Tip disp RK %s',cfg.label),'NumberTitle','off');
            plot(tFine/T1, uexTipF/uScaleRK, 'k--', 'LineWidth', 2); hold on;
            leg = {'Exact'};
        end

        % Pre-decimate to <= 2e5 step INDICES before simulating, and march
        % with viga_rk_stream_tip.m rather than mdof_rk.m: nRK can reach
        % into the millions (explicit-RK stability, not accuracy, sets
        % it - see viga_h_critico.m), and materialising the full
        % [n x nRK+1] state history for that many steps is what was
        % running this study out of memory. Streaming keeps only the
        % values at these pre-selected steps plus a running max, giving
        % the SAME decimated series decimate_hist would have produced
        % from the full array, without ever holding it.
        saveIdx = unique(round(linspace(0, nRK, min(nRK+1, cfg.nMaxSalvar))));

        % Calculada sob demanda na primeira iteracao e reaproveitada pelas
        % demais: a grade de instantes salvos e a mesma para os quatro RK.
        uexSalvo = [];

        for k = 1:numel(cfg.rkOrders)
            name = cfg.rkNames{k};
            [ts, us, uMax] = viga_rk_stream_tip(cfg.K, cfg.M, cfg.C, cfg.u0, cfg.v0, ...
                cfg.hRK, nRK, cfg.rkOrders(k), cfg.tipExtract, saveIdx);

            tresp.(name).t      = ts;
            tresp.(name).h      = cfg.hRK;
            tresp.(name).N      = nRK;
            tresp.(name).u_tip  = us;
            tresp.(name).u_tip_max = uMax;

            % Erro relativo, na MESMA definicao usada pela familia direta
            % (max|numerico - exato| normalizado), para que a tabela
            % N x erro do item 31 tenha as duas familias na mesma escala.
            % A solucao exata e avaliada nos instantes efetivamente
            % salvos - sao os mesmos passos, so decimados.
            if isempty(uexSalvo)
                UexS = viga_solucao_exata(ts, cfg.Mref, cfg.wAll, cfg.PhiAll, ...
                    cfg.xiAll, cfg.u0, cfg.v0);
                uexSalvo = cfg.tipExtract(UexS);
                clear UexS
            end
            tresp.(name).err = max(abs(us - uexSalvo))/max(uScaleRK, eps);

            if cfg.doPlots
                mask = (ts >= tIniRK - 0.5*cfg.hRK) & (ts <= tFimRK + 0.5*cfg.hRK);
                [tp,yp] = downsample_plot(ts(mask), us(mask)/uScaleRK, cfg.nMaxPlot);
                plot(tp/T1, yp, '-', 'Color', coresRK(k,:), 'LineWidth', 1.3);
                leg{end+1} = name; %#ok<AGROW>
            end
        end

        if cfg.doPlots
            grid on;
            xlabel('$t/T_1$', 'Interpreter','latex');
            ylabel('$w_{tip}/w_{ref}$', 'Interpreter','latex');
            title(sprintf('Tip displacement, Runge-Kutta (%s)', cfg.xiLabel), 'Interpreter','tex');
            legend(leg, 'Location','northeast', 'FontSize', 12);
            xHi = tFimRK/T1;
            if ~isempty(cfg.tMaxShow), xHi = min(xHi, cfg.tMaxShow); end
            xlim([tIniRK/T1 xHi]);  ylim([-cfg.YCLIP cfg.YCLIP]);
            set(fig,'Units','centimeters','Position',cfg.FIGPOS);
            set(gca,'FontName','Times New Roman','FontSize',cfg.FS);
            print(fig, fullfile(cfg.outDir, sprintf('viga_desloc_rk_%s.eps',cfg.label)), '-depsc2','-painters');
        end

        tresp.config.hRK = cfg.hRK;
        tresp.config.nStepsRK = nRK;
        tresp.config.N_RK = nRK;
    end
end

% ==========================================================================
%  LOCAL FUNCTIONS
% ==========================================================================
function v = pick(s, campo, padrao)
% Le s.(campo) se existir, senao devolve o padrao. Serve para que um
% chamador antigo, que ainda nao passe os campos do caso nao amortecido,
% continue funcionando com os valores do caso amortecido.
    if isfield(s, campo) && ~isempty(s.(campo))
        v = s.(campo);
    else
        v = padrao;
    end
end

function [tIni, tFim, info] = janela_plot(t, uex, T1, nCiclos, tol, amortecido)
%JANELA_PLOT  Decide a janela de tempo mostrada na figura.
%
% Com AMORTECIMENTO a resposta decai: passado certo instante ela e
% numericamente zero e plotar os ultimos ciclos da simulacao mostraria
% uma linha reta. Entao a janela e ancorada no DECAIMENTO - o ultimo
% instante em que a resposta exata ainda esta acima de tol vezes o seu
% proprio pico - e recua nCiclos ciclos de T1 a partir dali.
%
% SEM amortecimento nao ha decaimento a detectar, e a janela e
% simplesmente os ultimos nCiclos ciclos da simulacao.
%
% A deteccao usa sempre a solucao EXATA, nunca a resposta de um metodo:
% um esquema instavel diverge e arrastaria a janela para o lugar errado,
% justamente no caso em que a figura mais precisa estar certa.
%
% Como a resposta oscila, |uex| cruza qualquer limiar varias vezes por
% ciclo; o que interessa e o ULTIMO cruzamento, que e onde a envoltoria
% do decaimento atravessa o limiar - dai o find(..., 'last').
%
% OUTPUT
%   tIni, tFim : limites da janela [s]
%   info       : struct com o diagnostico (usado no console e salvo no
%                .mat), campos:
%                  .amortecido   true/false
%                  .tDecaimento  instante detectado [s] (NaN se nao houve)
%                  .decaiu       true se a resposta caiu abaixo de tol
%                                DENTRO da simulacao
%                  .tIni, .tFim  a janela, em segundos
%                  .cicloIni, .cicloFim  a janela, em ciclos de T1
%                  .nCiclosPedidos  nCiclos
%                  .truncada     true se a janela bateu em t=0 e ficou
%                                menor que nCiclos ciclos
%                  .tol          tolerancia usada

    info = struct('amortecido',amortecido, 'tDecaimento',NaN, 'decaiu',false, ...
                  'nCiclosPedidos',nCiclos, 'truncada',false, 'tol',tol);

    uPeak = max(abs(uex));

    if ~amortecido || uPeak == 0
        % Sem amortecimento (ou resposta identicamente nula): janela no fim.
        tFim = t(end);
    else
        idx = find(abs(uex) > tol*uPeak, 1, 'last');
        if isempty(idx) || idx == numel(t)
            % Nunca caiu abaixo da tolerancia dentro da simulacao: a
            % resposta ainda tem sinal no fim, entao a janela vai ate la.
            tFim = t(end);
        else
            tFim = t(idx);
            info.tDecaimento = tFim;
            info.decaiu = true;
        end
    end

    tIni = tFim - nCiclos*T1;
    if tIni < t(1)
        tIni = t(1);
        info.truncada = true;
    end

    info.tIni = tIni;      info.tFim = tFim;
    info.cicloIni = tIni/T1;  info.cicloFim = tFim/T1;
end

function print_janela(familia, info, T1, nPerCalculado)
% Diagnostico de uma linha no console: onde a janela caiu e, quando for o
% caso, quantos ciclos foram calculados a toa depois do decaimento.
    if info.decaiu
        fprintf(['   [janela %s] resposta decaiu a %.1e do pico em t/T1 = %.2f; ' ...
                 'plotando t/T1 em [%.2f, %.2f]\n'], ...
                familia, info.tol, info.tDecaimento/T1, info.cicloIni, info.cicloFim);
        sobra = nPerCalculado - info.cicloFim;
        if sobra > 1
            fprintf(['   [janela %s] AVISO: %.1f ciclos calculados depois do ' ...
                     'decaimento foram descartados - da para reduzir nPer para ~%d\n'], ...
                    familia, sobra, ceil(info.cicloFim) + 1);
        end
    else
        if info.amortecido
            fprintf(['   [janela %s] resposta NAO caiu abaixo de %.1e do pico dentro ' ...
                     'dos %d ciclos calculados; plotando o fim: t/T1 em [%.2f, %.2f]\n'], ...
                    familia, info.tol, nPerCalculado, info.cicloIni, info.cicloFim);
        else
            fprintf('   [janela %s] sem amortecimento; plotando t/T1 em [%.2f, %.2f]\n', ...
                    familia, info.cicloIni, info.cicloFim);
        end
    end
    if info.truncada
        fprintf(['   [janela %s] AVISO: janela truncada em t=0 - foram pedidos %g ciclos ' ...
                 'mas so ha %.2f antes do decaimento\n'], ...
                familia, info.nCiclosPedidos, info.cicloFim);
    end
end

function [t, U] = run_scheme_direct(name, K, M, C, u0, v0, h, nSteps)
    switch name
        case 'DC',      [t,U] = mdof_dc     (K,M,C,u0,v0,h,nSteps);
        case 'Houbolt', [t,U] = mdof_houbolt(K,M,C,u0,v0,h,nSteps);
        case 'Newmark', [t,U] = mdof_newmark(K,M,C,u0,v0,h,nSteps);
        case 'Wilson',  [t,U] = mdof_wilson (K,M,C,u0,v0,h,nSteps);
        case 'Bathe',   [t,U] = mdof_bathe  (K,M,C,u0,v0,h,nSteps);
        otherwise, error('viga_resposta_tempo:name','Unknown method "%s".', name);
    end
end

function [ts, ys] = downsample_plot(t, y, nMax)
    if nargin < 3, nMax = 5000; end
    if numel(t) <= nMax
        ts = t;  ys = y;  return
    end
    idx = unique([1, round(linspace(1, numel(t), nMax)), numel(t)]);
    ts = t(idx);  ys = y(idx);
end

function add_zoom_inset(axMain, tRef, uexTip, uScale, dirNames, tresp, T1, coresDir, cfg, skipDC)
% Small inset axes zooming on a late oscillation window (item 21).
% Window: cfg.zoomCycles = [c1 c2] in periods of T1, or, by default, the
% last 1.5 periods of the window already displayed on axMain. Both the
% exact curve and every method are normalised by the SAME uScale used on
% the main axes, so the inset is a literal magnification of it.
    xl = xlim(axMain);
    if isfield(cfg,'zoomCycles') && ~isempty(cfg.zoomCycles)
        zc = cfg.zoomCycles;
    else
        zc = [max(xl(2)-1.5, xl(1)), xl(2)];
    end

    pos = get(axMain,'Position');
    axZoom = axes('Position',[pos(1)+0.55*pos(3), pos(2)+0.12*pos(4), 0.40*pos(3), 0.38*pos(4)]);
    hold(axZoom,'on'); box(axZoom,'on');

    mZ = (tRef/T1 >= zc(1)) & (tRef/T1 <= zc(2));
    plot(axZoom, tRef(mZ)/T1, uexTip(mZ)/uScale, 'k--', 'LineWidth', 1.5);

    for k = 1:numel(dirNames)
        if strcmpi(dirNames{k},'DC') && skipDC, continue; end
        if ~isfield(tresp, dirNames{k}), continue; end
        t  = tresp.(dirNames{k}).t;
        u  = tresp.(dirNames{k}).u_tip;
        mZk = (t/T1 >= zc(1)) & (t/T1 <= zc(2));
        plot(axZoom, t(mZk)/T1, u(mZk)/uScale, '-', ...
            'Color', coresDir(k,:), 'LineWidth', 1.1);
    end
    xlim(axZoom, zc);
    grid(axZoom,'on');
    set(axZoom,'FontSize',max(cfg.FS-6,8));
    title(axZoom, sprintf('zoom: $t/T_1 \\in [%.1f,%.1f]$', zc(1), zc(2)), ...
        'Interpreter','latex','FontSize',max(cfg.FS-6,8));
end
