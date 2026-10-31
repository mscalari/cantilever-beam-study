function g = graficos()
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%
% Figuras e arquivos de saida: raio espectral, resposta no tempo, janela de
% plot, paleta de cores e salvamento em .eps e .png.
%
% As funcoes deste arquivo sao chamadas de fora como graficos().nome(...).
% So as listadas abaixo sao publicas; as demais sao de uso interno.

    g.aplica_config      = @aplica_config;
    g.caso_viga_completa = @caso_viga_completa;
    g.cores              = @cores;
    g.downsample_plot    = @downsample_plot;
    g.espectro           = @espectro;
    g.janela_plot        = @janela_plot;
    g.lista_obsoletos    = @lista_obsoletos;
    g.pasta_saida        = @pasta_saida;
    g.print_janela       = @print_janela;
    g.resposta_tempo     = @resposta_tempo;
    g.salva_figura       = @salva_figura;
end

%% ====================================================================
%  Figura de raio espectral
%  ====================================================================

function stab = espectro(w, xi, xiLabel, label, outDir, varargin)
% Figura de raio espectral dos nove metodos em um caso (w, xi): diretos em
% linha cheia e Runge-Kutta tracejados, contra o mesmo limite rho = 1.
%
% ENTRADA
%   w, xi   : [m x 1] frequencias e amortecimentos do caso
%   xiLabel : texto do titulo, ex. '\xi = 0.01'
%   label   : sufixo do arquivo, salvo como viga_rho_<label>
%   outDir  : pasta de saida
%   opcoes nome/valor: 'hSweep', 'methodsRK', 'methodsDir', 'dirLabels', 'FS',
%   'FIGPOS', 'doPlots', 'xScaleBy', 'xLabelStr'
%
% SAIDA
%   stab : .hSweep, .rho.direto.<metodo>, .rho.rk.<metodo>
    q = inputParser;
    q.addParameter('hSweep',     logspace(-9,-2,200));
    q.addParameter('methodsRK',  {'RK1','RK2','RK3','RK4'});
    q.addParameter('methodsDir', {'DC','Houbolt','Newmark','Wilson','Bathe'});
    q.addParameter('dirLabels',  {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'});
    q.addParameter('FS',      18);
    q.addParameter('FIGPOS',  [0 0 20 14]);
    q.addParameter('doPlots', true);
    q.addParameter('xScaleBy', 1);              % e.g. 1/T1 to plot h/T1
    q.addParameter('xLabelStr', '$h$ [s]');
    q.parse(varargin{:});
    o = q.Results;

    if o.doPlots && ~exist(outDir,'dir'), mkdir(outDir); end

    [coresDir, coresRK] = cores();

    stab.hSweep = o.hSweep;

    % 'tag' e o nome do campo dentro de stab.rho.
    fams = struct('tag',{'direto','rk'}, ...
                   'names',{o.methodsDir,o.methodsRK}, ...
                   'labels',{o.dirLabels,o.methodsRK}, ...
                   'cores',{coresDir,coresRK}, 'estilo',{'-','--'});

    if o.doPlots
        fig = figure('Name',sprintf('rho %s',label),'NumberTitle','off');
        hold on;
    end
    leg = {};
    for f = 1:numel(fams)
        names = fams(f).names;
        for k = 1:numel(names)
            r = arrayfun(@(hh) metodos().rho_beam(names{k}, w, xi, hh), o.hSweep);
            stab.rho.(fams(f).tag).(names{k}) = r;
            if o.doPlots
                plot(o.hSweep*o.xScaleBy, r, fams(f).estilo, ...
                    'Color', fams(f).cores(k,:), 'LineWidth', 1.5);
                leg{end+1} = fams(f).labels{k}; %#ok<AGROW>
            end
        end
    end

    if o.doPlots
        set(gca,'XScale','log','YScale','log');
        yline(1,'k--','HandleVisibility','off');
        grid on;
        xlabel(o.xLabelStr, 'Interpreter','latex');
        ylabel('$\rho(\mathbf{\hat{A}})$', 'Interpreter','latex');
        % 'tex' aceita texto comum misturado com '\xi' e '_'; 'latex' exigiria $...$.
        title(sprintf('Spectral radius, direct schemes and Runge-Kutta (%s)', xiLabel), ...
            'Interpreter','tex');
        legend(leg, 'Location','northwest', 'FontSize', 11, 'NumColumns', 2);
        ylim([1e-2 1e2]);
        set(fig,'Units','centimeters','Position',o.FIGPOS);
        set(gca,'FontName','Times New Roman','FontSize',o.FS);
        salva_figura(fig, outDir, sprintf('viga_rho_%s',label));
    end
end

function [coresDir, coresRK] = cores()
% Paleta unica dos nove metodos: cada metodo tem a mesma cor em toda figura.
% RK1 usa grena em vez do vermelho puro porque as duas familias aparecem
% juntas no raio espectral e RK1 ficaria identico a DC.
%
% SAIDA
%   coresDir : [5 x 3] RGB, ordem DC, Houbolt, Newmark, Wilson-Theta, Bathe
%   coresRK  : [4 x 3] RGB, ordem RK1, RK2, RK3, RK4
    coresDir = [1.0 0.0 0.0;    % DC           - vermelho
                0.0 0.2 1.0;    % Houbolt      - azul
                0.1 0.7 0.1;    % Newmark      - verde
                0.6 0.0 0.6;    % Wilson-Theta - roxo
                1.0 0.5 0.0];   % Bathe        - laranja

    coresRK  = [0.6 0.0 0.0;    % RK1 - grena, para nao repetir o vermelho da DC
                0.0 0.0 1.0;    % RK2 - azul   ('b')
                0.0 1.0 0.0;    % RK3 - verde  ('g')
                1.0 0.0 1.0];   % RK4 - magenta ('m')
end

%% ====================================================================
%  Resposta no tempo
%  ====================================================================

function tresp = resposta_tempo(cfg)
% Resposta da ponta no tempo para um caso qualquer, usada nos Estudos B-F:
% integra as duas familias, compara com a solucao exata, calcula o erro e
% desenha as figuras. A janela mostrada vem de janela_plot e a figura
% direta pode ganhar um inset com zoom (cfg.doZoom).
%
% O erro e medido na simulacao inteira e normalizado pelo pico GLOBAL da
% resposta exata, a mesma definicao do Estudo A.
%
% ENTRADA - cfg, campos obrigatorios
%   .K, .M, .C     matrizes do sistema (completo ou reduzido)
%   .u0, .v0       condicao inicial na mesma base
%   .Mref          massa usada na solucao exata (M ou eye(m))
%   .wAll, .PhiAll, .xiAll   dados modais na mesma base de Mref
%   .tipExtract    handle historico -> resposta da ponta
%   .T1            periodo de referencia
%   .hImpl, .hRK   passos das familias direta e Runge-Kutta
%   .nPerImpl, .nPerRK, .X_ciclos_final_plot   ciclos calculados e mostrados
%   .YCLIP         limite do eixo y normalizado
%   .label, .xiLabel, .outDir   sufixo do arquivo, titulo e pasta
%
% ENTRADA - cfg, campos opcionais
%   .nPerImpl_sem_xi, .nPerRK_sem_xi, .X_ciclos_final_plot_sem_xi   caso sem
%                  amortecimento (na falta, valem os do caso amortecido)
%   .tolDecaimento, .ciclosDecaimentoRapido, .verboseJanela
%   .dirNames, .dirLabel, .rkOrders, .rkNames
%   .doZoom, .zoomCycles   inset na figura direta
%   .tMaxShow      corta o eixo do tempo neste numero de periodos
%   .doPlots, .runDirect, .runRK
%   .nPontosExato, .nMaxSalvar, .nMaxPlot   amostragem (nao mudam a integracao)
%
% SAIDA
%   tresp : .exact, .exactRK, .<metodo> (t, h, N, u_tip, err) e .config
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
    if ~isfield(cfg,'nPontosExato'), cfg.nPontosExato = 4000;  end
    if ~isfield(cfg,'nMaxSalvar'),   cfg.nMaxSalvar   = 2e5;   end
    if ~isfield(cfg,'nMaxPlot'),     cfg.nMaxPlot     = 5000;  end
    if ~isfield(cfg,'tolDecaimento'),cfg.tolDecaimento = 1e-3; end
    if ~isfield(cfg,'ciclosDecaimentoRapido'), cfg.ciclosDecaimentoRapido = 6; end
    if ~isfield(cfg,'verboseJanela'),cfg.verboseJanela = true; end

    if cfg.doPlots && ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

    [coresDir, coresRK] = cores();

    T1 = cfg.T1;
    tresp = struct();

    % O xi do caso decide qual conjunto de contagens de ciclos vale.
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
    %  METODOS DIRETOS
    % ======================================================================
    if cfg.runDirect
        nS = round(nPerImpl*T1/cfg.hImpl);

        tRef = (0:nS)*cfg.hImpl;
        Uex  = modelo().solucao_exata(tRef, cfg.Mref, cfg.wAll, cfg.PhiAll, cfg.xiAll, cfg.u0, cfg.v0);
        uexTip = cfg.tipExtract(Uex);

        tresp.exact.t     = tRef;
        tresp.exact.u_tip = uexTip;

        % Janela mostrada: resposta inteira se o amortecimento a mata em
        % poucos ciclos; senao ancorada no decaimento; sem amortecimento,
        % o fim da simulacao (ver janela_plot).
        [tIni, tFim, infoJan] = janela_plot(tRef, uexTip, T1, ciclosPlot, ...
            cfg.tolDecaimento, amortecido, cfg.ciclosDecaimentoRapido);
        tresp.config.janelaDireto = infoJan;
        if cfg.verboseJanela
            print_janela('direct schemes', infoJan, T1, nPerImpl);
        end

        mask   = (tRef >= tIni - 0.5*cfg.hImpl) & (tRef <= tFim + 0.5*cfg.hImpl);
        uScale = max(abs(uexTip(mask)));
        if uScale == 0, uScale = 1; end

        % O ERRO e medido na simulacao INTEIRA e normalizado pelo pico
        % GLOBAL da resposta exata (mesma definicao do Estudo A). Nao pode
        % usar uScale: com amortecimento a janela cai no rabo da resposta,
        % onde o pico e uma fracao pequena do global, e o erro sairia
        % inflado por esse fator.
        uPeak = max(abs(uexTip));
        if uPeak == 0, uPeak = 1; end

        % DC fica fora da figura so quando esta instavel neste passo.
        hcrDC   = metodos().h_critico(cfg.K, cfg.M, cfg.C, 'DC');
        skipDC  = (hcrDC > 0) && (cfg.hImpl > hcrDC);

        if cfg.doPlots
            fig = figure('Name',sprintf('Tip disp direct %s',cfg.label),'NumberTitle','off');
            plot(tRef(mask)/T1, uexTip(mask)/uScale, 'k--', 'LineWidth', 2); hold on;
            leg = {'Exact'};
        end

        for k = 1:numel(cfg.dirNames)
            [t,U] = metodos().integra_metodo(cfg.dirNames{k}, cfg.K, cfg.M, cfg.C, cfg.u0, cfg.v0, cfg.hImpl, nS);
            uTip = cfg.tipExtract(U);

            tresp.(cfg.dirNames{k}).t     = t;
            tresp.(cfg.dirNames{k}).h     = cfg.hImpl;
            tresp.(cfg.dirNames{k}).N     = nS;
            tresp.(cfg.dirNames{k}).u_tip = uTip;
            tresp.(cfg.dirNames{k}).err   = max(abs(uTip - uexTip))/max(uPeak,eps);

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
            % 'tex' aceita texto comum misturado com '\xi' e '_' (ver espectro).
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

            salva_figura(fig, cfg.outDir, sprintf('viga_desloc_direto_%s',cfg.label));
        end

        tresp.config.hImpl = cfg.hImpl;
        tresp.config.nStepsImpl = nS;
        tresp.config.N_impl = nS;   % N = t_final/h
    end

    % ======================================================================
    %  METODOS RUNGE-KUTTA
    % ======================================================================
    if cfg.runRK
        nRK = round(nPerRK*T1/cfg.hRK);

        % A janela precisa ser decidida ANTES de montar a grade fina da
        % solucao exata (que so cobre a janela). Para isso, avalia-se a
        % solucao exata numa grade grosseira sobre a simulacao INTEIRA,
        % barata porque e forma fechada, so para localizar o decaimento.
        tVarre  = linspace(0, nPerRK*T1, max(cfg.nPontosExato, 2000));
        UexV    = modelo().solucao_exata(tVarre, cfg.Mref, cfg.wAll, cfg.PhiAll, cfg.xiAll, cfg.u0, cfg.v0);
        uexTipV = cfg.tipExtract(UexV);

        [tIniRK, tFimRK, infoJanRK] = janela_plot(tVarre, uexTipV, T1, ciclosPlot, ...
            cfg.tolDecaimento, amortecido, cfg.ciclosDecaimentoRapido);
        tresp.config.janelaRK = infoJanRK;
        if cfg.verboseJanela
            print_janela('Runge-Kutta', infoJanRK, T1, nPerRK);
        end

        tFine = linspace(tIniRK, tFimRK, cfg.nPontosExato);
        UexF  = modelo().solucao_exata(tFine, cfg.Mref, cfg.wAll, cfg.PhiAll, cfg.xiAll, cfg.u0, cfg.v0);
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

        % Os indices salvos sao escolhidos antes de simular, e metodos().rk_stream_tip
        % guarda so esses pontos: nRK pode passar de 1e6 passos.
        saveIdx = unique(round(linspace(0, nRK, min(nRK+1, cfg.nMaxSalvar))));

        % Calculada sob demanda na primeira iteracao e reaproveitada pelas
        % demais: a grade de instantes salvos e a mesma para os quatro RK.
        uexSalvo = [];

        for k = 1:numel(cfg.rkOrders)
            name = cfg.rkNames{k};
            [ts, us, uMax] = metodos().rk_stream_tip(cfg.K, cfg.M, cfg.C, cfg.u0, cfg.v0, ...
                cfg.hRK, nRK, cfg.rkOrders(k), cfg.tipExtract, saveIdx);

            tresp.(name).t      = ts;
            tresp.(name).h      = cfg.hRK;
            tresp.(name).N      = nRK;
            tresp.(name).u_tip  = us;
            tresp.(name).u_tip_max = uMax;

            % Erro relativo, na MESMA definicao usada pela familia direta
            % (max|numerico - exato| sobre a simulacao inteira, normalizado
            % pelo pico GLOBAL da exata), para que a tabela N x erro do
            % item 31 tenha as duas familias na mesma escala. A solucao
            % exata e avaliada nos instantes efetivamente salvos - sao os
            % mesmos passos, so decimados.
            if isempty(uexSalvo)
                UexS = modelo().solucao_exata(ts, cfg.Mref, cfg.wAll, cfg.PhiAll, ...
                    cfg.xiAll, cfg.u0, cfg.v0);
                uexSalvo = cfg.tipExtract(UexS);
                clear UexS
                uPeakRK = max(abs([uexTipV, uexSalvo]));
                if uPeakRK == 0, uPeakRK = 1; end
            end
            tresp.(name).err = max(abs(us - uexSalvo))/max(uPeakRK, eps);

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
            salva_figura(fig, cfg.outDir, sprintf('viga_desloc_rk_%s',cfg.label));
        end

        tresp.config.hRK = cfg.hRK;
        tresp.config.nStepsRK = nRK;
        tresp.config.N_RK = nRK;
    end
end

function cfgT = aplica_config(cfgT, cfg, varargin)
% Copia de cfg para cfgT (uma chamada de resposta_tempo) os campos que
% sao repassados direto: ciclos, janela, figuras e amostragem.
%
% ENTRADA
%   cfgT     : struct em construcao
%   cfg      : saida de modelo().config_estudos
%   varargin : pares nome/valor que sobrescrevem cfg (usado no Estudo D)
%
% SAIDA
%   cfgT : a mesma struct, com os campos preenchidos
    campos = { ...
        'nPerImpl', 'nPerRK', ...
        'nPerImpl_sem_xi', 'nPerRK_sem_xi', ...
        'X_ciclos_final_plot', 'X_ciclos_final_plot_sem_xi', ...
        'tolDecaimento', 'ciclosDecaimentoRapido', ...
        'FS', 'FIGPOS', 'YCLIP', ...
        'nPontosExato', 'nMaxSalvar', 'nMaxPlot'};

    for k = 1:numel(campos)
        if isfield(cfg, campos{k})
            cfgT.(campos{k}) = cfg.(campos{k});
        end
    end

    % --- Overrides opcionais -------------------------------------------
    if mod(numel(varargin), 2) ~= 0
        error('viga_aplica_config:pares', ...
            'Os overrides devem vir em pares nome/valor (recebidos %d argumentos).', ...
            numel(varargin));
    end
    for k = 1:2:numel(varargin)
        nome = varargin{k};
        if ~ischar(nome) && ~isstring(nome)
            error('viga_aplica_config:nome', ...
                'O nome do override %d nao e texto.', (k+1)/2);
        end
        if ~ismember(char(nome), campos)
            error('viga_aplica_config:desconhecido', ...
                ['"%s" nao e um campo repassavel. Campos validos: %s. ' ...
                 '(Campos que nao sao pass-through, como label ou outDir, ' ...
                 'sao atribuidos direto em cfgT pelo proprio estudo.)'], ...
                char(nome), strjoin(campos, ', '));
        end
        cfgT.(char(nome)) = varargin{k+1};
    end
end

function cfgT = caso_viga_completa(K, M, C, u0, v0, w, Phi, ximo, iTip)
% Campos de cfgT que descrevem a viga completa em resposta_tempo.
    cfgT = struct();
    cfgT.K = K; cfgT.M = M; cfgT.C = C;
    cfgT.u0 = u0; cfgT.v0 = v0;
    cfgT.Mref = M; cfgT.wAll = w; cfgT.PhiAll = Phi; cfgT.xiAll = ximo;
    cfgT.tipExtract = @(U) U(iTip,:);
end

function v = pick(s, campo, padrao)
% Le s.(campo) se existir; senao devolve o padrao.
    if isfield(s, campo) && ~isempty(s.(campo))
        v = s.(campo);
    else
        v = padrao;
    end
end

function [ts, ys] = downsample_plot(t, y, nMax)
% Reduz uma serie a no maximo nMax pontos (padrao 5000), mantendo o primeiro e
% o ultimo. Usada para desenhar e para guardar historicos longos.
    if nargin < 3, nMax = 5000; end
    if numel(t) <= nMax
        ts = t;  ys = y;  return
    end
    idx = unique([1, round(linspace(1, numel(t), nMax)), numel(t)]);
    ts = t(idx);  ys = y(idx);
end

function add_zoom_inset(axMain, tRef, uexTip, uScale, dirNames, tresp, T1, coresDir, cfg, skipDC)
% Inset com zoom em uma oscilacao tardia. Janela: cfg.zoomCycles = [c1 c2] em
% periodos de T1 ou, por padrao, o ultimo 1.5 periodo da janela mostrada. Usa
% o mesmo uScale do eixo principal, entao e uma ampliacao literal dele.
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

%% ====================================================================
%  Janela de plot
%  ====================================================================

function [tIni, tFim, info] = janela_plot(t, uex, T1, nCiclos, tol, amortecido, cicloMaxCompleta)
% Decide qual trecho da resposta aparece na figura. Sao tres regimes:
%
%   1) Amortecimento forte: se a resposta exata cai abaixo de tol*pico em ate
%      cicloMaxCompleta ciclos de T1, mostra a resposta inteira (de t = 0 ate
%      o decaimento). E onde se ve o transiente inicial.
%   2) Amortecimento fraco: mostra nCiclos antes do decaimento. Se a resposta
%      nao decair dentro da simulacao, a janela vai ate o fim.
%   3) Sem amortecimento: mostra os ultimos nCiclos da simulacao.
%
% A deteccao usa sempre a solucao EXATA: um metodo instavel diverge e levaria
% a janela para o lugar errado. So conta como decaida a resposta que fica
% abaixo da tolerancia por pelo menos um periodo T1 antes do fim; assim um
% cruzamento por zero no ultimo ponto nao e confundido com decaimento.
%
% ENTRADA
%   t, uex           : instantes e resposta exata da ponta na simulacao inteira
%   T1               : periodo de referencia [s]
%   nCiclos          : ciclos mostrados nos regimes 2 e 3
%   tol              : fracao do pico que conta como zero
%   amortecido       : true se xi ~= 0
%   cicloMaxCompleta : limite do regime 1 (opcional, 0 desliga)
%
% SAIDA
%   tIni, tFim : limites da janela [s]
%   info       : diagnostico (.tDecaimento, .decaiu, .completa, .truncada,
%                .cicloIni, .cicloFim, ...), impresso e salvo no .mat
    if nargin < 7 || isempty(cicloMaxCompleta), cicloMaxCompleta = 0; end

    info = struct('amortecido',amortecido, 'tDecaimento',NaN, 'decaiu',false, ...
                  'completa',false, 'nCiclosPedidos',nCiclos, ...
                  'cicloMaxCompleta',cicloMaxCompleta, 'truncada',false, 'tol',tol);

    uPeak = max(abs(uex));

    if ~amortecido || uPeak == 0
        % Sem amortecimento (ou resposta identicamente nula): janela no fim.
        tFim = t(end);
    else
        idx = find(abs(uex) > tol*uPeak, 1, 'last');
        % So conta como "decaiu" se sobrar pelo menos UM periodo T1 inteiro
        % abaixo da tolerancia no fim da simulacao. Sem essa folga, um
        % cruzamento por zero da oscilacao caindo justo no ultimo ponto
        % fazia o ultimo instante acima de tol ser o penultimo ponto, e
        % uma resposta que ainda esta longe de decair era dada como
        % decaida (com xi = 0.01 e tol = 1e-2, aos 50 ciclos ela esta a ~4%).
        if isempty(idx) || t(idx) > t(end) - T1
            % Nao caiu abaixo da tolerancia dentro da simulacao: a
            % resposta ainda tem sinal no fim, entao a janela vai ate la.
            tFim = t(end);
        else
            tFim = t(idx);
            info.tDecaimento = tFim;
            info.decaiu = true;
        end
    end

    if info.decaiu && tFim <= cicloMaxCompleta*T1
        % Regime 1: decaiu em poucos ciclos -> resposta inteira.
        tIni = t(1);
        info.completa = true;
    else
        tIni = tFim - nCiclos*T1;
        if tIni < t(1)
            tIni = t(1);
            info.truncada = true;
        end
    end

    info.tIni = tIni;      info.tFim = tFim;
    info.cicloIni = tIni/T1;  info.cicloFim = tFim/T1;
end

function print_janela(familia, info, T1, nPerCalculado)
% Imprime onde a janela caiu e quantos ciclos foram calculados a toa depois
% do decaimento.
%
% ENTRADA
%   familia       : rotulo impresso, ex. 'direct schemes'
%   info          : struct de janela_plot
%   T1            : periodo de referencia [s]
%   nPerCalculado : ciclos calculados nessa familia
    if info.completa
        fprintf(['   [janela %s] resposta decaiu a %.1e do pico em t/T1 = %.2f ' ...
                 '(<= %g ciclos): plotando a resposta INTEIRA, t/T1 em [%.2f, %.2f]\n'], ...
                familia, info.tol, info.tDecaimento/T1, info.cicloMaxCompleta, ...
                info.cicloIni, info.cicloFim);
    elseif info.decaiu
        fprintf(['   [janela %s] resposta decaiu a %.1e do pico em t/T1 = %.2f; ' ...
                 'plotando t/T1 em [%.2f, %.2f]\n'], ...
                familia, info.tol, info.tDecaimento/T1, info.cicloIni, info.cicloFim);
    end
    if info.decaiu
        % Ciclos calculados depois do decaimento sao jogados fora, em
        % qualquer um dos dois regimes acima.
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

%% ====================================================================
%  Arquivos de saida
%  ====================================================================

function ok = salva_figura(fig, outDir, nome)
% Saida unica de figuras: grava <nome>.eps (vetorial, para o paper) e
% <nome>.png (300 dpi) e confere se os dois existem, nao estao vazios e foram
% gravados agora. Um arquivo aberto em outro programa, por exemplo, nao passa
% em silencio. A falha vira warning, sem interromper a rodada.
%
% ENTRADA
%   fig    : handle da figura
%   outDir : pasta de saida (criada se preciso)
%   nome   : nome do arquivo sem extensao
%
% SAIDA
%   ok : true se os dois arquivos foram gravados e conferidos
    ok = false;
    tEntrada = now;                       %#ok<TNOW1>  (so p/ comparar datenum)

    outDir = char(outDir);
    if ~exist(outDir,'dir')
        try
            mkdir(outDir);
        catch ME
            aviso(outDir, nome, sprintf('nao consegui criar a pasta (%s)', ME.message));
            return
        end
    end
    % Caminho absoluto, para as mensagens nao dependerem do pwd.
    [~, attr] = fileattrib(outDir);
    if isstruct(attr), outDir = attr.Name; end

    if isempty(fig) || ~ishghandle(fig)
        aviso(outDir, nome, 'o handle da figura e invalido (a figura foi fechada?)');
        return
    end

    arqEps = fullfile(outDir, [nome '.eps']);
    arqPng = fullfile(outDir, [nome '.png']);

    try
        print(fig, arqEps, '-depsc2', '-painters');
    catch ME
        aviso(outDir, [nome '.eps'], ME.message);
        return
    end
    try
        print(fig, arqPng, '-dpng', '-r300');
    catch ME
        aviso(outDir, [nome '.png'], ME.message);
        return
    end

    % --- Conferencia: existe, nao esta vazio, e e DESTA execucao --------
    % (2 s de tolerancia: a resolucao do carimbo de data depende do sistema
    %  de arquivos / sincronizacao do OneDrive.)
    tol = 2/86400;
    for arq = {arqEps, arqPng}
        d = dir(arq{1});
        if isempty(d)
            aviso(outDir, arq{1}, 'o arquivo nao existe depois do print');
            return
        elseif d.bytes == 0
            aviso(outDir, arq{1}, 'o arquivo foi criado vazio');
            return
        elseif d.datenum < tEntrada - tol
            aviso(outDir, arq{1}, ['o arquivo NAO foi atualizado (data ' ...
                datestr(d.datenum,'yyyy-mm-dd HH:MM:SS') ') - provavelmente esta ' ...
                'aberto/bloqueado em outro programa']);
            return
        end
    end

    ok = true;
end

function aviso(outDir, arquivo, motivo)
    warning('viga_salva_figura:falhou', ...
        'FIGURA NAO SALVA: %s (pasta %s): %s', arquivo, outDir, motivo);
end

function d = pasta_saida(varargin)
% Pasta de saida: <raiz do projeto>/ImgViga, seja qual for o pwd do MATLAB.
%   pasta_saida()           -> <raiz>/ImgViga
%   pasta_saida('EstudoB')  -> <raiz>/ImgViga/EstudoB
% A pasta nao e criada aqui. A raiz e a pasta acima de src/.
    d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ImgViga', varargin{:});
end

function nomes = lista_obsoletos(pastas, tInicio)
% Lista (sem apagar) as figuras da pasta que nao foram regravadas desde
% tInicio, em geral sobras de nomes ou casos que mudaram.
%
% ENTRADA
%   pastas  : cellstr de pastas, sem recursao
%   tInicio : datenum do inicio da rodada
%
% SAIDA
%   nomes : caminhos das sobras, tambem impressos no console
    nomes = {};
    for p = 1:numel(pastas)
        for ext = {'*.eps','*.png'}
            d = dir(fullfile(pastas{p}, ext{1}));
            for k = 1:numel(d)
                if d(k).datenum < tInicio - 2/86400
                    nomes{end+1,1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
                end
            end
        end
    end

    if ~isempty(nomes)
        fprintf('\n AVISO: %d figura(s) antiga(s) NAO foram regravadas nesta execucao\n', numel(nomes));
        fprintf(' (sobras de execucoes anteriores - podem ser de um nome que mudou):\n');
        for k = 1:numel(nomes)
            d = dir(nomes{k});
            fprintf('   %s   (%s)\n', nomes{k}, datestr(d.datenum,'yyyy-mm-dd HH:MM'));
        end
        fprintf(' Nada foi apagado.\n');
    end
end
