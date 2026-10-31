function viga(acao)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%
% VIGA  Estudo de nove metodos de integracao no tempo aplicados a uma viga
% engastada de aco, discretizada por elementos finitos de Euler-Bernoulli.
%
% Este arquivo e o ponto de entrada. Cada etapa e chamada pelo nome, com a
% pasta src/ no path do MATLAB:
%
%   viga('validacao')  confere o modelo e os integradores (rode primeiro)
%   viga('janela')     teste rapido da janela de plot, so metodos diretos
%   viga('testes')     testes rapidos com assert
%   viga('A')          Estudo A: caso de referencia do paper
%   viga('B') ... viga('F')   Estudos B a F, um de cada vez
%   viga('estudos')    Estudos B a F, em sequencia
%   viga('tabelas')    consolida os resultados em tabelas CSV
%   viga('relatorio')  escreve o relatorio com todos os parametros
%
% As figuras (.eps e .png) e os .mat vao para ImgViga/, na raiz do projeto.
%
% ORGANIZACAO DO CODIGO
%   viga.m        ponto de entrada, validacao e Estudos A a F
%   modelo.m      viga, configuracao, elementos finitos, solucao exata
%   metodos.m     integradores, raio espectral, PE/AD, passo critico
%   graficos.m    figuras, janela de plot, arquivos de saida
%   resultados.m  tabelas, relatorio e testes
%
% Para mudar a viga, edite modelo.params. Para mudar o experimento numerico,
% edite modelo.config_estudos. O Estudo A tem a configuracao no proprio inicio.

    if nargin < 1
        help('viga');
        return
    end

    switch lower(char(acao))
        case 'validacao', validacao();
        case 'janela',    resultados().teste_janela();
        case 'testes',    resultados().testes();
        case 'a',         estudo_A();
        case 'b',         estudo_B();
        case 'c',         estudo_C();
        case 'd',         estudo_D();
        case 'e',         estudo_E();
        case 'f',         estudo_F();
        case 'estudos',   estudos_BF();
        case 'tabelas',   resultados().tabelas();
        case 'relatorio', resultados().relatorio();
        otherwise
            error('viga:acao', 'Acao desconhecida: "%s". Digite help viga.', char(acao));
    end
end

%% ====================================================================
%  Validacao
%  ====================================================================

function validacao()
% Validacao do modelo e dos integradores. Rode antes de confiar em qualquer
% resultado. Tres conferencias:
%
%   (a) ESPACIAL: as frequencias de EF convergem para as analiticas de
%       Euler-Bernoulli, entao o erro observado depois e so temporal.
%   (b) TEMPORAL: os integradores reproduzem a solucao de 1 GDL, cada um na
%       sua ordem de precisao.
%   (c) MODAL: o raio espectral do operador MDOF completo e igual ao maximo
%       dos raios escalares, o que justifica usar metodos().operador na viga.

    close all;


    p = modelo().params();

    fprintf('=====================================================\n');
    fprintf(' (a) SPATIAL DISCRETIZATION\n');
    fprintf('=====================================================\n');
    fprintf('  Steel, E = %.3e Pa, rho = %.0f kg/m^3\n', p.E, p.rho);
    fprintf('  L = %.3f m, b = %.3f m, h = %.4f m, L/h = %.0f\n\n', ...
        p.L, p.b, p.h, p.slenderness);

    nModes = 3;
    wAna   = modelo().modos_analiticos(p, nModes);

    fprintf('  %4s %6s %12s %12s %12s %14s\n', ...
        'N','ndof','err w1','err w2','err w3','w_max');
    for N = [2 5 10 20 40]
        [Kf, Mf, infoF] = modelo().matrizes_EF(p, N);
        wf = modelo().autovalores_EF(Kf, Mf);
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

    mdl = modelo().monta(p);
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
            [~, U] = metodos().integra_metodo(names{k}, Ks, Ms, Cs, 1, 0, hS, nsS);
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
    fprintf('  maximum built from metodos().operador.\n\n');

    fprintf('  %-6s %12s %18s %18s %10s\n','method','h [s]','full operator','modal maximum','diff');
    for meth = {'DC','RK2','RK4'}
        [~, rhoFull] = metodos().h_critico(mdl.K, mdl.M, mdl.C, meth{1});
        for hTest = [1e-6 1e-5 5e-5]
            rFull  = rhoFull(hTest);
            rModal = max(arrayfun(@(k) max(abs(eig( ...
                metodos().operador(meth{1}, mdl.w(k), mdl.damp.xi(k), hTest)))), ...
                1:numel(mdl.w)));
            fprintf('  %-6s %12.1e %18.10f %18.10f %10.1e\n', ...
                meth{1}, hTest, rFull, rModal, abs(rFull-rModal));
        end
    end
    fprintf('\n  Differences at round-off level confirm the modal decomposition.\n');
    fprintf('=====================================================\n');
end

%% ====================================================================
%  Estudo A
%  ====================================================================

function estudo_A()
% Estudo A: resposta ao impulso na ponta da viga, para xi = 0 e xi = 0.2.
% E o caso de referencia ja publicado no paper. Produz tres resultados:
%
%   1. Raio espectral de todos os metodos na viga de EF.
%   2. PE e AD (Secao 3.3) nas equacoes modais.
%   3. Deslocamento da ponta contra a resposta exata do sistema semidiscreto.
%
% Tudo fica na struct res, salva em ImgViga/viga_resultados.mat junto com as
% figuras. Os eixos sao centrados na resposta analitica, entao uma curva
% instavel sai do quadro em vez de achata-lo.

    close all;

    tInicio = now;                       %#ok<TNOW1>

    outDir = graficos().pasta_saida();
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fprintf(' Figuras (eps + png) e .mat em: %s\n\n', outDir);

    FS     = 18;
    FIGPOS = [0 0 20 14];
    YCLIP  = 1.6;                 % limite do eixo de deslocamento, normalizado

    % --- Ciclos simulados, em periodos fundamentais ----------------------
    % Os implicitos rodam em h/T_1 = 0.01: um periodo custa ~100 passos e
    % nPerImpl pode subir a vontade. Os Runge-Kutta ficam em h ~ 5e-7 s pela
    % estabilidade: CADA periodo custa ~2.2e5 passos por metodo e por caso
    % de amortecimento, entao aumentar nPerRK sai caro.
    nPerImpl = 12;                % periodos simulados, implicitos
    nPerRK   = 12;                 % periodos simulados, Runge-Kutta

    % --- Janela mostrada, em periodos fundamentais -----------------------
    % Mesma regra de tres regimes dos Estudos B-F (graficos().janela_plot).
    % n_amostras e o numero de ciclos mostrados nos regimes de amortecimento
    % fraco ou nulo. Cada figura e normalizada pelo pico da resposta exata
    % DENTRO da janela, para que uma cauda decaida nao colapse sobre o eixo.
    n_amostras = 3;
    n_amostras = max(n_amostras, eps);        % precisa ser positivo

    % Mesmo significado dos campos de mesmo nome em modelo().config_estudos.
    tolDecaimento = 1e-2;
    ciclosDecaimentoRapido = 6;

    % --- Zoom de um periodo nas quatro figuras de deslocamento ------------
    %   xi = 0   -> ultimo periodo simulado, onde o erro de fase acumulado e maior
    %   xi = 0.2 -> primeiro periodo, porque a resposta morre em ~3.8 ciclos
    % O periodo sai como inset na figura normal e como figura propria (_zoom).
    % O zoom so recorta e renormaliza series ja calculadas.

    % --- Salvamento ------------------------------------------------------
    % Os historicos sao guardados no GDL da ponta. As rodadas RK podem ter
    % milhoes de passos, entao sao decimadas para no maximo nSaveRKMax pontos
    % (Inf guarda todos os passos).
    nSaveRKMax   = 2e5;
    saveFullField = false;        % true guarda tambem todos os GDL dos implicitos

    %% ==================================================================
    %  MODELO
    %  ==================================================================
    p    = modelo().params();
    mdl  = modelo().monta(p);

    K = mdl.K;  M = mdl.M;  C = mdl.C;  Cz = mdl.Cz;
    u0 = mdl.u0;  v0 = mdl.v0;
    iTip = mdl.info.tipW;
    T1   = mdl.T1;

    xiZero = zeros(size(mdl.w));

    rkNames  = {'RK1','RK2','RK3','RK4'};
    dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'};
    dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'};
    caseTags = {'xi0','xi02'};

    [coresDir, coresRK] = graficos().cores();

    res = struct();
    res.model = mdl;

    fprintf('=====================================================\n');
    fprintf(' CANTILEVER BEAM - IMPULSE RESPONSE AT THE FREE END\n');
    fprintf('=====================================================\n');
    fprintf('  N = %d elements (%d DOF)\n', mdl.info.N, mdl.info.ndof);
    fprintf('  w_1 = %.4f rad/s,  T_1 = %.6f s\n', mdl.w(1), T1);
    fprintf('  w_max = %.4e rad/s  (w_max/w_1 = %.1f)\n', mdl.wmax, mdl.wmax/mdl.w(1));
    fprintf('  Rayleigh: a0 = %.4f, a1 = %.4e\n\n', mdl.damp.a0, mdl.damp.a1);

    %% ==================================================================
    %  1. RAIO ESPECTRAL DA VIGA DE EF
    %  ==================================================================
    % Com amortecimento proporcional, rho do operador MDOF e o maximo dos raios
    % escalares sobre os modos (ver validacao).
    hSweep = logspace(-9, -2, 200);
    res.stability.hSweep    = hSweep;
    res.stability.hSweep_T1 = hSweep/T1;

    for ic = 1:2
        tag = caseTags{ic};
        if ic == 1, ximo = xiZero;      xiL = 0.0;
        else,       ximo = mdl.damp.xi; xiL = p.xi_target;
        end
        res.stability.xi.(tag) = xiL;

        % Uma figura por caso: diretos em linha cheia, Runge-Kutta tracejados.
        fig = figure('Name',sprintf('rho %s',tag),'NumberTitle','off');
        hold on;
        leg = {};
        for fam = 1:2
            if fam == 1
                names = dirNames; cores = coresDir; labels = dirLabel; estilo = '-';
            else
                names = rkNames;  cores = coresRK;  labels = rkNames;  estilo = '--';
            end
            for k = 1:numel(names)
                r = arrayfun(@(hh) metodos().rho_beam(names{k}, mdl.w, ximo, hh), hSweep);
                res.stability.rho.(tag).(names{k}) = r;
                plot(hSweep/T1, r, estilo, 'Color', cores(k,:), 'LineWidth', 1.5);
                leg{end+1} = labels{k}; %#ok<SAGROW>
            end
        end
        set(gca,'XScale','log','YScale','log');
        yline(1,'k--','HandleVisibility','off');
        grid on;
        xlabel('$h/T_1$', 'Interpreter','latex');
        ylabel('$\rho(\mathbf{\hat{A}})$', 'Interpreter','latex');
        title(sprintf('Spectral radius, FE beam ($\\xi = %.1f$)', xiL), ...
            'Interpreter','latex');
        legend(leg, 'Location','northwest', 'FontSize', 11, 'NumColumns', 2);
        ylim([1e-2 1e2]);
        set(fig,'Units','centimeters','Position',FIGPOS);
        set(gca,'FontName','Times New Roman','FontSize',FS);
        graficos().salva_figura(fig, outDir, sprintf('viga_rho_%s',tag));
    end

    %% ==================================================================
    %  2. PASSOS CRITICOS
    %  ==================================================================
    condM = [rkNames, {'DC'}];
    hcr0 = zeros(numel(condM),1);
    hcrD = zeros(numel(condM),1);

    fprintf(' Critical time steps of the conditionally stable schemes:\n');
    fprintf('   %-5s %14s %12s %14s %12s\n','meth','h_cr (xi=0)','h_cr/T_1','h_cr (xi=0.2)','h_cr/T_1');
    for k = 1:numel(condM)
        hcr0(k) = metodos().h_critico(K, M, Cz, condM{k});
        hcrD(k) = metodos().h_critico(K, M, C , condM{k});
        fprintf('   %-5s %14.4e %12.3e %14.4e %12.3e\n', condM{k}, ...
            hcr0(k), hcr0(k)/T1, hcrD(k), hcrD(k)/T1);
    end
    fprintf('   (undamped RK1 and RK2 values are bisection-tolerance artifacts:\n');
    fprintf('    rho > 1 for every h > 0, as reported in Section 4.1)\n\n');

    res.stability.condMethods = condM;
    res.stability.hcr_xi0     = hcr0;
    res.stability.hcr_xi02    = hcrD;
    res.stability.hcr_table   = cell2table([condM(:), num2cell(hcr0), ...
        num2cell(hcr0/T1), num2cell(hcrD), num2cell(hcrD/T1)], ...
        'VariableNames', {'method','hcr_xi0','hcr_xi0_T1','hcr_xi02','hcr_xi02_T1'});

    % Passos de tempo usados daqui em diante
    hRat  = 0.01;                    % implicitos, como fracao de T1
    hImpl = hRat*T1;
    hRK   = 0.90*min(hcrD(1:4));     % explicitos, abaixo do limite amortecido

    %% ==================================================================
    %  3. ALONGAMENTO DE PERIODO E DECAIMENTO DE AMPLITUDE
    %  ==================================================================
    % Definicoes da Secao 3.3 aplicadas aos modos da viga. Sem amortecimento,
    % PE e AD em porcentagem; com amortecimento, razoes em que 1 reproduz
    % exatamente o periodo e o decaimento fisicos.
    %
    % Cada familia e avaliada no passo em que roda: implicitos em h/T_1 = 0.01
    % e Runge-Kutta em h = 0.90*h_cr. Avaliar RK no passo dos implicitos so
    % daria NaN, pois o operador fica fora da faixa oscilatoria.
    %
    % A tabela impressa mostra os tres primeiros modos; a struct guarda todos.
    allNames  = [dirNames, rkNames];
    allLabels = [dirLabel, rkNames];
    allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
    nModes    = numel(mdl.w);
    modesRep  = 1:3;

    res.pead.methods = allNames;
    res.pead.labels  = allLabels;
    res.pead.h       = allSteps;
    res.pead.modes   = (1:nModes).';
    res.pead.w       = mdl.w;

    for ic = 1:2
        tag = caseTags{ic};
        if ic == 1, ximo = xiZero;      xiL = 0.0;         unit = '[%]';
        else,       ximo = mdl.damp.xi; xiL = p.xi_target; unit = '[ratio]';
        end
        res.pead.xi.(tag)       = xiL;
        res.pead.xi_modal.(tag) = ximo;

        PEall  = nan(numel(allNames), nModes);
        ADall  = nan(numel(allNames), nModes);
        RHOall = nan(numel(allNames), nModes);

        fprintf(' PE and AD, xi = %.1f   %s\n', xiL, unit);
        fprintf('   implicit schemes at h/T_1 = %.3f,  RK schemes at h = %.3e s\n', ...
            hRat, hRK);
        fprintf('   %-20s %10s', 'method', 'h/T_1');
        for m = modesRep, fprintf('%10s%10s', sprintf('PE m%d',m), sprintf('AD m%d',m)); end
        fprintf('\n');

        for k = 1:numel(allNames)
            hk = allSteps(k);
            for m = 1:nModes
                [PEall(k,m), ADall(k,m), RHOall(k,m)] = ...
                    metodos().pe_ad(allNames{k}, mdl.w(m), ximo(m), hk);
            end
            fprintf('   %-20s %10.2e', allLabels{k}, hk/T1);
            for m = modesRep
                fprintf('%10.4f%10.4f', PEall(k,m), ADall(k,m));
            end
            fprintf('\n');
        end
        fprintf('\n');

        res.pead.PE.(tag)  = PEall;
        res.pead.AD.(tag)  = ADall;
        res.pead.rho.(tag) = RHOall;

        % Tabela dos tres primeiros modos
        res.pead.table.(tag) = cell2table( ...
            [allLabels(:), num2cell(allSteps(:)/T1), ...
             num2cell(PEall(:,1)), num2cell(ADall(:,1)), ...
             num2cell(PEall(:,2)), num2cell(ADall(:,2)), ...
             num2cell(PEall(:,3)), num2cell(ADall(:,3))], ...
            'VariableNames', {'method','h_T1','PE_m1','AD_m1','PE_m2','AD_m2','PE_m3','AD_m3'});
    end

    %% ==================================================================
    %  4. DESLOCAMENTO DA PONTA - METODOS DIRETOS
    %  ==================================================================
    nS = round(nPerImpl*T1/hImpl);

    for ic = 1:2
        tag = caseTags{ic};
        if ic == 1, Cm = Cz;  ximo = xiZero;      xiL = 0.0;         hcrDC = hcr0(end);
        else,       Cm = C;   ximo = mdl.damp.xi; xiL = p.xi_target; hcrDC = hcrD(end);
        end
        amortecido = xiL > 0;

        % Periodo do zoom (ver cabecalho): ultimo periodo sem amortecimento,
        % primeiro periodo com amortecimento.
        if amortecido
            zTIni = 0;                        zTFim = min(1,nPerImpl)*T1;
        else
            zTIni = max(nPerImpl-1,0)*T1;     zTFim = nPerImpl*T1;
        end

        tRef = (0:nS)*hImpl;
        [Uex, Vex, Aex] = modelo().solucao_exata(tRef, M, mdl.w, mdl.Phi, ximo, u0, v0);

        res.time.(tag).exact.t     = tRef;
        res.time.(tag).exact.u_tip = Uex(iTip,:);
        res.time.(tag).exact.v_tip = Vex(iTip,:);
        res.time.(tag).exact.a_tip = Aex(iTip,:);
        if saveFullField
            res.time.(tag).exact.U = Uex;
        end

        % Janela mostrada (so afeta o desenho, nao a simulacao).
        [tIni, tFim, infoJan] = graficos().janela_plot(tRef, Uex(iTip,:), T1, ...
            n_amostras, tolDecaimento, amortecido, ciclosDecaimentoRapido);
        res.window.(tag).direto = infoJan;
        graficos().print_janela('direct schemes', infoJan, T1, nPerImpl);

        mask   = (tRef >= tIni - 0.5*hImpl) & (tRef <= tFim + 0.5*hImpl);
        uScale = max(abs(Uex(iTip,mask)));
        if uScale == 0, uScale = 1; end

        % DC fica fora da figura so quando esta instavel neste passo.
        skipDC = (hcrDC > 0) && (hImpl > hcrDC);

        fig = figure('Name',sprintf('Tip disp direct %s',tag),'NumberTitle','off');
        plot(tRef(mask)/T1, Uex(iTip,mask)/uScale, 'k--', 'LineWidth', 2); hold on;
        leg = {'Exact'};

        for k = 1:numel(dirNames)
            [t,U,V,A] = metodos().integra_metodo(dirNames{k}, K, M, Cm, u0, v0, hImpl, nS);

            res.time.(tag).(dirNames{k}).t     = t;
            res.time.(tag).(dirNames{k}).h     = hImpl;
            res.time.(tag).(dirNames{k}).u_tip = U(iTip,:);
            res.time.(tag).(dirNames{k}).v_tip = V(iTip,:);
            res.time.(tag).(dirNames{k}).a_tip = A(iTip,:);
            res.time.(tag).(dirNames{k}).err   = ...
                max(abs(U(iTip,:) - Uex(iTip,:)))/max(abs(Uex(iTip,:)));
            if saveFullField
                res.time.(tag).(dirNames{k}).U = U;
            end

            if strcmpi(dirNames{k},'DC') && skipDC, continue; end

            plot(t(mask)/T1, U(iTip,mask)/uScale, '-', ...
                'Color', coresDir(k,:), 'LineWidth', 1.3);
            leg{end+1} = dirLabel{k}; %#ok<SAGROW>
        end

        grid on;
        xlabel('$t/T_1$', 'Interpreter','latex');
        ylabel('$w_{tip}/w_{ref}$', 'Interpreter','latex');
        title(sprintf('Tip displacement, implicit schemes ($h/T_1 = %.2f$, $\\xi = %.1f$)', ...
            hRat, xiL), 'Interpreter','latex');
        legend(leg, 'Location','northeast', 'FontSize', 12);
        xlim([tIni/T1 tFim/T1]);  ylim([-YCLIP YCLIP]);
        set(fig,'Units','centimeters','Position',FIGPOS);
        set(gca,'FontName','Times New Roman','FontSize',FS);
        axMain = gca;   % antes de o inset criar um eixo novo

        % --- Curvas que entram no zoom (as mesmas ja desenhadas acima) -----
        curvasZ = struct('t',{},'u',{},'cor',{},'nome',{});
        for k = 1:numel(dirNames)
            if strcmpi(dirNames{k},'DC') && skipDC, continue; end
            curvasZ(end+1) = struct('t', res.time.(tag).(dirNames{k}).t, ...
                'u', res.time.(tag).(dirNames{k}).u_tip, ...
                'cor', coresDir(k,:), 'nome', dirLabel{k}); %#ok<SAGROW>
        end

        % 1) Inset ampliado na figura de janela normal, antes de salva-la.
        add_zoom_periodo(axMain, tRef, Uex(iTip,:), curvasZ, zTIni, zTFim, T1, FS, YCLIP);
        graficos().salva_figura(fig, outDir, sprintf('viga_desloc_direto_%s',tag));

        % 2) Figura propria, com esse periodo como janela unica.
        figZ = figura_periodo(tRef, Uex(iTip,:), curvasZ, zTIni, zTFim, T1, ...
            sprintf('Tip displacement, implicit schemes, period %d ($h/T_1 = %.2f$, $\\xi = %.1f$)', ...
                round(zTFim/T1), hRat, xiL), ...
            sprintf('Tip disp direct %s zoom',tag), YCLIP, FIGPOS, FS);
        graficos().salva_figura(figZ, outDir, sprintf('viga_desloc_direto_%s_zoom',tag));
    end

    %% ==================================================================
    %  5. DESLOCAMENTO DA PONTA - RUNGE-KUTTA
    %  ==================================================================
    % Um unico passo para os quatro metodos: 90% do menor passo critico com
    % amortecimento, para que a comparacao isole o efeito do amortecimento
    % fisico. Com xi = 0 esse passo fica acima do limite (inexistente) de RK1
    % e RK2, que saem do quadro; com xi = 0.2 os quatro sao estaveis. E a
    % versao no tempo do resultado central da Secao 4.1.
    nRK = round(nPerRK*T1/hRK);

    fprintf(' Runge-Kutta simulations: h = %.4e s (%.3e T_1), %d steps\n', ...
        hRK, hRK/T1, nRK);
    fprintf('   this is the longest part of the run; please wait.\n\n');

    for ic = 1:2
        tag = caseTags{ic};
        if ic == 1, Cm = Cz;  ximo = xiZero;      xiL = 0.0;
        else,       Cm = C;   ximo = mdl.damp.xi; xiL = p.xi_target;
        end
        amortecido = xiL > 0;

        % Periodo do zoom (ver cabecalho), contado em periodos de T1 da
        % simulacao RK - que tem duracao propria (nPerRK), nao a dos diretos.
        if amortecido
            zTIni = 0;                      zTFim = min(1,nPerRK)*T1;
        else
            zTIni = max(nPerRK-1,0)*T1;     zTFim = nPerRK*T1;
        end
        % Referencia exata PROPRIA do zoom: tFine abaixo so cobre a janela
        % normal, e o periodo do zoom pode cair fora dela (depende de
        % n_amostras/tolDecaimento). Forma fechada, barata.
        tZoom   = linspace(zTIni, zTFim, 2000);
        UexZ    = modelo().solucao_exata(tZoom, M, mdl.w, mdl.Phi, ximo, u0, v0);
        uexZoom = UexZ(iTip,:);

        % Varredura grosseira na simulacao inteira (forma fechada) para achar a janela;
        % a grade fina de 4000 pontos so e montada depois, sobre ela.
        tVarre = linspace(0, nPerRK*T1, 2000);
        UexV   = modelo().solucao_exata(tVarre, M, mdl.w, mdl.Phi, ximo, u0, v0);
        [tIniRK, tFimRK, infoJanRK] = graficos().janela_plot(tVarre, UexV(iTip,:), T1, ...
            n_amostras, tolDecaimento, amortecido, ciclosDecaimentoRapido);
        res.window.(tag).rk = infoJanRK;
        graficos().print_janela('Runge-Kutta', infoJanRK, T1, nPerRK);

        tFine  = linspace(tIniRK, tFimRK, 4000);
        UexF   = modelo().solucao_exata(tFine, M, mdl.w, mdl.Phi, ximo, u0, v0);
        uScale = max(abs(UexF(iTip,:)));
        if uScale == 0, uScale = 1; end

        res.timeRK.(tag).exact.t     = tFine;
        res.timeRK.(tag).exact.u_tip = UexF(iTip,:);

        fig = figure('Name',sprintf('Tip disp RK %s',tag),'NumberTitle','off');
        plot(tFine/T1, UexF(iTip,:)/uScale, 'k--', 'LineWidth', 2); hold on;
        leg = {'Exact'};

        for k = 1:numel(rkNames)
            [t,U,V,A] = metodos().mdof_rk(K, M, Cm, u0, v0, hRK, nRK, k);

            % Historicos decimados para manter o .mat em tamanho razoavel
            [ts, us] = graficos().downsample_plot(t, U(iTip,:), nSaveRKMax);
            [~,  vs] = graficos().downsample_plot(t, V(iTip,:), nSaveRKMax);
            [~,  as] = graficos().downsample_plot(t, A(iTip,:), nSaveRKMax);
            res.timeRK.(tag).(rkNames{k}).t          = ts;
            res.timeRK.(tag).(rkNames{k}).h          = hRK;
            res.timeRK.(tag).(rkNames{k}).nSteps     = nRK;
            res.timeRK.(tag).(rkNames{k}).u_tip      = us;
            res.timeRK.(tag).(rkNames{k}).v_tip      = vs;
            res.timeRK.(tag).(rkNames{k}).a_tip      = as;
            res.timeRK.(tag).(rkNames{k}).decimated  = numel(ts) < numel(t);
            res.timeRK.(tag).(rkNames{k}).u_tip_max  = max(abs(U(iTip,:)));

            mask    = (t >= tIniRK - 0.5*hRK) & (t <= tFimRK + 0.5*hRK);
            [tp,yp] = graficos().downsample_plot(t(mask), U(iTip,mask)/uScale);
            plot(tp/T1, yp, '-', 'Color', coresRK(k,:), 'LineWidth', 1.3);
            leg{end+1} = rkNames{k}; %#ok<SAGROW>
        end

        grid on;
        xlabel('$t/T_1$', 'Interpreter','latex');
        ylabel('$w_{tip}/w_{ref}$', 'Interpreter','latex');
        title(sprintf('Tip displacement, Runge-Kutta ($h = %.2e$ s, $\\xi = %.1f$)', ...
            hRK, xiL), 'Interpreter','latex');
        legend(leg, 'Location','northeast', 'FontSize', 12);
        xlim([tIniRK/T1 tFimRK/T1]);  ylim([-YCLIP YCLIP]);
        set(fig,'Units','centimeters','Position',FIGPOS);
        set(gca,'FontName','Times New Roman','FontSize',FS);
        axMain = gca;   % antes de o inset criar um eixo novo

        % --- Curvas que entram no zoom (as series decimadas guardadas acima
        %     cobrem a simulacao INTEIRA, entao servem para qualquer periodo)
        curvasZ = struct('t',{},'u',{},'cor',{},'nome',{});
        for k = 1:numel(rkNames)
            curvasZ(end+1) = struct('t', res.timeRK.(tag).(rkNames{k}).t, ...
                'u', res.timeRK.(tag).(rkNames{k}).u_tip, ...
                'cor', coresRK(k,:), 'nome', rkNames{k}); %#ok<SAGROW>
        end

        % 1) Inset ampliado na figura de janela normal, antes de salva-la.
        add_zoom_periodo(axMain, tZoom, uexZoom, curvasZ, zTIni, zTFim, T1, FS, YCLIP);
        graficos().salva_figura(fig, outDir, sprintf('viga_desloc_rk_%s',tag));

        % 2) Figura propria, com esse periodo como janela unica.
        figZ = figura_periodo(tZoom, uexZoom, curvasZ, zTIni, zTFim, T1, ...
            sprintf('Tip displacement, Runge-Kutta, period %d ($h = %.2e$ s, $\\xi = %.1f$)', ...
                round(zTFim/T1), hRK, xiL), ...
            sprintf('Tip disp RK %s zoom',tag), YCLIP, FIGPOS, FS);
        graficos().salva_figura(figZ, outDir, sprintf('viga_desloc_rk_%s_zoom',tag));
    end

    %% ==================================================================
    %  SAVE
    %  ==================================================================
    % Estrutura completa de res: ver a secao "Dados salvos" do README.
    res.config = struct( ...
        'nPerImpl',     nPerImpl, ...
        'nPerRK',       nPerRK, ...
        'n_amostras',   n_amostras, ...
        'tolDecaimento', tolDecaimento, ...
        'ciclosDecaimentoRapido', ciclosDecaimentoRapido, ...
        'hRat',         hRat, ...
        'hImpl',        hImpl, ...
        'hRK',          hRK, ...
        'nStepsImpl',   nS, ...
        'nStepsRK',     nRK, ...
        'YCLIP',        YCLIP, ...
        'nSaveRKMax',   nSaveRKMax, ...
        'saveFullField',saveFullField, ...
        'tipDOF',       iTip, ...
        'T1',           T1, ...
        'runDate',      datestr(now));

    matFile = fullfile(outDir,'viga_resultados.mat');
    save(matFile, 'res', '-v7.3');

    d = dir(matFile);
    fprintf(' Figures saved to %s\n', outDir);
    fprintf(' All results saved to %s (%.1f MB)\n', matFile, d.bytes/1e6);
    fprintf('   load it with:  S = load(''%s'');  res = S.res;\n', 'viga_resultados.mat');
    graficos().lista_obsoletos({outDir}, tInicio);
    fprintf('=====================================================\n');
end

function uScaleZ = escala_periodo(tExact, uExact, zTIni, zTFim)
% Pico da resposta exata no periodo. A normalizacao do zoom e local: com
% xi = 0.2 a escala da figura grande deixaria o zoom achatado.
    mZ = (tExact >= zTIni) & (tExact <= zTFim);
    uScaleZ = max(abs(uExact(mZ)));
    if isempty(uScaleZ) || uScaleZ == 0, uScaleZ = 1; end
end

function add_zoom_periodo(axMain, tExact, uExact, curvas, zTIni, zTFim, T1, FS, YCLIP)
% Inset com um periodo dentro da figura normal, para diretos e Runge-Kutta.
%   axMain        eixo da figura principal
%   tExact,uExact referencia exata cobrindo ao menos [zTIni,zTFim]
%   curvas        struct array com .t, .u, .cor de cada metodo
%   zTIni,zTFim   limites do periodo [s]
%   T1, FS, YCLIP periodo de referencia, fonte e corte do eixo y
% YCLIP tambem vale aqui: sem ele, RK1/RK2 instaveis achatariam as outras curvas.
    uScaleZ = escala_periodo(tExact, uExact, zTIni, zTFim);

    pos = get(axMain,'Position');
    axZoom = axes('Parent', ancestor(axMain,'figure'), ...
        'Position',[pos(1)+0.55*pos(3), pos(2)+0.12*pos(4), ...
                    0.40*pos(3), 0.38*pos(4)]);
    hold(axZoom,'on'); box(axZoom,'on');

    mZ = (tExact >= zTIni) & (tExact <= zTFim);
    plot(axZoom, tExact(mZ)/T1, uExact(mZ)/uScaleZ, 'k--', 'LineWidth', 1.5);
    for k = 1:numel(curvas)
        mk = (curvas(k).t >= zTIni) & (curvas(k).t <= zTFim);
        [tp,yp] = graficos().downsample_plot(curvas(k).t(mk), curvas(k).u(mk)/uScaleZ);
        plot(axZoom, tp/T1, yp, '-', 'Color', curvas(k).cor, 'LineWidth', 1.1);
    end
    xlim(axZoom, [zTIni zTFim]/T1);  ylim(axZoom, [-YCLIP YCLIP]);
    grid(axZoom,'on');
    set(axZoom,'FontSize',max(FS-6,8));
    title(axZoom, sprintf('zoom: $t/T_1 \\in [%.0f,%.0f]$', zTIni/T1, zTFim/T1), ...
        'Interpreter','latex','FontSize',max(FS-6,8));
end

function figZ = figura_periodo(tExact, uExact, curvas, zTIni, zTFim, T1, ...
                               titulo, nomeFig, YCLIP, FIGPOS, FS)
% Figura propria (_zoom) em que a janela principal e so o periodo escolhido.
    uScaleZ = escala_periodo(tExact, uExact, zTIni, zTFim);

    figZ = figure('Name',nomeFig,'NumberTitle','off');
    mZ = (tExact >= zTIni) & (tExact <= zTFim);
    plot(tExact(mZ)/T1, uExact(mZ)/uScaleZ, 'k--', 'LineWidth', 2); hold on;
    leg = {'Exact'};
    for k = 1:numel(curvas)
        mk = (curvas(k).t >= zTIni) & (curvas(k).t <= zTFim);
        [tp,yp] = graficos().downsample_plot(curvas(k).t(mk), curvas(k).u(mk)/uScaleZ);
        plot(tp/T1, yp, '-', 'Color', curvas(k).cor, 'LineWidth', 1.3);
        leg{end+1} = curvas(k).nome; %#ok<AGROW>
    end

    grid on;
    xlabel('$t/T_1$', 'Interpreter','latex');
    ylabel('$w_{tip}/w_{ref}$', 'Interpreter','latex');
    title(titulo, 'Interpreter','latex');
    legend(leg, 'Location','northeast', 'FontSize', 12);
    xlim([zTIni zTFim]/T1);  ylim([-YCLIP YCLIP]);
    set(figZ,'Units','centimeters','Position',FIGPOS);
    set(gca,'FontName','Times New Roman','FontSize',FS);
end

%% ====================================================================
%  Estudos B a F
%  ====================================================================

function estudos_BF()
% Roda os Estudos B a F em sequencia. O Estudo A nao entra aqui porque e o
% caso de referencia e roda a parte, com viga('A').
    t0 = tic;
    estudos = {@estudo_B, @estudo_E, @estudo_C, @estudo_D, @estudo_F};
    for k = 1:numel(estudos)
        nome = func2str(estudos{k});
        fprintf('\n>>> Executando %s (%d/%d) ...\n\n', nome, k, numel(estudos));
        tk = tic;
        estudos{k}();
        fprintf('\n>>> %s concluido em %.1f s.\n', nome, toc(tk));
        close all;
    end
    fprintf('\n  ESTUDOS B-F CONCLUIDOS em %.1f min.\n', toc(t0)/60);
    fprintf('  Resultados em %s\n', graficos().pasta_saida());
end

function estudo_B()
% Estudo B: efeito do amortecimento fisico. Repete a resposta ao impulso do
% Estudo A para os valores de cfg.B.xi_list. Estruturas metalicas tem xi tipico
% entre 0.005 e 0.02, faixa em que o amortecimento numerico fica visivel.
%
% Um unico passo hRK, de 90% do h_cr do caso MAIS amortecido, vale para todos
% os xi. E isso que mostra RK1/RK2 divergindo com pouco amortecimento e
% estaveis com muito, usando o mesmo passo (Secao 4.1 do paper).
% Saida em ImgViga/EstudoB/.

    close all;

    tInicio = now;                       %#ok<TNOW1>

    outDir = graficos().pasta_saida('EstudoB');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

    cfg = modelo().config_estudos();
    FS = cfg.FS; FIGPOS = cfg.FIGPOS;
    hRat = cfg.hRat;

    p = modelo().params();
    [K, M, info] = modelo().matrizes_EF(p);
    [w, Phi]     = modelo().autovalores_EF(K, M);
    T1   = 2*pi/w(1);
    iTip = info.tipW;

    [dirNames, dirLabel, rkNames] = metodos().nomes_metodos();
    condM    = [rkNames, {'DC'}];

    xiList = cfg.B.xi_list;
    [xiRef, iRef] = max(xiList);   %#ok<ASGLU>  % o caso mais amortecido define o passo RK comum

    fprintf('=====================================================\n');
    fprintf(' ESTUDO B - PHYSICAL DAMPING SWEEP  (xi = %s)\n', mat2str(xiList));
    fprintf('=====================================================\n');

    % --- Matrizes de amortecimento, xi modal e h_cr de cada caso ---------
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
            [Cs{ix}, dtmp] = modelo().amortecimento(K, M, w, pp);
            ximoS{ix} = dtmp.xi;
            dampS{ix} = dtmp;
        end
        for k = 1:numel(condM)
            hcrS(k,ix) = metodos().h_critico(K, M, Cs{ix}, condM{k});
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
        [u0, v0] = modelo().ci_impulso(K, M, C, info, p);

        % --- Raio espectral ----------------------------------------------
        stab = graficos().espectro(w, ximo, sprintf('\\xi = %.3g', xiL), ['B_' tag], outDir, ...
            'hSweep', cfg.hSweep, 'FS', FS, 'FIGPOS', FIGPOS, ...
            'xScaleBy', 1/T1, 'xLabelStr', '$h/T_1$');
        res.stability.(tag) = stab;

        % --- PE/AD, nos mesmos passos usados na resposta no tempo --------
        allNames  = [dirNames, rkNames];
        allLabels = [dirLabel, rkNames];
        allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
        pead = metodos().pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
            'nModesReport', cfg.nModesReport, ...
            'printHeader', sprintf('Estudo B - PE/AD, xi = %.3g', xiL));
        res.pead.(tag) = pead;

        % --- Resposta no tempo, com zoom no caso de baixo amortecimento --
        cfgT = graficos().caso_viga_completa(K, M, C, u0, v0, w, Phi, ximo, iTip);
        cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = hRK;
        cfgT.label = ['B_' tag]; cfgT.xiLabel = sprintf('\\xi = %.3g', xiL);
        cfgT.outDir = outDir;
        cfgT = graficos().aplica_config(cfgT, cfg);
        % itens 20-21: zoom no caso escolhido em cfg.B.xiZoom
        cfgT.doZoom = ~isempty(cfg.B.xiZoom) && (xiL == cfg.B.xiZoom);

        tresp = graficos().resposta_tempo(cfgT);
        res.time.(tag) = tresp;
    end

    res.config = struct('hRat',hRat,'hImpl',hImpl,'hRK',hRK,'xiList',xiList, ...
        'runDate',datestr(now));
    res.cfg = cfg;   % guarda a configuracao usada nesta rodada

    save(fullfile(outDir,'estudoB_resultados.mat'), 'res', '-v7.3');
    fprintf('\nEstudo B concluido. Resultados e figuras em %s\n', outDir);
    graficos().lista_obsoletos({outDir}, tInicio);
end

function estudo_C()
% Estudo C: condicao inicial modal com passo fixo. Desloca a viga na forma de
% um modo por vez (cfg.C.modos) mantendo h = 0.01*T1 em termos ABSOLUTOS, para
% ver o efeito de h/T_modo crescer conforme o modo sobe. Com xi = 0.01 a
% dissipacao numerica aparece, em vez de ficar mascarada pela fisica.
%
% O limite do RK vem de w_max da malha, nao do modo excitado, entao hRK e
% calculado uma vez e vale para todos os casos.
% Saida em ImgViga/EstudoC/.

    close all;

    tInicio = now;                       %#ok<TNOW1>

    outDir = graficos().pasta_saida('EstudoC');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

    cfg = modelo().config_estudos();
    FS = cfg.FS; FIGPOS = cfg.FIGPOS;
    hRat = cfg.hRat;
    xiL  = cfg.C.xi;
    modos = cfg.C.modos;

    p = modelo().params();
    [K, M, info] = modelo().matrizes_EF(p);
    [w, Phi]     = modelo().autovalores_EF(K, M);
    T1   = 2*pi/w(1);
    iTip = info.tipW;

    pp = p; pp.xi_target = xiL;
    [C, damp] = modelo().amortecimento(K, M, w, pp);
    ximo = damp.xi;

    [dirNames, dirLabel, rkNames] = metodos().nomes_metodos();
    condM    = [rkNames, {'DC'}];

    hImpl = hRat*T1;
    hcr   = zeros(numel(condM),1);
    for k = 1:numel(condM)
        hcr(k) = metodos().h_critico(K, M, C, condM{k});
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

        [u0, v0] = modelo().ci_modal(K, M, C, Phi, modeIdx);

        % --- PE/AD no modo excitado e nos vizinhos -------------------------
        allNames  = [dirNames, rkNames];
        allLabels = [dirLabel, rkNames];
        allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
        % Relata pelo menos cfg.nModesReport modos, mas sempre o suficiente
        % para que o modo EXCITADO apareca na tabela.
        pead = metodos().pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
            'nModesReport', max(modeIdx, cfg.nModesReport), ...
            'printHeader', sprintf('Estudo C - PE/AD, modo %d excitado', modeIdx));
        res.pead.(tag) = pead;
        res.Tmodo.(tag) = Tmodo;
        res.hT_modo.(tag) = hImpl/Tmodo;

        % --- Resposta no tempo -------------------------------------------
        cfgT = graficos().caso_viga_completa(K, M, C, u0, v0, w, Phi, ximo, iTip);
        cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = hRK;
        cfgT.label = ['C_' tag];
        cfgT.xiLabel = sprintf('mode %d, \\xi = %.3f', modeIdx, xiL);
        cfgT.outDir = outDir;
        cfgT = graficos().aplica_config(cfgT, cfg);

        tresp = graficos().resposta_tempo(cfgT);
        res.time.(tag) = tresp;
    end

    res.config = struct('hRat',hRat,'hImpl',hImpl,'hRK',hRK,'xi',xiL, ...
        'modes',modos,'runDate',datestr(now));
    res.cfg = cfg;   % guarda a configuracao usada nesta rodada

    save(fullfile(outDir,'estudoC_resultados.mat'), 'res', '-v7.3');
    fprintf('\nEstudo C concluido. Resultados e figuras em %s\n', outDir);
    graficos().lista_obsoletos({outDir}, tInicio);
end

function estudo_D()
% Estudo D: condicao inicial modal com passo varrido. Para alguns modos
% (cfg.D.modos), varre h/T_modo (cfg.D.hRatios) e observa PE, AD e estabilidade
% conforme h cresce em relacao ao periodo do proprio modo.
%
% PE, AD e rho sao forma fechada e saem para toda combinacao. A simulacao no
% tempo fica nas razoes de cfg.D.ratiosSimuladas, com poucos ciclos, e cada
% rodada fica em try/catch: divergir aqui e resultado esperado, e um caso que
% falha e registrado sem parar o estudo.
% Saida em ImgViga/EstudoD/.

    close all;

    tInicio = now;                       %#ok<TNOW1>

    outDir = graficos().pasta_saida('EstudoD');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

    cfg = modelo().config_estudos();
    FS = cfg.FS; FIGPOS = cfg.FIGPOS;
    nPerD = cfg.D.nPer; ciclosPlotD = cfg.D.ciclosPlot;
    xiL  = cfg.D.xi;

    p = modelo().params();
    [K, M, info] = modelo().matrizes_EF(p);
    [w, Phi]     = modelo().autovalores_EF(K, M);
    T1   = 2*pi/w(1);
    iTip = info.tipW;

    pp = p; pp.xi_target = xiL;
    [C, damp] = modelo().amortecimento(K, M, w, pp);
    ximo = damp.xi;

    [dirNames, dirLabel, rkNames] = metodos().nomes_metodos();
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

    % --- PE/AD/rho em forma fechada para TODA combinacao (modo, razao) ---
    peadGrid = struct();
    for im = 1:numel(modes)
        modeIdx = modes(im);
        Tmodo = 2*pi/w(modeIdx);
        tagM = sprintf('modo%d', modeIdx);
        for ir = 1:numel(ratios)
            h = ratios(ir)*Tmodo;
            tagR = matlab.lang.makeValidName(sprintf('r_%.3g', ratios(ir)));
            allSteps = repmat(h, 1, numel(allNames));
            pead = metodos().pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
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
        [u0, v0] = modelo().ci_modal(K, M, C, Phi, modeIdx);

        for ir = ratiosPlot
            h = ratios(ir)*Tmodo;
            tagR = matlab.lang.makeValidName(sprintf('r_%.3g', ratios(ir)));
            fprintf('--- modo %d, h/T_modo = %.3g (h = %.3e s) ---\n', modeIdx, ratios(ir), h);

            cfgT = graficos().caso_viga_completa(K, M, C, u0, v0, w, Phi, ximo, iTip);
            cfgT.T1 = T1; cfgT.hImpl = h; cfgT.hRK = h;
            cfgT.label = sprintf('D_%s_%s', tagM, tagR);
            cfgT.xiLabel = sprintf('mode %d, h/T_{mode} = %.3g', modeIdx, ratios(ir));
            cfgT.outDir = outDir;
            % Estudo D usa contagens proprias de ciclos (cfg.D), bem menores
            % que as globais - ver o cabecalho e a secao 5 da configuracao.
            cfgT = graficos().aplica_config(cfgT, cfg, ...
                'nPerImpl', nPerD, 'nPerRK', nPerD, ...
                'nPerImpl_sem_xi', nPerD, 'nPerRK_sem_xi', nPerD, ...
                'X_ciclos_final_plot', ciclosPlotD, ...
                'X_ciclos_final_plot_sem_xi', ciclosPlotD);
    
            try
                tresp = graficos().resposta_tempo(cfgT);
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
    res.cfg = cfg;   % guarda a configuracao usada nesta rodada

    save(fullfile(outDir,'estudoD_resultados.mat'), 'res', '-v7.3');
    fprintf('\nEstudo D concluido. Resultados e figuras em %s\n', outDir);
    graficos().lista_obsoletos({outDir}, tInicio);
end

function estudo_E()
% Estudo E: Rayleigh sem xi alvo. Impoe a0 = 0 e a1 = 1e-4 ou 1e-5 direto
% (cfg.E.rayleigh). Sem o termo de massa, xi_k = a1*w_k/2 cresce linearmente
% com a frequencia e nao ha modos ancora; o xi de cada modo e calculado.
% Usa a mesma condicao de impulso do Estudo A, para comparacao direta.
% Saida em ImgViga/EstudoE/.

    close all;

    tInicio = now;                       %#ok<TNOW1>

    outDir = graficos().pasta_saida('EstudoE');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

    cfg = modelo().config_estudos();
    FS = cfg.FS; FIGPOS = cfg.FIGPOS;
    hRat = cfg.hRat;
    rayleigh = cfg.E.rayleigh;

    p = modelo().params();
    [K, M, info] = modelo().matrizes_EF(p);
    [w, Phi]     = modelo().autovalores_EF(K, M);
    T1   = 2*pi/w(1);
    iTip = info.tipW;

    [dirNames, dirLabel, rkNames] = metodos().nomes_metodos();
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

        [C, damp] = modelo().amortecimento(K, M, w, p, a0, a1);
        ximo = damp.xi;

        fprintf('   xi_1..xi_5 = %s\n', mat2str(round(ximo(1:5),4)));
        res.damp.(tag) = damp;

        [u0, v0] = modelo().ci_impulso(K, M, C, info, p);

        % --- Raio espectral -----------------------------------------------
        xiLabel = sprintf('a_0=0, a_1=%.0e', a1);
        stab = graficos().espectro(w, ximo, xiLabel, ['E_' tag], outDir, ...
            'hSweep', cfg.hSweep, 'FS', FS, 'FIGPOS', FIGPOS, ...
            'xScaleBy', 1/T1, 'xLabelStr', '$h/T_1$');
        res.stability.(tag) = stab;

        % --- Passos criticos --------------------------------------------------
        hcr = zeros(numel(condM),1);
        for k = 1:numel(condM)
            hcr(k) = metodos().h_critico(K, M, C, condM{k});
        end
        res.hcr.(tag) = cell2table([condM(:), num2cell(hcr), num2cell(hcr/T1)], ...
            'VariableNames', {'method','hcr','hcr_T1'});

        hImpl = hRat*T1;
        hRK   = cfg.fatorSegurancaRK*min(hcr(1:4));

        % --- PE/AD ----------------------------------------------------------
        allNames  = [dirNames, rkNames];
        allLabels = [dirLabel, rkNames];
        allSteps  = [repmat(hImpl,1,numel(dirNames)), repmat(hRK,1,numel(rkNames))];
        pead = metodos().pead_tabela(w, ximo, allNames, allLabels, allSteps, T1, ...
            'nModesReport', cfg.nModesReport, ...
            'printHeader', sprintf('Estudo E - PE/AD, a0=0, a1=%.0e', a1));
        res.pead.(tag) = pead;

        % --- Resposta no tempo ------------------------------------------------
        cfgT = graficos().caso_viga_completa(K, M, C, u0, v0, w, Phi, ximo, iTip);
        cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = max(hRK, eps);
        cfgT.label = ['E_' tag]; cfgT.xiLabel = xiLabel;
        cfgT.outDir = outDir;
        cfgT = graficos().aplica_config(cfgT, cfg);
        cfgT.runRK = hRK > 0;

        tresp = graficos().resposta_tempo(cfgT);
        res.time.(tag) = tresp;
    end

    res.config = struct('hRat',hRat,'rayleigh_alt',rayleigh,'runDate',datestr(now));
    res.cfg = cfg;   % guarda a configuracao usada nesta rodada

    save(fullfile(outDir,'estudoE_resultados.mat'), 'res', '-v7.3');
    fprintf('\nEstudo E concluido. Resultados e figuras em %s\n', outDir);
    graficos().lista_obsoletos({outDir}, tInicio);
end

function estudo_F()
% Estudo F: reducao modal, com xi = 0.2 (mesmo caso fisico do Estudo A).
%
% PARTE 1: raio espectral, PE/AD e resposta no tempo para o modelo completo
% e para os truncamentos de cfg.F.reduction_modes. Cada modelo usa o proprio
% passo critico. Responde se o sistema reduzido representa a resposta relevante.
%
% PARTE 2 (item 32): separa o efeito de um h maior do efeito da reducao, com
% RK4 no modelo de 4 modos:
%   Caso 1: modelo completo, h = 0.9*h_cr do completo
%   Caso 2: modelo reduzido, o MESMO h do caso 1
%   Caso 3: modelo reduzido, h = 0.9*h_cr do reduzido (maior)
% Os tres sao comparados com a resposta exata do modelo COMPLETO.
% Saida em ImgViga/EstudoF/.

    close all;

    tInicio = now;                       %#ok<TNOW1>

    outDir = graficos().pasta_saida('EstudoF');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fprintf(' Figuras (eps + png) e .mat em: %s\n', outDir);

    cfg = modelo().config_estudos();
    FS = cfg.FS; FIGPOS = cfg.FIGPOS; YCLIP = cfg.YCLIP;
    nPerImpl = cfg.nPerImpl; nPerRK = cfg.nPerRK;
    hRat = cfg.hRat;
    reductionModes = cfg.F.reduction_modes;

    p   = modelo().params();                 % xi_target = 0.2 por padrao
    mdl = modelo().monta(p);
    K = mdl.K; M = mdl.M; C = mdl.C; T1 = mdl.T1; iTip = mdl.info.tipW;

    [dirNames, dirLabel, rkNames] = metodos().nomes_metodos();
    condM    = [rkNames, {'DC'}];

    fprintf('=====================================================\n');
    fprintf(' ESTUDO F - MODAL REDUCTION  (xi = %.2f)\n', p.xi_target);
    fprintf('=====================================================\n');
    fprintf(' Full model: w_max = %.4e rad/s (mode %d)\n', mdl.wmax, mdl.info.ndof);

    res = struct();
    res.model_full = struct('K',K,'M',M,'w',mdl.w,'Phi',mdl.Phi,'xi',mdl.damp.xi, ...
        'T1',T1,'wmax',mdl.wmax);

    % ======================================================================
    %  PARTE 1 - comparacao: completo, 4 modos, 3 modos
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
        red = modelo().reducao_modal(mdl, m);
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

        stab = graficos().espectro(mo.w, mo.xi, sprintf('model: %s', mo.label), ['F_' mo.tag], outDir, ...
            'hSweep', cfg.hSweep, 'FS', FS, 'FIGPOS', FIGPOS, ...
            'xScaleBy', 1/mo.T1, 'xLabelStr', '$h/T_1$');
        res.stability.(mo.tag) = stab;

        hcr = zeros(numel(condM),1);
        for k = 1:numel(condM)
            hcr(k) = metodos().h_critico(mo.K, mo.M, mo.C, condM{k});
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
        pead = metodos().pead_tabela(mo.w, mo.xi, allNames, allLabels, allSteps, mo.T1, ...
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
        cfgT = graficos().aplica_config(cfgT, cfg);
        cfgT.runRK = hRK > 0;

        tresp = graficos().resposta_tempo(cfgT);
        res.time.(mo.tag) = tresp;
    end

    % --- Erro so de truncamento: solucao EXATA reduzida x solucao EXATA completa
    tCommon = (0:round(nPerImpl*T1/(hRat*T1)))*(hRat*T1);
    UexFull = modelo().solucao_exata(tCommon, M, mdl.w, mdl.Phi, mdl.damp.xi, mdl.u0, mdl.v0);
    uTipFull = UexFull(iTip,:);
    for i = 1:numel(reds)
        red = reds{i};
        UexRed = modelo().solucao_exata(tCommon, red.Mref, red.w, eye(red.m), red.xi, red.u0, red.v0);
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
    %  PARTE 2 - item 32: mesmo h x h maior, RK4, reducao de 4 modos
    % ======================================================================
    fprintf('\n=====================================================\n');
    fprintf(' PARTE 2 - Item 32: mesmo h x h maior (RK4, 4 modos)\n');
    fprintf('=====================================================\n');

    red4 = reds{1};                        % cfg.F.reduction_modes(1), hoje = 4
    rkP2 = sprintf('RK%d', cfg.F.item32_ordemRK);
    hcrFullRK4 = metodos().h_critico(K, M, C, rkP2);
    hcrRedRK4  = metodos().h_critico(red4.K, red4.M, red4.C, rkP2);
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
    UexFine  = modelo().solucao_exata(tRefFine, M, mdl.w, mdl.Phi, mdl.damp.xi, mdl.u0, mdl.v0);
    uexTipFine = UexFine(iTip,:);
    uPeakP2  = max(abs(uexTipFine));    % pico global: normaliza o ERRO

    % Mesma regra de janela dos demais estudos (graficos().janela_plot): com xi =
    % 0.2 a resposta morre em poucos ciclos, entao mostra-se a resposta
    % inteira em vez dos nPerP2 ciclos, quase todos achatados. So o que e
    % DESENHADO muda; a simulacao e o erro seguem cobrindo os nPerP2 ciclos.
    [tIniP2, tFimP2, infoP2] = graficos().janela_plot(tRefFine, uexTipFine, T1, ...
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
            % O caso 1 pode pedir milhoes de passos; a serie decimada basta para o interp1.
            saveIdx = unique(round(linspace(0, nSteps, min(nSteps+1, cfg.nMaxSalvar))));
            [t, uTip] = metodos().rk_stream_tip(cs.K, cs.M, cs.C, cs.u0, cs.v0, cs.h, ...
                nSteps, cfg.F.item32_ordemRK, cs.tipExtract, saveIdx);
            % Interpola na grade fina comum para comparar o erro de forma justa.
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
    graficos().salva_figura(fig, outDir, 'viga_item32_completo_vs_reduzido');

    res.config = struct('hRat',hRat,'reduction_modes',reductionModes, ...
        'xi',p.xi_target,'runDate',datestr(now));
    res.cfg = cfg;   % guarda a configuracao usada nesta rodada

    save(fullfile(outDir,'estudoF_resultados.mat'), 'res', '-v7.3');
    fprintf('\nEstudo F concluido. Resultados e figuras em %s\n', outDir);
    graficos().lista_obsoletos({outDir}, tInicio);
end
