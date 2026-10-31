function r = resultados()
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%
% Pos-processamento e testes: tabelas CSV, relatorio de parametros, teste
% rapido da janela de plot e testes com assert.
%
% As funcoes deste arquivo sao chamadas de fora como resultados().nome(...).
% So as listadas abaixo sao publicas; as demais sao de uso interno.

    r.relatorio    = @relatorio;
    r.tabelas      = @tabelas;
    r.teste_janela = @teste_janela;
    r.testes       = @testes;
end

%% ====================================================================
%  Tabelas CSV
%  ====================================================================

function tabelas()
% Consolida os .mat dos estudos em tabelas CSV em ImgViga/Tabelas/. Estudo
% ainda nao rodado e pulado com aviso, entao da para rodar a qualquer momento.
%
%   tabela_frequencias.csv      w_1..w_10 e T_1..T_10
%   tabela_amortecimento.csv    a0, a1 e xi de cada caso dos Estudos B e E
%   tabela_passos_criticos.csv  h_cr por metodo e por caso
%   tabela_N_erro.csv           N = t_final/h e erro relativo

    outDir  = graficos().pasta_saida();       % <pasta do codigo>/ImgViga, nao <pwd>/ImgViga
    tabDir  = fullfile(outDir,'Tabelas');
    if ~exist(tabDir,'dir'), mkdir(tabDir); end

    p = modelo().params();

    % ======================================================================
    %  1. FREQUENCIAS NATURAIS
    % ======================================================================
    [K, M, info] = modelo().matrizes_EF(p);
    w = modelo().autovalores_EF(K, M);
    nRep = min(10, numel(w));
    freqTable = table((1:nRep).', w(1:nRep), 2*pi./w(1:nRep), w(1:nRep)/w(1), ...
        'VariableNames', {'modo','w_rad_s','T_s','w_over_w1'});
    writetable(freqTable, fullfile(tabDir,'tabela_frequencias.csv'));
    fprintf('tabela_frequencias.csv escrita (%d modos, wmax = %.4e rad/s)\n', ...
        nRep, w(end));
    disp(freqTable);

    % ======================================================================
    %  2. AMORTECIMENTO
    % ======================================================================
    rows = {};
    fB = fullfile(outDir,'EstudoB','estudoB_resultados.mat');
    if exist(fB,'file')
        S = load(fB); resB = S.res;
        fns = fieldnames(resB.damp);
        for k = 1:numel(fns)
            d = resB.damp.(fns{k});
            rows(end+1,:) = {fns{k}, d.a0, d.a1, d.xi(1), d.xi(min(2,end)), d.xi(min(3,end))}; %#ok<SAGROW>
        end
    else
        warning('tabelas:missing', 'Estudo B nao encontrado (%s) - pulando xi=0/0.01/0.2 na tabela de amortecimento.', fB);
    end
    fE = fullfile(outDir,'EstudoE','estudoE_resultados.mat');
    if exist(fE,'file')
        S = load(fE); resE = S.res;
        fns = fieldnames(resE.damp);
        for k = 1:numel(fns)
            d = resE.damp.(fns{k});
            rows(end+1,:) = {['E_' fns{k}], d.a0, d.a1, d.xi(1), d.xi(min(2,end)), d.xi(min(3,end))}; %#ok<SAGROW>
        end
    else
        warning('tabelas:missing', 'Estudo E nao encontrado (%s) - pulando Rayleigh alternativo na tabela de amortecimento.', fE);
    end
    if ~isempty(rows)
        dampTable = cell2table(rows, 'VariableNames', {'caso','a0','a1','xi_1','xi_2','xi_3'});
        writetable(dampTable, fullfile(tabDir,'tabela_amortecimento.csv'));
        fprintf('\ntabela_amortecimento.csv escrita (%d casos)\n', height(dampTable));
        disp(dampTable);
    end

    % ======================================================================
    %  3. PASSOS CRITICOS
    % ======================================================================
    critRows = {};
    fA = fullfile(outDir,'viga_resultados.mat');
    if exist(fA,'file')
        S = load(fA); resA = S.res;
        hcrA = resA.stability.hcr_table;
        for r = 1:height(hcrA)
            critRows(end+1,:) = {'A_completo', hcrA.method{r}, hcrA.hcr_xi0(r), hcrA.hcr_xi02(r), NaN}; %#ok<SAGROW>
        end
    else
        warning('tabelas:missing', 'Estudo A nao encontrado (%s) - rode estudo_A primeiro.', fA);
    end
    if exist(fB,'file')
        S = load(fB); resB = S.res; %#ok<NASGU>
        hcrB = S.res.hcr;
        xiNames = hcrB.Properties.VariableNames;
        for r = 1:height(hcrB)
            vals = hcrB{r,:};
            for c = 1:numel(xiNames)
                critRows(end+1,:) = {sprintf('B_%s',xiNames{c}), hcrB.Properties.RowNames{r}, vals(c), NaN, NaN}; %#ok<SAGROW>
            end
        end
    end
    ffile = fullfile(outDir,'EstudoF','estudoF_resultados.mat');
    if exist(ffile,'file')
        S = load(ffile); resF = S.res;
        fns = fieldnames(resF.hcr);
        for k = 1:numel(fns)
            hcrF = resF.hcr.(fns{k});
            for r = 1:height(hcrF)
                critRows(end+1,:) = {['F_' fns{k}], hcrF.method{r}, hcrF.hcr(r), hcrF.hcr_T1(r), NaN}; %#ok<SAGROW>
            end
        end
    else
        warning('tabelas:missing', 'Estudo F nao encontrado (%s).', ffile);
    end
    if ~isempty(critRows)
        critTable = cell2table(critRows, 'VariableNames', {'caso','metodo','hcr_ou_valor1','valor2','valor3'});
        writetable(critTable, fullfile(tabDir,'tabela_passos_criticos.csv'));
        fprintf('\ntabela_passos_criticos.csv escrita (%d linhas)\n', height(critTable));
    end

    % ======================================================================
    %  4. N (numero de passos) x ERRO, para os metodos diretos
    % ======================================================================
    neRows = {};
    neRows = collect_N_err(fA, 'A', neRows);
    if exist(fB,'file')
        S = load(fB);
        fns = fieldnames(S.res.time);
        for k = 1:numel(fns)
            tresp = S.res.time.(fns{k});
            neRows = collect_N_err_from_tresp(tresp, sprintf('B_%s',fns{k}), neRows);
        end
    end
    if ~isempty(neRows)
        neTable = cell2table(neRows, 'VariableNames', {'caso','metodo','N_passos','h_s','erro_relativo'});
        writetable(neTable, fullfile(tabDir,'tabela_N_erro.csv'));
        fprintf('\ntabela_N_erro.csv escrita (%d linhas)\n', height(neTable));
    end

    fprintf('\nTodas as tabelas disponiveis foram escritas em %s\n', tabDir);
end

function rows = collect_N_err(matFile, caseTag, rows)
    if ~exist(matFile,'file'), return; end
    S = load(matFile);
    if ~isfield(S,'res') || ~isfield(S.res,'time'), return; end
    fns = fieldnames(S.res.time);
    for k = 1:numel(fns)
        tresp = S.res.time.(fns{k});
        rows = collect_N_err_from_tresp(tresp, sprintf('%s_%s',caseTag,fns{k}), rows);
    end
end

function rows = collect_N_err_from_tresp(tresp, tag, rows)
    if ~isstruct(tresp), return; end
    methodFields = setdiff(fieldnames(tresp), {'exact','exactRK','config'});
    for m = 1:numel(methodFields)
        mf = methodFields{m};
        entry = tresp.(mf);
        if isstruct(entry) && isfield(entry,'N') && isfield(entry,'h')
            errVal = NaN;
            if isfield(entry,'err'), errVal = entry.err; end
            rows(end+1,:) = {tag, mf, entry.N, entry.h, errVal}; %#ok<AGROW>
        end
    end
end

%% ====================================================================
%  Relatorio de parametros
%  ====================================================================

function relatorio()
% Escreve ImgViga/Tabelas/PARAMETROS_ESTUDOS.md com todos os parametros de cada
% estudo: modelo fisico, frequencias, configuracao numerica, casos varridos,
% passos usados, passos criticos, tabela N x erro e um glossario dos simbolos.
% Estudos ainda nao rodados sao pulados com aviso.

    outDir = graficos().pasta_saida();        % <pasta do codigo>/ImgViga, nao <pwd>/ImgViga
    tabDir = fullfile(outDir,'Tabelas');
    if ~exist(tabDir,'dir'), mkdir(tabDir); end
    mdFile = fullfile(tabDir,'PARAMETROS_ESTUDOS.md');

    cfg = modelo().config_estudos();
    p   = modelo().params();
    [K, M, info] = modelo().matrizes_EF(p);
    w = modelo().autovalores_EF(K, M);
    T1 = 2*pi/w(1);

    fid = fopen(mdFile,'w');
    if fid < 0
        error('relatorio:fopen','Nao consegui escrever em %s', mdFile);
    end
    closer = onCleanup(@() fclose(fid));

    % ======================================================================
    %  CABECALHO
    % ======================================================================
    w_(fid, '# Simulation parameters — cantilever beam time-integration study');
    w_(fid, '');
    w_(fid, 'Summary of every parameter used in the beam study.');
    w_(fid, 'Generated by `viga(''relatorio'')` on %s.', datestr(now,'yyyy-mm-dd HH:MM'));
    w_(fid, '');
    w_(fid, 'Nine time-integration schemes are applied to a damped cantilever beam');
    w_(fid, 'discretised by Euler-Bernoulli finite elements. Five schemes are *direct*');
    w_(fid, '(Central Differences, Houbolt, Newmark, Wilson-Theta, Bathe) and four are');
    w_(fid, 'explicit *Runge-Kutta* (RK1-RK4). Each study below varies one factor.');
    w_(fid, '');

    % ======================================================================
    %  MODELO FISICO
    % ======================================================================
    w_(fid, '## 1. Physical model (`modelo().params`)');
    w_(fid, '');
    w_(fid, 'Identical in every study — only the numerical experiment changes.');
    w_(fid, '');
    w_(fid, '| Parameter | Symbol | Value | Unit |');
    w_(fid, '|---|---|---|---|');
    w_(fid, '| Young''s modulus | E | %.4g | Pa |', p.E);
    w_(fid, '| Mass density | rho | %.4g | kg/m^3 |', p.rho);
    w_(fid, '| Length | L | %.4g | m |', p.L);
    w_(fid, '| Cross-section width | b | %.4g | m |', p.b);
    w_(fid, '| Cross-section height | h | %.4g | m |', p.h);
    w_(fid, '| Second moment of area | I | %.6g | m^4 |', p.I);
    w_(fid, '| Slenderness ratio | L/h | %.4g | - |', p.slenderness);
    w_(fid, '| Beam elements (2-node) | N | %d | - |', p.N);
    w_(fid, '| Degrees of freedom after clamping | ndof | %d | - |', info.ndof);
    w_(fid, '| Reference damping ratio (Estudo A) | xi_target | %.4g | - |', p.xi_target);
    w_(fid, '| Rayleigh anchor modes | i, j | %d, %d | - |', p.mode_i, p.mode_j);
    w_(fid, '| Initial tip velocity (impulse IC) | V0 | %.4g | m/s |', p.V0);
    w_(fid, '');
    w_(fid, '### Natural frequencies');
    w_(fid, '');
    w_(fid, '| Mode | omega [rad/s] | T [s] | omega/omega_1 |');
    w_(fid, '|---|---|---|---|');
    nRep = min(10, numel(w));
    for k = 1:nRep
        w_(fid, '| %d | %.6g | %.6g | %.6g |', k, w(k), 2*pi/w(k), w(k)/w(1));
    end
    w_(fid, '');
    w_(fid, 'Highest mesh frequency omega_max = %.6g rad/s (mode %d). This value, not', w(end), numel(w));
    w_(fid, 'omega_1, sets the stability limit of the conditionally stable schemes.');
    w_(fid, '');

    % ======================================================================
    %  CONFIGURACAO NUMERICA GLOBAL
    % ======================================================================
    w_(fid, '## 2. Numerical configuration (`modelo().config_estudos`)');
    w_(fid, '');
    w_(fid, 'Shared by Estudos B-F. Estudo A keeps its own inline configuration.');
    w_(fid, '');
    w_(fid, '### Time step');
    w_(fid, '');
    w_(fid, '| Key | Value | Meaning |');
    w_(fid, '|---|---|---|');
    w_(fid, '| `hRat` | %g | Direct-scheme step as a fraction of T_1: h_impl = hRat * T_1 = %.6g s |', cfg.hRat, cfg.hRat*T1);
    w_(fid, '| `fatorSegurancaRK` | %g | RK step as a fraction of the critical step: h_RK = %g * h_cr |', cfg.fatorSegurancaRK, cfg.fatorSegurancaRK);
    w_(fid, '');
    w_(fid, '### Simulated duration and plotted window');
    w_(fid, '');
    w_(fid, 'The response is integrated for `nPer*` cycles of T_1 but only a window of a');
    w_(fid, 'few cycles is plotted, because period (PE) and amplitude (AD) errors');
    w_(fid, 'accumulate and are only visible after several cycles.');
    w_(fid, '');
    w_(fid, 'Where that window sits depends on damping (three regimes):');
    w_(fid, '');
    w_(fid, '- **Strongly damped (xi != 0, decays within `ciclosDecaimentoRapido` cycles):**');
    w_(fid, '  the response dies out quickly, so the **whole response** is plotted, from');
    w_(fid, '  t = 0 to the decay instant (the last instant the *exact* response still');
    w_(fid, '  exceeds `tolDecaimento` times its own peak).');
    w_(fid, '- **Weakly damped (xi != 0, decays slower):** the code finds the same decay');
    w_(fid, '  instant and plots the `X_ciclos_final_plot` cycles before it. If the response');
    w_(fid, '  has not decayed by the end of the simulation, the last cycles are plotted.');
    w_(fid, '- **Undamped (xi == 0):** nothing decays, so the window is simply the last');
    w_(fid, '  `X_ciclos_final_plot_sem_xi` cycles of the simulation.');
    w_(fid, '');
    w_(fid, 'A response only counts as decayed if it stays below the tolerance for at least');
    w_(fid, 'one full period T_1 before the end of the simulation. Decay is always measured');
    w_(fid, 'on the exact solution, never on a scheme''s own response.');
    w_(fid, '');
    w_(fid, '| Key | Value | Meaning |');
    w_(fid, '|---|---|---|');
    w_(fid, '| `nPerImpl` | %g | Cycles computed, direct schemes, damped case |', cfg.nPerImpl);
    w_(fid, '| `nPerRK` | %g | Cycles computed, RK schemes, damped case |', cfg.nPerRK);
    w_(fid, '| `X_ciclos_final_plot` | %g | Cycles plotted, damped case (back from the decay instant) |', cfg.X_ciclos_final_plot);
    w_(fid, '| `nPerImpl_sem_xi` | %g | Cycles computed, direct schemes, undamped case |', cfg.nPerImpl_sem_xi);
    w_(fid, '| `nPerRK_sem_xi` | %g | Cycles computed, RK schemes, undamped case |', cfg.nPerRK_sem_xi);
    w_(fid, '| `X_ciclos_final_plot_sem_xi` | %g | Cycles plotted, undamped case (back from the end) |', cfg.X_ciclos_final_plot_sem_xi);
    w_(fid, '| `tolDecaimento` | %.1e | Fraction of its own peak below which the response counts as zero |', cfg.tolDecaimento);
    w_(fid, '| `ciclosDecaimentoRapido` | %g | If the response decays within this many cycles, the whole response is plotted |', cfg.ciclosDecaimentoRapido);
    w_(fid, '');
    w_(fid, '### Figures and sampling');
    w_(fid, '');
    w_(fid, '| Key | Value | Meaning |');
    w_(fid, '|---|---|---|');
    w_(fid, '| `FS` | %g | Axis font size (pt) |', cfg.FS);
    w_(fid, '| `FIGPOS` | [%g %g %g %g] | Exported figure size (cm) |', cfg.FIGPOS);
    w_(fid, '| `YCLIP` | %g | y-axis clip of the normalised response |', cfg.YCLIP);
    w_(fid, '| `nModesReport` | %g | Modes printed in the console PE/AD tables |', cfg.nModesReport);
    w_(fid, '| `hSweep` | %d points, %.1e to %.1e s | Step sweep of the spectral-radius figures (one figure per case, all nine schemes together: direct solid, Runge-Kutta dashed) |', numel(cfg.hSweep), cfg.hSweep(1), cfg.hSweep(end));
    w_(fid, '| `nPontosExato` | %g | Points on the exact reference curve |', cfg.nPontosExato);
    w_(fid, '| `nMaxSalvar` | %g | Max samples stored per RK method |', cfg.nMaxSalvar);
    w_(fid, '| `nMaxPlot` | %g | Max points drawn per curve |', cfg.nMaxPlot);
    w_(fid, '');
    w_(fid, 'The last three are memory/plot limits only: every time step is still');
    w_(fid, 'integrated, they just control how much is stored and drawn.');
    w_(fid, '');

    % ======================================================================
    %  UM CAPITULO POR ESTUDO
    % ======================================================================
    w_(fid, '## 3. Studies');
    w_(fid, '');

    % --- Estudo A ---------------------------------------------------------
    fA = fullfile(outDir,'viga_resultados.mat');
    w_(fid, '### Estudo A — reference case (already in the paper)');
    w_(fid, '');
    w_(fid, 'Source: `estudo_A`. Keeps its own configuration, as the reference case;');
    w_(fid, 'it does NOT read `modelo().config_estudos`.');
    w_(fid, '');
    if exist(fA,'file')
        S = load(fA);
        w_(fid, '- Impulse initial condition, damping ratios xi = 0 and 0.2');
        w_(fid, '- Status: **run**, results in `ImgViga/viga_resultados.mat`');
        if isfield(S,'res') && isfield(S.res,'stability') && isfield(S.res.stability,'hcr_table')
            hcrA = S.res.stability.hcr_table;
            w_(fid, '');
            w_(fid, '| Scheme | h_cr at xi=0 [s] | h_cr at xi=0.2 [s] |');
            w_(fid, '|---|---|---|');
            for r = 1:height(hcrA)
                w_(fid, '| %s | %.6g | %.6g |', hcrA.method{r}, hcrA.hcr_xi0(r), hcrA.hcr_xi02(r));
            end
        end
    else
        w_(fid, '- Status: **not run** (`%s` missing)', fA);
    end
    w_(fid, '');

    % --- Estudos B a F ----------------------------------------------------
    estudos = struct( ...
      'letra', {'B','C','D','E','F'}, ...
      'titulo', { ...
         'physical damping sweep', ...
         'modal initial conditions, fixed time step', ...
         'modal initial conditions, swept time step', ...
         'alternative Rayleigh damping', ...
         'modal reduction (truncation)'}, ...
      'script', { ...
         'estudo_B', ...
         'estudo_C', ...
         'estudo_D', ...
         'estudo_E', ...
         'estudo_F'}, ...
      'pergunta', { ...
         'How do the schemes behave as the physical damping changes, with the time step held the same?', ...
         'How do they behave when a single, progressively higher mode is excited, with h held fixed?', ...
         'What happens as h grows relative to the period of the excited mode, up to the stability limit?', ...
         'What changes when Rayleigh coefficients are imposed directly instead of anchored to a target xi?', ...
         'Can the relevant response be reproduced by a truncated model, and does that relax the time step?'});

    for ie = 1:numel(estudos)
        es = estudos(ie);
        mf = fullfile(outDir, ['Estudo' es.letra], sprintf('estudo%s_resultados.mat', es.letra));

        w_(fid, '### Estudo %s — %s', es.letra, es.titulo);
        w_(fid, '');
        w_(fid, 'Source: `%s`', es.script);
        w_(fid, '');
        w_(fid, '**Question:** %s', es.pergunta);
        w_(fid, '');

        % Casos varridos, direto da configuracao
        switch es.letra
            case 'B'
                w_(fid, '**Cases swept:** xi = %s', mat2str(cfg.B.xi_list));
                w_(fid, '');
                w_(fid, 'A single RK step, sized from the most damped case, is shared across all');
                w_(fid, 'xi values on purpose: that is what lets RK1/RK2 be shown diverging at low');
                w_(fid, 'damping and stable at high damping with the *same* step.');
                if ~isempty(cfg.B.xiZoom)
                    w_(fid, '');
                    w_(fid, 'Zoom inset on the xi = %g case.', cfg.B.xiZoom);
                end
            case 'C'
                w_(fid, '**Cases swept:** excited mode = %s, at xi = %g', mat2str(cfg.C.modos), cfg.C.xi);
                w_(fid, '');
                w_(fid, 'h is held fixed in ABSOLUTE terms (%.6g s) across every mode, so h/T_mode', cfg.hRat*T1);
                w_(fid, 'grows as the excited mode gets higher — that ratio is the point of the study.');
            case 'D'
                w_(fid, '**Cases swept:** excited mode = %s, h/T_mode = %s, at xi = %g', ...
                    mat2str(cfg.D.modos), mat2str(cfg.D.hRatios), cfg.D.xi);
                w_(fid, '');
                w_(fid, '| Key | Value | Meaning |');
                w_(fid, '|---|---|---|');
                w_(fid, '| `D.nPer` | %g | Cycles computed (counted in T_1, the fundamental period - NOT in the excited mode''s period) |', cfg.D.nPer);
                w_(fid, '| `D.ciclosPlot` | %g | Cycles plotted (in T_1) |', cfg.D.ciclosPlot);
                w_(fid, '| `D.ratiosSimuladas` | `%s` | Which h ratios get a time-domain run |', cfg.D.ratiosSimuladas);
                w_(fid, '');
                w_(fid, 'Closed-form PE/AD is computed for EVERY (mode, ratio) pair — that part is');
                w_(fid, 'cheap. Time-domain simulation runs only for the ratios selected above, with');
                w_(fid, 'fewer cycles, because the high ratios are expected to diverge (that');
                w_(fid, 'divergence IS the result) and each run is otherwise expensive.');
            case 'E'
                w_(fid, '**Cases swept:**');
                w_(fid, '');
                w_(fid, '| Case | a0 | a1 |');
                w_(fid, '|---|---|---|');
                for k = 1:numel(cfg.E.rayleigh)
                    w_(fid, '| %d | %g | %g |', k, cfg.E.rayleigh(k).a0, cfg.E.rayleigh(k).a1);
                end
                w_(fid, '');
                w_(fid, 'With a0 = 0 only the stiffness-proportional term survives, so');
                w_(fid, 'xi_k = a1*omega_k/2 grows LINEARLY with frequency: there is no anchor');
                w_(fid, 'mode and no two modes forced to share a xi. The resulting per-mode xi is');
                w_(fid, 'computed, not imposed.');
            case 'F'
                w_(fid, '**Cases swept:** full model vs. modal truncations keeping %s modes, at xi = %g', ...
                    mat2str(cfg.F.reduction_modes), p.xi_target);
                w_(fid, '');
                w_(fid, 'Part 2 isolates the effect of a larger step from the effect of the');
                w_(fid, 'reduction itself, using RK%d on the %d-mode reduction:', ...
                    cfg.F.item32_ordemRK, cfg.F.reduction_modes(1));
                w_(fid, '');
                w_(fid, '- Case 1: full model, h = %g * h_cr,full', cfg.fatorSegurancaRK);
                w_(fid, '- Case 2: reduced model, the SAME h as Case 1');
                w_(fid, '- Case 3: reduced model, h = %g * h_cr,reduced (larger)', cfg.fatorSegurancaRK);
                w_(fid, '');
                w_(fid, 'All three are compared against the FULL model''s exact response.');
        end
        w_(fid, '');

        % Resultados que existem em disco
        if ~exist(mf,'file')
            w_(fid, 'Status: **not run** (`%s` missing).', mf);
            w_(fid, '');
            warning('relatorio:missing', ...
                'Estudo %s nao encontrado (%s) - secao incompleta no report.', es.letra, mf);
            continue
        end

        S = load(mf);
        res = S.res;
        if isfield(res,'config') && isfield(res.config,'runDate')
            w_(fid, 'Status: **run** on %s', res.config.runDate);
        else
            w_(fid, 'Status: **run**');
        end
        w_(fid, '');

        % Passos derivados
        if isfield(res,'config')
            c = res.config;
            linhas = {};
            if isfield(c,'hImpl'), linhas{end+1} = sprintf('| h_impl | %.6g s |', c.hImpl); end %#ok<SAGROW>
            if isfield(c,'hRK'),   linhas{end+1} = sprintf('| h_RK | %.6g s |', c.hRK); end %#ok<SAGROW>
            if ~isempty(linhas)
                w_(fid, '**Time steps actually used**');
                w_(fid, '');
                w_(fid, '| Quantity | Value |');
                w_(fid, '|---|---|');
                for k = 1:numel(linhas), w_(fid, '%s', linhas{k}); end
                w_(fid, '');
            end
        end

        % Passos criticos
        if isfield(res,'hcr')
            w_(fid, '**Critical time steps**');
            w_(fid, '');
            escreve_hcr(fid, res.hcr);
            w_(fid, '');
        end

        % N e erro por caso
        escreve_N_erro(fid, res);

        % Estudo F, parte 2
        if strcmp(es.letra,'F') && isfield(res,'item32')
            w_(fid, '**Part 2 — full vs. reduced model**');
            w_(fid, '');
            w_(fid, '| Case | h [s] | Steps N | Max relative error |');
            w_(fid, '|---|---|---|---|');
            fns = fieldnames(res.item32);
            for k = 1:numel(fns)
                it = res.item32.(fns{k});
                if isfield(it,'failed') && it.failed
                    w_(fid, '| %s | %.6g | - | FAILED: %s |', strrep(fns{k},'_',' '), it.h, it.message);
                else
                    w_(fid, '| %s | %.6g | %d | %.4e |', strrep(fns{k},'_',' '), it.h, it.N, it.err);
                end
            end
            w_(fid, '');
        end
        if strcmp(es.letra,'F') && isfield(res,'truncation')
            w_(fid, '**Modal truncation error** (exact solutions only, no time integration)');
            w_(fid, '');
            w_(fid, '| Reduced model | Max relative tip error |');
            w_(fid, '|---|---|');
            fns = fieldnames(res.truncation);
            for k = 1:numel(fns)
                w_(fid, '| %s | %.4e |', fns{k}, res.truncation.(fns{k}).err);
            end
            w_(fid, '');
        end
    end

    % ======================================================================
    %  GLOSSARIO
    % ======================================================================
    w_(fid, '## 4. Glossary');
    w_(fid, '');
    w_(fid, '| Symbol | Meaning |');
    w_(fid, '|---|---|');
    w_(fid, '| T_1 | Period of the first (fundamental) mode |');
    w_(fid, '| omega_k | Natural frequency of mode k |');
    w_(fid, '| omega_max | Highest frequency of the FE mesh; sets the stability limit |');
    w_(fid, '| xi | Modal damping ratio |');
    w_(fid, '| a0, a1 | Rayleigh damping coefficients, C = a0*M + a1*K |');
    w_(fid, '| h | Time step |');
    w_(fid, '| h_cr | Critical (largest stable) time step of a conditionally stable scheme |');
    w_(fid, '| h_impl | Step used by the direct schemes, = hRat * T_1 |');
    w_(fid, '| h_RK | Step used by the Runge-Kutta schemes, = fatorSegurancaRK * h_cr |');
    w_(fid, '| N | Number of time steps, = t_final / h |');
    w_(fid, '| PE | Period Elongation — ratio of numerical to exact period |');
    w_(fid, '| AD | Amplitude Decay — ratio of numerical to exact amplitude per cycle |');
    w_(fid, '| rho | Spectral radius of the amplification operator; rho > 1 means unstable |');
    w_(fid, '| w_tip | Transverse displacement of the free end |');
    w_(fid, '');
    w_(fid, '**Schemes.** Direct: Central Differences (conditionally stable), Houbolt,');
    w_(fid, 'Newmark, Wilson-Theta, Bathe. Explicit Runge-Kutta: RK1 (Euler), RK2');
    w_(fid, '(modified Euler), RK3 (Heun), RK4 (classical) — all conditionally stable.');
    w_(fid, '');

    clear closer;   % fecha o arquivo
    fprintf('Report escrito em %s\n', mdFile);
    fprintf('(%d KB)\n', round(dir(mdFile).bytes/1024));
end

function w_(fid, fmt, varargin)
% Escreve uma linha no markdown. Sem argumentos extras, o texto vai
% literal (evita que um '%' ou '\' no conteudo seja interpretado).
    if nargin < 3
        fprintf(fid, '%s\n', fmt);
    else
        fprintf(fid, [fmt '\n'], varargin{:});
    end
end

function escreve_hcr(fid, hcr)
% hcr pode ser uma table (Estudo B) ou uma struct de tables (E, F).
    if istable(hcr)
        vn = hcr.Properties.VariableNames;
        rn = hcr.Properties.RowNames;
        % 'xi_0_01' -> 'xi = 0.01' (makeValidName trocou o ponto por '_')
        cab = regexprep(vn, '^xi_', 'xi = ');
        cab = regexprep(cab, '(?<=\d)_(?=\d)', '.');
        fprintf(fid, '| Scheme | %s |\n', strjoin(cab, ' | '));
        fprintf(fid, '|---%s|\n', repmat('|---', 1, numel(vn)));
        for r = 1:height(hcr)
            vals = hcr{r,:};
            celulas = arrayfun(@(v) sprintf('%.6g', v), vals, 'UniformOutput', false);
            fprintf(fid, '| %s | %s |\n', rn{r}, strjoin(celulas, ' | '));
        end
    elseif isstruct(hcr)
        fns = fieldnames(hcr);
        for k = 1:numel(fns)
            t = hcr.(fns{k});
            fprintf(fid, '\n*Case `%s`:*\n\n', fns{k});
            fprintf(fid, '| Scheme | h_cr [s] | h_cr/T_1 |\n|---|---|---|\n');
            for r = 1:height(t)
                fprintf(fid, '| %s | %.6g | %.6g |\n', t.method{r}, t.hcr(r), t.hcr_T1(r));
            end
        end
    end
end

function escreve_N_erro(fid, res)
% Numero de passos e erro relativo por metodo, para cada caso do estudo.
    if ~isfield(res,'time'), return; end
    casos = fieldnames(res.time);
    escreveu = false;
    for c = 1:numel(casos)
        tresp = res.time.(casos{c});
        if ~isstruct(tresp), continue; end
        % Estudo D aninha mais um nivel (modo -> razao)
        if ~isfield(tresp,'config')
            sub = fieldnames(tresp);
            for s = 1:numel(sub)
                escreveu = bloco_N_erro(fid, tresp.(sub{s}), ...
                    sprintf('%s / %s', casos{c}, sub{s}), escreveu) || escreveu;
            end
        else
            escreveu = bloco_N_erro(fid, tresp, casos{c}, escreveu) || escreveu;
        end
    end
    if escreveu, fprintf(fid, '\n'); end
end

function escreveu = bloco_N_erro(fid, tresp, tag, jaEscreveu)
    escreveu = false;
    if ~isstruct(tresp) || isfield(tresp,'failed'), return; end
    nomesM = setdiff(fieldnames(tresp), {'exact','exactRK','config','cfg'});
    linhas = {};
    for m = 1:numel(nomesM)
        e = tresp.(nomesM{m});
        if isstruct(e) && isfield(e,'N') && isfield(e,'h')
            if isfield(e,'err')
                linhas{end+1} = sprintf('| %s | %s | %.6g | %d | %.4e |', ...
                    tag, nomesM{m}, e.h, e.N, e.err); %#ok<AGROW>
            else
                linhas{end+1} = sprintf('| %s | %s | %.6g | %d | - |', ...
                    tag, nomesM{m}, e.h, e.N); %#ok<AGROW>
            end
        end
    end
    if isempty(linhas), return; end
    if ~jaEscreveu
        fprintf(fid, '**Steps and error per scheme**\n\n');
        fprintf(fid, '| Case | Scheme | h [s] | Steps N | Relative error |\n');
        fprintf(fid, '|---|---|---|---|---|\n');
    end
    for k = 1:numel(linhas), fprintf(fid, '%s\n', linhas{k}); end
    escreveu = true;
end

%% ====================================================================
%  Testes
%  ====================================================================

function teste_janela()
% Teste rapido da janela de plot: confere em segundos onde a janela cai, sem
% pagar uma bateria com Runge-Kutta. Roda so os metodos diretos em tres casos,
% um por regime de viga_janela_plot: muito amortecido, pouco amortecido e sem
% amortecimento. Para cada um imprime a janela e a estimativa de ciclos ate o
% decaimento, n = -ln(tolDecaimento)/(2*pi*xi).
%
% Ajuste cfg (secao 2 de modelo().config_estudos), rode viga('janela') e olhe as
% figuras em ImgViga/TesteJanela/ antes de rodar os estudos.

    close all;


    tInicio = now;                       %#ok<TNOW1>
    outDir = graficos().pasta_saida('TesteJanela');   % <pasta do codigo>/ImgViga/TesteJanela
    if ~exist(outDir,'dir'), mkdir(outDir); end

    cfg = modelo().config_estudos();
    p   = modelo().params();

    [K, M, info] = modelo().matrizes_EF(p);
    [w, Phi]     = modelo().autovalores_EF(K, M);
    T1   = 2*pi/w(1);
    iTip = info.tipW;

    fprintf('#############################################################\n');
    fprintf('  TESTE DA JANELA DE PLOT (so metodos diretos - rapido)\n');
    fprintf('#############################################################\n');
    fprintf(' T1 = %.6f s\n', T1);
    fprintf(' Configuracao atual (modelo().config_estudos, secao 2):\n');
    fprintf('   AMORTECIDO      : calcula %d ciclos, plota os ultimos %g antes do decaimento\n', ...
        cfg.nPerImpl, cfg.X_ciclos_final_plot);
    fprintf('   NAO AMORTECIDO  : calcula %d ciclos, plota os ultimos %g\n', ...
        cfg.nPerImpl_sem_xi, cfg.X_ciclos_final_plot_sem_xi);
    fprintf('   tolDecaimento   : %.1e do pico\n', cfg.tolDecaimento);
    fprintf('   ciclosDecaimentoRapido : %g (decaindo em ate isso, plota a resposta inteira)\n\n', ...
        cfg.ciclosDecaimentoRapido);

    hImpl = cfg.hRat*T1;

    % ======================================================================
    %  OS TRES CASOS
    % ======================================================================
    casos = struct( ...
        'tag',   {'muito_amortecido','pouco_amortecido','sem_amortecimento'}, ...
        'xi',    {p.xi_target, cfg.C.xi, 0}, ...
        'label', {sprintf('\\xi = %.3g', p.xi_target), ...
                  sprintf('\\xi = %.3g', cfg.C.xi), '\xi = 0'});

    res = struct();

    for ic = 1:numel(casos)
        cs = casos(ic);
        fprintf('--- Caso: %s (xi = %.3g) ---\n', cs.tag, cs.xi);

        % Estimativa analitica de quantos ciclos a resposta leva para cair a
        % tolDecaimento, para comparar com quantos ciclos serao calculados.
        if cs.xi > 0
            nMorte = -log(cfg.tolDecaimento)/(2*pi*cs.xi);
            nCalc  = cfg.nPerImpl;
            fprintf('   estimativa: decai a %.0e do pico em ~%.1f ciclos (serao calculados %d)\n', ...
                cfg.tolDecaimento, nMorte, nCalc);
            if nMorte <= cfg.ciclosDecaimentoRapido
                fprintf('   -> decaimento RAPIDO (<= %g ciclos): sera plotada a resposta inteira\n', ...
                    cfg.ciclosDecaimentoRapido);
            end
            if nMorte > nCalc
                fprintf(['   -> nPerImpl e PEQUENO demais para a deteccao disparar neste xi; ' ...
                         'seriam precisos ~%d ciclos\n'], ceil(nMorte)+1);
            end
        end

        if cs.xi == 0
            C    = zeros(size(K));
            ximo = zeros(size(w));
        else
            pp = p; pp.xi_target = cs.xi;
            [C, damp] = modelo().amortecimento(K, M, w, pp);
            ximo = damp.xi;
        end

        [u0, v0] = modelo().ci_impulso(K, M, C, info, p);

        cfgT = graficos().caso_viga_completa(K, M, C, u0, v0, w, Phi, ximo, iTip);
        cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = hImpl;   % hRK nao e usado
        cfgT.label = ['teste_' cs.tag];
        cfgT.xiLabel = cs.label;
        cfgT.outDir = outDir;
        cfgT = graficos().aplica_config(cfgT, cfg);
        cfgT.runRK = false;    % sem Runge-Kutta o teste roda em segundos

        tresp = graficos().resposta_tempo(cfgT);
        res.(cs.tag) = tresp;

        close all;

        j = tresp.config.janelaDireto;
        txtRegime = '';
        if j.completa, txtRegime = ', resposta inteira'; end
        fprintf('   -> janela final: t/T1 em [%.2f, %.2f]  (%.2f ciclos%s)\n\n', ...
            j.cicloIni, j.cicloFim, j.cicloFim - j.cicloIni, txtRegime);
    end

    fprintf('#############################################################\n');
    fprintf('  Figuras (eps + png) em %s\n', outDir);
    graficos().lista_obsoletos({outDir}, tInicio);
    fprintf('  Se as janelas estiverem boas, rode os estudos.\n');
    fprintf('#############################################################\n');
end

function testes()
% Testes rapidos (asserts) das funcoes principais. Roda em cerca de um minuto
% e grava figuras so em ImgViga/TesteFumaca/. Para na primeira falha.

    close all;


    p = modelo().params();
    [K, M, info] = modelo().matrizes_EF(p);
    [w, Phi] = modelo().autovalores_EF(K, M);
    T1 = 2*pi/w(1);
    fprintf('OK modelo().params/matrizes/autovalores: w1=%.4f, T1=%.6f, ndof=%d\n', w(1), T1, info.ndof);

    % --- viga_amortecimento: modo ancorado e modo com coeficientes diretos ---
    [Ca, dampA] = modelo().amortecimento(K, M, w, p);
    assert(isempty(setdiff({'a0','a1','xi','anchors','xi_target'}, fieldnames(dampA))), 'campos faltando (modo ancorado)');
    [Cd, dampD] = modelo().amortecimento(K, M, w, p, 0, 1e-4);
    assert(isequal(dampD.a0,0) && isequal(dampD.a1,1e-4), 'coeficientes diretos nao aplicados');
    assert(isempty(dampD.anchors), 'anchors deveria ser vazio no modo direto');
    fprintf('OK modelo().amortecimento (2 modos): xi1_ancorado=%.4f, xi1_direto=%.4f\n', dampA.xi(1), dampD.xi(1));

    % --- modelo().ci_modal ---
    [u0m, v0m] = modelo().ci_modal(K, M, Ca, Phi, 3);
    assert(norm(v0m) == 0, 'v0 deveria ser zero');
    q0check = Phi.' * (M * u0m);
    [~, idxmax] = max(abs(q0check));
    assert(idxmax == 3, 'CI modal nao projeta no modo esperado');
    fprintf('OK viga_ci_modal: q0 concentrado no modo %d (valor %.4f, resto <= %.2e)\n', ...
        idxmax, q0check(idxmax), max(abs(q0check([1 2 4:end]))));

    % --- modelo().reducao_modal ---
    mdl = modelo().monta(p);
    red4 = modelo().reducao_modal(mdl, 4);
    assert(isequal(size(red4.K),[4 4]) && isequal(size(red4.M),[4 4]), 'tamanho errado do modelo reduzido');
    assert(all(abs(diag(red4.M)-1) < 1e-10), 'M reduzido deveria ser identidade');
    assert(red4.wmax == mdl.w(4), 'wmax reduzido incorreto');
    fprintf('OK viga_reducao_modal: wmax_full=%.3e, wmax_red4=%.3e (%.1fx menor)\n', ...
        mdl.wmax, red4.wmax, mdl.wmax/red4.wmax);

    % --- metodos().rho_beam / graficos().espectro (sem plot, sweep pequeno) ---
    r = metodos().rho_beam('RK4', mdl.w, mdl.damp.xi, 1e-6);
    assert(isfinite(r) && r >= 0, 'rho invalido');
    stab = graficos().espectro(mdl.w, mdl.damp.xi, '\xi=0.2', 'smoketest', pwd, ...
        'hSweep', logspace(-7,-4,5), 'doPlots', false);
    assert(isfield(stab.rho.rk,'RK4') && numel(stab.rho.rk.RK4)==5, 'graficos().espectro nao retornou o esperado');
    fprintf('OK metodos().rho_beam/viga_espectro: rho(RK4,h=1e-6)=%.6f\n', r);

    % --- graficos().espectro com plot: UMA figura com os 9 metodos (diretos + RK) ---
    spDir = graficos().pasta_saida('TesteFumaca');
    graficos().espectro(mdl.w, mdl.damp.xi, '\xi=0.2', 'smoketest', spDir, ...
        'hSweep', logspace(-7,-4,5));
    assert(exist(fullfile(spDir,'viga_rho_smoketest.eps'),'file')==2, 'figura combinada nao foi gerada (.eps)');
    assert(exist(fullfile(spDir,'viga_rho_smoketest.png'),'file')==2, 'copia .png nao foi gerada');
    figRho = findobj('Type','figure','Name','rho smoketest');
    assert(numel(figRho)==1, 'deveria haver UMA figura de raio espectral');
    assert(numel(findobj(figRho,'Type','line','LineStyle','-')) >= 5 && ...
           numel(findobj(figRho,'Type','line','LineStyle','--')) >= 4, ...
           'figura deveria ter 5 diretos (cheio) e 4 RK (tracejado)');
    close(figRho);
    fprintf('OK viga_espectro: 1 figura, 5 diretos + 4 RK\n');

    % --- viga_janela_plot: os tres regimes ---
    tJ = linspace(0, 50*T1, 4001);
    [~,~,iJ] = graficos().janela_plot(tJ, exp(-2*pi*0.2*tJ/T1).*sin(2*pi*tJ/T1), T1, 2, 1e-2, true, 6);
    assert(iJ.completa && iJ.cicloIni == 0, 'xi=0.2 deveria plotar a resposta inteira');
    [~,~,iJ] = graficos().janela_plot(tJ, exp(-2*pi*0.1*tJ/T1).*sin(2*pi*tJ/T1), T1, 2, 1e-2, true, 6);
    assert(~iJ.completa && iJ.decaiu, 'xi=0.1 deveria ser janela ancorada no decaimento');
    [~,~,iJ] = graficos().janela_plot(tJ, exp(-2*pi*0.01*tJ/T1).*sin(2*pi*tJ/T1), T1, 2, 1e-2, true, 6);
    assert(~iJ.decaiu, 'xi=0.01 nao decai em 50 ciclos (ultimo ponto e um zero: nao pode dar falso decaimento)');
    fprintf('OK viga_janela_plot: resposta inteira / janela / sem decaimento\n');

    % --- viga_salva_figura: grava eps + png, confere, e AVISA quando falha ---
    fTmp = figure('Visible','off'); plot(1:3);
    dTmp = fullfile(spDir,'salva');
    assert(graficos().salva_figura(fTmp, dTmp, 'fig_teste'), 'graficos().salva_figura devia retornar true');
    for ext = {'.eps','.png'}
        dd = dir(fullfile(dTmp, ['fig_teste' ext{1}]));
        assert(~isempty(dd) && dd.bytes > 0, ['faltou fig_teste' ext{1}]);
    end
    % regravar por cima (arquivo ja existente) tambem tem que dar certo e atualizar a data
    d0 = dir(fullfile(dTmp,'fig_teste.png')).datenum; pause(2.1);
    assert(graficos().salva_figura(fTmp, dTmp, 'fig_teste'), 'regravar por cima falhou');
    assert(dir(fullfile(dTmp,'fig_teste.png')).datenum > d0, 'o .png antigo nao foi atualizado');
    % handle invalido -> warning + false (nao derruba a rodada)
    close(fTmp);
    lastwarn('');
    okInv = graficos().salva_figura(fTmp, dTmp, 'fig_invalida');
    [~, wid] = lastwarn;
    assert(~okInv && strcmp(wid,'viga_salva_figura:falhou'), 'handle invalido deveria dar warning e false');
    fprintf('OK viga_salva_figura: eps+png, sobrescrita e aviso de falha\n');

    % --- viga_pasta_saida: independe do pwd ---
    pastaA = graficos().pasta_saida('EstudoB');
    pwd0 = pwd; cd(tempdir); pastaB = graficos().pasta_saida('EstudoB'); cd(pwd0);
    assert(strcmp(pastaA, pastaB), 'graficos().pasta_saida depende do pwd');
    assert(strcmp(fileparts(fileparts(pastaA)), fileparts(fileparts(mfilename('fullpath')))), 'saida fora da raiz do projeto');
    fprintf('OK viga_pasta_saida: %s (igual de qualquer pwd)\n', pastaA);

    % --- metodos().pead_tabela ---
    allNames = {'DC','Houbolt','Newmark','Wilson','Bathe','RK1','RK2','RK3','RK4'};
    allLabels = allNames;
    allSteps = repmat(0.01*T1, 1, numel(allNames));
    pead = metodos().pead_tabela(mdl.w, mdl.damp.xi, allNames, allLabels, allSteps, T1, 'nModesReport', 3);
    assert(isequal(size(pead.PE), [numel(allNames), numel(mdl.w)]), 'tamanho errado de PE');
    fprintf('OK viga_pead_tabela: PE(Newmark,modo1)=%.4f, AD(Houbolt,modo3)=%.4f\n', ...
        pead.PE(3,1), pead.AD(2,3));

    % --- graficos().resposta_tempo (poucos periodos, com zoom) ---
    cfgT = struct();
    cfgT.K = mdl.K; cfgT.M = mdl.M; cfgT.C = mdl.C;
    cfgT.u0 = mdl.u0; cfgT.v0 = mdl.v0;
    cfgT.Mref = mdl.M; cfgT.wAll = mdl.w; cfgT.PhiAll = mdl.Phi; cfgT.xiAll = mdl.damp.xi;
    cfgT.tipExtract = @(U) U(mdl.info.tipW,:);
    cfgT.T1 = T1; cfgT.hImpl = 0.01*T1; cfgT.hRK = 0.9*metodos().h_critico(mdl.K,mdl.M,mdl.C,'RK4');
    cfgT.nPerImpl = 8; cfgT.nPerRK = 2; cfgT.X_ciclos_final_plot = 2;
    cfgT.tolDecaimento = 1e-2; cfgT.ciclosDecaimentoRapido = 6;
    cfgT.YCLIP = 1.6; cfgT.label = 'smoketest'; cfgT.xiLabel = 'smoke test';
    cfgT.outDir = spDir; cfgT.doZoom = true;
    cfgT.zoomCycles = [0.5 2];

    tresp = graficos().resposta_tempo(cfgT);
    assert(isfield(tresp,'DC') && isfield(tresp,'RK4'), 'graficos().resposta_tempo faltando metodos');
    assert(isfield(tresp.DC,'err'), 'campo err faltando');
    assert(tresp.config.janelaDireto.completa, 'xi=0.2 (decai em ~3.8 ciclos) deveria plotar a resposta inteira');
    assert(~tresp.config.janelaRK.completa, 'RK com 2 ciclos calculados nao chega a decair');
    assert(all(isfinite([tresp.Newmark.err tresp.RK4.err])), 'err invalido');
    for arq = {'viga_desloc_direto_smoketest','viga_desloc_rk_smoketest'}
        for ext = {'.eps','.png'}
            assert(exist(fullfile(cfgT.outDir,[arq{1} ext{1}]),'file')==2, ['faltou ' arq{1} ext{1}]);
        end
    end
    fprintf('OK viga_resposta_tempo: err(Newmark)=%.4e, N_impl=%d, N_RK=%d\n', ...
        tresp.Newmark.err, tresp.config.N_impl, tresp.config.N_RK);

    % --- Estudo F: um integrador roda no modelo reduzido ---
    [t,Q] = metodos().mdof_bathe(red4.K, red4.M, red4.C, red4.u0, red4.v0, 1e-4, 50);
    u_tip_red = red4.tipRow * Q;
    assert(numel(u_tip_red) == 51, 'tamanho errado da resposta reduzida');
    fprintf('OK integrador no modelo reduzido: max|u_tip_red| = %.4e\n', max(abs(u_tip_red)));

    fprintf('\n===== TODOS OS TESTES DE FUMACA PASSARAM =====\n');
end
