%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  Gera um report markdown com TODOS   %
%  os parametros de cada estudo.       %
%--------------------------------------%
%
% Le os .mat salvos por viga_main_simulacao.m (Estudo A) e pelos
% viga_estudo_{B,C,D,E,F}_*.m e escreve
%
%     ImgViga/Tabelas/PARAMETROS_ESTUDOS.md
%
% O formato e pensado para ser COLADO NUMA IA como contexto: tudo
% explicito e auto-contido (sem "ver arquivo tal"), uma secao por estudo,
% pares chave: valor com unidade, e um glossario no fim explicando o que
% cada simbolo significa. Escrito em ingles, como o paper.
%
% Estudos ainda nao rodados sao pulados com um aviso, entao da para rodar
% isto a qualquer momento enquanto as baterias ainda estao saindo.
%
% Rode DEPOIS de viga_main_estudos.m (e de viga_main_simulacao.m, se
% quiser o Estudo A no report).

clc; clear;

outDir = viga_pasta_saida();        % <pasta do codigo>/ImgViga, nao <pwd>/ImgViga
tabDir = fullfile(outDir,'Tabelas');
if ~exist(tabDir,'dir'), mkdir(tabDir); end
mdFile = fullfile(tabDir,'PARAMETROS_ESTUDOS.md');

cfg = viga_config_estudos();
p   = viga_params();
[K, M, info] = viga_matrizes_EF(p);
w = viga_autovalores_EF(K, M);
T1 = 2*pi/w(1);

fid = fopen(mdFile,'w');
if fid < 0
    error('viga_report_parametros:fopen','Nao consegui escrever em %s', mdFile);
end
closer = onCleanup(@() fclose(fid));

% ======================================================================
%  CABECALHO
% ======================================================================
w_(fid, '# Simulation parameters — cantilever beam time-integration study');
w_(fid, '');
w_(fid, 'Machine-readable summary of every parameter used in the beam study.');
w_(fid, 'Generated automatically by `viga_report_parametros.m` on %s.', datestr(now,'yyyy-mm-dd HH:MM'));
w_(fid, '');
w_(fid, 'Nine time-integration schemes are applied to a damped cantilever beam');
w_(fid, 'discretised by Euler-Bernoulli finite elements. Five schemes are *direct*');
w_(fid, '(Central Differences, Houbolt, Newmark, Wilson-Theta, Bathe) and four are');
w_(fid, 'explicit *Runge-Kutta* (RK1-RK4). Each study below varies one factor.');
w_(fid, '');

% ======================================================================
%  MODELO FISICO
% ======================================================================
w_(fid, '## 1. Physical model (`viga_params.m`)');
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
w_(fid, '## 2. Numerical configuration (`viga_config_estudos.m`)');
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
w_(fid, 'Source: `viga_main_simulacao.m`. Kept untouched as a regression reference;');
w_(fid, 'it does NOT read `viga_config_estudos.m`.');
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
     'viga_estudo_B_amortecimento.m', ...
     'viga_estudo_C_modal_h_fixo.m', ...
     'viga_estudo_D_modal_h_variavel.m', ...
     'viga_estudo_E_rayleigh_alt.m', ...
     'viga_estudo_F_reducao_modal.m'}, ...
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

    % Cases swept, straight from the configuration
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

    % Results actually on disk
    if ~exist(mf,'file')
        w_(fid, 'Status: **not run** (`%s` missing).', mf);
        w_(fid, '');
        warning('viga_report_parametros:missing', ...
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

    % Derived steps
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

    % Critical steps
    if isfield(res,'hcr')
        w_(fid, '**Critical time steps**');
        w_(fid, '');
        escreve_hcr(fid, res.hcr);
        w_(fid, '');
    end

    % Per-case N and error
    escreve_N_erro(fid, res);

    % Estudo F part 2
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

% ==========================================================================
%  LOCAL FUNCTIONS
% ==========================================================================
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
    metodos = setdiff(fieldnames(tresp), {'exact','exactRK','config','cfg'});
    linhas = {};
    for m = 1:numel(metodos)
        e = tresp.(metodos{m});
        if isstruct(e) && isfield(e,'N') && isfield(e,'h')
            if isfield(e,'err')
                linhas{end+1} = sprintf('| %s | %s | %.6g | %d | %.4e |', ...
                    tag, metodos{m}, e.h, e.N, e.err); %#ok<AGROW>
            else
                linhas{end+1} = sprintf('| %s | %s | %.6g | %d | - |', ...
                    tag, metodos{m}, e.h, e.N); %#ok<AGROW>
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
