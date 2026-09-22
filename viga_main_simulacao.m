%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  STEP 2 - Results of the cantilever  %
%  beam impulse response example.      %
%--------------------------------------%
%
% Model problem: uniform steel cantilever, clamped at x = 0, subjected to
% an impulsive transverse velocity at the free end.  The response is read
% at the free end, that is, as far as possible from the support.
%
% Three sets of results only:
%
%   1. Spectral radius of every scheme applied to the FE beam (direct and
%      Runge-Kutta schemes together in ONE figure per damping case).
%   2. Period elongation and amplitude decay, using the definitions of
%      Section 3.3 applied to the decoupled modal equations.
%   3. Tip displacement histories against the exact response of the
%      semi-discrete system.
%
% All displacement axes are centred on the analytical response, so an
% unstable curve simply leaves the frame instead of flattening it.
%
% EVERYTHING computed here is stored in the struct "res" and saved to
% ImgViga/viga_resultados.mat, alongside the .eps figures.  See the
% "SAVE" section at the bottom for the layout of the struct.
%
% Run after viga_main_validacao.m.

clc; clear; close all;
tInicio = now;                       %#ok<TNOW1>

% Saida SEMPRE em <pasta do codigo>/ImgViga, nao em <pwd>/ImgViga. Cada
% figura sai em .eps E .png (viga_salva_figura.m).
outDir = viga_pasta_saida();
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf(' Figuras (eps + png) e .mat em: %s\n\n', outDir);

FS     = 18;
FIGPOS = [0 0 20 14];
YCLIP  = 1.6;                 % displacement axis limit, normalised

% --- Simulated window, in fundamental periods -----------------------
% The implicit schemes run at h/T_1 = 0.01, so a period costs ~100 steps
% and nPerImpl can be raised freely.  The Runge-Kutta schemes are forced
% to h ~ 5e-7 s by their stability limit, so EACH period costs about
% 2.2e5 steps per method per damping case: raising nPerRK is expensive.
nPerImpl = 12;                % periods simulated, implicit schemes
nPerRK   = 12;                 % periods simulated, Runge-Kutta schemes

% --- Plotted window, in fundamental periods -------------------------
% Same three-regime rule as the Estudos B-F (viga_janela_plot.m), so every
% figure in the project - A included - is read the same way:
%
%   1) STRONG damping (exact response decays below tolDecaimento within
%      ciclosDecaimentoRapido cycles of T1): the WHOLE response is shown,
%      from t = 0 to the decay instant - this is xi = 0.2 here, which
%      decays in ~3.7 cycles (< 6), so the figure now shows the initial
%      transient (e.g. the Wilson-Theta overshoot) that the old "last 10
%      cycles" window cut off.
%   2) WEAK damping: n_amostras cycles are shown, anchored on the decay
%      instant (last cycles of the simulation if it never decays).
%   3) NO damping (xi = 0): nothing decays, so n_amostras is simply the
%      last cycles of the run - unchanged from before.
%
% n_amostras plays the role of X_ciclos_final_plot in viga_config_estudos.m
% for regimes 2-3. Example: nPerImpl = 100 with n_amostras = 10 runs 100
% periods and, in regime 2/3, plots periods 90 to 100.
% The normalisation of each figure uses the peak of the exact response
% INSIDE the displayed window, so a strongly decayed tail stays centred
% on the analytical curve instead of collapsing onto the axis.
n_amostras = 10;
n_amostras = max(n_amostras, eps);        % must be positive

% tolDecaimento / ciclosDecaimentoRapido: same meaning and same DEFAULT
% values as cfg.tolDecaimento/cfg.ciclosDecaimentoRapido in
% viga_config_estudos.m (section 2), duplicated here - not read from that
% file - so Estudo A stays self-contained, per the file's own header note
% and section 6 of viga_config_estudos.m. Keep the two in sync by hand if
% you change the defaults there.
tolDecaimento = 1e-2;
ciclosDecaimentoRapido = 6;

% --- Saving ---------------------------------------------------------
% Time histories are stored at the tip degree of freedom.  The RK runs
% can reach millions of steps, so their saved histories are decimated to
% at most nSaveRKMax points (set it to Inf to store every step).
nSaveRKMax   = 2e5;
saveFullField = false;        % true also stores all DOFs of the implicit runs

%% ==================================================================
%  MODEL
%  ==================================================================
p    = viga_params();
mdl  = viga_modelo(p);

K = mdl.K;  M = mdl.M;  C = mdl.C;  Cz = mdl.Cz;
u0 = mdl.u0;  v0 = mdl.v0;
iTip = mdl.info.tipW;
T1   = mdl.T1;

xiZero = zeros(size(mdl.w));

rkNames  = {'RK1','RK2','RK3','RK4'};
dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'};
dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'};
caseTags = {'xi0','xi02'};

[coresDir, coresRK] = viga_cores();   % paleta unica - ver viga_cores.m

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
%  1. SPECTRAL RADIUS OF THE FE BEAM
%  ==================================================================
% MDOF counterpart of Figures 2 and 3.  Because the damping is
% proportional, the spectral radius of the MDOF operator is the maximum
% of the scalar spectral radii over the modes of the mesh (verified in
% viga_main_validacao.m), so the operators of viga_operador.m are used
% directly.
hSweep = logspace(-9, -2, 200);
res.stability.hSweep    = hSweep;
res.stability.hSweep_T1 = hSweep/T1;

for ic = 1:2
    tag = caseTags{ic};
    if ic == 1, ximo = xiZero;      xiL = 0.0;
    else,       ximo = mdl.damp.xi; xiL = p.xi_target;
    end
    res.stability.xi.(tag) = xiL;

    % ONE figure per damping case with all nine schemes: direct schemes as
    % solid lines, Runge-Kutta as dashed lines (same coresDir/coresRK used
    % everywhere else - viga_cores.m already keeps RK1 distinguishable
    % from DC for exactly this combined-figure case).
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
            r = arrayfun(@(hh) rhoBeam(names{k}, mdl.w, ximo, hh), hSweep);
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
    viga_salva_figura(fig, outDir, sprintf('viga_rho_%s',tag));
end

%% ==================================================================
%  2. CRITICAL TIME STEPS
%  ==================================================================
condM = [rkNames, {'DC'}];
hcr0 = zeros(numel(condM),1);
hcrD = zeros(numel(condM),1);

fprintf(' Critical time steps of the conditionally stable schemes:\n');
fprintf('   %-5s %14s %12s %14s %12s\n','meth','h_cr (xi=0)','h_cr/T_1','h_cr (xi=0.2)','h_cr/T_1');
for k = 1:numel(condM)
    hcr0(k) = viga_h_critico(K, M, Cz, condM{k});
    hcrD(k) = viga_h_critico(K, M, C , condM{k});
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

% Time steps used from here on
hRat  = 0.01;                    % implicit schemes, as a fraction of T1
hImpl = hRat*T1;
hRK   = 0.90*min(hcrD(1:4));     % explicit schemes, below the damped limit

%% ==================================================================
%  3. PERIOD ELONGATION AND AMPLITUDE DECAY
%  ==================================================================
% Section 3.3 definitions applied to the modes of the beam.  Undamped:
% PE and AD in percent.  Damped: PE and AD as ratios, with 1 meaning the
% physical period / physical decay is reproduced exactly.
%
% Each family is evaluated at the step it is actually run with: the
% implicit schemes at h/T_1 = 0.01, the Runge-Kutta schemes at
% h = 0.90*h_cr.  Evaluating RK at the implicit step would only return
% NaN, since the operator is then far outside its oscillatory range.
%
% The printed table covers the first three modes; the stored struct
% covers EVERY mode of the mesh.
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
                viga_pe_ad(allNames{k}, mdl.w(m), ximo(m), hk);
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

    % Readable table of the first three modes
    res.pead.table.(tag) = cell2table( ...
        [allLabels(:), num2cell(allSteps(:)/T1), ...
         num2cell(PEall(:,1)), num2cell(ADall(:,1)), ...
         num2cell(PEall(:,2)), num2cell(ADall(:,2)), ...
         num2cell(PEall(:,3)), num2cell(ADall(:,3))], ...
        'VariableNames', {'method','h_T1','PE_m1','AD_m1','PE_m2','AD_m2','PE_m3','AD_m3'});
end

%% ==================================================================
%  4. TIP DISPLACEMENT - UNCONDITIONALLY STABLE SCHEMES
%  ==================================================================
nS = round(nPerImpl*T1/hImpl);

for ic = 1:2
    tag = caseTags{ic};
    if ic == 1, Cm = Cz;  ximo = xiZero;      xiL = 0.0;         hcrDC = hcr0(end);
    else,       Cm = C;   ximo = mdl.damp.xi; xiL = p.xi_target; hcrDC = hcrD(end);
    end
    amortecido = xiL > 0;

    tRef = (0:nS)*hImpl;
    [Uex, Vex, Aex] = viga_solucao_exata(tRef, M, mdl.w, mdl.Phi, ximo, u0, v0);

    res.time.(tag).exact.t     = tRef;
    res.time.(tag).exact.u_tip = Uex(iTip,:);
    res.time.(tag).exact.v_tip = Vex(iTip,:);
    res.time.(tag).exact.a_tip = Aex(iTip,:);
    if saveFullField
        res.time.(tag).exact.U = Uex;
    end

    % Janela mostrada: as mesmas TRES regras dos Estudos B-F (resposta
    % inteira / janela ancorada no decaimento / fim da simulacao, ver
    % viga_janela_plot.m). So o que e DESENHADO muda - nS e a simulacao em
    % si continuam exatamente os de sempre.
    [tIni, tFim, infoJan] = viga_janela_plot(tRef, Uex(iTip,:), T1, ...
        n_amostras, tolDecaimento, amortecido, ciclosDecaimentoRapido);
    res.window.(tag).direto = infoJan;
    viga_print_janela('direct schemes', infoJan, T1, nPerImpl);

    mask   = (tRef >= tIni - 0.5*hImpl) & (tRef <= tFim + 0.5*hImpl);
    uScale = max(abs(Uex(iTip,mask)));
    if uScale == 0, uScale = 1; end

    % Central Differences fica fora da figura so quando de fato instavel
    % neste passo (generaliza o skip fixo de antes, que so valia porque DC
    % sempre excedia o proprio limite a h/T1 = 0.01 - mesma regra usada
    % pelos Estudos B-F em viga_resposta_tempo.m).
    skipDC = (hcrDC > 0) && (hImpl > hcrDC);

    fig = figure('Name',sprintf('Tip disp direct %s',tag),'NumberTitle','off');
    plot(tRef(mask)/T1, Uex(iTip,mask)/uScale, 'k--', 'LineWidth', 2); hold on;
    leg = {'Exact'};

    for k = 1:numel(dirNames)
        [t,U,V,A] = run_scheme(dirNames{k}, K, M, Cm, u0, v0, hImpl, nS);

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
    viga_salva_figura(fig, outDir, sprintf('viga_desloc_direto_%s',tag));
end

%% ==================================================================
%  5. TIP DISPLACEMENT - RUNGE-KUTTA SCHEMES
%  ==================================================================
% A single time step is used for all four schemes: 90 percent of the
% smallest damped critical step, so that the comparison isolates the
% effect of physical damping.  With xi = 0 the same step is above the
% (nonexistent) limit of RK1 and RK2 and they leave the frame; with
% xi = 0.2 all four are stable.  This is the time-domain counterpart of
% the central finding of Section 4.1.
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

    % Janela mostrada: mesma regra de tres regimes usada acima para os
    % metodos diretos (viga_janela_plot.m). A deteccao usa uma varredura
    % grosseira sobre a simulacao INTEIRA (barata, forma fechada); so
    % DEPOIS a grade fina de 4000 pontos e montada sobre a janela achada -
    % mesmo padrao de viga_resposta_tempo.m.
    tVarre = linspace(0, nPerRK*T1, 2000);
    UexV   = viga_solucao_exata(tVarre, M, mdl.w, mdl.Phi, ximo, u0, v0);
    [tIniRK, tFimRK, infoJanRK] = viga_janela_plot(tVarre, UexV(iTip,:), T1, ...
        n_amostras, tolDecaimento, amortecido, ciclosDecaimentoRapido);
    res.window.(tag).rk = infoJanRK;
    viga_print_janela('Runge-Kutta', infoJanRK, T1, nPerRK);

    tFine  = linspace(tIniRK, tFimRK, 4000);
    UexF   = viga_solucao_exata(tFine, M, mdl.w, mdl.Phi, ximo, u0, v0);
    uScale = max(abs(UexF(iTip,:)));
    if uScale == 0, uScale = 1; end

    res.timeRK.(tag).exact.t     = tFine;
    res.timeRK.(tag).exact.u_tip = UexF(iTip,:);

    fig = figure('Name',sprintf('Tip disp RK %s',tag),'NumberTitle','off');
    plot(tFine/T1, UexF(iTip,:)/uScale, 'k--', 'LineWidth', 2); hold on;
    leg = {'Exact'};

    for k = 1:numel(rkNames)
        [t,U,V,A] = mdof_rk(K, M, Cm, u0, v0, hRK, nRK, k);

        % Stored histories are decimated to keep the .mat manageable
        [ts, us] = decimate_hist(t, U(iTip,:), nSaveRKMax);
        [~,  vs] = decimate_hist(t, V(iTip,:), nSaveRKMax);
        [~,  as] = decimate_hist(t, A(iTip,:), nSaveRKMax);
        res.timeRK.(tag).(rkNames{k}).t          = ts;
        res.timeRK.(tag).(rkNames{k}).h          = hRK;
        res.timeRK.(tag).(rkNames{k}).nSteps     = nRK;
        res.timeRK.(tag).(rkNames{k}).u_tip      = us;
        res.timeRK.(tag).(rkNames{k}).v_tip      = vs;
        res.timeRK.(tag).(rkNames{k}).a_tip      = as;
        res.timeRK.(tag).(rkNames{k}).decimated  = numel(ts) < numel(t);
        res.timeRK.(tag).(rkNames{k}).u_tip_max  = max(abs(U(iTip,:)));

        mask    = (t >= tIniRK - 0.5*hRK) & (t <= tFimRK + 0.5*hRK);
        [tp,yp] = downsample_plot(t(mask), U(iTip,mask)/uScale);
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
    viga_salva_figura(fig, outDir, sprintf('viga_desloc_rk_%s',tag));
end

%% ==================================================================
%  SAVE
%  ==================================================================
% Layout of "res":
%
%   res.model                 output of viga_modelo(): K, M, C, w, Phi,
%                             damping data, initial conditions, T1, wmax
%   res.config                every parameter used in this run
%   res.stability.hSweep      time steps of the spectral radius sweep
%   res.stability.rho.<case>.<method>   spectral radius curves
%   res.stability.hcr_xi0 / .hcr_xi02   critical steps
%   res.stability.hcr_table   the same, as a MATLAB table
%   res.pead.PE.<case>        [nMethods x nModes] period elongation
%   res.pead.AD.<case>        [nMethods x nModes] amplitude decay
%   res.pead.rho.<case>       [nMethods x nModes] scalar spectral radii
%   res.pead.table.<case>     first three modes, as a MATLAB table
%   res.time.<case>.<method>  implicit runs: t, h, u_tip, v_tip, a_tip
%   res.time.<case>.exact     exact reference on the same time grid
%   res.timeRK.<case>.<method>  Runge-Kutta runs (decimated histories)
%   res.window.<case>.direto/.rk  window diagnostic from viga_janela_plot.m
%                             (which regime was used, decay instant, etc.)
%                             - purely informative, nothing else reads it
%
% <case> is 'xi0' or 'xi02'; <method> is DC, Houbolt, Newmark, Wilson,
% Bathe, RK1, RK2, RK3 or RK4.

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
viga_lista_obsoletos({outDir}, tInicio);
fprintf('=====================================================\n');

%% ==================================================================
%  LOCAL FUNCTIONS
%  ==================================================================
function r = rhoBeam(method, w, xi, h)
% Spectral radius of the MDOF operator: maximum of the scalar spectral
% radii over the modes of the mesh (exact for proportional damping).
    r = 0;
    for k = 1:numel(w)
        rk = max(abs(eig(viga_operador(method, w(k), xi(k), h))));
        if rk > r, r = rk; end
    end
end

function [t, U, V, A] = run_scheme(name, K, M, C, u0, v0, h, nSteps)
    switch name
        case 'DC',      [t,U,V,A] = mdof_dc     (K,M,C,u0,v0,h,nSteps);
        case 'Houbolt', [t,U,V,A] = mdof_houbolt(K,M,C,u0,v0,h,nSteps);
        case 'Newmark', [t,U,V,A] = mdof_newmark(K,M,C,u0,v0,h,nSteps);
        case 'Wilson',  [t,U,V,A] = mdof_wilson (K,M,C,u0,v0,h,nSteps);
        case 'Bathe',   [t,U,V,A] = mdof_bathe  (K,M,C,u0,v0,h,nSteps);
        otherwise, error('run_scheme:name','Unknown method "%s".', name);
    end
end

function [ts, ys] = decimate_hist(t, y, nMax)
% Uniform decimation used when STORING long histories.  Keeps the first
% and last samples so the run length is preserved.
    if isinf(nMax) || numel(t) <= nMax
        ts = t;  ys = y;  return
    end
    idx = unique([1, round(linspace(1, numel(t), nMax)), numel(t)]);
    ts = t(idx);  ys = y(idx);
end

function [ts, ys] = downsample_plot(t, y, nMax)
% Keeps the FIGURES light when a conditionally stable scheme needs
% hundreds of thousands of steps.  Preserves the first and last points.
    if nargin < 3, nMax = 5000; end
    if numel(t) <= nMax
        ts = t;  ys = y;  return
    end
    idx = unique([1, round(linspace(1, numel(t), nMax)), numel(t)]);
    ts = t(idx);  ys = y(idx);
end
