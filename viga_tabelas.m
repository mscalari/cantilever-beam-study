%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  Consolidates the results of every   %
%  Estudo into export-ready tables     %
%  (item 36 of the advisor's list).    %
%--------------------------------------%
%
% Reads the .mat files saved by viga_main_simulacao.m (Estudo A) and by
% viga_estudo_{B,C,D,E,F}_*.m (run viga_main_estudos.m first, or run the
% scripts individually) and writes a set of CSV tables to
% ImgViga/Tabelas/. Any study not yet run is skipped with a warning
% instead of stopping the script, so this can be re-run at any point
% while the batteries are still being completed.
%
% Produces:
%   tabela_frequencias.csv     w_1..w_10, T_1..T_10 (full model)
%   tabela_amortecimento.csv   xi por caso (Estudo B) + Rayleigh alt (E)
%   tabela_passos_criticos.csv h_cr por metodo, por caso de xi/modelo
%   tabela_N_erro.csv          N = t_final/h e erro relativo, por caso

outDir  = viga_pasta_saida();       % <pasta do codigo>/ImgViga, nao <pwd>/ImgViga
tabDir  = fullfile(outDir,'Tabelas');
if ~exist(tabDir,'dir'), mkdir(tabDir); end

p = viga_params();

% ======================================================================
%  1. FREQUENCIAS NATURAIS
% ======================================================================
[K, M, info] = viga_matrizes_EF(p);
w = viga_autovalores_EF(K, M);
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
    warning('viga_tabelas:missing', 'Estudo B nao encontrado (%s) - pulando xi=0/0.01/0.2 na tabela de amortecimento.', fB);
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
    warning('viga_tabelas:missing', 'Estudo E nao encontrado (%s) - pulando Rayleigh alternativo na tabela de amortecimento.', fE);
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
    warning('viga_tabelas:missing', 'Estudo A nao encontrado (%s) - rode viga_main_simulacao.m primeiro.', fA);
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
    warning('viga_tabelas:missing', 'Estudo F nao encontrado (%s).', ffile);
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

% ==========================================================================
%  LOCAL FUNCTIONS
% ==========================================================================
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
