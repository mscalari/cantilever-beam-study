%% Smoke test - novas funcoes do plano do orientador (rapido, nao gera .eps)
clc; clear; close all;

p = viga_params();
[K, M, info] = viga_matrizes_EF(p);
[w, Phi] = viga_autovalores_EF(K, M);
T1 = 2*pi/w(1);
fprintf('OK viga_params/matrizes/autovalores: w1=%.4f, T1=%.6f, ndof=%d\n', w(1), T1, info.ndof);

% --- viga_amortecimento: modo anexado (original) e modo direto (novo) ---
[Ca, dampA] = viga_amortecimento(K, M, w, p);
assert(isempty(setdiff({'a0','a1','xi','anchors','xi_target'}, fieldnames(dampA))), 'campos faltando (modo anexado)');
[Cd, dampD] = viga_amortecimento(K, M, w, p, 0, 1e-4);
assert(isequal(dampD.a0,0) && isequal(dampD.a1,1e-4), 'coeficientes diretos nao aplicados');
assert(isempty(dampD.anchors), 'anchors deveria ser vazio no modo direto');
fprintf('OK viga_amortecimento (2 modos): xi1_anexado=%.4f, xi1_direto=%.4f\n', dampA.xi(1), dampD.xi(1));

% --- viga_ci_modal ---
[u0m, v0m] = viga_ci_modal(K, M, Ca, Phi, 3);
assert(norm(v0m) == 0, 'v0 deveria ser zero');
q0check = Phi.' * (M * u0m);
[~, idxmax] = max(abs(q0check));
assert(idxmax == 3, 'CI modal nao projeta no modo esperado');
fprintf('OK viga_ci_modal: q0 concentrado no modo %d (valor %.4f, resto <= %.2e)\n', ...
    idxmax, q0check(idxmax), max(abs(q0check([1 2 4:end]))));

% --- viga_reducao_modal ---
mdl = viga_modelo(p);
red4 = viga_reducao_modal(mdl, 4);
assert(isequal(size(red4.K),[4 4]) && isequal(size(red4.M),[4 4]), 'tamanho errado do modelo reduzido');
assert(all(abs(diag(red4.M)-1) < 1e-10), 'M reduzido deveria ser identidade');
assert(red4.wmax == mdl.w(4), 'wmax reduzido incorreto');
fprintf('OK viga_reducao_modal: wmax_full=%.3e, wmax_red4=%.3e (%.1fx menor)\n', ...
    mdl.wmax, red4.wmax, mdl.wmax/red4.wmax);

% --- viga_rho_beam / viga_espectro (sem plot, sweep pequeno) ---
r = viga_rho_beam('RK4', mdl.w, mdl.damp.xi, 1e-6);
assert(isfinite(r) && r >= 0, 'rho invalido');
stab = viga_espectro(mdl.w, mdl.damp.xi, '\xi=0.2', 'smoketest', pwd, ...
    'hSweep', logspace(-7,-4,5), 'doPlots', false);
assert(isfield(stab.rho.rk,'RK4') && numel(stab.rho.rk.RK4)==5, 'viga_espectro nao retornou o esperado');
fprintf('OK viga_rho_beam/viga_espectro: rho(RK4,h=1e-6)=%.6f\n', r);

% --- viga_espectro com plot: UMA figura com os 9 metodos (diretos + RK) ---
spDir = fullfile(pwd,'smoketest_out');
viga_espectro(mdl.w, mdl.damp.xi, '\xi=0.2', 'smoketest', spDir, ...
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
[~,~,iJ] = viga_janela_plot(tJ, exp(-2*pi*0.2*tJ/T1).*sin(2*pi*tJ/T1), T1, 2, 1e-2, true, 6);
assert(iJ.completa && iJ.cicloIni == 0, 'xi=0.2 deveria plotar a resposta inteira');
[~,~,iJ] = viga_janela_plot(tJ, exp(-2*pi*0.1*tJ/T1).*sin(2*pi*tJ/T1), T1, 2, 1e-2, true, 6);
assert(~iJ.completa && iJ.decaiu, 'xi=0.1 deveria ser janela ancorada no decaimento');
[~,~,iJ] = viga_janela_plot(tJ, exp(-2*pi*0.01*tJ/T1).*sin(2*pi*tJ/T1), T1, 2, 1e-2, true, 6);
assert(~iJ.decaiu, 'xi=0.01 nao decai em 50 ciclos (ultimo ponto e um zero: nao pode dar falso decaimento)');
fprintf('OK viga_janela_plot: resposta inteira / janela / sem decaimento\n');

% --- viga_salva_figura: grava eps + png, confere, e AVISA quando falha ---
fTmp = figure('Visible','off'); plot(1:3);
dTmp = fullfile(pwd,'smoketest_out','salva');
assert(viga_salva_figura(fTmp, dTmp, 'fig_teste'), 'viga_salva_figura devia retornar true');
for ext = {'.eps','.png'}
    dd = dir(fullfile(dTmp, ['fig_teste' ext{1}]));
    assert(~isempty(dd) && dd.bytes > 0, ['faltou fig_teste' ext{1}]);
end
% regravar por cima (arquivo ja existente) tambem tem que dar certo e atualizar a data
d0 = dir(fullfile(dTmp,'fig_teste.png')).datenum; pause(2.1);
assert(viga_salva_figura(fTmp, dTmp, 'fig_teste'), 'regravar por cima falhou');
assert(dir(fullfile(dTmp,'fig_teste.png')).datenum > d0, 'o .png antigo nao foi atualizado');
% handle invalido -> warning + false (nao derruba a rodada)
close(fTmp);
lastwarn('');
okInv = viga_salva_figura(fTmp, dTmp, 'fig_invalida');
[~, wid] = lastwarn;
assert(~okInv && strcmp(wid,'viga_salva_figura:falhou'), 'handle invalido deveria dar warning e false');
fprintf('OK viga_salva_figura: eps+png, sobrescrita e aviso de falha\n');

% --- viga_pasta_saida: independe do pwd ---
pastaA = viga_pasta_saida('EstudoB');
pwd0 = pwd; cd(tempdir); pastaB = viga_pasta_saida('EstudoB'); cd(pwd0);
assert(strcmp(pastaA, pastaB), 'viga_pasta_saida depende do pwd');
assert(strcmp(fileparts(fileparts(pastaA)), fileparts(which('viga_params'))), 'saida fora da pasta do codigo');
fprintf('OK viga_pasta_saida: %s (igual de qualquer pwd)\n', pastaA);

% --- viga_pead_tabela ---
allNames = {'DC','Houbolt','Newmark','Wilson','Bathe','RK1','RK2','RK3','RK4'};
allLabels = allNames;
allSteps = repmat(0.01*T1, 1, numel(allNames));
pead = viga_pead_tabela(mdl.w, mdl.damp.xi, allNames, allLabels, allSteps, T1, 'nModesReport', 3);
assert(isequal(size(pead.PE), [numel(allNames), numel(mdl.w)]), 'tamanho errado de PE');
fprintf('OK viga_pead_tabela: PE(Newmark,modo1)=%.4f, AD(Houbolt,modo3)=%.4f\n', ...
    pead.PE(3,1), pead.AD(2,3));

% --- viga_resposta_tempo (poucos periodos, com zoom, sem tocar Estudo A) ---
cfgT = struct();
cfgT.K = mdl.K; cfgT.M = mdl.M; cfgT.C = mdl.C;
cfgT.u0 = mdl.u0; cfgT.v0 = mdl.v0;
cfgT.Mref = mdl.M; cfgT.wAll = mdl.w; cfgT.PhiAll = mdl.Phi; cfgT.xiAll = mdl.damp.xi;
cfgT.tipExtract = @(U) U(mdl.info.tipW,:);
cfgT.T1 = T1; cfgT.hImpl = 0.01*T1; cfgT.hRK = 0.9*viga_h_critico(mdl.K,mdl.M,mdl.C,'RK4');
cfgT.nPerImpl = 8; cfgT.nPerRK = 2; cfgT.X_ciclos_final_plot = 2;
cfgT.tolDecaimento = 1e-2; cfgT.ciclosDecaimentoRapido = 6;
cfgT.YCLIP = 1.6; cfgT.label = 'smoketest'; cfgT.xiLabel = 'smoke test';
cfgT.outDir = fullfile(pwd,'smoketest_out'); cfgT.doZoom = true;
cfgT.zoomCycles = [0.5 2];

tresp = viga_resposta_tempo(cfgT);
assert(isfield(tresp,'DC') && isfield(tresp,'RK4'), 'viga_resposta_tempo faltando metodos');
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

% --- Estudo F: modelo reduzido de fato roda um integrador ---
[t,Q] = mdof_bathe(red4.K, red4.M, red4.C, red4.u0, red4.v0, 1e-4, 50);
u_tip_red = red4.tipRow * Q;
assert(numel(u_tip_red) == 51, 'tamanho errado da resposta reduzida');
fprintf('OK integrador no modelo reduzido: max|u_tip_red| = %.4e\n', max(abs(u_tip_red)));

fprintf('\n===== TODOS OS TESTES DE FUMACA PASSARAM =====\n');
