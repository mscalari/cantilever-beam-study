%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  TESTE RAPIDO da janela de plot -    %
%  rode ISTO antes dos Estudos B-F.    %
%--------------------------------------%
%
% PARA QUE SERVE
% Conferir, em segundos, se a janela de plot esta caindo onde voce quer,
% SEM pagar o custo de uma bateria completa. Roda so a familia de metodos
% DIRETOS (implicitos), que a hRat = 0.01 custa ~100 passos por ciclo -
% isto e, o script inteiro roda em alguns segundos, contra os minutos de
% um estudo com Runge-Kutta.
%
% O QUE ELE FAZ
% Monta tres casos com a mesma viga, um para cada ramo da logica:
%   1. MUITO AMORTECIDO (xi = p.xi_target, tipicamente 0.2) -> a resposta
%      decai em poucos ciclos (<= cfg.ciclosDecaimentoRapido), entao a
%      figura mostra a resposta INTEIRA, de t = 0 ate o decaimento.
%   2. POUCO AMORTECIDO (xi = cfg.C.xi, tipicamente 0.01)   -> decai
%      devagar. Se morre dentro dos ciclos calculados, plota os
%      X_ciclos_final_plot ciclos anteriores ao decaimento; se nao morre,
%      cai no fallback (plota o fim) e o console avisa.
%   3. NAO AMORTECIDO   (xi = 0)                            -> nao ha o
%      que detectar; plota os X_ciclos_final_plot_sem_xi ciclos finais.
% Para cada um, imprime a janela detectada e gera a figura.
%
% QUANTOS CICLOS PRECISO PARA A RESPOSTA MORRER?
% A envoltoria do modo 1 cai como exp(-xi*2*pi*n) depois de n ciclos,
% entao para chegar a tolDecaimento:
%
%       n = -ln(tolDecaimento) / (2*pi*xi)
%
% Com tol = 1e-2: xi = 0.2 morre em ~3.7 ciclos, xi = 0.05 em ~15, e
% xi = 0.01 so em ~73. (Com tol = 1e-3: ~5.5, ~22 e ~110.) O script
% imprime essa estimativa para cada caso, e diz em qual regime ele cai,
% para voce ver se nPerImpl esta grande o bastante para a deteccao valer
% - se nao estiver, o codigo avisa e plota o fim da simulacao.
%
% COMO USAR
%   1. Edite os valores em viga_config_estudos.m (secao 2).
%   2. Rode este script.
%   3. Olhe as figuras em ImgViga/TesteJanela/ e as linhas "[janela ...]"
%      no console.
%   4. Repita ate gostar. So entao rode viga_main_estudos.m.
%
% Para varrer varios valores de uma vez sem editar a configuracao a cada
% tentativa, use o bloco VARREDURA OPCIONAL no fim deste arquivo.
%
% Este script NAO escreve nada em ImgViga/EstudoX/ nem mexe nos .mat dos
% estudos - a saida dele fica isolada em ImgViga/TesteJanela/.

clc; clear; close all;

tInicio = now;                       %#ok<TNOW1>
outDir = viga_pasta_saida('TesteJanela');   % <pasta do codigo>/ImgViga/TesteJanela
if ~exist(outDir,'dir'), mkdir(outDir); end

cfg = viga_config_estudos();
p   = viga_params();

[K, M, info] = viga_matrizes_EF(p);
[w, Phi]     = viga_autovalores_EF(K, M);
T1   = 2*pi/w(1);
iTip = info.tipW;

fprintf('#############################################################\n');
fprintf('  TESTE DA JANELA DE PLOT (so metodos diretos - rapido)\n');
fprintf('#############################################################\n');
fprintf(' T1 = %.6f s\n', T1);
fprintf(' Configuracao atual (viga_config_estudos.m, secao 2):\n');
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
        [C, damp] = viga_amortecimento(K, M, w, pp);
        ximo = damp.xi;
    end

    [u0, v0] = viga_ci_impulso(K, M, C, info, p);

    cfgT = struct();
    cfgT.K = K; cfgT.M = M; cfgT.C = C;
    cfgT.u0 = u0; cfgT.v0 = v0;
    cfgT.Mref = M; cfgT.wAll = w; cfgT.PhiAll = Phi; cfgT.xiAll = ximo;
    cfgT.tipExtract = @(U) U(iTip,:);
    cfgT.T1 = T1; cfgT.hImpl = hImpl; cfgT.hRK = hImpl;   % hRK nao e usado
    cfgT.label = ['teste_' cs.tag];
    cfgT.xiLabel = cs.label;
    cfgT.outDir = outDir;
    cfgT = viga_aplica_config(cfgT, cfg);
    cfgT.runRK = false;    % <<< e isto que torna o teste barato

    tresp = viga_resposta_tempo(cfgT);
    res.(cs.tag) = tresp;

    % viga_resposta_tempo ja gravou viga_desloc_direto_teste_<caso>.eps e
    % .png (300 dpi) em outDir - o PNG e a copia para conferir na tela.
    % (Antes havia aqui um laco que reimprimia TODAS as figuras abertas
    % (findobj) num unico preview_<caso>.png; como as figuras nao eram
    % fechadas entre os casos, os tres previews saiam identicos, com a
    % figura do primeiro caso.) Fecha tudo para nada vazar para o proximo.
    close all;

    j = tresp.config.janelaDireto;
    txtRegime = '';
    if j.completa, txtRegime = ', resposta inteira'; end
    fprintf('   -> janela final: t/T1 em [%.2f, %.2f]  (%.2f ciclos%s)\n\n', ...
        j.cicloIni, j.cicloFim, j.cicloFim - j.cicloIni, txtRegime);
end

fprintf('#############################################################\n');
fprintf('  Figuras (eps + png) em %s\n', outDir);
viga_lista_obsoletos({outDir}, tInicio);
fprintf('  Se as janelas estiverem boas, rode viga_main_estudos.m.\n');
fprintf('#############################################################\n');

% ======================================================================
%  VARREDURA OPCIONAL - varios X_ciclos_final_plot de uma vez
% ======================================================================
% Descomente o bloco abaixo para gerar, numa so rodada, a figura do caso
% amortecido para varios numeros de ciclos plotados, e escolher olhando.
% Nao precisa editar viga_config_estudos.m entre as tentativas.
%
% ciclosTeste = [2 5 10];
% pp = p; pp.xi_target = cfg.C.xi;
% [Cv, dv] = viga_amortecimento(K, M, w, pp);
% [u0v, v0v] = viga_ci_impulso(K, M, Cv, info, p);
% for nc = ciclosTeste
%     cfgV = struct();
%     cfgV.K = K; cfgV.M = M; cfgV.C = Cv;
%     cfgV.u0 = u0v; cfgV.v0 = v0v;
%     cfgV.Mref = M; cfgV.wAll = w; cfgV.PhiAll = Phi; cfgV.xiAll = dv.xi;
%     cfgV.tipExtract = @(U) U(iTip,:);
%     cfgV.T1 = T1; cfgV.hImpl = hImpl; cfgV.hRK = hImpl;
%     cfgV.label = sprintf('teste_ciclos%g', nc);
%     cfgV.xiLabel = sprintf('\\xi = %.3g, %g cycles', cfg.C.xi, nc);
%     cfgV.outDir = outDir;
%     cfgV = viga_aplica_config(cfgV, cfg, 'X_ciclos_final_plot', nc);
%     cfgV.runRK = false;
%     fprintf('--- Varredura: X_ciclos_final_plot = %g ---\n', nc);
%     viga_resposta_tempo(cfgV);
% end
