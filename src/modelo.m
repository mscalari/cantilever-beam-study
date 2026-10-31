function m = modelo()
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%
% O modelo fisico: parametros da viga, configuracao dos estudos, elementos
% finitos, modos, amortecimento, condicoes iniciais, solucao exata e reducao
% modal.
%
% As funcoes deste arquivo sao chamadas de fora como modelo().nome(...).
% So as listadas abaixo sao publicas; as demais sao de uso interno.

    m.amortecimento    = @amortecimento;
    m.autovalores_EF   = @autovalores_EF;
    m.ci_impulso       = @ci_impulso;
    m.ci_modal         = @ci_modal;
    m.config_estudos   = @config_estudos;
    m.matrizes_EF      = @matrizes_EF;
    m.modos_analiticos = @modos_analiticos;
    m.monta            = @monta;
    m.params           = @params;
    m.reducao_modal    = @reducao_modal;
    m.solucao_exata    = @solucao_exata;
end

%% ====================================================================
%  Parametros da viga e configuracao dos estudos
%  ====================================================================

function p = params()
% Propriedades geometricas e de material da viga, malha e amortecimento.
%
% Viga de aco retangular, engastada em x = 0 e livre em x = L. Com L/h = 100
% o modelo fica bem dentro do regime de Euler-Bernoulli (sem cisalhamento).
%
% Aqui fica O QUE A VIGA E: mudar algo neste bloco muda o problema fisico.
% COMO a simulacao e rodada fica em config_estudos.
%
% SAIDA
%   p : struct com as propriedades da viga
    p.E   = 210e9;      % Modulo de Young            [Pa]
    p.rho = 7850;       % Densidade                  [kg/m^3]

    % --- Geometria ---
    p.L = 1.00;         % Comprimento                [m]
    p.b = 0.050;        % Largura da secao           [m]
    p.h = 0.010;        % Altura da secao            [m]

    % --- Propriedades derivadas da secao ---
    p.A = p.b*p.h;              % Area                        [m^2]
    p.I = p.b*p.h^3/12;         % Momento de inercia          [m^4]
    p.slenderness = p.L/p.h;    % Esbeltez                    [-]

    % --- Discretizacao ---
    p.N = 10;           % Numero de elementos de 2 nos

    % --- Amortecimento ---
    p.xi_target = 0.2;  % Amortecimento modal alvo (o mesmo do estudo de 1 GDL)
    p.mode_i    = 1;    % Primeiro modo ancora de Rayleigh
    p.mode_j    = 3;    % Segundo modo ancora de Rayleigh

    % --- Condicao inicial de impulso ---
    p.V0 = 1.0;         % Velocidade transversal inicial na ponta [m/s]

    % --- Conferencia ---
    if p.slenderness < 10
        warning('viga_params:slenderness', ...
            'L/h = %.1f < 10: Euler-Bernoulli theory may be inaccurate.', ...
            p.slenderness);
    end
end

function cfg = config_estudos()
% Configuracao numerica dos Estudos B-F: passo de tempo, numero de ciclos,
% casos varridos e aparencia das figuras. Cada estudo chama esta funcao uma
% vez e le tudo daqui; nenhum deles tem numero de simulacao fixo no corpo.
%
% O Estudo A e o caso de referencia ja publicado no paper e tem a propria
% configuracao no inicio de estudo_A, para servir de ponto fixo de comparacao.
%
% SAIDA
%   cfg : struct com os campos comentados abaixo
cfg.hRat = 0.01;

% fatorSegurancaRK: fracao do passo critico usada pelos metodos RK
%   (hRK = fatorSegurancaRK * min(h_cr dos RK)).
%   Deve ficar ABAIXO de 1.0 - em 1.0 o esquema esta exatamente no limite
%   de estabilidade e o resultado fica na corda bamba (raio espectral = 1).
%   0.90 da 10% de folga, suficiente para as figuras do paper.
%   ATENCAO AO CUSTO: h_cr e ditado por w_max da MALHA (modo 20, ~9e4
%   rad/s), nao pelo modo 1. Na malha completa isso da hRK ~ 5e-7 s, ou
%   seja ~2e5 passos POR PERIODO POR METODO. Reduzir este fator encarece
%   proporcionalmente.
cfg.fatorSegurancaRK = 0.90;

% ======================================================================
%  2. DURACAO DA SIMULACAO E JANELA MOSTRADA
% ======================================================================
% A simulacao roda nPer* ciclos, mas a figura mostra so uma JANELA de
% poucos ciclos. Isso e proposital: os erros de periodo (PE) e de
% amplitude (AD) se ACUMULAM, entao a diferenca entre os metodos, e entre
% eles e a solucao exata, so fica visivel depois de varios ciclos.
% Mostrar desde t=0 esconderia justamente o efeito estudado.
%
% ONDE a janela fica depende do amortecimento - sao TRES regimes:
%
%   xi != 0, amortecimento FORTE (decai em poucos ciclos): a resposta
%             morre logo, entao nao ha o que recortar - a figura mostra a
%             resposta INTEIRA, de t = 0 ate o instante do decaimento.
%             "Poucos" e definido por cfg.ciclosDecaimentoRapido (abaixo).
%
%   xi != 0, amortecimento FRACO: a resposta DECAI, mas devagar. Depois de
%             certo ponto ela e numericamente zero e a figura vira uma
%             linha reta - mostrar os ultimos ciclos da SIMULACAO nao
%             mostraria nada. Entao o codigo detecta o instante em que a
%             resposta exata cai abaixo de cfg.tolDecaimento e plota os
%             ultimos X_ciclos_final_plot ciclos ANTES desse instante, que
%             e onde ainda ha sinal.
%
%   xi == 0 : a resposta nao decai nunca, nao ha o que detectar - a
%             janela sao os ultimos X_ciclos_final_plot_sem_xi ciclos do
%             fim da simulacao.
%
% Por isso os dois casos tem contagens INDEPENDENTES: com amortecimento
% costuma valer a pena calcular muitos ciclos (para conter o decaimento
% inteiro) e plotar poucos; sem amortecimento nao ha razao para calcular
% uma margem grande.
%
% A escolha entre os dois blocos e AUTOMATICA: graficos().resposta_tempo olha
% o xi do caso que recebeu e usa o conjunto correspondente.

% --- CASO AMORTECIDO (xi != 0) ----------------------------------------
% nPerImpl / nPerRK: quantos ciclos de T1 sao CALCULADOS.
%   nPerImpl (metodos diretos) e barato: a hRat=0.01 sao ~100 passos por
%   ciclo, entao 50 ciclos custam ~5000 passos.
%   nPerRK e CARO e domina o tempo total de execucao: cada ciclo custa
%   ~2e5 passos POR METODO (ver fatorSegurancaRK), ou seja, 50 ciclos sao
%   ~1e7 passos por metodo, por caso. Aumente com parcimonia.
%   Se a resposta decair antes de nPer ciclos, os ciclos restantes sao
%   calculados e jogados fora - o console avisa quando isso acontece,
%   para voce poder reduzir.
cfg.nPerImpl = 50;
cfg.nPerRK   = 50;

% X_ciclos_final_plot: quantos ciclos aparecem na figura, contados PARA
%   TRAS a partir do instante de decaimento detectado.
%   Exemplo do uso pretendido: calcular 30 ciclos e plotar os 5 ultimos
%   antes de a resposta morrer -> nPerImpl = 30, X_ciclos_final_plot = 5.
cfg.X_ciclos_final_plot = 2;

% --- CASO NAO AMORTECIDO (xi == 0) ------------------------------------
% Contagens proprias: sem decaimento para conter, nao ha motivo para
% calcular a mesma margem do caso amortecido.
cfg.nPerImpl_sem_xi = 50;
cfg.nPerRK_sem_xi   = 50;

% X_ciclos_final_plot_sem_xi: ciclos mostrados, contados para tras a
%   partir do FIM da simulacao.
cfg.X_ciclos_final_plot_sem_xi = 2;

% --- DETECCAO DO DECAIMENTO (so usada quando xi != 0) -----------------
% tolDecaimento: a partir de que fracao do pico a resposta conta como
%   "zero". O codigo acha o ULTIMO instante em que
%   |w_tip_exato| > tolDecaimento * max|w_tip_exato| e ancora a janela
%   ali. Sempre medido na solucao EXATA, nunca na resposta de um metodo -
%   um metodo instavel diverge e arrastaria a janela para o lugar errado.
%   1e-3  -> a resposta caiu a 0.1% do pico
%   1e-2  -> a resposta caiu a 1% do pico (valor atual): janela mais cedo,
%            ainda com bastante sinal
%   1e-6  -> janela mais tarde, ja perto do ruido numerico
%   Tambem define QUANDO o decaimento e "rapido" (ver abaixo): a tolerancia
%   e o que fixa o instante em que a resposta conta como morta.
cfg.tolDecaimento = 1e-2;

% ciclosDecaimentoRapido: se a resposta exata cair abaixo de tolDecaimento
%   em ate este numero de ciclos de T1, o amortecimento e "forte" e a
%   figura mostra a resposta INTEIRA (t = 0 ate o decaimento) em vez de uma
%   janela de X_ciclos_final_plot ciclos. Acima disso vale a janela
%   ancorada no decaimento.
%   Ciclos ate decair a tol:   n = -ln(tol) / (2*pi*xi)
%   Com tol = 1e-2:  xi = 0.2 -> 3.7 ciclos (resposta inteira)
%                    xi = 0.1 -> 7.3 ciclos (janela)
%   6     -> "decai entre 1 e 6 ciclos" = amortecimento forte
%   0     -> desliga: sempre usa a janela
cfg.ciclosDecaimentoRapido = 6;

% ======================================================================
%  3. FIGURAS
% ======================================================================
% Toda figura e gravada em <pasta do codigo>/ImgViga/... como .eps (vetorial,
% para o paper) E .png (300 dpi, copia para conferir), por
% graficos().salva_figura - que tambem avisa quando um arquivo nao foi gravado
% (ex.: aberto em outro programa). A pasta nao depende do pwd do MATLAB.

% FS: tamanho da fonte dos eixos das figuras (pt).
cfg.FS = 18;

% FIGPOS: tamanho da figura exportada, [x y largura altura] em CENTIMETROS.
%   Usado no .eps final; 20x14 cm e o padrao do resto do paper.
cfg.FIGPOS = [0 0 20 14];

% YCLIP: limite do eixo y nas figuras de resposta no tempo, em unidades da
%   resposta normalizada (1.0 = amplitude da solucao exata na janela).
%   Existe porque metodos instaveis (RK1 com pouco amortecimento, DC acima
%   do h_cr) divergem para valores enormes e, sem corte, achatariam todas
%   as outras curvas contra o eixo. 1.6 deixa 60% de folga acima da
%   solucao exata: o suficiente para ver a curva divergindo sair da
%   figura, sem perder a resolucao das que se comportam.
cfg.YCLIP = 1.6;

% nModesReport: quantos modos aparecem por padrao nas tabelas de PE/AD
%   impressas no console (metodos().pead_tabela). Nao afeta o que e SALVO -
%   o .mat guarda PE/AD de todos os modos, isso e so a impressao.
%   Estudos C e D aumentam este valor quando excitam um modo mais alto,
%   para que o modo excitado apareca na tabela.
cfg.nModesReport = 3;

% hSweep: varredura de passos usada nas figuras de raio espectral
%   (graficos().espectro), em segundos. Precisa cobrir varias ordens de
%   grandeza porque os h_cr dos nove metodos vao de ~1e-9 s (RK1 sem
%   amortecimento) a ~1e-2 s. Calculo em forma fechada, barato: aumentar
%   o numero de pontos so deixa a curva mais lisa.
cfg.hSweep = logspace(-9, -2, 200);

% ======================================================================
%  4. AMOSTRAGEM E ARMAZENAMENTO DAS RESPOSTAS
% ======================================================================
% Estes tres parametros nao mudam a FISICA nem a precisao da integracao -
% todo passo continua sendo calculado. Eles controlam quantos pontos sao
% GUARDADOS e DESENHADOS, e existem por memoria: uma corrida de RK pode
% ter milhoes de passos, e guardar o historico inteiro estoura a RAM
% (o caso completo pode passar de 1e6 passos por metodo).

% nPontosExato: pontos da solucao exata de referencia (forma modal
%   fechada, solucao_exata) usados para desenhar a curva "Exact"
%   e normalizar as demais. 4000 pontos ja deixam a curva lisa.
cfg.nPontosExato = 4000;

% nMaxSalvar: maximo de pontos por metodo guardados no .mat e usados por
%   metodos().rk_stream_tip como grade de amostragem. Os passos entre esses
%   pontos sao calculados normalmente, so nao sao armazenados.
cfg.nMaxSalvar = 2e5;

% nMaxPlot: maximo de pontos efetivamente desenhados por curva. Acima
%   disso o .eps fica gigante e lento de abrir, sem ganho visual.
cfg.nMaxPlot = 5000;

% ======================================================================
%  5. CASOS VARRIDOS POR CADA ESTUDO
% ======================================================================
% Cada sub-struct e lido so pelo estudo correspondente.

% --- ESTUDO B: varredura de amortecimento fisico -----------------------
% xi_list: razoes de amortecimento modal varridas.
%   0     -> caso nao amortecido (a dissipacao que sobrar e 100% numerica)
%   0.002, 0.005, 0.01 -> tipico de estruturas metalicas (faixa de
%            0.005-0.02); 0.01 e o caso onde o amortecimento numerico mais
%            se destaca. Com tol = 1e-2 estes tres levam ~366, ~147 e ~73
%            ciclos para decair (n = -ln(tol)/(2*pi*xi)) - se nPer for
%            menor que isso, a deteccao nao dispara e a janela vai ate o
%            fim da simulacao (o console avisa).
%   0.20  -> mesmo valor do Estudo A, mantido para comparacao direta
% Custo: cada valor adicionado roda a bateria inteira de novo.
cfg.B.xi_list = [0, 0.005, 0.002, 0.01, 0.20];

% xiZoom: qual caso de xi ganha o "zoom" na oscilacao tardia (item 21 do
%   orientador). Deve ser um dos valores de xi_list; use [] para nenhum.
cfg.B.xiZoom = 0.01;

% --- ESTUDO C: CI modal, passo FIXO ------------------------------------
% xi: amortecimento usado (fixo, o caso de baixo amortecimento).
cfg.C.xi = 0.01;

% modos: quais modos sao excitados, um de cada vez, deslocando a viga na
%   forma daquele modo (ci_modal). O passo h NAO e reajustado
%   entre eles - esse e o ponto do estudo: ver h/T_modo crescer conforme o
%   modo sobe (modo 10 tem T ~260x menor que o modo 1).
cfg.C.modos = [1, 2, 3, 5, 10];

% --- ESTUDO D: CI modal, passo VARRIDO ---------------------------------
% xi: mesmo amortecimento do Estudo C, para os dois serem comparaveis.
cfg.D.xi = 0.01;

% modos: subconjunto representativo (baixo / medio / alto) de cfg.C.modos.
%   E um subconjunto porque aqui cada modo e combinado com cada razao de
%   passo abaixo, e varrer os 5 modos sairia caro no Runge-Kutta.
cfg.D.modos = [1, 3, 10];

% hRatios: razoes h/T_modo varridas, relativas ao periodo do modo
%   EXCITADO (nao ao T1).
%   0.01 -> regime atual, bem resolvido
%   0.10 -> intermediario, erros ja visiveis
%   0.50 -> meio caminho entre o intermediario e o limite
%   0.90 -> perto/acima do limite de estabilidade de varios esquemas.
%           A divergencia aqui e ESPERADA e faz parte do resultado; os
%           casos que explodem sao capturados por try/catch e reportados.
cfg.D.hRatios = [0.01, 0.10, 0.5, 0.90];

% nPer / ciclosPlot: ciclos calculados e mostrados SO no Estudo D.
%   Menores que os globais porque aqui cada combinacao (modo x razao) e
%   uma simulacao inteira, e as razoes altas divergem rapido - mais ciclos
%   nao acrescentariam informacao.
%   ATENCAO: os ciclos aqui sao contados em T1 (periodo do modo 1), o mesmo
%   eixo t/T1 dos demais estudos - NAO em T do modo excitado. Para o modo
%   10 (T ~250x menor que T1), 1 ciclo de T1 sao ~250 oscilacoes.
%   Como o Estudo D roda sempre amortecido (cfg.D.xi = 0.01), a deteccao
%   de decaimento da secao 2 vale aqui tambem (inclusive a regra da
%   resposta inteira, cfg.ciclosDecaimentoRapido), com estes valores no
%   lugar dos globais.
cfg.D.nPer = 7;
cfg.D.ciclosPlot = 2;

% ratiosSimuladas: quais razoes de hRatios ganham simulacao no tempo.
%   'extremos' -> so a menor e a maior (padrao, contem o custo; PE/AD em
%                 forma fechada continuam sendo calculados para TODAS)
%   'todas'    -> simula tambem as intermediarias
cfg.D.ratiosSimuladas = 'todas';

% --- ESTUDO E: amortecimento de Rayleigh alternativo -------------------
% rayleigh: pares (a0, a1) impostos DIRETAMENTE, em vez de resolver os
%   coeficientes para atingir um xi alvo em dois modos ancora.
%   Com a0 = 0 (termo proporcional a massa desligado), sobra so o termo
%   proporcional a rigidez e xi_k = a1*w_k/2 cresce LINEARMENTE com a
%   frequencia - os modos altos ficam fortemente amortecidos. O xi
%   resultante de cada modo e calculado, nao imposto.
%   Valores sugeridos pelo orientador (item 18): a1 = 1e-4 e 1e-5.
cfg.E.rayleigh = struct('a0', {0, 0}, 'a1', {1e-4, 1e-5});

% --- ESTUDO F: reducao modal (truncamento) -----------------------------
% reduction_modes: quantos modos manter em cada modelo reduzido, sempre
%   comparados contra o modelo completo (2N = 20 GDL).
%   4 modos -> w_max ~50x menor que o completo (passo critico ~50x maior)
%   3 modos -> w_max ~97x menor, mas ja com erro de truncamento visivel
%   O primeiro da lista e o usado no experimento do item 32 (Parte 2).
cfg.F.reduction_modes = [4, 3];

% item32_ordemRK: ordem do Runge-Kutta usada no experimento do item 32
%   (completo com h pequeno x reduzido com o mesmo h x reduzido com h
%   maior). RK4 por ser o mais preciso da familia, isolando o efeito da
%   reducao do efeito do erro de integracao.
cfg.F.item32_ordemRK = 4;

end

%% ====================================================================
%  Elementos finitos e modos
%  ====================================================================

function [K, M, info] = matrizes_EF(p, N)
% Matrizes K e M da viga de Euler-Bernoulli por elementos finitos.
%
% Elemento de 2 nos com 4 GDL (w1, theta1, w2, theta2), interpolacao cubica
% de Hermite e massa consistente. O engaste w(0) = w'(0) = 0 e imposto
% removendo as duas primeiras linhas/colunas do sistema montado.
%
% ENTRADA
%   p : struct de params
%   N : numero de elementos (opcional, padrao p.N)
%
% SAIDA
%   K, M : [n x n] matrizes globais, sem os GDL engastados
%   info : .N, .Le, .ndof (n = 2*N), .tipW (GDL transversal da ponta),
%          .tipTh (GDL de rotacao da ponta), .xNodes
    if nargin < 2 || isempty(N), N = p.N; end

    Le   = p.L/N;
    ndof = 2*(N+1);

    K = zeros(ndof);
    M = zeros(ndof);

    % --- Matrizes do elemento (constantes na malha uniforme) ---
    Ke = (p.E*p.I/Le^3) * [ 12,      6*Le,   -12,      6*Le;
                             6*Le,   4*Le^2, -6*Le,    2*Le^2;
                           -12,     -6*Le,    12,     -6*Le;
                             6*Le,   2*Le^2, -6*Le,    4*Le^2 ];

    Me = (p.rho*p.A*Le/420) * [ 156,     22*Le,   54,    -13*Le;
                                 22*Le,   4*Le^2, 13*Le,  -3*Le^2;
                                 54,     13*Le,  156,    -22*Le;
                                -13*Le,  -3*Le^2,-22*Le,   4*Le^2 ];

    % --- Montagem ---
    for e = 1:N
        d = [2*e-1, 2*e, 2*e+1, 2*e+2];
        K(d,d) = K(d,d) + Ke;
        M(d,d) = M(d,d) + Me;
    end

    % --- Engaste em x = 0: remove w1 e theta1 ---
    K = K(3:end, 3:end);
    M = M(3:end, 3:end);

    % --- Impoe simetria exata (elimina assimetria de arredondamento) ---
    K = 0.5*(K + K.');
    M = 0.5*(M + M.');

    info.N      = N;
    info.Le     = Le;
    info.ndof   = size(K,1);          % = 2*N
    info.tipW   = size(K,1) - 1;      % GDL transversal da ponta livre
    info.tipTh  = size(K,1);          % GDL de rotacao da ponta livre
    info.xNodes = (0:N).'*Le;
end

function [w, Phi] = autovalores_EF(K, M)
% Resolve K*phi = w^2*M*phi para o sistema discreto.
%
% SAIDA
%   w   : [n x 1] frequencias naturais em ordem crescente [rad/s]
%   Phi : [n x n] modos normalizados pela massa (Phi'*M*Phi = I)
    [V, D] = eig(full(K), full(M), 'chol');

    w2 = real(diag(D));
    w2(w2 < 0) = 0;                     % protege contra arredondamento
    [w2, idx] = sort(w2, 'ascend');
    w   = sqrt(w2);
    Phi = real(V(:, idx));

    % --- Normalizacao pela massa: phi' * M * phi = 1 ---
    for k = 1:size(Phi,2)
        Phi(:,k) = Phi(:,k) / sqrt(Phi(:,k).' * M * Phi(:,k));
    end
end

function mdl = monta(p, N)
% Monta de uma vez o modelo completo, na ordem correta.
%
% ENTRADA
%   p : struct de params (opcional)
%   N : numero de elementos (opcional, padrao p.N)
%
% SAIDA
%   mdl : .p, .info, .K, .M, .C (xi = p.xi_target), .Cz (xi = 0), .w, .Phi,
%         .damp, .u0, .v0, .a0v, .T1 (periodo fundamental), .wmax
    if nargin < 1 || isempty(p), p = params(); end
    if nargin < 2 || isempty(N), N = p.N;           end

    [K, M, info]  = matrizes_EF(p, N);
    [w, Phi]      = autovalores_EF(K, M);
    [C, damp]     = amortecimento(K, M, w, p);
    [u0, v0, a0v] = ci_impulso(K, M, C, info, p);

    mdl.p    = p;
    mdl.info = info;
    mdl.K    = K;
    mdl.M    = M;
    mdl.C    = C;
    mdl.Cz   = zeros(size(K));
    mdl.w    = w;
    mdl.Phi  = Phi;
    mdl.damp = damp;
    mdl.u0   = u0;
    mdl.v0   = v0;
    mdl.a0v  = a0v;
    mdl.T1   = 2*pi/w(1);
    mdl.wmax = w(end);
end

function [wn, betaL, phi] = modos_analiticos(p, nModes)
% Modos analiticos da viga continua, usados so para validar a malha.
%
% Equacao caracteristica da viga engastada-livre: cos(bL)*cosh(bL) + 1 = 0,
% com wn = (bL)^2 * sqrt(E*I/(rho*A*L^4)).
%
% SAIDA
%   wn    : [nModes x 1] frequencias naturais [rad/s]
%   betaL : [nModes x 1] raizes da equacao caracteristica
%   phi   : handle phi(x,n) com a forma do modo n em x
    if nargin < 2, nModes = 6; end

    f = @(x) cos(x).*cosh(x) + 1;

    betaL = zeros(nModes,1);
    for k = 1:nModes
        % A estimativa assintotica (2k-1)*pi/2 e otima para k >= 2 e
        % serve como ponto de partida para k = 1.
        guess = (2*k-1)*pi/2;
        lo = max(guess - 0.9, 1e-3);
        hi = guess + 0.9;
        % Garante troca de sinal antes de chamar fzero
        while f(lo)*f(hi) > 0
            lo = max(lo - 0.2, 1e-3);
            hi = hi + 0.2;
        end
        betaL(k) = fzero(f, [lo hi]);
    end

    wn = betaL.^2 * sqrt(p.E*p.I/(p.rho*p.A*p.L^4));

    % Funcoes de forma dos modos, sem normalizacao na ponta: para a
    % validacao so a forma importa.
    phi = @(x,k) modeShape(x, betaL(k), p.L);
end

function y = modeShape(x, bL, L)
    beta  = bL/L;
    sigma = (sinh(bL) + sin(bL)) / (cosh(bL) + cos(bL));
    y = cosh(beta*x) - cos(beta*x) - sigma*(sinh(beta*x) - sin(beta*x));
end

%% ====================================================================
%  Amortecimento e condicoes iniciais
%  ====================================================================

function [C, damp] = amortecimento(K, M, w, p, a0in, a1in)
% Amortecimento de Rayleigh C = a0*M + a1*K, com xi_k = a0/(2*w_k) + a1*w_k/2.
%
% Duas formas de uso:
%   1) amortecimento(K,M,w,p): escolhe a0 e a1 para que xi = p.xi_target
%      nos modos ancora p.mode_i e p.mode_j.
%   2) amortecimento(K,M,w,p,a0,a1): usa os coeficientes dados (Estudo E)
%      e o xi de cada modo e resultado, nao imposicao.
%
% O termo proporcional a rigidez faz xi crescer linearmente com a frequencia,
% entao os modos altos da malha ficam superamortecidos (xi >> 1). Esse e
% justamente o regime em que o amortecimento fisico interfere na estabilidade
% dos integradores.
%
% SAIDA
%   C    : [n x n] matriz de amortecimento
%   damp : .a0, .a1, .xi (todos os modos), .anchors, .xi_target
    if nargin < 4 || isempty(p), p = params(); end

    if nargin >= 6 && ~isempty(a0in) && ~isempty(a1in)
        % ---- Modo 2: coeficientes diretos ----
        a0 = a0in;
        a1 = a1in;
        damp.anchors    = [];
        damp.xi_target  = [];
    else
        % ---- Modo 1: ancorado em dois modos ----
        i = p.mode_i;  j = p.mode_j;
        wi = w(i);     wj = w(j);
        xi = p.xi_target;

        a1 = 2*xi/(wi + wj);
        a0 = a1*wi*wj;

        damp.anchors    = [i j];
        damp.xi_target  = xi;
    end

    C = a0*M + a1*K;
    C = 0.5*(C + C.');

    damp.a0 = a0;
    damp.a1 = a1;
    damp.xi = a0./(2*w) + a1*w/2;      % todos os modos da malha
end

function [u0, v0, a0] = ci_impulso(K, M, C, info, p)
% Condicao inicial de impulso: velocidade transversal V0 no GDL da ponta.
%
% O problema vira vibracao livre M*udd + C*ud + K*u = 0, com u(0) = 0 e
% ud(0) = v0. O sistema e linear, entao V0 so escala a resposta e nao altera
% nenhuma conclusao de estabilidade, PE ou AD.
%
% SAIDA
%   u0, v0 : [n x 1] deslocamento (nulo) e velocidade iniciais
%   a0     : [n x 1] aceleracao inicial, do equilibrio em t = 0
    if nargin < 5, p = params(); end

    n  = size(K,1);
    u0 = zeros(n,1);
    v0 = zeros(n,1);
    v0(info.tipW) = p.V0;

    a0 = M \ (-C*v0 - K*u0);
end

function [u0, v0, a0] = ci_modal(K, M, C, Phi, modeIdx, u0scale)
% Condicao inicial modal: u(0) = u0scale*phi_j e ud(0) = 0.
%
% Como Phi e normalizado pela massa e o amortecimento e proporcional, essa
% condicao projeta apenas no modo j. A resposta e uma oscilacao pura em w(j),
% o que permite isolar o comportamento dos integradores em uma frequencia
% (Estudos C e D), ao contrario do impulso, que excita varios modos juntos.
%
% ENTRADA
%   K, M, C : matrizes do sistema
%   Phi     : modos normalizados pela massa
%   modeIdx : modo excitado (coluna de Phi)
%   u0scale : amplitude modal (opcional, padrao 1)
    if nargin < 6 || isempty(u0scale), u0scale = 1; end

    n = size(K,1);
    if modeIdx < 1 || modeIdx > size(Phi,2)
        error('viga_ci_modal:modeIdx', ...
            'modeIdx = %d is out of range for a model with %d retained mode(s).', ...
            modeIdx, size(Phi,2));
    end

    u0 = u0scale * Phi(:,modeIdx);
    v0 = zeros(n,1);
    a0 = M \ (-C*v0 - K*u0);
end

%% ====================================================================
%  Solucao exata e reducao modal
%  ====================================================================

function [U, V, A] = solucao_exata(t, M, w, Phi, xi, u0, v0)
% Resposta exata do sistema semidiscreto, a referencia de todos os metodos.
%
% Com amortecimento proporcional o sistema desacopla na base modal e cada
% coordenada obedece a equacao de 1 GDL:
%
%     eta_k'' + 2*xi_k*w_k*eta_k' + w_k^2*eta_k = 0
%
% A resposta nao tem erro de integracao no tempo. Trata os tres regimes:
% subamortecido, critico e superamortecido. O ramo superamortecido usa
% exponenciais decrescentes em vez de sinh para evitar overflow nos modos altos.
%
% ENTRADA
%   t      : [1 x nt] instantes [s]
%   M      : matriz de massa (projeta as condicoes iniciais)
%   w, Phi : frequencias e modos normalizados pela massa
%   xi     : [n x 1] amortecimento modal
%   u0, v0 : condicoes iniciais
%
% SAIDA
%   U, V, A : [n x nt] deslocamento, velocidade e aceleracao
    t  = t(:).';
    n  = numel(w);
    nt = numel(t);

    U = zeros(n, nt);
    V = zeros(n, nt);
    A = zeros(n, nt);

    % --- Projecao modal das condicoes iniciais ---
    % Com Phi'*M*Phi = I a massa generalizada e unitaria.
    q0  = Phi.' * (M * u0);
    dq0 = Phi.' * (M * v0);

    tolc = 1e-10;

    for k = 1:n
        wk = w(k);
        zk = xi(k);
        e0 = q0(k);
        d0 = dq0(k);

        if wk < eps
            % Modo de corpo rigido (nao existe na viga engastada, fica por seguranca)
            eta  = e0 + d0*t;
            deta = d0*ones(1,nt);
            dde  = zeros(1,nt);

        elseif zk < 1 - tolc
            % ---- Subamortecido ----
            wd  = wk*sqrt(1 - zk^2);
            s   = zk*wk;
            ex  = exp(-s*t);
            cwd = cos(wd*t);
            swd = sin(wd*t);
            Ac  = e0;
            Bc  = (d0 + s*e0)/wd;
            eta  = ex.*(Ac*cwd + Bc*swd);
            deta = ex.*((-s*Ac + wd*Bc)*cwd + (-s*Bc - wd*Ac)*swd);
            dde  = -2*zk*wk*deta - wk^2*eta;

        elseif zk < 1 + tolc
            % ---- Amortecimento critico ----
            ex   = exp(-wk*t);
            eta  = ex.*(e0 + (d0 + wk*e0)*t);
            deta = ex.*((d0 + wk*e0) - wk*(e0 + (d0 + wk*e0)*t));
            dde  = -2*wk*deta - wk^2*eta;

        else
            % ---- Superamortecido (forma exponencial, sem overflow) ----
            sg = sqrt(zk^2 - 1);
            r1 = -(zk - sg)*wk;     % raiz de decaimento lento  (r1 < 0)
            r2 = -(zk + sg)*wk;     % raiz de decaimento rapido (r2 < r1 < 0)
            c1 = (d0 - r2*e0)/(r1 - r2);
            c2 = e0 - c1;
            E1 = exp(r1*t);
            E2 = exp(r2*t);
            eta  = c1*E1 + c2*E2;
            deta = c1*r1*E1 + c2*r2*E2;
            dde  = c1*r1^2*E1 + c2*r2^2*E2;
        end

        U = U + Phi(:,k)*eta;
        V = V + Phi(:,k)*deta;
        A = A + Phi(:,k)*dde;
    end
end

function red = reducao_modal(mdl, m)
% Reducao modal: mantem os m primeiros modos, em coordenadas modais q.
%
%     I*qdd + diag(2*xi_k*w_k)*qd + diag(w_k^2)*q = 0,   k = 1..m
%
% Com Phi normalizado pela massa e amortecimento proporcional, retirar os
% modos m+1..n nao introduz aproximacao numerica; a unica hipotese e fisica
% (esses modos nao importam para a resposta), e o Estudo F existe para
% testa-la. O sistema reduzido continua sendo um sistema linear, entao
% integradores e analises funcionam nele sem nenhuma alteracao.
%
% ENTRADA
%   mdl : modelo completo (monta)
%   m   : numero de modos mantidos (1 <= m <= n)
%
% SAIDA
%   red : .m, .K, .M, .C (diagonais [m x m]), .Phi, .w, .xi, .wmax, .T1,
%         .u0, .v0 (condicao inicial projetada nos modos mantidos),
%         .tipRow (u_tip = tipRow*q recupera a ponta), .Mref (eye(m), a
%         massa usada em solucao_exata junto com Phi = eye(m))
    n = size(mdl.K,1);
    if m < 1 || m > n
        error('viga_reducao_modal:m', ...
            'm = %d is out of range for a model with %d DOF.', m, n);
    end

    w  = mdl.w(1:m);
    xi = mdl.damp.xi(1:m);
    Phi_r = mdl.Phi(:,1:m);

    red.m    = m;
    red.M    = eye(m);
    red.K    = diag(w.^2);
    red.C    = diag(2*xi(:).*w(:));
    red.Phi  = Phi_r;
    red.w    = w;
    red.xi   = xi;
    red.wmax = w(end);
    red.T1   = 2*pi/w(1);

    red.u0 = Phi_r.' * (mdl.M * mdl.u0);
    red.v0 = Phi_r.' * (mdl.M * mdl.v0);

    red.tipRow = Phi_r(mdl.info.tipW, :);
    red.Mref   = eye(m);
end
