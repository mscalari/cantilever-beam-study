function cfg = viga_config_estudos()
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  CONFIGURACAO NUMERICA CENTRAL dos   %
%  Estudos B-F. Mexa AQUI, em mais     %
%  lugar nenhum.                       %
%--------------------------------------%
%
% ESTE e o unico arquivo que voce precisa editar para mudar QUALQUER
% parametro numerico de simulacao dos Estudos B-F: quanto tempo simular,
% qual passo usar, quais casos varrer, como as figuras saem. Cada script
% viga_estudo_*.m chama esta funcao uma vez, no inicio, e le tudo daqui -
% nenhum deles tem mais numero "chumbado" no corpo do codigo.
%
% DIVISAO DE RESPONSABILIDADE (importante nao confundir os dois arquivos):
%
%   viga_params.m          -> o que a VIGA e: material (E, rho), geometria
%                             (L, b, h), malha (N), o amortecimento de
%                             referencia do Estudo A (xi_target, modos
%                             ancora) e a CI de impulso (V0). Mudar algo
%                             la muda o PROBLEMA FISICO sendo resolvido.
%
%   viga_config_estudos.m  -> COMO a simulacao e rodada e mostrada: passo
%       (este arquivo)        de tempo, numero de periodos, quais casos
%                             varrer, tamanho de fonte das figuras. Mudar
%                             algo aqui nao muda a viga, muda o
%                             experimento numerico feito sobre ela.
%
% ESTUDO A NAO LE DAQUI, DE PROPOSITO. viga_main_simulacao.m e o caso de
% referencia ja validado e publicado no paper, e foi mantido intocado
% (com a configuracao propria dele embutida) justamente para servir de
% regressao: se algo aqui for mexido e um resultado mudar, Estudo A
% continua sendo o ponto fixo de comparacao. Se um dia voce quiser que
% ele tambem leia daqui, os valores equivalentes estao documentados na
% secao 6 no fim deste arquivo.
%
% OUTPUT
%   cfg : struct com todos os campos descritos abaixo.

% ======================================================================
%  1. PASSO DE TEMPO
% ======================================================================
% Os nove metodos se dividem em duas familias com criterios de passo
% completamente diferentes, e por isso ha dois passos por simulacao:
%
%   - Metodos DIRETOS (Houbolt, Newmark, Wilson-Theta, Bathe): sao
%     incondicionalmente estaveis, entao o passo e escolhido por PRECISAO.
%     Convencao do paper: uma fracao do periodo fundamental T1.
%     (Diferencas Centrais e a excecao: e direto mas condicionalmente
%     estavel, e por isso viga_resposta_tempo.m checa o h_cr dele e o
%     retira da figura quando hImpl o ultrapassa.)
%
%   - Metodos RUNGE-KUTTA (RK1 a RK4): sao explicitos e condicionalmente
%     estaveis, entao o passo e escolhido pela ESTABILIDADE - um pouco
%     abaixo do passo critico h_cr calculado por viga_h_critico.m.

% hRat: passo dos metodos diretos, como fracao de T1 (hImpl = hRat*T1).
%   Menor  -> mais preciso e mais caro (custo cresce como 1/hRat).
%   Tipico -> 0.01 (o mesmo valor do Estudo A e das figuras do paper).
%   Barato de aumentar: a 0.01, um periodo custa so ~100 passos.
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
% A escolha entre os dois blocos e AUTOMATICA: viga_resposta_tempo.m olha
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
%   ancorada no decaimento, como antes.
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
% viga_salva_figura.m - que tambem avisa quando um arquivo nao foi gravado
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
%   impressas no console (viga_pead_tabela.m). Nao afeta o que e SALVO -
%   o .mat guarda PE/AD de todos os modos, isso e so a impressao.
%   Estudos C e D aumentam este valor quando excitam um modo mais alto,
%   para que o modo excitado apareca na tabela.
cfg.nModesReport = 3;

% hSweep: varredura de passos usada nas figuras de raio espectral
%   (viga_espectro.m), em segundos. Precisa cobrir varias ordens de
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
% (foi exatamente o que aconteceu na primeira execucao completa - ver
% README_viga.md, secao "Execucao completa").

% nPontosExato: pontos da solucao exata de referencia (forma modal
%   fechada, viga_solucao_exata.m) usados para desenhar a curva "Exact"
%   e normalizar as demais. 4000 pontos ja deixam a curva lisa.
cfg.nPontosExato = 4000;

% nMaxSalvar: maximo de pontos por metodo guardados no .mat e usados por
%   viga_rk_stream_tip.m como grade de amostragem. Os passos entre esses
%   pontos sao calculados normalmente, so nao sao armazenados.
cfg.nMaxSalvar = 2e5;

% nMaxPlot: maximo de pontos efetivamente desenhados por curva. Acima
%   disso o .eps fica gigante e lento de abrir, sem ganho visual.
cfg.nMaxPlot = 5000;

% ======================================================================
%  5. CASOS VARRIDOS POR CADA ESTUDO
% ======================================================================
% Estes campos vinham de viga_params.m (secao "Extended study
% parameters"); foram movidos para ca porque definem o EXPERIMENTO, nao a
% viga. Cada sub-struct e lido so pelo estudo correspondente.

% --- ESTUDO B: varredura de amortecimento fisico -----------------------
% xi_list: razoes de amortecimento modal varridas.
%   0     -> caso nao amortecido (a dissipacao que sobrar e 100% numerica)
%   0.002, 0.005, 0.01 -> tipico de estruturas metalicas (o orientador cita
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
%   forma daquele modo (viga_ci_modal.m). O passo h NAO e reajustado
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

% ======================================================================
%  6. EQUIVALENTES NO ESTUDO A (referencia, NAO lidos daqui)
% ======================================================================
% Para comparacao, estes sao os valores correspondentes embutidos em
% viga_main_simulacao.m. Eles NAO sao usados por este arquivo - estao
% aqui so para voce saber o que teria de casar se um dia quiser alinhar
% os dois, ou entender por que uma figura do Estudo A tem menos ciclos
% que a equivalente do Estudo B.
%
%   nPerImpl  = 12      (aqui: cfg.nPerImpl / cfg.nPerImpl_sem_xi)
%   nPerRK    = 12      (aqui: cfg.nPerRK / cfg.nPerRK_sem_xi)
%   n_amostras = 10     (aqui: X_ciclos_final_plot - ciclos plotados nos
%                        regimes 2/3, ver abaixo)
%   tolDecaimento = 1e-2, ciclosDecaimentoRapido = 6   (iguais aos
%                        defaults abaixo - duplicados no proprio arquivo
%                        de A, nao lidos daqui)
%   YCLIP     = 1.6     (igual)
%   hRat      = 0.01    (igual)
%   xi        = [0, 0.2]  -> Estudo B varre cfg.B.xi_list
%
% JANELA DE PLOT: desde 21/09/2026 o Estudo A usa a MESMA regra de tres
% regimes dos Estudos B-F (viga_janela_plot.m) - resposta inteira quando o
% amortecimento e forte (xi=0.2 decai em ~3.7 ciclos, dentro do limiar de
% 6), janela ancorada no decaimento quando e fraco, fim da simulacao sem
% amortecimento. So o que e DESENHADO muda; nPerImpl/nPerRK e todo o resto
% numerico de A continuam os de sempre.
%
% Nas figuras de RAIO ESPECTRAL o Estudo A ja segue o mesmo layout dos
% demais: metodos diretos e Runge-Kutta juntos, uma figura por xi.
%
% CORES: todo mundo (Estudo A e B-F) usa a mesma paleta, viga_cores.m.
%
% (Estudo A tambem chama mdof_rk.m com os quatro argumentos de saida
% [t,U,V,A], enquanto os Estudos B-F usam viga_rk_stream_tip.m, que nao
% guarda o historico completo. Mesmo resultado numerico, memoria muito
% menor - ver README_viga.md.)

end
