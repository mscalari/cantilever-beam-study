function [tIni, tFim, info] = viga_janela_plot(t, uex, T1, nCiclos, tol, amortecido, cicloMaxCompleta)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Janela de tempo mostrada nas figuras%
%  de resposta (Estudos B-F).          %
%--------------------------------------%
%
% Decide QUE TRECHO da resposta aparece na figura. Sao tres regimes:
%
% 1) AMORTECIMENTO FORTE - resposta inteira.
%    Se a resposta exata cai abaixo de tol vezes o proprio pico em NO MAXIMO
%    cicloMaxCompleta ciclos de T1 (ex.: 6), ela morre logo, e nao ha o que
%    escolher: a figura mostra a resposta TODA, de t = 0 ate o instante do
%    decaimento. E o unico jeito de ver o transiente inicial, que e onde os
%    metodos mais diferem entre si nesse caso.
%
% 2) AMORTECIMENTO FRACO - janela ancorada no decaimento.
%    A resposta decai, mas devagar: plotar o fim da SIMULACAO mostraria uma
%    linha reta. Entao a janela e ancorada no DECAIMENTO - o ultimo instante
%    em que a resposta exata ainda esta acima de tol vezes o seu proprio
%    pico - e recua nCiclos ciclos de T1 a partir dali. (Se a resposta nem
%    chegou a decair dentro da simulacao, a janela vai ate o fim dela.)
%
% 3) SEM AMORTECIMENTO - fim da simulacao.
%    Nao ha decaimento a detectar; a janela e simplesmente os ultimos
%    nCiclos ciclos da simulacao.
%
% A deteccao usa sempre a solucao EXATA, nunca a resposta de um metodo:
% um esquema instavel diverge e arrastaria a janela para o lugar errado,
% justamente no caso em que a figura mais precisa estar certa.
%
% Como a resposta oscila, |uex| cruza qualquer limiar varias vezes por
% ciclo; o que interessa e o ULTIMO cruzamento, que e onde a envoltoria do
% decaimento atravessa o limiar - dai o find(..., 'last'). Para nao
% confundir um cruzamento por zero no fim da simulacao com decaimento, so
% se considera decaida a resposta que fica abaixo da tolerancia por pelo
% menos um periodo T1 inteiro antes do fim.
%
% INPUT
%   t, uex            : instantes [s] e resposta EXATA da ponta, sobre a
%                       simulacao inteira (a grade so precisa resolver o
%                       decaimento, nao cada oscilacao)
%   T1                : periodo de referencia [s] (contagem de ciclos)
%   nCiclos           : ciclos mostrados no regime 2 (ou 3)
%   tol               : fracao do pico que conta como "zero"
%   amortecido        : true se o caso tem xi ~= 0
%   cicloMaxCompleta  : (opcional, default 0 = regime 1 desligado) maximo
%                       de ciclos de T1 ate o decaimento para que a
%                       resposta inteira seja mostrada
%
% OUTPUT
%   tIni, tFim : limites da janela [s]
%   info       : struct com o diagnostico (usado no console e salvo no
%                .mat), campos:
%                  .amortecido   true/false
%                  .tDecaimento  instante detectado [s] (NaN se nao houve)
%                  .decaiu       true se a resposta caiu abaixo de tol
%                                DENTRO da simulacao
%                  .completa     true se caiu no regime 1 (resposta inteira)
%                  .tIni, .tFim  a janela, em segundos
%                  .cicloIni, .cicloFim  a janela, em ciclos de T1
%                  .nCiclosPedidos  nCiclos
%                  .cicloMaxCompleta  o limiar do regime 1
%                  .truncada     true se a janela bateu em t=0 e ficou
%                                menor que nCiclos ciclos (regime 2)
%                  .tol          tolerancia usada

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
