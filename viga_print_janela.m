function viga_print_janela(familia, info, T1, nPerCalculado)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Diagnostico de console da janela de %
%  plot - usado pelo Estudo A e B-F.   %
%--------------------------------------%
%
% Extraida da funcao local de mesmo nome em viga_resposta_tempo.m (mesmo
% padrao de viga_rho_beam.m/viga_print_janela.m) para que viga_main_
% simulacao.m tambem possa imprimir o diagnostico da janela, sem duplicar
% o texto.
%
% Uma linha dizendo onde a janela caiu e, quando for o caso, quantos
% ciclos foram calculados a toa depois do decaimento.
%
% INPUT
%   familia       : rotulo impresso, ex. 'direct schemes' ou 'Runge-Kutta'
%   info          : struct devolvida por viga_janela_plot.m
%   T1            : periodo de referencia [s]
%   nPerCalculado : ciclos efetivamente calculados nessa familia

    if info.completa
        fprintf(['   [janela %s] resposta decaiu a %.1e do pico em t/T1 = %.2f ' ...
                 '(<= %g ciclos): plotando a resposta INTEIRA, t/T1 em [%.2f, %.2f]\n'], ...
                familia, info.tol, info.tDecaimento/T1, info.cicloMaxCompleta, ...
                info.cicloIni, info.cicloFim);
    elseif info.decaiu
        fprintf(['   [janela %s] resposta decaiu a %.1e do pico em t/T1 = %.2f; ' ...
                 'plotando t/T1 em [%.2f, %.2f]\n'], ...
                familia, info.tol, info.tDecaimento/T1, info.cicloIni, info.cicloFim);
    end
    if info.decaiu
        % Ciclos calculados depois do decaimento sao jogados fora, em
        % qualquer um dos dois regimes acima.
        sobra = nPerCalculado - info.cicloFim;
        if sobra > 1
            fprintf(['   [janela %s] AVISO: %.1f ciclos calculados depois do ' ...
                     'decaimento foram descartados - da para reduzir nPer para ~%d\n'], ...
                    familia, sobra, ceil(info.cicloFim) + 1);
        end
    else
        if info.amortecido
            fprintf(['   [janela %s] resposta NAO caiu abaixo de %.1e do pico dentro ' ...
                     'dos %d ciclos calculados; plotando o fim: t/T1 em [%.2f, %.2f]\n'], ...
                    familia, info.tol, nPerCalculado, info.cicloIni, info.cicloFim);
        else
            fprintf('   [janela %s] sem amortecimento; plotando t/T1 em [%.2f, %.2f]\n', ...
                    familia, info.cicloIni, info.cicloFim);
        end
    end
    if info.truncada
        fprintf(['   [janela %s] AVISO: janela truncada em t=0 - foram pedidos %g ciclos ' ...
                 'mas so ha %.2f antes do decaimento\n'], ...
                familia, info.nCiclosPedidos, info.cicloFim);
    end
end
