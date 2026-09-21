function cfgT = viga_aplica_config(cfgT, cfg, varargin)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Repassa a configuracao central para %
%  uma chamada de viga_resposta_tempo. %
%--------------------------------------%
%
% Copia de cfg (viga_config_estudos.m) para cfgT (a struct de UMA chamada
% de viga_resposta_tempo.m) todos os campos que sao simples pass-through:
% contagens de ciclos, janela de plot, cosmetica das figuras e limites de
% amostragem. Existe para que os cinco viga_estudo_*.m nao repitam a
% mesma lista de atribuicoes - e para que um campo novo em
% viga_config_estudos.m chegue automaticamente em todos eles.
%
% Campos repassados:
%   nPerImpl, nPerRK                    ciclos calculados (amortecido)
%   nPerImpl_sem_xi, nPerRK_sem_xi      idem, caso nao amortecido
%   X_ciclos_final_plot                 ciclos plotados (amortecido)
%   X_ciclos_final_plot_sem_xi          idem, caso nao amortecido
%   tolDecaimento                       limiar de deteccao do decaimento
%   FS, FIGPOS, YCLIP                   cosmetica das figuras
%   nPontosExato, nMaxSalvar, nMaxPlot  amostragem/armazenamento
%
% Qual conjunto de ciclos vale (amortecido ou nao) NAO e decidido aqui:
% viga_resposta_tempo.m decide isso sozinho, olhando o xi do caso.
%
% INPUT
%   cfgT     : struct em construcao (pode ja ter K, M, C, u0, ... )
%   cfg      : saida de viga_config_estudos()
%   varargin : pares nome/valor que SOBRESCREVEM o que veio de cfg, para
%              um estudo que precise de valores proprios. Ex. (Estudo D):
%                 cfgT = viga_aplica_config(cfgT, cfg, ...
%                     'nPerImpl', cfg.D.nPer, 'nPerRK', cfg.D.nPer, ...
%                     'X_ciclos_final_plot', cfg.D.ciclosPlot);
%
% OUTPUT
%   cfgT : a mesma struct, com os campos acima preenchidos.

    campos = { ...
        'nPerImpl', 'nPerRK', ...
        'nPerImpl_sem_xi', 'nPerRK_sem_xi', ...
        'X_ciclos_final_plot', 'X_ciclos_final_plot_sem_xi', ...
        'tolDecaimento', ...
        'FS', 'FIGPOS', 'YCLIP', ...
        'nPontosExato', 'nMaxSalvar', 'nMaxPlot'};

    for k = 1:numel(campos)
        if isfield(cfg, campos{k})
            cfgT.(campos{k}) = cfg.(campos{k});
        end
    end

    % --- Overrides opcionais -------------------------------------------
    if mod(numel(varargin), 2) ~= 0
        error('viga_aplica_config:pares', ...
            'Os overrides devem vir em pares nome/valor (recebidos %d argumentos).', ...
            numel(varargin));
    end
    for k = 1:2:numel(varargin)
        nome = varargin{k};
        if ~ischar(nome) && ~isstring(nome)
            error('viga_aplica_config:nome', ...
                'O nome do override %d nao e texto.', (k+1)/2);
        end
        if ~ismember(char(nome), campos)
            error('viga_aplica_config:desconhecido', ...
                ['"%s" nao e um campo repassavel. Campos validos: %s. ' ...
                 '(Campos que nao sao pass-through, como label ou outDir, ' ...
                 'sao atribuidos direto em cfgT pelo proprio estudo.)'], ...
                char(nome), strjoin(campos, ', '));
        end
        cfgT.(char(nome)) = varargin{k+1};
    end
end
