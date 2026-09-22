function [coresDir, coresRK] = viga_cores()
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Paleta de cores dos 9 metodos -     %
%  FONTE UNICA, usada em toda figura   %
%  que mostra mais de um metodo.       %
%--------------------------------------%
%
% Cada metodo tem SEMPRE a mesma cor, em toda figura do projeto (Estudo A
% ou B-F, resposta no tempo ou raio espectral) - para mudar uma cor, mude
% AQUI, nao em cada arquivo que desenha.
%
% coresDir (metodos diretos, nesta ordem: DC, Houbolt, Newmark,
%   Wilson-Theta, Bathe): vermelho, azul, verde, roxo, laranja.
%
% coresRK (RK1..RK4): as cores primarias do MATLAB 'r','b','g','m',
%   EXCETO RK1, que usa um vermelho mais escuro (grena, [0.6 0 0]) em vez
%   do 'r' puro. E a UNICA mudanca em relacao ao pedido original (DC =
%   vermelho puro [1 0 0], RK1 = 'r' = o MESMO vermelho puro): nas figuras
%   que mostram diretos e RK juntos (raio espectral), as duas curvas
%   ficariam literalmente identicas em cor. Os outros pares (RK2/Houbolt,
%   RK3/Newmark, RK4/Wilson-Theta) sao PARECIDOS mas nao iguais em RGB, e
%   continuam distinguiveis pelo estilo de linha (RK tracejado, direto
%   continuo) - nenhum outro ajuste foi feito.
%
% OUTPUT
%   coresDir : [5 x 3] RGB, ordem DC, Houbolt, Newmark, Wilson-Theta, Bathe
%   coresRK  : [4 x 3] RGB, ordem RK1, RK2, RK3, RK4

    coresDir = [1.0 0.0 0.0;    % DC           - vermelho
                0.0 0.2 1.0;    % Houbolt      - azul
                0.1 0.7 0.1;    % Newmark      - verde
                0.6 0.0 0.6;    % Wilson-Theta - roxo
                1.0 0.5 0.0];   % Bathe        - laranja

    coresRK  = [0.6 0.0 0.0;    % RK1 - vermelho escuro/grena (nao 'r' puro - ver acima)
                0.0 0.0 1.0;    % RK2 - azul   ('b')
                0.0 1.0 0.0;    % RK3 - verde  ('g')
                1.0 0.0 1.0];   % RK4 - magenta ('m')
end
