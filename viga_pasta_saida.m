function d = viga_pasta_saida(varargin)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Pasta de saida das figuras e .mat,  %
%  FIXA no local do codigo.            %
%--------------------------------------%
%
% Antes, cada script montava a saida como fullfile(pwd,'ImgViga',...): se a
% pasta atual do MATLAB nao fosse a do codigo, as figuras iam parar em
% outro lugar sem nenhum aviso (e a pasta ImgViga "de sempre" ficava com os
% arquivos antigos, dando a impressao de que nada foi salvo).
%
% Agora a saida e sempre <pasta onde esta viga_params.m>/ImgViga,
% independente do pwd.
%
% USO
%   viga_pasta_saida()              -> <codigo>/ImgViga
%   viga_pasta_saida('EstudoB')     -> <codigo>/ImgViga/EstudoB
%
% A pasta NAO e criada aqui (quem grava cria).

    raiz = fileparts(which('viga_params'));
    if isempty(raiz)
        raiz = pwd;
        warning('viga_pasta_saida:semCodigo', ...
            ['viga_params.m nao esta no path; usando a pasta atual (%s) como raiz. ' ...
             'Adicione a pasta do codigo ao path.'], raiz);
    end
    d = fullfile(raiz, 'ImgViga', varargin{:});
end
