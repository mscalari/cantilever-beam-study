function ok = viga_salva_figura(fig, outDir, nome)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  UNICA porta de saida das figuras:   %
%  grava .eps E .png, e CONFERE.       %
%--------------------------------------%
%
% Toda figura do projeto passa por aqui. Para cada figura grava
%
%     <outDir>/<nome>.eps    vetorial, para o paper (-depsc2 -painters)
%     <outDir>/<nome>.png    copia raster, 300 dpi, para conferir/colar
%
% e em seguida CONFERE que os dois arquivos realmente existem, nao estao
% vazios e foram gravados agora (e nao sao um arquivo antigo que ficou na
% pasta). Se algo falhar - figura invalida, arquivo aberto em outro
% programa, pasta sem permissao - avisa com o caminho ABSOLUTO em vez de
% deixar a rodada seguir como se tivesse dado certo.
%
% Uma falha NAO interrompe a rodada (uma bateria pode levar minutos):
% vira warning('viga_salva_figura:falhou', ...) e ok = false.
%
% INPUT
%   fig    : handle da figura
%   outDir : pasta de saida (criada se nao existir)
%   nome   : nome do arquivo SEM extensao, ex. 'viga_rho_B_xi_0_2'
%
% OUTPUT
%   ok     : true se .eps e .png foram gravados e conferidos

    ok = false;
    tEntrada = now;                       %#ok<TNOW1>  (so p/ comparar datenum)

    outDir = char(outDir);
    if ~exist(outDir,'dir')
        try
            mkdir(outDir);
        catch ME
            aviso(outDir, nome, sprintf('nao consegui criar a pasta (%s)', ME.message));
            return
        end
    end
    % Caminho absoluto, para as mensagens nao dependerem do pwd.
    [~, attr] = fileattrib(outDir);
    if isstruct(attr), outDir = attr.Name; end

    if isempty(fig) || ~ishghandle(fig)
        aviso(outDir, nome, 'o handle da figura e invalido (a figura foi fechada?)');
        return
    end

    arqEps = fullfile(outDir, [nome '.eps']);
    arqPng = fullfile(outDir, [nome '.png']);

    try
        print(fig, arqEps, '-depsc2', '-painters');
    catch ME
        aviso(outDir, [nome '.eps'], ME.message);
        return
    end
    try
        print(fig, arqPng, '-dpng', '-r300');
    catch ME
        aviso(outDir, [nome '.png'], ME.message);
        return
    end

    % --- Conferencia: existe, nao esta vazio, e e DESTA execucao --------
    % (2 s de tolerancia: a resolucao do carimbo de data depende do sistema
    %  de arquivos / sincronizacao do OneDrive.)
    tol = 2/86400;
    for arq = {arqEps, arqPng}
        d = dir(arq{1});
        if isempty(d)
            aviso(outDir, arq{1}, 'o arquivo nao existe depois do print');
            return
        elseif d.bytes == 0
            aviso(outDir, arq{1}, 'o arquivo foi criado vazio');
            return
        elseif d.datenum < tEntrada - tol
            aviso(outDir, arq{1}, ['o arquivo NAO foi atualizado (data ' ...
                datestr(d.datenum,'yyyy-mm-dd HH:MM:SS') ') - provavelmente esta ' ...
                'aberto/bloqueado em outro programa']);
            return
        end
    end

    ok = true;
end

function aviso(outDir, arquivo, motivo)
    warning('viga_salva_figura:falhou', ...
        'FIGURA NAO SALVA: %s (pasta %s): %s', arquivo, outDir, motivo);
end
