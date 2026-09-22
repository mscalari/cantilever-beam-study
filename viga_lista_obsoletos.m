function nomes = viga_lista_obsoletos(pastas, tInicio)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Lista figuras que NAO foram gravadas%
%  nesta execucao (sobras antigas).    %
%--------------------------------------%
%
% Figuras .eps/.png que ja estavam na pasta e NAO foram regravadas desde
% tInicio: sobras de execucoes anteriores, tipicamente de um nome que mudou
% (ex.: viga_rho_rk_*.eps virou viga_rho_*.eps) ou de um caso que saiu da
% configuracao. Sao facilmente confundidas com resultado atual.
%
% So LISTA - nunca apaga nada.
%
% INPUT
%   pastas  : cellstr de pastas a examinar (sem recursao)
%   tInicio : datenum do inicio da execucao (now, no comeco do script)
%
% OUTPUT
%   nomes   : cellstr com os caminhos das sobras (tambem impressa no console)

    nomes = {};
    for p = 1:numel(pastas)
        for ext = {'*.eps','*.png'}
            d = dir(fullfile(pastas{p}, ext{1}));
            for k = 1:numel(d)
                if d(k).datenum < tInicio - 2/86400
                    nomes{end+1,1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
                end
            end
        end
    end

    if ~isempty(nomes)
        fprintf('\n AVISO: %d figura(s) antiga(s) NAO foram regravadas nesta execucao\n', numel(nomes));
        fprintf(' (sobras de execucoes anteriores - podem ser de um nome que mudou):\n');
        for k = 1:numel(nomes)
            d = dir(nomes{k});
            fprintf('   %s   (%s)\n', nomes{k}, datestr(d.datenum,'yyyy-mm-dd HH:MM'));
        end
        fprintf(' Nada foi apagado.\n');
    end
end
