%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%  MASTER DRIVER - runs Estudos B-F    %
%  in the order recommended by the     %
%  plan (cheapest first).              %
%--------------------------------------%
%
% Does NOT run Estudo A: viga_main_simulacao.m is the untouched
% reference case and should be (re-)run on its own, as a regression
% check, before trusting anything here (see the plan's Etapa 2).
%
% Each Estudo script is self-contained (clears the workspace, builds its
% own model) and saves its own .mat + figures under ImgViga/EstudoX/, so
% this driver is only a convenience: run it, or run the scripts one at a
% time to inspect intermediate results.
%
% Estudo D is the most expensive per combination it actually simulates
% in the time domain (see its own header), but it caps periods and only
% simulates the two extreme h-ratios, so it stays bounded.

clc;
t0 = tic;

fprintf('############################################################\n');
fprintf('  ESTUDOS B-F - revisao do orientador (%s)\n', datestr(now));
fprintf('############################################################\n\n');

scripts = { ...
    'viga_estudo_B_amortecimento.m', ...
    'viga_estudo_E_rayleigh_alt.m', ...
    'viga_estudo_C_modal_h_fixo.m', ...
    'viga_estudo_D_modal_h_variavel.m', ...
    'viga_estudo_F_reducao_modal.m'};

for k = 1:numel(scripts)
    fprintf('\n>>> Executando %s (%d/%d) ...\n\n', scripts{k}, k, numel(scripts));
    tk = tic;
    run_estudo_isolado(scripts{k});
    fprintf('\n>>> %s concluido em %.1f s.\n', scripts{k}, toc(tk));
    close all;   % keep the figure count from a script from leaking into the next
end

fprintf('\n############################################################\n');
fprintf('  TODOS OS ESTUDOS CONCLUIDOS em %.1f min.\n', toc(t0)/60);
fprintf('  Resultados em %s\n', fullfile(pwd,'ImgViga'));
fprintf('############################################################\n');

% ==========================================================================
%  LOCAL FUNCTION
% ==========================================================================
% Each Estudo script starts with "clc; clear; close all;" so it can also
% be run standalone. Calling run() directly from this driver SCRIPT would
% execute it in the same base workspace as this driver, so that internal
% `clear` would wipe out this loop's own variables (scripts, k, t0) and
% abort the "for" loop after the first iteration. Wrapping the call in a
% local FUNCTION gives run() a fresh function workspace to clear instead,
% leaving the driver's loop state untouched.
function run_estudo_isolado(scriptname)
    run(scriptname);
end
