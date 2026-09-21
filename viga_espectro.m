function stab = viga_espectro(w, xi, xiLabel, label, outDir, varargin)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Spectral radius of every scheme,    %
%  generalised for the Estudos B-F.    %
%--------------------------------------%
%
% Generalises Section 1 of viga_main_simulacao.m (Figs. 9-10 of the
% paper) to an ARBITRARY case, identified only by its modal set (w, xi).
% This works identically for the full beam model or for a reduced
% (truncated) model from viga_reducao_modal.m, since viga_rho_beam.m
% evaluates the operator mode by mode.
%
% viga_main_simulacao.m is left untouched - Estudo A keeps producing its
% own figures exactly as before. This function is called once per case
% by the new Estudo B-F scripts.
%
% INPUT (required)
%   w, xi    : [m x 1] modal frequencies/damping ratios of the case
%   xiLabel  : title string, e.g. '\xi = 0.01' or 'a_1 = 10^{-4}'
%   label    : filename tag, e.g. 'B_xi001' or 'F_red4'
%   outDir   : folder for the .eps figures (created if missing)
%
% INPUT (name-value, optional)
%   'hSweep'     : time-step sweep, default logspace(-9,-2,200)
%   'methodsRK'  : default {'RK1','RK2','RK3','RK4'}
%   'methodsDir' : default {'DC','Houbolt','Newmark','Wilson','Bathe'}
%   'dirLabels'  : default {'Central Differences','Houbolt','Newmark', ...
%                            'Wilson-Theta','Bathe'}
%   'FS'         : font size, default 18
%   'FIGPOS'     : figure size [x y w h] cm, default [0 0 20 14]
%   'doPlots'    : true/false, default true (set false to only compute
%                  the curves, e.g. for a quick smoke test)
%
% OUTPUT
%   stab : struct with
%            .hSweep
%            .rho.rk.<method>, .rho.direto.<method>   [1 x numel(hSweep)]

    q = inputParser;
    q.addParameter('hSweep',     logspace(-9,-2,200));
    q.addParameter('methodsRK',  {'RK1','RK2','RK3','RK4'});
    q.addParameter('methodsDir', {'DC','Houbolt','Newmark','Wilson','Bathe'});
    q.addParameter('dirLabels',  {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'});
    q.addParameter('FS',      18);
    q.addParameter('FIGPOS',  [0 0 20 14]);
    q.addParameter('doPlots', true);
    q.addParameter('xScaleBy', 1);              % e.g. 1/T1 to plot h/T1
    q.addParameter('xLabelStr', '$h$ [s]');
    q.parse(varargin{:});
    o = q.Results;

    if o.doPlots && ~exist(outDir,'dir'), mkdir(outDir); end

    coresRK  = [1.00 0.00 0.00; 0.00 0.20 1.00; 0.10 0.70 0.10; 0.60 0.00 0.60];
    coresDir = [1.00 0.00 0.00; 0.00 0.20 1.00; 0.10 0.70 0.10; 0.60 0.00 0.60; 1.00 0.50 0.00];

    stab.hSweep = o.hSweep;

    % 'tag' vai para o NOME DO ARQUIVO (viga_rho_<tag>_<label>.eps) e
    % 'titulo' vai para o titulo da figura, em ingles.
    fams = struct('tag',{'rk','direto'}, 'titulo',{'Runge-Kutta','direct schemes'}, ...
                   'names',{o.methodsRK,o.methodsDir}, ...
                   'labels',{o.methodsRK,o.dirLabels}, 'cores',{coresRK,coresDir});

    for f = 1:numel(fams)
        names = fams(f).names;  labels = fams(f).labels;  cores = fams(f).cores;
        if o.doPlots
            fig = figure('Name',sprintf('rho %s %s',fams(f).tag,label),'NumberTitle','off');
            hold on;
        end
        for k = 1:numel(names)
            r = arrayfun(@(hh) viga_rho_beam(names{k}, w, xi, hh), o.hSweep);
            stab.rho.(fams(f).tag).(names{k}) = r;
            if o.doPlots
                plot(o.hSweep*o.xScaleBy, r, '-', 'Color', cores(k,:), 'LineWidth', 1.5);
            end
        end
        if o.doPlots
            set(gca,'XScale','log','YScale','log');
            yline(1,'k--','HandleVisibility','off');
            grid on;
            xlabel(o.xLabelStr, 'Interpreter','latex');
            ylabel('$\rho(\mathbf{\hat{A}})$', 'Interpreter','latex');
            % Interpreter 'tex' (not 'latex') for the title: xiLabel comes
            % from the caller in many different shapes (plain text, '\xi',
            % 'a_0=..., a_1=...', 'modo 3, \xi = 0.01', ...) and 'tex' is
            % happy to mix plain words with backslash/underscore markup
            % inline, unlike 'latex', which demands the WHOLE string be
            % valid LaTeX (raw '_' or a bare '\xi' outside '$...$' is a
            % syntax error there, and silently breaks the exported title).
            title(sprintf('Spectral radius, %s (%s)', fams(f).titulo, xiLabel), ...
                'Interpreter','tex');
            legend(labels, 'Location','northwest', 'FontSize', 12);
            ylim([1e-2 1e2]);
            set(fig,'Units','centimeters','Position',o.FIGPOS);
            set(gca,'FontName','Times New Roman','FontSize',o.FS);
            print(fig, fullfile(outDir, sprintf('viga_rho_%s_%s.eps',fams(f).tag,label)), ...
                '-depsc2','-painters');
        end
    end
end
