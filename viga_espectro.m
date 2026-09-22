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
% ONE figure per case, with ALL nine schemes together: the five direct
% schemes as solid lines and the four Runge-Kutta schemes as dashed lines,
% so both families are read against the same rho = 1 stability limit.
% (They used to be drawn in two separate figures.) Colours come from
% viga_cores.m, the single palette shared by every figure in the project.
%
% viga_main_simulacao.m draws its own figures in the same layout. This
% function is called once per case by the Estudo B-F scripts.
%
% INPUT (required)
%   w, xi    : [m x 1] modal frequencies/damping ratios of the case
%   xiLabel  : title string, e.g. '\xi = 0.01' or 'a_1 = 10^{-4}'
%   label    : filename tag, e.g. 'B_xi001' or 'F_red4'; the figure is
%              saved as viga_rho_<label>.eps and .png (viga_salva_figura.m)
%   outDir   : folder for the .eps figure (created if missing)
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

    [coresDir, coresRK] = viga_cores();   % paleta unica - ver viga_cores.m

    stab.hSweep = o.hSweep;

    % 'tag' is the field name inside stab.rho (unchanged from the old
    % two-figure layout, so stored results keep the same structure).
    fams = struct('tag',{'direto','rk'}, ...
                   'names',{o.methodsDir,o.methodsRK}, ...
                   'labels',{o.dirLabels,o.methodsRK}, ...
                   'cores',{coresDir,coresRK}, 'estilo',{'-','--'});

    if o.doPlots
        fig = figure('Name',sprintf('rho %s',label),'NumberTitle','off');
        hold on;
    end
    leg = {};
    for f = 1:numel(fams)
        names = fams(f).names;
        for k = 1:numel(names)
            r = arrayfun(@(hh) viga_rho_beam(names{k}, w, xi, hh), o.hSweep);
            stab.rho.(fams(f).tag).(names{k}) = r;
            if o.doPlots
                plot(o.hSweep*o.xScaleBy, r, fams(f).estilo, ...
                    'Color', fams(f).cores(k,:), 'LineWidth', 1.5);
                leg{end+1} = fams(f).labels{k}; %#ok<AGROW>
            end
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
        title(sprintf('Spectral radius, direct schemes and Runge-Kutta (%s)', xiLabel), ...
            'Interpreter','tex');
        legend(leg, 'Location','northwest', 'FontSize', 11, 'NumColumns', 2);
        ylim([1e-2 1e2]);
        set(fig,'Units','centimeters','Position',o.FIGPOS);
        set(gca,'FontName','Times New Roman','FontSize',o.FS);
        viga_salva_figura(fig, outDir, sprintf('viga_rho_%s',label));
    end
end
