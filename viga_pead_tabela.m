function pead = viga_pead_tabela(w, xi, allNames, allLabels, allSteps, T1, varargin)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  PE and AD table, generalised for    %
%  the Estudos B-F.                    %
%--------------------------------------%
%
% Generalises Section 3 of viga_main_simulacao.m (Section 3.3 of the
% paper) to an arbitrary modal set (w, xi) and an arbitrary list of
% methods/steps - one call per case, instead of the fixed xi0/xi02 loop.
% viga_pe_ad.m itself already decides per mode whether PE/AD are
% percent (xi_k = 0) or ratios (0 < xi_k < 1) or NaN (xi_k >= 1), so no
% single "unit" label is imposed here - the xi column tells the reader
% which regime each mode is in.
%
% viga_main_simulacao.m is left untouched.
%
% INPUT (required)
%   w, xi     : [m x 1] modal frequencies/damping ratios of the case
%   allNames  : cellstr of methods, e.g. [dirNames, rkNames]
%   allLabels : cellstr of display labels (same length as allNames)
%   allSteps  : [1 x numel(allNames)] time step used by each method
%   T1        : reference period for the printed h/T1 column
%
% INPUT (name-value, optional)
%   'nModesReport' : number of modes kept in the summary table (default 3,
%                    capped at numel(w))
%   'printHeader'  : string printed above the table (default '')
%
% OUTPUT
%   pead : struct with .methods, .labels, .h, .modes, .w, .xi,
%          .PE, .AD, .rho  ([nMethods x nModes], EVERY mode), and
%          .table (the first nModesReport modes, as a MATLAB table)

    q = inputParser;
    q.addParameter('nModesReport', 3);
    q.addParameter('printHeader', '');
    q.parse(varargin{:});
    o = q.Results;

    nModes = numel(w);
    nMeth  = numel(allNames);
    nR     = min(o.nModesReport, nModes);

    PEall  = nan(nMeth, nModes);
    ADall  = nan(nMeth, nModes);
    RHOall = nan(nMeth, nModes);

    if ~isempty(o.printHeader)
        fprintf(' %s\n', o.printHeader);
    end
    fprintf('   %-20s %10s', 'method', 'h/T_1');
    for m = 1:nR, fprintf('%10s%10s', sprintf('PE m%d',m), sprintf('AD m%d',m)); end
    fprintf('\n');

    for k = 1:nMeth
        hk = allSteps(k);
        for m = 1:nModes
            [PEall(k,m), ADall(k,m), RHOall(k,m)] = ...
                viga_pe_ad(allNames{k}, w(m), xi(m), hk);
        end
        fprintf('   %-20s %10.2e', allLabels{k}, hk/T1);
        for m = 1:nR
            fprintf('%10.4f%10.4f', PEall(k,m), ADall(k,m));
        end
        fprintf('\n');
    end

    pead.methods = allNames;
    pead.labels  = allLabels;
    pead.h       = allSteps;
    pead.modes   = (1:nModes).';
    pead.w       = w;
    pead.xi      = xi;
    pead.PE      = PEall;
    pead.AD      = ADall;
    pead.rho     = RHOall;

    data     = [allLabels(:), num2cell(allSteps(:)/T1)];
    varnames = {'method','h_T1'};
    for m = 1:nR
        data     = [data, num2cell(PEall(:,m)), num2cell(ADall(:,m))]; %#ok<AGROW>
        varnames = [varnames, sprintf('PE_m%d',m), sprintf('AD_m%d',m)]; %#ok<AGROW>
    end
    pead.table = cell2table(data, 'VariableNames', varnames);
end
