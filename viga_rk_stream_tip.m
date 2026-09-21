function [tSave, uSave, uMax] = viga_rk_stream_tip(K, M, C, u0, v0, h, nSteps, order, tipExtract, saveIdx)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Memory-light RK time march, TIP     %
%  response only (Estudos B-F).        %
%--------------------------------------%
%
% RK step counts are set by explicit-scheme STABILITY, not accuracy, so
% nSteps can reach into the millions at this mesh's w_max (see
% viga_h_critico.m) - large enough that mdof_rk.m's [n x nSteps+1] U (let
% alone V and A) can exceed available memory on a modest machine, even
% after mdof_rk.m stopped allocating V/A when they are not requested
% (nargout-gated there). This function runs the EXACT SAME explicit-RK
% recursion (identical Ahat truncated-Taylor operator, so numerically
% identical to mdof_rk.m step for step) but never materialises the full
% state history: only the scalar tip value at the caller-chosen SAVED
% steps, and the running max(|tip|), are kept.
%
% INPUT
%   K, M, C, u0, v0, h, nSteps, order : same as mdof_rk.m
%   tipExtract : function handle mapping a state column ([n x 1], the
%                displacement part y(1:n)) to the scalar tip response -
%                the SAME handles used elsewhere in the Estudos, e.g.
%                @(U) U(tipIdx,:) for the full model or @(Q) red.tipRow*Q
%                for a reduced model (viga_reducao_modal.m); both work
%                unchanged on a single column.
%   saveIdx    : sorted, 0-based step indices (0..nSteps) at which to
%                save the tip value, e.g. from
%                unique(round(linspace(0, nSteps, min(nSteps+1, 2e5))))
%
% OUTPUT
%   tSave, uSave : [1 x numel(saveIdx)] time / tip value at saveIdx
%   uMax         : max(|tip value|) over ALL nSteps+1 steps (not just the
%                  saved ones - matches what max(abs(uTip)) gave when
%                  uTip was the full mdof_rk-based array)

    n = size(K,1);

    MiK = M \ K;
    MiC = M \ C;
    Ast = [ zeros(n), eye(n) ;
            -MiK    , -MiC   ];

    Ahat = eye(2*n);
    term = eye(2*n);
    for i = 1:order
        term = (term*(h*Ast))/i;
        Ahat = Ahat + term;
    end

    saveIdx = unique(saveIdx(:).');
    nSave   = numel(saveIdx);
    tSave   = saveIdx*h;
    uSave   = zeros(1, nSave);

    y     = [u0(:); v0(:)];
    uTip0 = tipExtract(y(1:n));
    uMax  = abs(uTip0);

    si = 1;
    if si <= nSave && saveIdx(si) == 0
        uSave(si) = uTip0;
        si = si + 1;
    end

    for k = 1:nSteps
        y    = Ahat*y;
        uTip = tipExtract(y(1:n));
        uMax = max(uMax, abs(uTip));
        if si <= nSave && saveIdx(si) == k
            uSave(si) = uTip;
            si = si + 1;
        end
    end
end
