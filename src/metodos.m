function m = metodos()
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%       Prof. Thales Peixoto           %
%--------------------------------------%
%
% Os nove metodos de integracao no tempo e as ferramentas de analise deles:
% operador de aproximacao, raio espectral, PE/AD e passo critico.
%
% As funcoes deste arquivo sao chamadas de fora como metodos().nome(...).
% So as listadas abaixo sao publicas; as demais sao de uso interno.

    m.h_critico      = @h_critico;
    m.integra_metodo = @integra_metodo;
    m.mdof_bathe     = @mdof_bathe;
    m.mdof_rk        = @mdof_rk;
    m.nomes_metodos  = @nomes_metodos;
    m.operador       = @operador;
    m.pe_ad          = @pe_ad;
    m.pead_tabela    = @pead_tabela;
    m.rho_beam       = @rho_beam;
    m.rk_stream_tip  = @rk_stream_tip;
end

%% ====================================================================
%  Metodos diretos
%  ====================================================================

function [t, U, V, A] = mdof_dc(K, M, C, u0, v0, h, nSteps)
% Diferencas Centrais, condicionalmente estavel: h_cr = 2/w_max, em que w_max
% e a MAIOR frequencia da malha, nao a fundamental.
%
% Velocidade e aceleracao sao estimativas centrais e dependem dos vizinhos;
% a primeira coluna vem do estado inicial e a ultima de estimativa unilateral.
%
% ENTRADA
%   K, M, C : [n x n] matrizes do sistema
%   u0, v0  : [n x 1] condicoes iniciais
%   h       : passo de tempo [s]
%   nSteps  : numero de passos
%
% SAIDA
%   t       : [1 x nSteps+1] instantes
%   U, V, A : [n x nSteps+1] deslocamento, velocidade e aceleracao
    n = size(K,1);

    % --- Constantes de integracao ---
    a0 = 1/h^2;
    a1 = 1/(2*h);
    a2 = 2*a0;
    a3 = 1/a2;

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);

    % --- Deslocamento ficticio em t = -h (expansao de Taylor) ---
    uPrev = u - h*v + a3*a;

    % --- Massa efetiva, fatorada uma vez ---
    Meff = a0*M + a1*C;
    dMe  = decomposition(Meff);

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);
    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    uCurr = u;
    for k = 1:nSteps
        Reff  = -(K - a2*M)*uCurr - (a0*M - a1*C)*uPrev;
        uNext = dMe \ Reff;

        % Estimativas centrais no passo ATUAL
        if k > 1
            V(:,k) = a1*(uNext - uPrev);
            A(:,k) = a0*(uNext - 2*uCurr + uPrev);
        end

        uPrev = uCurr;
        uCurr = uNext;
        U(:,k+1) = uCurr;
    end

    % Estimativa unilateral da ultima coluna
    if nSteps >= 1
        V(:,end) = (U(:,end) - U(:,max(end-1,1)))/h;
        A(:,end) = M \ (-C*V(:,end) - K*U(:,end));
    end

    t = (0:nSteps)*h;
end

function [t, U, V, A] = mdof_houbolt(K, M, C, u0, v0, h, nSteps)
% Houbolt, metodo de tres passos e portanto nao auto-inicializavel: os
% deslocamentos em t = h e t = 2h vem de dois passos de Newmark de mesmo h.
%
% ENTRADA
%   K, M, C : [n x n] matrizes do sistema
%   u0, v0  : [n x 1] condicoes iniciais
%   h       : passo de tempo [s]
%   nSteps  : numero de passos
%
% SAIDA
%   t       : [1 x nSteps+1] instantes
%   U, V, A : [n x nSteps+1] deslocamento, velocidade e aceleracao
    n = size(K,1);

    % --- Constantes de integracao ---
    a0 = 2/h^2;
    a1 = 11/(6*h);
    a2 = 5/h^2;
    a3 = 3/h;
    a4 = -2*a0;
    a5 = -a3/2;
    a6 = a0/2;
    a7 = a3/9;

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    % --- Inicializacao: dois passos de Newmark ---
    nStart = min(2, nSteps);
    [~, Us, Vs, As] = mdof_newmark(K, M, C, u0, v0, h, nStart);
    U(:,1:nStart+1) = Us;
    V(:,1:nStart+1) = Vs;
    A(:,1:nStart+1) = As;

    if nSteps < 3
        t = (0:nSteps)*h;
        return
    end

    % --- Rigidez efetiva, fatorada uma vez ---
    Keff = K + a0*M + a1*C;
    dKe  = decomposition(Keff);

    % --- Laco no tempo (o indice k corresponde a t = (k-1)*h) ---
    for k = 3:nSteps
        uc = U(:,k);      % t
        u1 = U(:,k-1);    % t - h
        u2 = U(:,k-2);    % t - 2h

        Reff = M*(a2*uc + a4*u1 + a6*u2) + C*(a3*uc + a5*u1 + a7*u2);
        uNew = dKe \ Reff;

        U(:,k+1) = uNew;
        A(:,k+1) = a0*uNew - a2*uc - a4*u1 - a6*u2;
        V(:,k+1) = a1*uNew - a3*uc - a5*u1 - a7*u2;
    end

    t = (0:nSteps)*h;
end

function [t, U, V, A] = mdof_newmark(K, M, C, u0, v0, h, nSteps, alpha, delta)
% Newmark. Com alpha = 1/4 e delta = 1/2 (padrao) e a regra trapezoidal de
% aceleracao media constante, incondicionalmente estavel. A rigidez efetiva e
% constante em analise linear, entao e fatorada uma unica vez.
%
% ENTRADA
%   K, M, C : [n x n] matrizes do sistema
%   u0, v0  : [n x 1] condicoes iniciais
%   h       : passo de tempo [s]
%   nSteps  : numero de passos
%   alpha, delta : parametros de Newmark (opcionais)
%
% SAIDA
%   t       : [1 x nSteps+1] instantes
%   U, V, A : [n x nSteps+1] deslocamento, velocidade e aceleracao
    if nargin < 8 || isempty(alpha), alpha = 1/4; end
    if nargin < 9 || isempty(delta), delta = 1/2; end

    n = size(K,1);

    % --- Constantes de integracao (Bathe, Finite Element Procedures) ---
    a0 = 1/(alpha*h^2);
    a1 = delta/(alpha*h);
    a2 = 1/(alpha*h);
    a3 = 1/(2*alpha) - 1;
    a4 = delta/alpha - 1;
    a5 = (h/2)*(delta/alpha - 2);
    a6 = h*(1 - delta);
    a7 = delta*h;

    % --- Rigidez efetiva, fatorada uma vez ---
    Keff = K + a0*M + a1*C;
    dKe  = decomposition(Keff);

    % --- Alocacao e estado inicial ---
    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);

    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    % --- Laco no tempo ---
    for k = 1:nSteps
        Reff = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a);
        uNew = dKe \ Reff;
        aNew = a0*(uNew - u) - a2*v - a3*a;
        vNew = v + a6*a + a7*aNew;

        u = uNew;  v = vNew;  a = aNew;
        U(:,k+1) = u;  V(:,k+1) = v;  A(:,k+1) = a;
    end

    t = (0:nSteps)*h;
end

function [t, U, V, A] = mdof_wilson(K, M, C, u0, v0, h, nSteps, theta)
% Wilson-Theta, com theta = 1.4 por padrao. Em vibracao livre a carga
% extrapolada em t + theta*h e nula.
%
% O metodo apresenta overshoot nos primeiros passos. Na malha, ele e puxado
% pelos modos com maior w*h, entao fica forte na condicao de impulso, que
% excita todo o espectro.
%
% ENTRADA
%   K, M, C : [n x n] matrizes do sistema
%   u0, v0  : [n x 1] condicoes iniciais
%   h       : passo de tempo [s]
%   nSteps  : numero de passos
%   theta   : parametro do metodo (opcional)
%
% SAIDA
%   t       : [1 x nSteps+1] instantes
%   U, V, A : [n x nSteps+1] deslocamento, velocidade e aceleracao
    if nargin < 8 || isempty(theta), theta = 1.4; end

    n = size(K,1);

    % --- Constantes de integracao ---
    a0 = 6/(theta*h)^2;
    a1 = 3/(theta*h);
    a2 = 2*a1;
    a3 = (theta*h)/2;
    a4 = a0/theta;
    a5 = -a2/theta;
    a6 = 1 - 3/theta;
    a7 = h/2;
    a8 = h^2/6;

    % --- Rigidez efetiva, fatorada uma vez ---
    Keff = K + a0*M + a1*C;
    dKe  = decomposition(Keff);

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);
    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    % --- Laco no tempo ---
    for k = 1:nSteps
        Reff  = M*(a0*u + a2*v + 2*a) + C*(a1*u + 2*v + a3*a);
        uTh   = dKe \ Reff;                 % deslocamento em t + theta*h

        aNew = a4*(uTh - u) + a5*v + a6*a;
        vNew = v + a7*(aNew + a);
        uNew = u + h*v + a8*(aNew + 2*a);

        u = uNew;  v = vNew;  a = aNew;
        U(:,k+1) = u;  V(:,k+1) = v;  A(:,k+1) = a;
    end

    t = (0:nSteps)*h;
end

function [t, U, V, A] = mdof_bathe(K, M, C, u0, v0, h, nSteps)
% Bathe, com dois subpassos iguais (gamma = 1/2):
%   subpasso 1: regra trapezoidal em h/2, equilibrio em t + h/2
%   subpasso 2: Euler regressivo de tres pontos, equilibrio em t + h
%
% As duas rigidezes efetivas sao diferentes e cada uma e fatorada uma vez;
% esse e o custo extra do Bathe em relacao a um metodo de passo unico.
%
% ENTRADA
%   K, M, C : [n x n] matrizes do sistema
%   u0, v0  : [n x 1] condicoes iniciais
%   h       : passo COMPLETO (os dois subpassos) [s]
%   nSteps  : numero de passos
%
% SAIDA
%   t       : [1 x nSteps+1] instantes
%   U, V, A : [n x nSteps+1] deslocamento, velocidade e aceleracao
    n = size(K,1);

    % --- Constantes do subpasso 1 (regra trapezoidal em h/2) ---
    b0 = 16/h^2;
    b1 = 4/h;
    b2 = 8/h;

    % --- Constantes do subpasso 2 (Euler regressivo de tres pontos) ---
    c1 =  1/h;
    c2 = -4/h;
    c3 =  3/h;

    Keff1 = K + b0*M + b1*C;
    Keff2 = K + c3^2*M + c3*C;
    dK1   = decomposition(Keff1);
    dK2   = decomposition(Keff2);

    U = zeros(n, nSteps+1);
    V = zeros(n, nSteps+1);
    A = zeros(n, nSteps+1);

    u = u0(:);  v = v0(:);
    a = M \ (-C*v - K*u);
    U(:,1) = u;  V(:,1) = v;  A(:,1) = a;

    % --- Laco no tempo ---
    for k = 1:nSteps
        % ----- Primeiro subpasso: regra trapezoidal, t -> t + h/2 -----
        R1 = M*(b0*u + b2*v + a) + C*(b1*u + v);
        uH = dK1 \ R1;
        vH = b1*(uH - u) - v;
        aH = b1*(vH - v) - a;

        % ----- Segundo subpasso: Euler regressivo de 3 pontos, t + h/2 -> t + h --
        R2 = -M*(c1*v + c2*vH + c3*c1*u + c3*c2*uH) - C*(c1*u + c2*uH);
        uNew = dK2 \ R2;
        vNew = c1*u + c2*uH + c3*uNew;
        aNew = c1*v + c2*vH + c3*vNew;

        u = uNew;  v = vNew;  a = aNew;
        U(:,k+1) = u;  V(:,k+1) = v;  A(:,k+1) = a;
    end

    t = (0:nSteps)*h;
end

%% ====================================================================
%  Runge-Kutta
%  ====================================================================

function [t, U, V, A] = mdof_rk(K, M, C, u0, v0, h, nSteps, order)
% Runge-Kutta explicito de ordem 1 a 4 (RK1 = Euler, RK2 = Euler modificado,
% RK3 = Heun, RK4 = classico). O sistema e escrito em espaco de estados,
%
%     yd = Ast*y,   y = {u; ud},   Ast = [0 I; -M\K -M\C]
%
% e o operador e a serie de Taylor truncada de exp(h*Ast), Eq. (14) do paper.
% O operador e montado uma vez e aplicado em cada passo.
%
% ENTRADA
%   K, M, C : [n x n] matrizes do sistema
%   u0, v0  : [n x 1] condicoes iniciais
%   h       : passo de tempo [s]
%   nSteps  : numero de passos
%   order   : 1, 2, 3 ou 4 (padrao 4)
%
% SAIDA
%   t       : [1 x nSteps+1] instantes
%   U, V, A : [n x nSteps+1] deslocamento, velocidade e aceleracao
    if nargin < 8 || isempty(order), order = 4; end
    if ~ismember(order, 1:4)
        error('mdof_rk:order','order must be 1, 2, 3 or 4.');
    end

    n = size(K,1);

    % --- Matriz de estado ---
    MiK = M \ K;
    MiC = M \ C;
    Ast = [ zeros(n), eye(n) ;
            -MiK    , -MiC   ];

    Ahat = rk_operador(Ast, h, order);

    % O passo do RK e limitado pela estabilidade e nSteps pode passar de 1e6;
    % V e A so sao alocados quando pedidos, para economizar memoria.
    computeV = nargout >= 3;
    computeA = nargout >= 4;

    U = zeros(n, nSteps+1);
    if computeV, V = zeros(n, nSteps+1); end
    if computeA, A = zeros(n, nSteps+1); end

    y = [u0(:); v0(:)];
    U(:,1) = y(1:n);
    if computeV, V(:,1) = y(n+1:end); end
    if computeA, A(:,1) = -MiK*y(1:n) - MiC*y(n+1:end); end

    % --- Laco no tempo ---
    for k = 1:nSteps
        y = Ahat*y;
        U(:,k+1) = y(1:n);
        if computeV, V(:,k+1) = y(n+1:end); end
        if computeA, A(:,k+1) = -MiK*y(1:n) - MiC*y(n+1:end); end
    end

    t = (0:nSteps)*h;
end

function Ahat = rk_operador(Ast, h, order)
% Operador do RK explicito: I + soma de (h*Ast)^i/i! ate a ordem do metodo.
    m    = size(Ast,1);
    Ahat = eye(m);
    term = eye(m);
    for i = 1:order
        term = (term*(h*Ast))/i;
        Ahat = Ahat + term;
    end
end

function [tSave, uSave, uMax] = rk_stream_tip(K, M, C, u0, v0, h, nSteps, order, tipExtract, saveIdx)
% Mesma recursao de mdof_rk, mas sem guardar o historico completo: so o valor
% da ponta nos passos escolhidos e o maximo de |ponta| em todos os passos.
% Com milhoes de passos, a matriz U de mdof_rk nao caberia na memoria.
%
% ENTRADA
%   K, M, C, u0, v0, h, nSteps, order : como em mdof_rk
%   tipExtract : handle coluna de estado -> valor na ponta
%   saveIdx    : indices de passo (0..nSteps) a guardar
%
% SAIDA
%   tSave, uSave : instantes e valores da ponta em saveIdx
%   uMax         : max(|ponta|) em todos os passos
    n = size(K,1);

    MiK = M \ K;
    MiC = M \ C;
    Ast = [ zeros(n), eye(n) ;
            -MiK    , -MiC   ];

    Ahat = rk_operador(Ast, h, order);

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

%% ====================================================================
%  Selecao do metodo
%  ====================================================================

function varargout = integra_metodo(name, K, M, C, u0, v0, h, nSteps)
% Chama o integrador pelo nome. Aceita a sigla ('DC') ou o nome por extenso.
% SAIDA: [t, U, V, A], com tantas saidas quantas forem pedidas.
    args = {K, M, C, u0, v0, h, nSteps};
    switch name
        case {'DC','Central Differences'}, [varargout{1:nargout}] = mdof_dc(args{:});
        case 'Houbolt',                    [varargout{1:nargout}] = mdof_houbolt(args{:});
        case 'Newmark',                    [varargout{1:nargout}] = mdof_newmark(args{:});
        case {'Wilson','Wilson-Theta'},    [varargout{1:nargout}] = mdof_wilson(args{:});
        case 'Bathe',                      [varargout{1:nargout}] = mdof_bathe(args{:});
        case {'RK1','RK2','RK3','RK4'},    [varargout{1:nargout}] = mdof_rk(args{:}, str2double(name(3)));
        otherwise, error('integra_metodo:name','Unknown method "%s".', name);
    end
end

function [dirNames, dirLabel, rkNames] = nomes_metodos()
% Siglas e rotulos dos nove metodos, na ordem usada em todas as figuras e tabelas.
    dirNames = {'DC','Houbolt','Newmark','Wilson','Bathe'};
    dirLabel = {'Central Differences','Houbolt','Newmark','Wilson-Theta','Bathe'};
    rkNames  = {'RK1','RK2','RK3','RK4'};
end

%% ====================================================================
%  Operador, raio espectral, PE/AD e passo critico
%  ====================================================================

function Ahat = operador(method, wn, xi, h)
% Operador de aproximacao Ahat de cada metodo para a equacao modal
%
%     xdd + 2*xi*wn*xd + wn^2*x = 0      (Secao 2.2 e Apendice A do paper)
%
% Isso basta para a viga porque, com Rayleigh, M, C e K diagonalizam juntas
% nos modos: o operador MDOF e bloco-diagonal na base modal e seu raio
% espectral e o maximo dos raios escalares. validacao confere isso.
%
% ENTRADA
%   method : 'DC','Houbolt','Newmark','Wilson','Bathe','RK1'..'RK4'
%   wn, xi : frequencia [rad/s] e amortecimento do modo
%   h      : passo de tempo [s]
%
% SAIDA
%   Ahat   : operador 2x2 ou 3x3
    switch upper(method)

        case 'DC'       % ---- Central Differences ----
            Ahat = [ (2 - wn^2*h^2)/(1 + xi*wn*h), -(1 - xi*wn*h)/(1 + xi*wn*h);
                      1,                            0 ];

        case 'HOUBOLT'  % ---- Houbolt, Eqs. (A.1)-(A.2) ----
            bt = (2/(wn^2*h^2) + 11*xi/(3*wn*h) + 1)^(-1);
            ki = (xi*bt)/(wn*h);
            Ahat = [ 5*bt/(wn^2*h^2) + 6*ki, -(4*bt/(wn^2*h^2) + 3*ki), bt/(wn^2*h^2) + 2*ki/3;
                     1, 0, 0;
                     0, 1, 0 ];

        case 'NEWMARK'  % ---- Newmark, alpha = 1/4, delta = 1/2 ----
            a  = 1 + xi*wn*h + (wn^2*h^2)/4;
            c1 = -wn^2/a;
            c2 = -(2*xi*wn + wn^2*h)/a;
            c3 = -(xi*wn*h + (wn^2*h^2)/4)/a;
            Ahat = [ 1 + (h^2/4)*c1, h + (h^2/4)*c2, (h^2/4)*(1 + c3);
                         (h/2)*c1,   1 + (h/2)*c2,   (h/2)*(1 + c3);
                              c1,          c2,             c3        ];

        case 'WILSON'   % ---- Wilson-theta, theta = 1.4, Eqs. (A.3)-(A.4) ----
            th = 1.4;
            b3 = (th/(wn^2*h^2) + xi*th^2/(wn*h) + th^3/6)^(-1);
            k3 = (xi*b3)/(wn*h);
            Ahat = [ 1 - b3*th^2/3 - 1/th - k3*th,        (-b3*th - 2*k3)/h,      -b3/h^2;
                     h*(1 - 1/(2*th) - b3*th^2/6 - k3*th/2), 1 - b3*th/2 - k3,   -b3/(2*h);
                     h^2*(1/2 - 1/(6*th) - b3*th^2/18 - k3*th/6), h*(1 - b3*th/6 - k3/3), 1 - b3/6 ];

        case 'BATHE'    % ---- Bathe, Eqs. (A.5)-(A.6) ----
            b1 = 16 + 8*xi*wn*h + wn^2*h^2;
            b2 =  9 + 6*xi*wn*h + wn^2*h^2;
            Ahat = (1/(b1*b2)) * [ ...
                -4*wn*h*(24*xi + 7*wn*h), ...
                 wn*(-288*xi + 14*xi*wn^2*h^2 - 144*wn*h + 5*wn^3*h^3 + 48*xi^2*wn*h), ...
                 wn^2*(24*xi*wn*h + 19*wn^2*h^2 - 144); ...
                -4*h*(-12 + wn^2*h^2), ...
                 144 - 47*wn^2*h^2 - 8*xi*wn^3*h^3 - 24*xi*wn*h, ...
                 wn^2*h*(-96 - 24*xi*wn*h + wn^2*h^2); ...
                 4*h^2*(7 + 2*xi*wn*h), ...
                 h*(144 - 5*wn^2*h^2 + 80*xi*wn*h + 16*xi^2*wn^2*h^2), ...
                -19*wn^2*h^2 + 144 + 168*xi*wn*h + 48*xi^2*wn^2*h^2 - 2*xi*wn^3*h^3 ];

        case {'RK1','RK2','RK3','RK4'}   % ---- Explicit RK, Eq. (14) ----
            j = str2double(method(3));
            A = [0 1; -wn^2 -2*xi*wn];
            Ahat = eye(2);
            term = eye(2);
            for i = 1:j
                term = (term*(h*A))/i;
                Ahat = Ahat + term;
            end

        otherwise
            error('viga_operador:method','Unknown method "%s".', method);
    end
end

function r = rho_beam(method, w, xi, h)
% Raio espectral do operador MDOF: maximo dos raios escalares sobre os modos
% (exato para amortecimento proporcional).
%
% ENTRADA
%   method : nome do metodo, como em operador
%   w, xi  : frequencias e amortecimentos dos modos do caso
%   h      : passo de tempo [s]
    r = 0;
    for k = 1:numel(w)
        rk = max(abs(eig(operador(method, w(k), xi(k), h))));
        if rk > r, r = rk; end
    end
end

function [PE, AD, rho, lambda] = pe_ad(method, wn, xi, h)
% Alongamento de periodo (PE) e decaimento de amplitude (AD), Secao 3.3.
%
% Sem amortecimento (xi = 0), Eqs. (19)-(21), PE e AD em PORCENTAGEM:
%     xibar      = -ln(rho)/Omega_n
%     Omegabar_d = Omega_n*sqrt(1 - xibar^2)
%     PE = (Omega_n/Omegabar_d - 1)*100
%     AD = (1 - exp(-2*pi*xibar*Omega_n/Omegabar_d))*100
%
% Com amortecimento (0 < xi < 1), Eqs. (23)-(26), PE e AD como RAZOES, em que
% 1 significa reproduzir exatamente o periodo e o decaimento fisicos:
%     PE = Omega_d/angle(lambda)
%     AD = exp(-2*pi*xibar*Omega_n/angle(lambda)) / exp(-xi*wn*T_d)
%
% As formulacoes do paper para RK (Eqs. 23-24) e para os diretos (Eqs. 25-26)
% sao algebricamente iguais, entao um unico ramo atende as duas familias.
%
% PE e AD voltam NaN quando nao estao definidos, e rho continua valido:
%   (i)   xi >= 1: o modo fisico nao oscila e nao ha periodo para comparar;
%         na viga isso acontece do modo 7 em diante (xi chega a 18.4);
%   (ii)  angle(lambda) ~ 0: a resposta numerica nao oscila;
%   (iii) |xibar| >= 1: o metodo esta muito fora da faixa estavel.
%
% SAIDA
%   PE, AD : conforme acima
%   rho    : raio espectral de Ahat
%   lambda : autovalor dominante de Ahat
    Ahat = operador(method, wn, xi, h);

    e = eig(Ahat);
    [rho, idx] = max(abs(e));
    lambda = e(idx);

    Omega_n = wn*h;
    th      = angle(lambda);

    % Amortecimento numerico introduzido pelo metodo, Eq. (19)
    xibar = -log(rho)/Omega_n;

    PE = NaN;  AD = NaN;

    % --- (i) modo fisico superamortecido ou com amortecimento critico --
    % Testado PRIMEIRO, antes de qualquer uso de sqrt(1-xi^2).
    if xi >= 1
        return
    end

    % --- (ii) a resposta numerica nao oscila ---------------------------
    if th <= eps
        return
    end

    % --- (iii) amortecimento numerico fora da faixa subamortecida ------
    if abs(xibar) >= 1
        return
    end

    % --- Definicoes -----------------------------------------------------
    if xi == 0
        % Sem amortecimento, Eqs. (20)-(21), em porcentagem
        Omegabar_d = Omega_n*sqrt(1 - xibar^2);
        PE = (Omega_n/Omegabar_d - 1)*100;
        AD = (1 - exp(-2*pi*xibar*Omega_n/Omegabar_d))*100;
    else
        % Com amortecimento, Eqs. (23)-(26), como razoes. Chegar aqui
        % garante 0 < xi < 1, entao wd e Td sao reais.
        wd      = wn*sqrt(1 - xi^2);
        Td      = 2*pi/wd;
        Omega_d = wd*h;
        PE = Omega_d/th;
        AD = exp(-2*pi*xibar*Omega_n/th) / exp(-xi*wn*Td);
    end
end

function pead = pead_tabela(w, xi, allNames, allLabels, allSteps, T1, varargin)
% Calcula e imprime a tabela de PE e AD de um caso (w, xi) qualquer.
%
% ENTRADA
%   w, xi     : [m x 1] frequencias e amortecimentos do caso
%   allNames  : metodos, ex. [dirNames, rkNames]
%   allLabels : rotulos impressos
%   allSteps  : passo de cada metodo
%   T1        : periodo de referencia da coluna h/T1
%   'nModesReport' : modos na tabela resumida (padrao 3)
%   'printHeader'  : texto acima da tabela
%
% SAIDA
%   pead : .methods, .labels, .h, .modes, .w, .xi, .PE, .AD, .rho
%          ([nMetodos x nModos], todos os modos) e .table (modos do resumo)
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
                pe_ad(allNames{k}, w(m), xi(m), hk);
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

function [hcr, rhoFun] = h_critico(K, M, C, method, hHi, tolRho)
% Passo critico de um metodo condicionalmente estavel, por bissecao em
% rho(Ahat(h)) <= 1, usando o operador de bloco completo de ordem 2n.
%
% Na malha, o limite e ditado pela MAIOR frequencia do sistema discreto, por
% isso os valores ficam ordens de grandeza abaixo do limite h/Tn de 1 GDL.
% Referencias sem amortecimento: DC = 2/w_max, RK3 = sqrt(3)/w_max,
% RK4 = sqrt(8)/w_max.
%
% ENTRADA
%   K, M, C : matrizes do sistema
%   method  : 'RK1','RK2','RK3','RK4' ou 'DC'
%   hHi     : limite superior da busca (opcional, padrao 1e-1)
%   tolRho  : tolerancia em rho = 1 (opcional, padrao 1e-8; tolerancia zero
%             nao funciona porque o operador fica mal condicionado com h -> 0)
%
% SAIDA
%   hcr    : passo critico [s]; 0 se nenhum passo estavel for encontrado
%   rhoFun : handle h -> rho(Ahat(h))
    if nargin < 5 || isempty(hHi),    hHi    = 1e-1; end
    if nargin < 6 || isempty(tolRho), tolRho = 1e-8; end

    n = size(K,1);

    switch upper(method)
        case {'RK1','RK2','RK3','RK4'}
            order = str2double(method(3));
            Ast = [zeros(n), eye(n); -(M\K), -(M\C)];
            rhoFun = @(h) rhoRK(Ast, h, order);

        case 'DC'
            rhoFun = @(h) rhoCD(K, M, C, h);

        otherwise
            error('viga_h_critico:method', ...
                'Unsupported method "%s". Use RK1..RK4 or DC.', method);
    end

    % --- Ja e estavel no limite superior? ---
    if rhoFun(hHi) <= 1 + tolRho
        hcr = hHi;
        warning('viga_h_critico:bracket', ...
            'Method %s still stable at h = %g; enlarge hHi.', method, hHi);
        return
    end

    % --- Busca para baixo um limite inferior estavel ---
    hLo   = hHi;
    found = false;
    for it = 1:200
        hLo = hLo/2;
        if rhoFun(hLo) <= 1 + tolRho
            found = true;
            break
        end
    end
    if ~found
        hcr = 0;                     % instavel para todo passo pratico
        return
    end

    % --- Bissecao ---
    for it = 1:200
        hMid = 0.5*(hLo + hHi);
        if rhoFun(hMid) <= 1 + tolRho
            hLo = hMid;
        else
            hHi = hMid;
        end
    end
    hcr = hLo;
end

function r = rhoRK(Ast, h, order)
% Raio espectral do RK explicito.
    r = max(abs(eig(rk_operador(Ast, h, order))));
end

function r = rhoCD(K, M, C, h)
% Recorrencia de dois passos das Diferencas Centrais
%   (a0*M + a1*C) u_{n+1} = (a2*M - K) u_n + (a1*C - a0*M) u_{n-1}
% escrita como operador de bloco de primeira ordem sobre {u_{n+1}; u_n}.
    n  = size(K,1);
    a0 = 1/h^2;  a1 = 1/(2*h);  a2 = 2*a0;
    Meff = a0*M + a1*C;
    B1 = Meff \ (a2*M - K);
    B2 = Meff \ (a1*C - a0*M);
    Ahat = [B1, B2; eye(n), zeros(n)];
    r = max(abs(eig(Ahat)));
end
