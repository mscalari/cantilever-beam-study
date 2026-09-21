function [w, Phi] = viga_autovalores_EF(K, M)
%% Murilo Barros Scalari - RA 256830 %%
%   PAPER - Time Integration Schemes   %
%--------------------------------------%
%  Modal properties of the discrete    %
%  (semi-discrete) FE beam system.     %
%--------------------------------------%
%
% Solves the generalised eigenvalue problem  K*phi = w^2 * M*phi.
%
% INPUT
%   K, M : symmetric positive definite FE matrices
%
% OUTPUT
%   w    : [n x 1] natural circular frequencies, ascending   [rad/s]
%   Phi  : [n x n] mass-normalised eigenvectors (Phi'*M*Phi = I)
%
% Mass normalisation makes the modal decomposition used by
% viga_solucao_exata.m straightforward and numerically well behaved.

    [V, D] = eig(full(K), full(M), 'chol');

    w2 = real(diag(D));
    w2(w2 < 0) = 0;                     % guard against round-off
    [w2, idx] = sort(w2, 'ascend');
    w   = sqrt(w2);
    Phi = real(V(:, idx));

    % --- Mass normalisation: phi' * M * phi = 1 ---
    for k = 1:size(Phi,2)
        Phi(:,k) = Phi(:,k) / sqrt(Phi(:,k).' * M * Phi(:,k));
    end
end
