function [ksd, w, wksd] = comp_ksd(X, G, Gam)
% [ksd, w, wksd] = comp_ksd(X, G, Gam)
% Inputs:
% X    - n x d array, each row a sample from MCMC.
% G    - n x d array, each row the gradient of the log target.
% Gam  - d x d positive definite matrix, the preconditioner in the kernel.
%
% Outputs:
% ksd  - scalar, the KSD of the unweighted empirical distribution based on X.
% w    - n x 1 vector, containing the optimal weights for the point set X.
% wksd - scalar, the KSD of the optimally weighted empirical distribution
%        based on X.

% dimensions
[n,~] = size(X);

% vectorised computation of Stein kernel matrix
tmp0 = trace(inv(Gam));
tmp1 = dot(repmat(G,n,1),repelem(G,n,1),2);
tmp2 = repmat(X,n,1) - repelem(X,n,1);
tmp3 = repmat(G,n,1) - repelem(G,n,1);
tmp4 = (Gam \ (tmp2'))';
tmp5 = - 3 * dot(tmp4,tmp4,2) ./ ((1 + dot(tmp2,tmp4,2)).^(5/2)) ...
       + (tmp0 + dot(tmp3,tmp4,2)) ./ ((1 + dot(tmp2,tmp4,2)).^(3/2)) ...
       + tmp1 ./ ((1 + dot(tmp2,tmp4,2)).^(1/2));
K = reshape(tmp5,n,n);

% compute outputs
ksd = (1/n) * sqrt(ones(1,n) * K * ones(n,1));
tmp6 = K \ ones(n,1);
w = tmp6 / sum(tmp6);
wksd = 1 / sqrt(sum(tmp6));

end
