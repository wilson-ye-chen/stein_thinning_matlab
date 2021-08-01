function [pi, ksd] = thin(X, G, m, std, precon)
% [pi, ksd] = thin(X, G, m, std, precon)
%
% Version 1.0 of Stein Thinning.
%
% Inputs:
% X      - n x d array, each row a sample from MCMC.
% G      - n x d array, each row the gradient of the log target.
% m      - desired number of points.
% std    - optional logical, either 'true' (default) or 'false', indicating
%          whether or not to standardise the colums of X.
% precon - optional string, either 'id' (default), 'med', 'sclmed', or
%          'smpcov', specifying the preconditioner to be used.
%
% Outputs:
% pi     - m x 1 vector, containing the row indices in X of the selected
%          points.
% ksd    - m x 1 vector, with jth entry giving the KSD of the unweighted
%          empirical distribution based on the first j points selected.

% dimensions
[n,d] = size(X);

% error checking
if (n == 0) || (d == 0)
    error('X is empty')
end
if (size(G,1) ~= n) || (size(G,2) ~= d)
    error('Dimensions of X and G are inconsistent')
end
if (sum(isnan(X)) + sum(isnan(G))) > 0
    error('One of X or G contains NaNs')
end
if (sum(isinf(X)) + sum(isinf(G))) > 0
    error('One of X or G contains infs')
end

% defaults
if nargin == 3
    std = true;
    precon = 'id';
elseif nargin == 4
    precon = 'id';
end

% standardisation
if std == true
    scl = mad(X,0,1); % mean absolute deviation
    if min(scl) == 0
        error('Too few unique samples in X')
    end
    X = X ./ scl;
    G = G .* scl; % using the chain rule
end

% preconditioner
if strcmp(precon,'id')
    Gam = eye(d);
elseif strcmp(precon,'med')
    n0 = 1000;
    ix = unique(round(linspace(1, n, n0)));
    ell = median(pdist(X(ix,:)));
    if ell == 0
        error('Too few unique samples in X')
    end
    Gam = ell^2 * eye(d);
elseif strcmp(precon,'sclmed')
    n0 = 1000;
    ix = unique(round(linspace(1, n, n0)));
    ell = median(pdist(X(ix,:)));
    if ell == 0
        error('Too few unique samples in X')
    end
    Gam = ell^2 * eye(d) / log(min(n, n0));
elseif strcmp(precon,'smpcov')
    C = cov(X);
    evals = eig(C);
    if ~all(evals > 0)
        error('Too few unique samples in X')
    end
    Gam = C;
else
    error('Incorrect preconditioner type.');
end

% Stein kernel sub-matrix K(X,pi) (i.e., just store entries that we need)
K_X_pi = zeros(n,m);

% also compute the diagonal elements of the Stein kernel matrix
K_diag = d + dot(G,G,2);

% main loop
pi = zeros(m,1);
ksd = zeros(m,1);
tmp0 = X * inv(Gam);
for j = 1:m
    % monitor
    disp(['Selecting point ',num2str(j,'%u'),' of ',num2str(m,'%u')])

    % select next point
    if j == 1
        [~,pi(1)] = min((1/2) * K_diag);
    else
        [~,pi(j)] = min( (1/2) * K_diag + sum( K_X_pi(:,1:(j-1)) , 2) );
    end

    % populate row and column of kernel matrix associate with new point
    tmp1 = dot( tmp0 - repmat(X(pi(j),:) * inv(Gam),n,1)  , ...
                tmp0 - repmat(X(pi(j),:) * inv(Gam),n,1) , 2);
    tmp2 = dot( G - repmat(G(pi(j),:),n,1) , ...
                tmp0 - repmat(X(pi(j),:) * inv(Gam),n,1) , 2);
    tmp3 = dot( (X - repmat(X(pi(j),:),n,1)) , ...
                tmp0 - repmat(X(pi(j),:)* inv(Gam),n,1) , 2);
    tmp4 = dot( G , repmat(G(pi(j),:),n,1) , 2 );
    K_pi_j = - 3 * tmp1 ./ ((1 + tmp3).^(5/2)) ...
             + (trace(inv(Gam)) + tmp2) ./ ((1 + tmp3).^(3/2)) ...
             + tmp4 ./ ((1 + tmp3).^(1/2));
    K_X_pi(:,j) = K_pi_j;

    % compute kernel Stein discrepancy
    ksd(j) = (1/j) * sqrt( ones(1,j) * (K_X_pi(pi(1:j),1:j) * ones(j,1)) );
end

end
