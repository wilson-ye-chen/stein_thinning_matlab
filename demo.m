
X = csvread('smp.csv');
G = csvread('scr.csv');

m = 40;
[idx, ksd] = thin(X, G, m);

% Plot point-set over trace
plot(X(:, 1), X(:, 2), '-')
hold on;
plot(X(idx, 1), X(idx, 2), '.', 'markersize', 30)

% Plot KSD path
figure()
semilogy(1:m, ksd)
