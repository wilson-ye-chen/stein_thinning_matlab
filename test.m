
X = csvread('smp.csv');
G = csvread('scr.csv');

m = 40;
[pi, ksd] = thin(X, G, m);

% Plot point-set over trace
plot(X(:, 1), X(:, 2), '-')
hold on;
plot(X(pi, 1), X(pi, 2), '.', 'markersize', 30)

% Plot KSD
figure()
semilogy(1:m, ksd)
