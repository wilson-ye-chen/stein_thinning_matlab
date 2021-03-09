% Post-process output from Stan.


% Simple bivariate Gaussian model
mc = {
    'parameters {vector[2] x;}'
    'model {x ~ multi_normal([0, 0], [[1, 0.8], [0.8, 1]]);}'
};

fit = stan('model_code', mc, 'iter', 1000);

% Extract sampled points and gradients
% smp = fit['x']
% scr = np.apply_along_axis(fit.grad_log_prob, 1, smp)

% Obtain a subset of 40 points
% idx = thin(smp, scr, 40)

% Plot point-set over trace
% figure()
% plot(smp[:,0], smp[:,1], color=(0.4, 0.4, 0.4), linewidth=1)
% plot(smp[idx, 0], smp[idx, 1], 'r.', markersize=10)
