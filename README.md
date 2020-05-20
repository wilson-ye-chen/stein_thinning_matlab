# Stein Thinning
This MATLAB code implements an algorithm for optimally compressing
sampling algorithm outputs by minimising a kernel Stein discrepancy.
Please see the accompanying paper "Optimal Thinning of MCMC Output"
([arXiv](https://arxiv.org/pdf/2005.03952.pdf)) for details of the
algorithm.

# Getting Started
For example, correlated samples from a posterior distribution are
obtained using a MCMC algorithm and stored in the matrix `smpl`,
and the corresponding gradients of the log-posterior are stored in
another matrix `grad`. One can then perform Stein Thinning to obtain
a subset of 40 sample points by running the following code:
```matlab
idx = thin(smpl, grad, 40)
```
The `thin` function returns a vector containing the row indices in
`smpl` (and `grad`) of the selected points. Please refer to `demo.m`
as a starting example.

The default usage requires no additional user input and is based on
the `sclmed` heuristic. Alternatively, the user can choose to specify
which heuristic to use for computing the preconditioning matrix by
setting the option string `pre` to either `med`,  `sclmed`, `smpcov`,
`bayesian`, or `avehess`. For example, the default setting corresponds
to:
```matlab
idx = thin(smpl, grad, 40, pre='sclmed')
```
The details for each of the heuristics are documented in Section 3.4 of
the accompanying paper.
