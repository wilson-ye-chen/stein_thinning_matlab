# Stein Thinning MATLAB
This MATLAB code implements an algorithm for optimally compressing
sampling algorithm outputs by minimising a kernel Stein discrepancy.
Please see the accompanying paper "Optimal Thinning of MCMC Output"
([arXiv](https://arxiv.org/pdf/2005.03952.pdf)) for details of the
algorithm.

# Installation
To install the toolbox, clone the repository `stein_thinning_matlab`
to a suitable directory and then add MATLAB path to that directory.
```
git clone https://github.com/wilson-ye-chen/stein_thinning_matlab.git
```

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
the identity (`id`) preconditioning matrix and standardised sample.
Alternatively, the user can choose to specify which heuristic to use
for computing the preconditioning matrix by setting the option string
to either `id`, `med`,  `sclmed`, or `smpcov`. Standardisation can be
disabled by setting the fourth argument to `false`. For example, the
default setting corresponds to:
```matlab
idx = thin(smpl, grad, 40, true, 'id')
```
The details for each of the heuristics are documented in Section 2.3 of
the accompanying paper.
