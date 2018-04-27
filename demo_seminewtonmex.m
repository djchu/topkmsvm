%% Demo script for our semismooth Newton method vs. sorting-based method

clc, clear;
warning off;
seed = rng;

dim = 10000;
k = 3;  % k<dim
rho = 1;
r = 1;

a = randn(dim, 1);

tic
[x_newtonMEX, t_newtonMex] = proj_seminewtonmex(a, k, rho);
toc;

opts = [];
opts.prox = 'topk_cone_biased';
opts.k = k;
opts.rhs = r;
opts.rho = rho;
tic
x_NIPS = matsdca_prox(a, opts);
toc

% [x_newtonMEX x_NIPS]

