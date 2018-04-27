% opts.prox = topk_cone_biased
clc;

% load a_cdwrong.mat;
seed = rng;
% rng(seed);
a = randn(10,1);
dim = length(a);
% load t_isnotunique.mat

rho = 1.2;
k = 2;
r=0.93;

save debug.mat a rho k r;

addpath('/home/freda/Mycodes/topksvm/src/');
opts = [];
opts.prox = 'topk_cone_biased';
opts.k = k;
opts.rhs = r;
opts.rho = rho;
B_topk = matsdca_prox(a, opts);