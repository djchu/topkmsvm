% opts.prox = topk_cone_biased
clc;

% a = load('randomdata');
a = load('topkprojection\build\randomdata');

dim = length(a);

% opts = [];
% opts.prox = 'topk_cone_biased';
% opts.k = k;
% opts.rhs = r;
% opts.rho = rho;
% B_topk = matsdca_prox(a, opts)

addpath('../');
cvx_begin
    variable x(dim);
    minimize((a-x)'*(a-x)+ rho * sum(x)*sum(x))
    subject to
        sum(x) <= r
        0 <= x <= sum(x)/k
cvx_end

[x B_topk]

xCPP_newton = proj_generalmex(a, k, rho);

figure;
plot_fg(a, rho, k, r);
a_sort = sort(a, 'descend');
s_min = sum(a_sort(1:k))/(1+k*rho);
t_max = a_sort(k) - s_min/k;
t_min = a_sort(k+1) - s_min/k;

s0 = s_min;
t0 = t_max;
[t_newton, s_newton, iter_newton] = newtonraphson(a, k, rho, s0, t0);
[t_cd, s_cd] = coordinatedescent(a, k, rho);
x_newton = min(max(0, a-t_newton), s_newton/k);
x_cd = min(max(0, a-t_cd), s_cd/k);
[a x x_newton B_topk]