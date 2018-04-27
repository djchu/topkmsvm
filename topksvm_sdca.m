function model = topksvm_sdca(xtrn, ytrn, lambda, k, epsilon, Tmax, Method)
% topksvm_sdca: SDCA solver for solving top-k multiclass SVM problem
% INPUT
% xtrn: dim-by-num data matrix.
% ytrn: 1-by-num label vector
% lambda: the regularization parameter
% k: top-k
%
% OUTPUT
% model.A: dual variable
% model.W: primal variable
% model.pobj: primal objective function value
% model.dobj: dual objective function value
% model.accuracy: accuracy

global TEST;

if nargin < 7
    Method = 1;
end
if nargin < 6
    Tmax = 1000;
end
if nargin < 5
    epsilon = 1e-3;
end
if nargin < 4
    k = 1;
end

[dim, num] = size(xtrn);
nY = length(unique(ytrn));
Dim = dim*nY;

W = zeros(Dim, 1);
Alpha = zeros(nY, num);

accuracy = 0;
if TEST
    global ytst;
    global Xtst;
    [~,accuracy] = mypredict(Xtst, ytst, reshape(W, dim, nY));
end

fprintf('SDCA for top-k multiclass SVM\n' );

c = ones(nY, num);  % Iyi = 1 - e_yi ==> c_i
rhos = zeros(1, num);  % rho(i) = xi'*xi/(num*lambda)

for i=1:num
    c(ytrn(i),i) = 0;
    rhos(i) = xtrn(:,i)'*xtrn(:,i)/(num*lambda);
end

%{
  The whole optimization: biased projection problem + knapsack probelm.
  timeVec(1,:): save time for biased projection problem, i.e. Semismooth
                Newton or sorting based method;
  timeVec(2,:): save time for knapsack problem, i.e. Newton method or fixing
                variable method
  total_times(1): save total time for biased projection problem;
  total_times(2): save total time for knapsack probelm.
%}
total_times = zeros(2,1);
pobj = zeros(1, Tmax);
dobj = zeros(1, Tmax);
accuracy = zeros(1, Tmax);

rng('shuffle');
for iter = 1: Tmax
    InxRand = randperm(num);
    timeVec = zeros(2,num);
    
    for i=1:num
        inx = InxRand(i);
        
        ci = c(:,inx);
        eyi = zeros(nY,1);
        eyi(ytrn(inx)) = 1;
        Xiai_old = xtrn(:,inx) * (Alpha(:,inx) - sum(Alpha(:,inx)).*eyi)';
        
        What = reshape(W, dim, nY) - Xiai_old./(num*lambda);
        XiWhat = What' * xtrn(:,inx);
        XiWhat = XiWhat - XiWhat(ytrn(inx));
        
        rho = rhos(inx);
        a = (ci + XiWhat)./rho;
        yi = ytrn(inx);
        a(yi) = [];
        
        r = 1;
        switch Method
            case 1
                %% Method I: my mex implementation
                % Our Semi-smooth Newton + Newton, warm start and without
                % fixing variable strategy in knapsack problem.
                tic
                [z_topk1, t_newton1] = proj_seminewtonmex(a, k, r);
                timeVec(1,i) = toc;
                
                tic
                if sum(z_topk1)>r
                    t0 = t_newton1; % warm start in knapsack problem
                    FIX = 0;
                    z_topk1 = proj_knapsackmex(a, k, r, t0, FIX);
                end
                timeVec(2,i) = toc;
                
                z = z_topk1;
                
            case 2
                %% Method II: my mex implementation
                % Our Semi-smooth Newton + NIPS15's Variable Fixing Algorithm
                tic;
                z_topk2 = proj_seminewtonmex(a, k, r);
                timeVec(1,i) = toc;
                
                tic
                if sum(z_topk2)>r
                    opts = [];
                    opts.prox = 'knapsack';
                    opts.k = k;
                    opts.hi = 1/k;
                    opts.rhs = 1;
                    
                    z_topk2 = matsdca_prox(a, opts);
                end
                timeVec(2,i) = toc;
                z = z_topk2;
                
            case 3
                %% Method III: NIPS'15 implementation
                % Variable Fixing Algorithm + Sorting based Algorithm
                tic
                opts = [];
                opts.prox = 'knapsack';
                opts.lo = 0;
                opts.hi = r/k;
                opts.rhs = r;
                z_topk3 = matsdca_prox(a, opts);
                timeVec(2,i) = toc;
                
                [z_topktmp, tstar] = proj_knapsackmex(a, k, r, 0, 0);
                if(norm(z_topktmp-z_topk3)>1e-6)
                    disp('Something is wrong in Method 3.');
                end
                
                tic
                newi = max(0, a-tstar-r/k);
                check_lam = tstar+sum(newi)/k-1*r;
                if check_lam <0
                    opts = [];
                    opts.prox = 'topk_simplex_biased';
                    opts.k = k;
                    opts.rho = 1;
                    opts.rhs = r;
                    
                    z_topk3 = matsdca_prox(a, opts);
                end
                timeVec(1,i) = toc;
                z = z_topk3;
        end
        
        alphai = zeros(nY, 1);
        alphai(1:yi-1) = -z(1:yi-1);
        alphai(yi+1:end) = -z(yi:end);
        
        %% update dual variable and primal variable
        Alpha(:,inx) = alphai;
        Xiai = xtrn(:,inx) * (Alpha(:,inx) - sum(Alpha(:,inx)).*eyi)';
        W = W + (Xiai(:)-Xiai_old(:)) / (num*lambda);
    end
    
    total_times(1) = total_times(1) + sum(timeVec(1,:));
    total_times(2) = total_times(2) + sum(timeVec(2,:));
    
    pobj(iter) = topksvm_pfv(xtrn, ytrn, lambda, k, reshape(W, dim, nY));
    dobj(iter)   = topksvm_dualfvmex(xtrn(:), ytrn, lambda, Alpha(:), nY, dim, num);
    
    if TEST
        [~, accuracy(:, iter)] = mypredict(Xtst, ytst, reshape(W, dim, nY));
    end

    if( (pobj(iter)-dobj(iter))/pobj(iter) < epsilon...
            || pobj(iter)-dobj(iter) <epsilon )
        disp('The stopping criterion has been reached.');
        pobj(iter+1:end) = [];
        dobj(iter+1:end) = [];
        accuracy(iter+1:end) = [];
        break;
    end
end

wdual_sdca = reshape(W, dim, nY);

model.iter = iter;
model.A = Alpha;
model.W = wdual_sdca;
model.accuracy = accuracy;
model.pobj = pobj;
model.dobj = dobj;
model.times = total_times;

end
