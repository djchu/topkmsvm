%
% Demo script for testing matsdca mex files
%
clear;
rng(0);

% Test prox
if 0
  fprintf('Demo: prox\n');
  d = 100;
  n = 5;

  A = randn(d,n);

  opts = [];
  opts.prox = 'knapsack';
  B = matsdca_prox(A, opts);

  fprintf('All columns should sum up to 1.\n');
  disp(sum(B));

  % Further details:
  % matsdca_prox('help')

end


% Test solver
if 1
  fprintf('Demo: solver\n');
%   num_trn = 200;
%   num_val = 200;
%   num_tst = 200000;

%   [Xtrn,Ytrn,Xval,Yval,Xtst,Ytst] = getdata(num_trn,num_val,num_tst);
  [ytrn, xtrn] = libsvmread('../example4_train.light');
  [ytst, xtst] = libsvmread( '../example4_test.light' );
  [dim, num] = size(xtrn');
  lambda = 1e-2;

  opts = [];
  opts.objective = 'l2_hinge_topk';
%   opts.objective = 'softmax'; % or 'msvm_smooth'
%   opts.objective = 'l2_multiclass_hinge';
%   opts.objective = 'l2_topk_hinge';
  
  opts.c = 1/(num*lambda);
  opts.k = 2;
  opts.epsilon = 1e-8;
  opts.max_epoch = 1000;
  opts.eval_epoch = 1; % only for demo, in practice would be wasteful
%   opts.log_level = 'verbose'; % or 'none'

%   model = matsdca_fit({Xtrn, Xval, Xtst}, {Ytrn, Yval, Ytst}, opts);
  model = matsdca_fit({full(xtrn'), full(xtst')}, {ytrn', ytst'}, opts);
  disp(model);
  fprintf('trn accuracy: %.2f\n', 100*model.train(end).accuracy);
%   fprintf('val accuracy: %.2f\n', 100*model.test(end,1).accuracy);
%   fprintf('tst accuracy: %.2f\n', 100*model.test(end,2).accuracy);
  fprintf('tst accuracy: %.2f\n', 100*model.test(end,1).accuracy);
  
  pfv = [];
  dfv = [];
  itermax = 1000;
  for i=1:itermax
      trn = model.train;
      pfv(i) = trn(i).primal*lambda;
      dfv(i) = trn(i).dual*lambda;
  end
  figure;
  plot(1:itermax, pfv, 'b--');
  hold on;
  plot(1:itermax, dfv, 'k--');

  % Further details:
  % matsdca_fit('help')

end
