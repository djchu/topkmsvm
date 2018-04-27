function pfv = topksvm_pfv(xtrn, ytrn, lambda, k, W)
% Return the primal objective function value of Crammer and Singer's multiclass SVM
%
% Iy: nY x num
% x: dim x num
% lambda: regularization parameter
% k: top-k
% W: dim x nY

assert(1<=k && k<size(W, 2));

[dim, num] = size(xtrn);
nY = length(unique(ytrn));

if nargin < 5
    W = zeros(dim, nY);
end
if nargin < 4
    k = 1;
end

scores = W'*xtrn;
for i=1:numel(ytrn)
   scores(:,i) = 1+scores(:,i) - scores(ytrn(i),i);
   if k==1
      scores(ytrn(i),i) = 0;
   else	
      scores(ytrn(i),i) = -inf;    %% Lapin's setting
   end
end
scores_sort = sort(scores, 'descend');

loss = max(0, sum(scores_sort(1:k,:))/k);
pfv = lambda/2*W(:)'*W(:) + sum(loss)/num;

end
