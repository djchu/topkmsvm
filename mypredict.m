function [ypred, top1_acc, topk_acc] = mypredict(testX, testy, W, k)
% testing function for multiclass svm
% 
% testX: dim x num
% testy: 1 x num
% W: dim x k
% k: top-k accuracy

if nargin<4
    k = 1;
end
dim_tst = size(testX, 1);
dim = min(size(W,1), dim_tst);
W = W(1:dim,:);

dec_values = W'*testX;
[~, ypred] = max(dec_values);
top1_acc = 100*sum(ypred(:)==testy(:))/length(testy);

if nargout>2
    [~, ypred_topk] = sort(dec_values, 1, 'descend');
    topk_acc = bsxfun(@eq, ypred_topk(1:k,:), testy);
    topk_acc =  100*sum(sum(topk_acc))/length(testy);
end


