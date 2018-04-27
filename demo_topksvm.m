%% demo for top-k svm optimization on example4_train dataset

clc,clear;
close all;
warning off;

global TEST;
TEST = 1;

load('./data/example4_train.mat');
xtrn = full(X);
ytrn = y';

[dim, num] = size(xtrn);

if TEST
    global ytst;
    global Xtst;
    load('./data/example4_test.mat' );
    Xtst = full(X);
    ytst = y';
end
clear X y;

lambda = 0.2;
k = 3;
epsilon = 1e-3;
Tmax = 1000;

c = 1/(num*lambda);

%% SDCA for top-k multiclass SVM
Method = 1;
model_ours1 = topksvm_sdca(xtrn, ytrn, lambda, k, epsilon, Tmax, Method);
Method = 2;
model_ours2 = topksvm_sdca(xtrn, ytrn, lambda, k, epsilon, Tmax, Method);
Method = 3;
model_ours3 = topksvm_sdca(xtrn, ytrn, lambda, k, epsilon, Tmax, Method);

timesMatrix1 = model_ours1.times;
timesMatrix2 = model_ours2.times;
timesMatrix3 = model_ours3.times;

[timesMatrix1 timesMatrix2 timesMatrix3]

