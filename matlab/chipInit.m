function [Sigma, beta, mu]=chipInit(data,X);
%CHIPINIT initialises the chip parameters
g=cov(data');
Sigma=X'*diag(diag(g))*diag(diag(X*X').^(-1))*X;
%mu=pdinv(Sigma)*X'*mean(data,2);
mu=zeros(size(X,2),1);
beta=100;