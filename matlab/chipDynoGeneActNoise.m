function [list,maxActivity,maxActivityError]=chipDynoGeneActNoise(data, ...
                                                  X,Sigma,beta,precs,gamma,mu, ...
                                                  TransNames, ...
                                                  annotation,geneName);
% CHIPDYNOGENEACTNOISE given a gene, lists activators in decreasing order

% CHIPDYNO

I=find(strcmp(geneName,annotation));
activeNames=TransNames(find(X(I,:)));
nTransFact=sum(X(I,:));
maxActivity=[];
maxActivityError=[];
for i=1:nTransFact
  [tf,tfError]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation, ...
                                         activeNames(i),geneName);
  [act,index]=max(tf-mu(find(strcmp(activeNames(i),TransNames)))*ones(size(tf)));
  maxActivity=[maxActivity,act];
  maxActivityError=[maxActivityError,tfError(index)];
end
[maxActivity,index]=sort(maxActivity,'descend');
maxActivityError=maxActivityError(index);
list=activeNames(index);