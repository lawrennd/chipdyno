function [list,maxActivity,maxActivityError]=chipDynoGeneAct(data, ...
                                                  X,Sigma,beta,gamma,mu, ...
                                                  transNames, ...
                                                  annotation,geneName);
%CHIPDYNOGENEACT given a gene, lists activators in decreasing order

%CHIPDYNO

I=find(strcmp(geneName,annotation));
activeNames=transNames(find(X(I,:)));
nTransFact=sum(X(I,:));
maxActivity=[];
maxActivityError=[];
for i=1:nTransFact
  [tf,tfError]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotation, ...
                                         activeNames(i),geneName);
  ind=find(strcmp(activeNames(i),transNames));
  tf=tf-mu(ind)*ones(size(tf));
  [act,index]=max(tf);
  maxActivity=[maxActivity,act];
  maxActivityError=[maxActivityError,tfError(index)];
end
[maxActivity,index]=sort(maxActivity,'descend');
maxActivityError=maxActivityError(index);
list=activeNames(index);