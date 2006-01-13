function [list,maxActivity,maxActivityError]=chipDynoGeneAct(data, ...
                                                  X,Sigma,beta,tau,mu, ...
                                                  transNames, ...
                                                  annotation,geneName);
%CHIPDYNOGENEACT given a gene, lists activators in decreasing order

%CHIPDYNO

I=find(strcmp(geneName,annotation));
activeNames=transNames(find(X(I,:));
nTransFact=sum(X(I,:));
maxActivity=[];
maxActivityError=[];
for i=1:nTransFact
  [tf,tfError]=chipDynoExpectationsFast(data,X,Sigma,beta,tau,mu, ...
                                         transNames, annotations, ...
                                         activeNames(i),geneName);
  [act,index]=max(tf);
  maxActivity=[maxActivity,act];
  maxActivityError=[maxActivityError,tfError(index)];
end
[maxActivity,index]=sort(maxAcitvity,'descend');
maxActivityError=maxActivityError(index);
list=activeNames(index);