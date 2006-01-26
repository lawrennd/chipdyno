function [tf,tfErrors]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         transNames, annotations, ...
                                         transName,geneName);
% CHIPDYNOEXPECTATIONSFASTNOISE computes posterior expectations TFA.

% CHIPDYNO

npts=size(data,2);
nTrans=size(X,2);

c=class(geneName);
if c(1)=='c'
    x=X(find(strcmp(geneName,annotations)),:)';
    data=data(find(strcmp(geneName,annotations)),:);
    precs=precs(find(strcmp(geneName,annotations)),:);
elseif c(1)=='d'
    x=X(geneName,:)';
    data=data(geneName,:);
    precs=precs(geneName,:);
else
    error('Genes can be identified either by number or name\n')
end
expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);
index=find(strcmp(transName,transNames));
if x(index)==0
  error('The gene selected is not a target of the transcription factor \n')
end
tf=expectations.b(:,index);
tfErrors=sqrt(expectations.bTb(:,index));