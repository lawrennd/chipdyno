function [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation, ...
                                         transName,geneName);

% CHIPDYNOEXPECTATIONSFASTNOISE computes posterior expectations TFA.
%
%	Description:
%	[tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
%                                         TransNames, annotation, ...
%                                         transName,geneName);
%% 	chipDynoExpectationsFastNoise.m version 1.4


npts=size(data,2);
nTrans=size(X,2);

c=class(geneName);
if c(1)=='c'
    x=X(find(strcmp(geneName,annotation)),:)';
    data=data(find(strcmp(geneName,annotation)),:);
    precs=precs(find(strcmp(geneName,annotation)),:);
elseif c(1)=='d'
    x=X(geneName,:)';
    data=data(geneName,:);
    precs=precs(geneName,:);
else
    error('Genes can be identified either by number or name\n')
end
expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);
index=find(strcmp(transName,TransNames));
if x(index)==0
  error('The gene selected is not a target of the transcription factor \n')
end
tf=expectations.b(:,index);
ind=find(strcmp(transName,TransNames(find(x))));
tfErrors=expectations.tfError(:,ind);
tfErrorsDiffs=[expectations.tfErrorDiffs(:,:,ind)];