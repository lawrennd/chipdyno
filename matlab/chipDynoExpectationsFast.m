function [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         transName,geneName);
%CHIPDYNOEXPECTATIONS computes posterior expectations of TFA.

%CHIPDYNO
npts=size(data,2);
nTrans=size(X,2);

c=class(geneName);
if c(1)=='c'
    x=X(find(strcmp(geneName,annotations)),:)';
    data=data(find(strcmp(geneName,annotations)),:);
elseif c(1)=='d'
    x=X(geneName,:)';
    data=data(geneName,:);
else
    error('Genes can be identified either by number or name\n')
end
expectations=chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);
index=find(strcmp(transName,transNames));
if x(index)==0
  error('The gene selected is not a target of the transcription factor \n')
end
tf=expectations.b(:,index);
ind=find(strcmp(transName,transNames(find(x))));
tfErrors=expectations.tfError(:,ind);
tfErrorsDiffs=[expectations.tfErrorDiffs(:,:,ind)];