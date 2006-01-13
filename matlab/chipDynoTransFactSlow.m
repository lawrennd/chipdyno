function [TF,TFError]=chipDynoTransFactSlow(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                        name);
%CHIPDYNOTRANSFACT given a transcription factor, provides
%gene-specific TFAs with errorbars.
index=find(strcmp(name,transNames));
genesIn=find(X(:,index));
anno=annotations(find(X(:,index)));
nTargets=size(anno,1);
npts=size(data,2);
TF=zeros(nTargets,npts);
TFError=zeros(nTargets,npts);
for i=1:nTargets
  [TF(i,:),TFError(i,:)]=chipDynoExpectationsSlow(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         name,genesIn(i));
end
