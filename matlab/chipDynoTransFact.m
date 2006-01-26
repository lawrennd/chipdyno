function [TF,TFError]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                        name);
% CHIPDYNOTRANSFACT provides gene-specific TFAs with errorbars.

% CHIPDYNO

index=find(strcmp(name,transNames));
genesIn=find(X(:,index));
anno=annotations(find(X(:,index)));
nTargets=size(anno,1);
npts=size(data,2);
TF=zeros(nTargets,npts);
TFError=zeros(nTargets,npts);
for i=1:nTargets
  [TF(i,:),TFError(i,:)]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         name,genesIn(i));
end
