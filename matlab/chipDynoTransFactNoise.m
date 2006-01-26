function [TF,TFError]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         transNames, annotations, ...
                                        name);
% CHIPDYNOTRANSFACTNOISE given a transcription factor, provides TFAs.

% CHIPDYNO

index=find(strcmp(name,transNames));
genesIn=find(X(:,index));
anno=annotations(find(X(:,index)));
nTargets=size(anno,1);
npts=size(data,2);
TF=zeros(nTargets,npts);
TFError=zeros(nTargets,npts);
for i=1:nTargets
  [TF(i,:),TFError(i,:)]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         transNames, annotations, ...
                                         name,genesIn(i));
end
