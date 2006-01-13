%DEMSPELLMAN demonstrates dynamical chipCHIP on Spellman data.

%CHIPDYNO
clear all
[data,X,annotation,TransNames]=chipDynoLoadData();
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
options=foptions;
options(1)=1;
options(14)=10000;
g=cov(data');
muIn=zeros(nTrans,1);

[R,C,V,nEffectGenes]=chipReduceVariables(X);
diagonal=0.5*ones(1,nTrans);
beta=3;
gamma=pi/4;
params=[beta,gamma,muIn',0.1*V', diagonal];
params=scg('chipDynoLikeStat',params,options,'chipDynoLikeStatGrad',data, ...
           X, nEffectGenes,R, C);
V=params(nTrans+3:end-nTrans)';
preSigma=sparse(R,C,V,nEffectGenes,nTrans);
diagonal=params(end-nTrans+1:end);
Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
beta=params(1);
gamma=params(2);
mu=params(3:2+nTrans);
save ResultsSpellmanStatSCG2 params