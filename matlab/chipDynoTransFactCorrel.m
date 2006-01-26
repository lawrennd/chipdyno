function [names, correls]=chipDynoTransFactCorrel(Sigma, TransNames, name)

% CHIPDYNOTRANSFACTCORREL correlations between transcription factors

% CHIPDYNO

auxSigma=sqrt(diag(Sigma))*sqrt(diag(Sigma))';
Corrs=Sigma./auxSigma;
index=find(strcmp(name, TransNames));
corres=Corrs(index,:);
[correls, indices]=sort(corres,'descend');
names=TransNames(indices);