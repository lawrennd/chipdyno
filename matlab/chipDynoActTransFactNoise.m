function [list,newX, newXVals]=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation);
% CHIPDYNOACTTRANSFACTNOISE identifies significantly varying TFs.

% CHIPDYNO
nTrans=size(TransNames,1);
list=[];
newX=zeros(size(X));
newXVals=zeros(size(X));
for i=1:nTrans
    [TF,TFError,TFErrorDiff]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation, ...
                                        TransNames(i));
    maxVars=chipDynoMaxDiff(TF,TFErrorDiff);
    pvals=cumGaussian(-maxVars));
    sigVars=pvals(find(pvals<0.02));
    list=[list, size(sigVars,2)];
    index=find(X(:,i));
    newX(index(find(pvals<0.02)),i)=1;
     newXVals(index(find(pvals<0.02)),i)=pvals((find(pvals<0.02)));
end