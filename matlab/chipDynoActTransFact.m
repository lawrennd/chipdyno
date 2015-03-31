function [list,newX, newXVals]=chipDynoActTransFact(data,X,Sigma,beta,gamma,mu, ...
                                         TransNames, annotation,sigLev);
% CHIPDYNOACTTRANSFACT identifies significantly varying TFs.

% CHIPDYNO

nTrans=size(TransNames,1);
list=[];
newX=zeros(size(X));
newXVals=zeros(size(X));
for i=1:nTrans
    [TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
                                         TransNames, annotation, ...
                                        TransNames(i));
    %vars=max(abs((TF-mu(i)*ones(size(TF)))'./TFError'));
    maxVars=chipDynoMaxDiff(TF,TFErrorDiff);
    sigVars=maxVars(find(maxVars>sigLev));
    list=[list, size(sigVars,2)];
    index=find(X(:,i));
    newX(index(find(maxVars>sigLev)),i)=1;
    newXVals(index(find(maxVars>sigLev)),i)=maxVars((find(maxVars>sigLev)));
end