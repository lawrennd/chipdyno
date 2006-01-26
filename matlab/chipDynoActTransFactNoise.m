function list=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation);

% CHIPDYNOACTTRANSFACT identifies significantly varying TFs.

% CHIPDYNO

nTrans=size(TransNames,1);
list=[];
for i=1:nTrans
    [TF,TFError]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation, ...
                                        TransNames(i));
    vars=sqrt(var(TF'))./max(TFError');
    sigVars=vars(find(vars>2));
    list=[list, size(sigVars,2)];
end