function [tf,tfErrors]=chipDynoStatPostEstSlow(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         transName,geneName);
%CHIPDYNOEXPECTATIONS computes posterior expectations of gene-specific
%TFA. 
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
invSigma=pdinv(Sigma);
postCov=chipStatMatrixInverterSlow(Sigma,gamma,beta,x,npts);
% [invC,auxMat,lambda,Y,nAct]=chipMatrixInverter(invSigma, tau,beta,x,npts);
% YYT=Y*Y';
% XXT=x*x';
% XYT=x*Y';
% YXT=Y*x';

index=find(strcmp(transName,transNames));
if x(index)==0
  error('The gene selected is not a target of the transcription factor \n')
end
Z=invSigma*mu';
factor=cos(gamma);
preMean=[beta^2*data(1)*x+(1+factor)^-1*Z];
for i=1:npts-1
    preMean=[preMean;beta^2*data(i)*x+(1-factor)*(1+factor)^-1*Z];
end
Mean=[Mean;beta^2*data(end)*x+(1+factor)*Z];
m=postCov*Mean;
tf=zeros(1,npts);
for i=1:npts
    tf(i)=m((i-1)*nTrans+index);
end
errors=zeros(1,npts);
for i=1:npts
  errors(i)=postCov((i-1)*nTrans+index,(i-1)*nTrans+index);
end
tfErrors=sqrt(errors);
     
     