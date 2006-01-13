function f=chipDynoLikeStatNoise(params,data,precs,X,nEffectGenes,R,C);
%CHIPDYNOLIKESTATNOISE marginal likelihood for chipChip dynamical model

%CHIPDYNO
nGenes=size(data,1);
nTrans=size(X,2);
npts=size(data,2);
beta=params(1);
gamma=params(2);
mu=params(3:2+nTrans)';
V=params(3+nTrans:end-nTrans)';
preSigma=sparse(R,C,V,nEffectGenes,nTrans);
diagonal=params(end-nTrans+1:end);
Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
factor=cos(gamma)^2;
preLike=zeros(1, nGenes);
k=zeros(1,npts);
alpha=zeros(1,npts);
for i=1:nGenes
  coeff=(1-factor^2)*X(i,:)*Sigma*X(i,:)';
  number=X(i,:)*mu;
  alpha(1)=beta^2+precs(i,npts);
  k(1)=(data(i,npts)-(1-factor)*number)/factor;
  sigma=(beta^2+precs(i,npts-1))^-1+(alpha(1)^-1+coeff)/factor^2;
  preLike(i)=preLike(i)+ 0.5*log(sigma)+...
       0.5*((k(1)-data(i,npts-1))^2)/(sigma);
  for j=2:npts-1
    alpha(j)=(beta^2+precs(i,npts-j+1)+factor^2*(alpha(j-1)^-1+coeff)^-1);
    k(j)=alpha(j)^-1*(data(i,npts-j+1)*beta^2+factor^2*k(j-1)*(alpha(j-1)^-1+coeff)^-1);
    sigma=(beta^2+precs(i,npts-j))^-1+(alpha(j)^-1+coeff)/factor^2;
    preLike(i)=preLike(i)+ 0.5*log(sigma)+...
      0.5*((k(j)-data(i,npts-j))^2)/(sigma);
  end
  alpha(end)=(beta^2+precs(i,1)+factor^2*(alpha(end-1)^-1+coeff)^-1);
  k(end)=alpha(end)^-1*(data(i,1)*beta^2+factor^2*k(end-1)*(alpha(end-1)^-1+coeff)^-1);
  preLike(i)=preLike(i)+0.5*log(alpha(end)^-1+X(i,:)*Sigma*X(i,:)')+...
      0.5*((k(end)-number)^2)/(alpha(end)^-1+X(i,:)*Sigma*X(i,:)');
end
f=sum(preLike);