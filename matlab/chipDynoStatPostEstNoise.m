function expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);

% CHIPDYNOPOSTEST computes posterior expectations.

% CHIPDYNO


npts=size(data,2);
nTrans=size(x,1);
invC=chipStatMatrixInverterNoise(Sigma, gamma, beta,precs,x,npts);
Y=Sigma*x;
lambda=Y'*x;
factor=cos(gamma)^2;
coeff=x'*mu';
%YYT=Y*Y';
%invSigma=pdinv(Sigma);
%Z=invSigma*mu';
Mean=[];

for i=1:npts

  Mean=[Mean;sum((beta^-2*ones(1,npts)+...
      precs.^-1).^-1.*data.*(invC.Sigma(i,:)+lambda*invC.YYT(i,:)))*Y'+ ...
        (1+factor)^-1*(invC.Sigma(i,1)*mu+coeff*invC.YYT(i,1)*Y')+ ...
     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(i,2:end-1))*mu+coeff* ...
        sum(invC.YYT(i,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(i, ...
                                                  end)*mu+coeff* ...
     invC.YYT(i,end)*Y')];
end

       
expectations.b=Mean;
expectations.bTb=[];
for i=1:npts
  expectations.bTb=[expectations.bTb;invC.Sigma(i,i)*diag(Sigma)'+ ...
                    invC.YYT(i,i)*Y'.*Y'];
end


  
     
     