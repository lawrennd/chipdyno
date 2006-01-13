function [tf,tfErrors]=chipDynoExpectations(data,X,Sigma,beta,tau,mu, ...
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
postCov=chipMatrixInverterSlow(invSigma,tau,beta,x,npts);
% [invC,auxMat,lambda,Y,nAct]=chipMatrixInverter(invSigma, tau,beta,x,npts);
% YYT=Y*Y';
% XXT=x*x';
% XYT=x*Y';
% YXT=Y*x';
for i=1:npts
  adjData(i)=data(npts-i+1);
end
index=find(strcmp(transName,transNames));
if x(index)==0
  error('The gene selected is not a target of the transcription factor \n')
end
Z=invSigma*mu';
% W=auxMat*Z;
% Mean=[];
% for i=1:npts-1
% 
%   Mean=[Mean;beta^2*(sum((invC.AM(i,:)+invC.YXT(i,:)*nAct+ ...
%                           invC.YYT(i,:)*lambda).*adjData)*Y(index)+ ...
%                      sum((invC.Id(i,:)+invC.XXT(i,:)*nAct+ ...
%                           invC.XYT(i,:)*lambda).*adjData)*x(index))];
% end
% Mean=[Mean;sum(invC.AM(npts,:))*W(index)+sum(invC.YYT(npts,:))*(Y'*Z)*Y(index)+...
%     sum(invC.XYT(npts,:))*(Y'*Z)*x(index)+sum(invC.YXT(npts,:))*(x'*Z)*Y(index)+...
%     sum(invC.XXT(npts,:))*(x'*Z)*x(index)+...
%     beta^2*(sum((invC.AM(npts,:)+invC.YXT(npts,:)*nAct+ ...
%                           invC.YYT(npts,:)*lambda).*adjData)*Y(index)+ ...
%                      sum((invC.Id(npts,:)+invC.XXT(npts,:)*nAct+ ...
%                           invC.XYT(npts,:)*lambda).*adjData)*x(index))];  
% 
% tf=zeros(1,npts);
% for i=1:npts
%   tf(i)=Mean(npts-i+1);
% end
% 
% expectations.bbT=zeros(nTrans,nTrans,npts);
% for i=1:npts
%   expectations.bbT(:,:,i)=invC.AM(i,i)*auxMat+ ...
%       invC.YYT(i,i)*YYT+ invC.XXT(i,i)*XXT+ ...
%       invC.YXT(i,i)*YXT+invC.XYT(i,i)*XYT;
% end
Mean=[];
for i=1:npts-1
    Mean=[Mean;beta^2*adjData(i)*x];
end
Mean=[Mean;beta^2*adjData(end)*x+Z];
m=postCov*Mean;
tf=zeros(1,npts);
for i=1:npts
    tf(i)=m((npts-i)*nTrans+index);
end
errors=zeros(1,npts);
for i=1:npts
  errors(i)=postCov((npts-i)*nTrans+index,(npts-i)*nTrans+index);
end
tfErrors=sqrt(errors);
     
     