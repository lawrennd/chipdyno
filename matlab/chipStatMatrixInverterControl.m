function [C,f]=chipStatMatrixInverterControl(Sigma, gamma, beta, x, npts);
%CHIPMATRIXINVERTERCONTROL inverts chip matrix in the inefficient (but
%safe)way.
invSigma=pdinv(Sigma);
factor=cos(gamma);
nTrans=size(x,1);
row1=[beta^2*x*x'+(1-factor^2)^-1*invSigma,-factor*(1-factor^2)^-1*invSigma,zeros(nTrans,nTrans*(npts-2))];
C=row1;
for i=2:npts-1
    row=[zeros(nTrans,(i-2)*nTrans), -factor*(1-factor^2)^-1*invSigma,...
        beta^2*x*x'+(1+factor^2)*(1-factor^2)^-1*invSigma,-factor*(1- ...
                                                      factor^2)^-1* ...
         invSigma,zeros(nTrans,nTrans*(npts-1-i))];
    
    
    C=[C;row];
end
row=[zeros(nTrans,(npts-2)*nTrans),-factor*(1-factor^2)^-1*invSigma,...
    beta^2*x*x'+(1-factor^2)^-1*invSigma];
C=[C;row];
f=inv(C);
