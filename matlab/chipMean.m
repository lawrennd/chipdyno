function chipMeanACE2(data,X,Sigma,beta,gamma,mu,transNames, ...
                      annotations)
%CHIPMEANACE2 computes the mean over the significantly regulated
%genes of ACEs TFA

[tf1,tfError1]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         'ACE2','CTS1');
[tf2,tfError2]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         'ACE2','SCW11');
[tf3,tfError3]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         'ACE2','YER124C');
[tf4,tfError4]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         'ACE2','YHR143W');
TF=[tf1;tf2;tf3;tf4];
TFerror=[tfError1;tfError2;tfError3;tfError4];
meanTF=mean(TF.*TFerror,1)./sum(TFerror,1);
meanTFerror=mean(TFerror,1);
chipPlotter(meanTF,meanTFerror);