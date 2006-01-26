% CHIPREGRESSION calculates TFAs by simple regression for comparison

% CHIPDYNO

[probeName, data] = chipTextRead('SpellmanMicro.txt');
[probeName2, annotation, dataChip] = chipChipTextRead( ...
    'Yeast_Connectivity.txt', 'Connectivity_Matrix.txt');
TransNames=textread('Trans_Names.txt','%q', ...
                    'headerlines',1,'whitespace','','delimiter','\t');
TransNames=TransNames(4:end-2);
index=zeros(size(dataChip,1),1);
for i=1:size(dataChip,1)
    index(i)=sum(strcmp(probeName2(i),probeName));
end
dataChip=dataChip(find(index),:);
index=zeros(size(data,1),1);
for i=1:size(data,1)
    index(i)=sum(strcmp(probeName(i),probeName2));
end
data=data(find(index),:);
X=zeros(size(dataChip,1),size(dataChip,2));
I=find(dataChip<1e-3);
X(I)=1;
npts = size(data,2);
nTrans = size(X,2);
Cov = (data'*data)/size(data,1);
B = zeros(nTrans, npts);
totalVar=sum(eig(Cov));
sigma2=(totalVar-sum(eig(X'*X)))/(npts-nTrans);
if sigma2 < 0
   sigma2=1e-6;
end
invM = inv(X'*X+sigma2*eye(nTrans));
f = invM*X'*data;