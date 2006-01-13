function [data,vars,X,annotation,TransNames]=chipDynoTuLoadData();
%CHIPDYNOTULOADDATA loads Tu Data with Lee et al ChIP data.

%CHIPDYNO
[probeName, data, vars] = chipTuTextRead;
dataChip=load('../../data/Connectivity2.txt');
[probeName2, annotation] = textread('../../data/annotations2.txt','%q%q');
TransNames=textread('../../data/Trans_Names2.txt','%q');
TransNames=TransNames(2:end);

index=zeros(size(dataChip,1),1);
for i=1:size(dataChip,1)
    index(i)=sum(strcmp(probeName2(i),probeName));
end
dataChip=dataChip(find(index),:);
annotation=annotation(find(index));
index=zeros(size(data,1),1);
for i=1:size(data,1)
    index(i)=sum(strcmp(probeName(i),probeName2));
end
data=data(find(index),:);
vars=vars(find(index),:);
probeName=probeName(find(index));
X=zeros(size(dataChip,1),size(dataChip,2));
I=find(dataChip<1e-3);
X(I)=1;

fakeX=sum(X,2);
X=X(find(fakeX),:);
annotation=annotation(find(fakeX));
effectX=sum(X,1);
TransNames=TransNames(find(effectX));
X=X(:,find(effectX));
data=data(find(fakeX),:);
vars=vars(find(fakeX),:);