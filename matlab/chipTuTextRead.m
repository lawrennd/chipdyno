function [ORF, data, vars]=chipTuTextRead()
%CHIPTUTEXTREAD assigns common names to probe IDs

%CHIPDYNO
[IDs,crap,orfs,geneCN]=textread('../../data/MetabolData/dictionary.txt','%q%q%q%q');
IDSnew=textread('../../data/MetabolData/probeIDTu.txt','%q');
data=load('../../data/MetabolData/YeastMetabolism_exprs.txt');
vars=load('../../data/MetabolData/YeastMetabolism_se.txt');
ORF=[];
for i=1:size(IDSnew,1)
    index=strmatch(IDSnew(i), IDs);
    ORF=[ORF;orfs(index)];
end