function list=chipDynoPerGenes(anno)

% CHIPDYNOPERGENES genes periodic according to Tu et al

% CHIPDYNO

TuPerGenes=textread('./data/MetabolData/PerGenes.txt', ...
                     '%q');
index=[];
for i=1:size(anno,1)
  index=[index,size(find(strcmp(anno(i),TuPerGenes)),1)];
end
list=anno(find(index));