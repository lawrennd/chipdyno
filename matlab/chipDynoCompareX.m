function indicator=chipDynoCompareX(X1,probeNames1,TransNames1,X2,probeNames2,TransNames2);

% CHIPDYNOCOMPAREX compares connectivity for different pvals

% CHIPDYNO

[row,col]=find(X1);
indicator=[];
for i=1:size(row,1)
    indicator=[indicator,X2(find(strcmp(probeNames1(row(i)),probeNames2)), find(strcmp(TransNames1(col(i)),TransNames2)))];
end
