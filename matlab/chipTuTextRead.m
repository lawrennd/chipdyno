function [ORF, data, vars]=chipTuTextRead()

% CHIPTUTEXTREAD assigns common names to probe IDs

% CHIPDYNO

[IDs,crap,ORF,geneCN]=textread('dictionary.txt','%q%q%q%q');
IDSnew=textread('probeIDTu.txt','%q');
data=load('YeastMetabolism_exprs.txt');
vars=load('YeastMetabolism_se.txt');
