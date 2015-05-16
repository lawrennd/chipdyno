CHIPDYNO software
Version 0.12		Thursday 04 May 2006 at 15:20
Copyright (c) 2006 Guido Sanguinetti


This code contains the scripts for recreating the results in G.Sanguinetti, M.
Rattray and Neil D. Lawrence, "A probabilistic dynamical model for quantitative inference of the regulatory mechanism of transcription".

Version 0.11
------------

Release of code for final version of Bioinformatics paper.

Version 0.1
-----------

Initial release.

MATLAB Files
------------

Matlab files associated with the toolbox are:

%CHIPDYNOGENEACT given a gene, lists activators in decreasing order
%CHIPDYNOEXPECTATIONSFASTNOISE computes posterior expectations TFA.
%CHIPDYNOACTTRANSFACT identifies significantly varying TFs.
%CHIPDYNOLIKESTAT marginal likelihood for chipChip dynamical model
chipDynoCompareTrans.m: compares list of active tfs at different p
%CHIPDYNOLIKESTATGRAD gradient of CHIPDYNOLIKESTAT.
%CHIPDYNOLIKESTATNOISE marginal likelihood for chipChip dynamical model
%CHIPDYNOLIKESTATNOISEGRAD gradient of CHIPDYNOLIKESTATNOISE.
%CHIPDYNOLOADDATA loads Spellman Data with Lee et al ChIP data.
chipDynoCompareX.m: compares connectivity for different pvals
%CHIPDYNONEWTRANSFACTNOISE tfs active for us and not for Tu et al.
%CHIPDYNOPERGENES genes periodic according to Tu et al
chipDynoMaxDiff.m: computes most significant changes in TFAs
%CHIPDYNOTRANSFACT provides gene-specific TFAs with errorbars.
%CHIPDYNOTRANSFACTCORREL correlations between transcription factors
%CHIPDYNOTRANSFACTNOISE given a transcription factor, provides TFAs.
%CHIPDYNOTULOADDATA loads Tu Data with Lee et al ChIP data.
chipDynoSignRels.m: finds the regulatory strength conserved rels
%CHIPTEXTREAD reads TXT file for the Spellman data files.
%DEMTU demonstrates dynamical chipCHIP on Tu data.
%CHIPREDUCEVARIABLES reduce  number of variables in chipDyno model
%CHIPTUTEXTREAD assigns common names to probe IDs
