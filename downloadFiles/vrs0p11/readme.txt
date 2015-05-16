CHIPDYNO software
Version 0.11		Thursday 26 Jan 2006 at 20:41
Copyright (c) 2006 Guido Sanguinetti


This code contains the scripts for recreating the results in G.Sanguinetti, M.
Rattray and Neil D. Lawrence, "A probabilistic dynamical model for quantitative inference of the regulatory mechanism of transcription".

Version 0.1
-----------

Initial release.

MATLAB Files
------------

Matlab files associated with the toolbox are:

chipChipTextRead.m: reads TXT file for the Lee ChIP data files.
chipDynoActTransFactNoise.m: identifies significantly varying TFs.
chipDynoGeneActNoise.m: given a gene, lists activators in decreasing order
chipDynoGeneAct.m: given a gene, lists activators in decreasing order
chipDynoExpectationsFast.m: computes posterior expectations of TFA.
chipDynoExpectationsFastNoise.m: computes posterior expectations TFA.
chipDynoLikeStat.m: marginal likelihood for chipChip dynamical model
chipRegression.m: calculates TFAs by simple regression for comparison
chipDynoLikeStatGrad.m: gradient of CHIPDYNOLIKESTAT.
chipDynoLikeStatNoise.m: marginal likelihood for chipChip dynamical model
chipDynoLikeStatNoiseGrad.m: gradient of CHIPDYNOLIKESTATNOISE.
chipDynoLoadData.m: loads Spellman Data with Lee et al ChIP data.
chipDynoNewTransFactNoise.m: tfs active for us and not for Tu et al.
chipDynoPerGenes.m: genes periodic according to Tu et al
chipDynoStatPostEst.m: computes posterior expectations 
chipDynoStatPostEstNoise.m: computes posterior expectations.
chipDynoTransFact.m: provides gene-specific TFAs with errorbars.
chipDynoTransFactCorrel.m: correlations between transcription factors
chipDynoTransFactNoise.m: given a transcription factor, provides TFAs.
chipDynoTuLoadData.m: loads Tu Data with Lee et al ChIP data.
chipTextRead.m: reads TXT file for the Spellman data files.
demTu.m: demonstrates dynamical chipCHIP on Tu data.
chipReduceVariables.m: reduce  number of variables in chipDyno model
chipStatMatrixInverter.m: inverts block tridiagonal matrices for chipChip
chipStatMatrixInverterNoise.m: inverts block tridiagonal matrices for chipChip
chipTuTextRead.m: assigns common names to probe IDs
demSpellmanDynoStat.m: demonstrates dynamical chipCHIP on Spellman data.
