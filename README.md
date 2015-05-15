|----------------------------------------------|-----------------------------------------------------|----------------------------------------------------------|
| [home](http://www.cs.manchester.ac.uk/)  |   | [research](http://www.cs.manchester.ac.uk/ai/)  |   | [people](http://www.cs.manchester.ac.uk/ai/people/)  |   |

|-----|---------------------------------------------------------------------------------------------------------|
|     | ![University Logo](http://www.cs.manchester.ac.uk/ai/pictures/icons/est1824.gif "Manchester Est. 1824") |

n
Dynamical Model for Transcription
=================================

This page describes examples of how to use the Dynamical Model for Transcription detailed in [this paper](http://www.cs.man.ac.uk/neill-bin/publications/bibpage.cgi?keyName=Sanguinetti:chipdyno06&printAbstract=1). The code is available for [download here](http://www.cs.man.ac.uk/neill-bin/software/downloadForm.cgi?toolbox=chipdyno).
### Release Information

Current release is 0.11.

Release 0.11 is the second release associated with the paper, containing scripts for recreating the results given.

The toolboxes required to run the code are listed below.

The CHIPDYNO software can be downloaded [here](http://www.cs.man.ac.uk/neill-bin/software/downloadForm.cgi?toolbox=chipdyno).

Release Information
-------------------

Current release is 0.121.

As well as downloading the CHIPDYNO software you need to obtain the toolboxes specified below. These can be downloaded using the *same* password you get from registering for the CHIPDYNO software.

|--------------------------------------------------------------------------|-------------|
| **Toolbox**                                                              | **Version** |
| [NDLUTIL](http://www.cs.man.ac.uk/~neill/ndlutil/downloadFiles/vrs0p159) | 0.159       |
| [NETLAB](http://www.cs.man.ac.uk/~neill/netlab/downloadFiles/vrs3p3)     | 3.3         |

#### Version 0.12

Release version with some files missing (blank files for several demSpellman scripts).
Finally you will also need the [NETLAB toolbox](http://www.ncrg.aston.ac.uk/netlab/) in your path.

Results in the Paper
--------------------

The main scripts are `demTu` and `demSpellman` which run the code on the metabolic and cell cycle data sets. The function `chipdynoTransFact` computes the posterior, gene specific, transcription factor activities for a given transcription factor. The function `chipdynoTransFactNoise` does the same thing for expression levels with associated noise (such as computed by mmgMOS). The function `chipdynoGeneAct` sorts the transcription factor activities for a particular gene. The function `chipdynoTransFactCorrel` normalises the correlations and presents the most correlated TFs. Some scripts have been added to recreate the further validatory experiments detailed in the [Supplementary material](supplement.pdf). Specifically,`demSpellmanp1`, `demSpellmanp2` and `demSpellmanp5` run the code on the cell cycle data set with connectivity matrix obtained using different cut-offs at p=0.001, p=0.0002 and p=0.005 respectively. `demSpellmanp1p2` and `demSpellmanp5p2` do the same thing but only use the genes and transcription factors active at p=0.0002, so that the systems have the same dimensions. The script `chipDynoCompareX` can then be used to compare the matrices of effective regulations in the different cases.

n
Page updated on Fri Mar 14 17:11:24 2008
| [Disclaimer](http://www.manchester.ac.uk/aboutus/documents/disclaimer/ "Disclaimer") | [Privacy](http://www.manchester.ac.uk/aboutus/documents/privacy/ "Privacy") | [Copyright notice](http://www.manchester.ac.uk/aboutus/documents/copyright/ "Copyright Notice") | [Accessibility](http://www.manchester.ac.uk/aboutus/documents/accessibility/ "Accessibility") | [Freedom of information](http://www.manchester.ac.uk/aboutus/documents/foi/ "Freedom of information") | [Feedback](http://www.manchester.ac.uk/aboutus/contact/feedback/ "Feedback") |

Please contact <webmaster.cs@manchester.ac.uk> with comments and suggestions


