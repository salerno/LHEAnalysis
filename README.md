LHEAnalysis
===========

Simple code to make analysis on the LHE files in CMSSW:

```
scram p -n CMSSW_5_3_11_p6_LHE CMSSW CMSSW_5_3_11_patch6
cd CMSSW_5_3_11_p6_LHE/src
cmsenv
git clone https://github.com/salerno/LHEAnalysis LHEAnalysis
scram b -j 8
```

To run the code to perform the generator level analysis: 
```
cd  LHEAnalysis/LHEAna/test 
cmsRun generatorana.py
```
and it will produce a flat ROOT tree
