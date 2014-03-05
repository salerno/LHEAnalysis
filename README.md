LHEAnalysis
===========

Simple code to make analysis on the LHE files or read the generator level information in CMSSW files:

```
scram p -n CMSSW_5_3_11_p6_LHE CMSSW CMSSW_5_3_11_patch6
cd CMSSW_5_3_11_p6_LHE/src
cmsenv
git clone https://github.com/salerno/LHEAnalysis LHEAnalysis
scram b -j 8
```

* LHEAna: class to read a LHE file and produce ROOT histograms (in the case of HZZ4l event). 

* GeneratorAna: class to read generator information from a CMSSW file and produce a flat ROOT tree, to run the code
  ```
  cd  LHEAnalysis/LHEAna/test 
  cmsRun generatorana.py
  ```
