import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/1CA79834-4BF9-E211-B6CF-1CC1DE040FE8.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/265D8272-54F9-E211-8B81-00266CFE79A4.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/309B6433-6EF9-E211-8496-1CC1DE1CE026.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/38B03E96-64F9-E211-BA12-1CC1DE046F18.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/500325EE-60F9-E211-AC36-0017A4770C1C.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/58BB0B39-60F9-E211-9EB2-1CC1DE048FD0.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/9652011F-69F9-E211-A385-1CC1DE055158.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/96DE1A13-66F9-E211-9C84-00266CFFCD00.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/9AA81DA0-6CF9-E211-8CFE-0017A4770430.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/AA0E9D93-6AF9-E211-94B6-1CC1DE047F98.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/BA51A643-5EF9-E211-A4A6-1CC1DE1CED1C.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/D4C00D7C-67F9-E211-8957-1CC1DE046F18.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/E2AF1AE8-64F9-E211-BD0D-0017A4770C1C.root',
                                      #'/store/mc/Summer12_DR53X/GluGluToHHTo2B2Tau_mH-125_8TeV-madgraph-pythia6-tauola/AODSIM/PU_S10_START53_V19-v1/10000/F824ED06-70F9-E211-B526-00266CFF0044.root'
                                      'file:/home/llr/cms/salerno/CMSSW/DoubleHiggs/1CA79834-4BF9-E211-B6CF-1CC1DE040FE8.root'
                                                                             )


)

process.demo = cms.EDAnalyzer('GeneratorAna')


process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root") )

process.p = cms.Path(process.demo)
