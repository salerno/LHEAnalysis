#ifndef GeneratorAna_H
#define GeneratorAna_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//MC
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ROOT include
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//
// class declaration
//

using namespace std;

class GeneratorAna : public edm::EDAnalyzer {
public:
  explicit GeneratorAna(const edm::ParameterSet&);
  ~GeneratorAna();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
 // ----------member data ---------------------------
    TTree * mytree_;
    TLorentzVector myvector ;
    TClonesArray * gen_particle_;
    std::vector<Int_t> gen_particle_ID_;
    std::vector<Int_t> gen_particle_status_;
    std::vector<Int_t> gen_particle_mother_;
    
    typedef math::XYZTLorentzVector LorentzVector ;
    void setMomentum (TLorentzVector & myvector , const LorentzVector & mom) ;
    
};
#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//
