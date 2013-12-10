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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
// ROOT include
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
    TClonesArray * gen_particle_;	 // array with all interesting particles' 4vectors
    TClonesArray * gen_jet_; 		// array with all interesting jets produced in event
	 										// that is, that are associated with something
											
    std::vector<Int_t> gen_particle_ID_; // vector with all particles pdg ID
    std::vector<Int_t> gen_particle_status_; // " " " status
    std::vector<Int_t> gen_particle_mother_; // " " " mother
    std::vector<bool> gen_particle_issignal_; // " " " flag that indicates if 
                           					// particle is signal (true) or originating from radiation (false)
	 std::vector<Int_t> gen_particle_assoc_jet_; // >= 0: Indicates the jet (number in the vector)
	                                             // was originated by the particle
    					                     // -1: particle did not produce a jet
													// -2: particle produced a jet, but match could not find it
    
	 std::vector<Int_t> gen_jets_constituents_; // vector with all jets number of constituents
   
    
    typedef math::XYZTLorentzVector LorentzVector ;
    void setMomentum (TLorentzVector & myvector , const LorentzVector & mom) ;
    bool IsSignal (const reco::Candidate* particle);
    bool IsLast (const reco::Candidate* particle); 
    bool IsSonOfSignalTau(const reco::Candidate* particle);
    void PrintParticleInfo(const reco::Candidate* p, int counter);
    void PrintDaughters(const reco::Candidate* p);
    const reco::GenJet* GetAssociatedJet (const reco::Candidate* p, int& flag, edm::Handle<reco::GenJetCollection>* genJets);

	 // just the sign function
	 template <typename T> int sgn(T val) {
			return (T(0) < val) - (val < T(0));
			}
			
    int bgoodmatch;			
    int bnomatch;
    int bmultmatch;
    int taugoodmatch;
    int taunomatch;
    int taumultmatch;
};
#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//
