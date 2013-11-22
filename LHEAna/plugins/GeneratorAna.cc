#include "GeneratorAna.h"

// CMSSW include
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT include
#include "TLorentzVector.h"


using namespace std;
using namespace edm;
using namespace reco;

// ------------ method called for each event  ------------
// ===============================================================================
void GeneratorAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ===============================================================================
{
    //clear variables
    gen_particle_->Clear();
    gen_particle_ID_.clear();
    gen_particle_status_.clear();
	gen_particle_mother_.clear();
	
    edm::Handle<View<Candidate> > genCandidatesCollection;
	iEvent.getByLabel("genParticles", genCandidatesCollection);
	
	TClonesArray &gen_particle = *gen_particle_;
    
	int counter = 0;
	
	// ----------------------------
	//      Loop on particles
	// ----------------------------
	for( View<Candidate>::const_iterator p = genCandidatesCollection->begin();
        p != genCandidatesCollection->end();
        ++p ) {
		
		if ( p->pdgId() == 25 ||
             abs(p->pdgId()) == 15 ||
             abs(p->pdgId()) == 5 ||
            (abs(p->pdgId()) == 11  && abs(p->mother(0)->pdgId()) == 15) ||
            (abs(p->pdgId()) == 13  && abs(p->mother(0)->pdgId()) == 15) ) { // start if Higgs/b/tau/e/mu
			
            setMomentum(myvector,p->p4());
			new (gen_particle[counter]) TLorentzVector(myvector);
			counter++;
            gen_particle_ID_.push_back(p->pdgId());
			gen_particle_status_.push_back(p->status());
			gen_particle_mother_.push_back(p->mother(0)->pdgId());
		
            //cout << "PdgId=" << p->pdgId() << " status=" << p->status() << " Mass=" << myvector.M()  << " X,Y,Z=" << p->px() << "," << p->py() << "," << p->pz() << " Mother=" << p->mother(0)->pdgId() << endl;
	
        } // end if Higgs
    }

  mytree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
GeneratorAna::beginJob()
{

  edm::Service<TFileService> fs;
  mytree_  = fs->make <TTree>("simpleRootTree","simpleRootTree");

  gen_particle_ = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("gen_particle", "TClonesArray", &gen_particle_, 256000,0);
  mytree_->Branch ("gen_particle_ID",&gen_particle_ID_);
  mytree_->Branch ("gen_particle_status",&gen_particle_status_);
  mytree_->Branch ("gen_particle_mother",&gen_particle_mother_);
    
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GeneratorAna::endJob() 
{
    delete gen_particle_;

}

// ------------ method called when starting to processes a run  ------------
void 
GeneratorAna::beginRun(edm::Run const&, edm::EventSetup const&)
{

 

}

// ------------ method called when ending the processing of a run  ------------
void 
GeneratorAna::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GeneratorAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GeneratorAna::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GeneratorAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//
// constructors and destructor
//
GeneratorAna::GeneratorAna(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
}


GeneratorAna::~GeneratorAna()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
void GeneratorAna::setMomentum(TLorentzVector & myvector, const LorentzVector & mom)
{
	myvector.SetPx (mom.Px());
	myvector.SetPy (mom.Py());
	myvector.SetPz (mom.Pz());
	myvector.SetE (mom.E());
}


//define this as a plug-in
DEFINE_FWK_MODULE(GeneratorAna);
