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
    gen_jet_->Clear();
    gen_particle_ID_.clear();
    gen_particle_status_.clear();
    gen_particle_mother_.clear();
    gen_jets_constituents_.clear();
    gen_particle_issignal_.clear();
    gen_particle_assoc_jet_.clear();
		
    edm::Handle<View<Candidate> > genCandidatesCollection;
    iEvent.getByLabel("genParticles", genCandidatesCollection);
	
	
    edm::Handle<GenJetCollection> genJets;
    iEvent.getByLabel("ak5GenJets", genJets);
    //iEvent.getByLabel("kt4GenJets", genJets);
    
	 		
    TClonesArray &gen_particle = *gen_particle_;
    TClonesArray &gen_jet = *gen_jet_;
    
    int counterp = 0;
    int counterj = 0;

    int nconstituents;
    const GenJet* bufJet;
    int jmatchflag;
    bool flSignal;
    int tausons;
/*
    // stores the pointers to Candidate for those particle that will produce jet
	 // and the match jet will loop on them 
    std::vector<const reco::Candidate*> JetMother;
*/
	 
// ----------------------------
//      Loop on particles
// ----------------------------

//	cout << "*************\n";
//	cout << "* PARTICLES *\n";
//	cout << "*************\n";
	

	for( View<Candidate>::const_iterator p = genCandidatesCollection->begin();
  	   p != genCandidatesCollection->end();
 	   ++p )
	{

	/* PARTICLES (antiparticles have -)
	**
	** 11 e- 
	** 12 nu_e
	** 13 mu-
	** 14 nu_mu
	** 15 tau-
	** 16 nu_tau
	**
	** 25 H
	** 5 b
	**
	*/
	 	
		// save Higgs, b and taus (only if last)
		if (
				IsLast(&(*p)) == true &&
				(
					(p->pdgId() == 25) ||
					(abs(p->pdgId()) == 15) ||
					(abs(p->pdgId()) == 5)  
				)
			 )

		{
			// check if it is signal; NB: Higgs is always signal (set by hand)
			flSignal = IsSignal(&(*p));
			if (p->pdgId() == 25)
			{flSignal = true;}
		
			// import particle 4vector and store in the TClonesArray
			setMomentum(myvector,p->p4());			
			new (gen_particle[counterp]) TLorentzVector(myvector);
			counterp++;
			
			// fill std::vectors with interesting data (id, mother, status, signal flag)
			gen_particle_ID_.push_back(p->pdgId());
			gen_particle_status_.push_back(p->status());
			gen_particle_mother_.push_back(p->mother(0)->pdgId());
			gen_particle_issignal_.push_back(flSignal);
			
			PrintParticleInfo(&(*p), counterp);	
			PrintDaughters(&(*p));
	
			// Now find associated jets
			// Higgs won't produce jet
			// b will always produce jet
			// for tau depends from the decay
			
			// is a b or a tau, not a higgs
			if ( p->pdgId() != 25)
			{
				if (abs(p->pdgId()) == 5)
				{
					bufJet = GetAssociatedJet ( &(*p), jmatchflag, &genJets );
					
						if (jmatchflag != -2)
						{
							// import jet 4vector and store in the TCLonesArray
							setMomentum (myvector, bufJet->p4());
							new (gen_jet[counterj]) TLorentzVector(myvector);
							counterj++;
							// store num of constituents for each jet and relative particle
							nconstituents = (bufJet->getGenConstituents()).size();
							gen_jets_constituents_.push_back(nconstituents);
							gen_particle_assoc_jet_.push_back(counterj-1);
							
							cout << "MATCHED B jet - Pt: " << bufJet->pt() << " - Comp: " << nconstituents << endl;
							cout << "Jet is the " << counterj-1 << " in the list" << endl;
						}
					
						else
						{
							gen_particle_assoc_jet_.push_back(-2);
							cout << "NO B JET MATCHED!" << endl;
						}
					
						// update number of jet match
						switch (jmatchflag)
						{
							case(0):
								bgoodmatch++;
								break;
							case(-1):
								bmultmatch++;
								break;
							case(-2):
								bnomatch++;
								break;
						}
					
				}
			
				// analyze tau decay
				if (abs(p->pdgId()) == 15)
				{
					tausons = 0;
					for (unsigned int son = 0; son < p->numberOfDaughters(); son++)
					{
						// if tau has sons different from leptonic products (including radiative decays)
						if (
								abs(p->daughter(son)->pdgId()) != 11 &&
								abs(p->daughter(son)->pdgId()) != 12 &&
								abs(p->daughter(son)->pdgId()) != 13 &&
								abs(p->daughter(son)->pdgId()) != 14 && 
								abs(p->daughter(son)->pdgId()) != 16 &&
								abs(p->daughter(son)->pdgId()) != 22 						
							)
						{tausons++;}
					}
				
					// hadronic decay
					if (tausons > 0)
					{
						bufJet = GetAssociatedJet ( &(*p), jmatchflag, &genJets );
									
						if (jmatchflag != -2)
						{
							// import jet 4vector and store in the TCLonesArray
							setMomentum (myvector, bufJet->p4());
							new (gen_jet[counterj]) TLorentzVector(myvector);
							counterj++;
							// store num of constituents for each jet and relative particle
							nconstituents = (bufJet->getGenConstituents()).size();
							gen_jets_constituents_.push_back(nconstituents);
							gen_particle_assoc_jet_.push_back(counterj-1);
							
							cout << "MATCHED TAU jet - Pt: " << bufJet->pt() << " - Comp: " << nconstituents << endl;
							cout << "Jet is the " << counterj-1 << " in the list" << endl;
						}
					
						else
						{
							gen_particle_assoc_jet_.push_back(-2);
							cout << "NO TAU JET MATCHED!" << endl;
						}
					
						// update number of jet match
						switch (jmatchflag)
						{
							case(0):
								taugoodmatch++;
								break;
							case(-1):
								taumultmatch++;
								break;
							case(-2):
								taunomatch++;
								break;
						}
					}
				
					// leptonic decay --> no jet
					else
					{gen_particle_assoc_jet_.push_back(-1);}
				}
			}
			
			// particle did not satisfy previous requirements --> did not produce a jet
			else
			{gen_particle_assoc_jet_.push_back(-1);}
			
		/*	
			// this is to better analyze taus
			int y;
			
			if (abs(p->pdgId()) == 15 )
			{
				int nnn = 0;
				for (unsigned int i = 0; i < p->numberOfDaughters(); i++)
				{
					if (
							abs(p->daughter(i)->pdgId()) != 11 &&
							abs(p->daughter(i)->pdgId()) != 12 &&
							abs(p->daughter(i)->pdgId()) != 13 &&
							abs(p->daughter(i)->pdgId()) != 14 && 
							abs(p->daughter(i)->pdgId()) != 16 &&
							abs(p->daughter(i)->pdgId()) != 20213 &&
							abs(p->daughter(i)->pdgId()) != 213 &&
							abs(p->daughter(i)->pdgId()) != 24 &&
							abs(p->daughter(i)->pdgId()) != 211 &&
							abs(p->daughter(i)->pdgId()) != 321 &&
							abs(p->daughter(i)->pdgId()) != 323
					 )
						{
							//nnn++;
							y = i;
							nnn = 1000;
						}
				}
					
				if (nnn == 1000)
				{
					//PrintDaughters(&(*p));
					
					for (unsigned int k = 0; k < p->daughter(y)->numberOfDaughters(); k++)
					{
						cout << "    " << k << " Id: " << p->daughter(y)->daughter(k)->pdgId() << " status: "
								 << p->daughter(y)->daughter(k)->status() << endl;
					}
					
					//PrintDaughters(&(*p));
					char help;
					cin >>help;
				}
				
			}
		*/			
			
		}
		
		
		
		// now store all leptons and neutrinos
		// flag IsSignal in this case is "true" only if they come
		// from a tau decay (not b or other)
		if (
				IsLast(&(*p)) == true &&
				(
					(abs(p->pdgId()) == 11) ||
					(abs(p->pdgId()) == 12) ||
					(abs(p->pdgId()) == 13) ||
					(abs(p->pdgId()) == 14) ||
					(abs(p->pdgId()) == 16) 
				)
			 )
		{
			// are nu/leptons from signal tau?
			flSignal = IsSonOfSignalTau (&(*p));
			
			// import particle 4vector and store in the TClonesArray
			setMomentum(myvector,p->p4());			
			new (gen_particle[counterp]) TLorentzVector(myvector);
			counterp++;
			
			// fill std::vectors with interesting data (id, mother, status, signal flag)
			gen_particle_ID_.push_back(p->pdgId());
			gen_particle_status_.push_back(p->status());
			gen_particle_mother_.push_back(p->mother(0)->pdgId());
			gen_particle_issignal_.push_back(flSignal);
			gen_particle_assoc_jet_.push_back(-1); // leptons do not produce jet
			//PrintParticleInfo(&(*p), counterp);
	
		}
		
	}


// -------------------------------------
//      Loop on jets -- just to print
// -------------------------------------

//	cout << "*************\n";
//	cout << "***  JETS ***\n";
//	cout << "*************\n";

	// print only stored jets;
	for (int jj = 0; jj < gen_jet_ -> GetEntries(); jj++)
	{
		cout << jj << " - pt:" << ((TLorentzVector*)gen_jet_->At(jj) )->Perp() << " - nconst: " << gen_jets_constituents_.at(jj) << endl;
	}

	cout << "\n\n all jets\n";

	// print all jets
	for( GenJetCollection::const_iterator jets = genJets->begin();
        jets != genJets->end();
        ++jets )
	{

   		// print on screen pt, eta, phi and # constituents of jet	
		cout << "  " << counterj << " pt=" << jets->pt() << " eta=" << jets->eta() <<" phi=" << jets->phi()
		     << " costituents=" << (jets->getGenConstituents()).size()  << endl;
   		
		// loop on jet constituents -- print on screen ID, pt and eta of any jet constituent
		for (unsigned int i = 0 ; i < (jets->getGenConstituents()).size() ; ++i)
		{cout << "      genConstID="  << (jets->getGenConstituent(i))->pdgId() << " pt=" <<  (jets->getGenConstituent(i))->pt() << " eta=" << (jets->getGenConstituent(i))->eta() <<  " phi=" << (jets->getGenConstituent(i))->phi() << endl;}
	
	
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
  gen_jet_ = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("gen_particle", "TClonesArray", &gen_particle_, 256000,0);
  mytree_->Branch ("gen_jet", "TClonesArray", &gen_jet_, 256000,0);
  mytree_->Branch ("gen_particle_ID",&gen_particle_ID_);
  mytree_->Branch ("gen_particle_status",&gen_particle_status_);
  mytree_->Branch ("gen_particle_mother",&gen_particle_mother_);
  mytree_->Branch ("gen_particle_issignal", &gen_particle_issignal_);
  mytree_->Branch ("gen_particle_assoc_jet", &gen_particle_assoc_jet_);
  mytree_->Branch ("gen_jets_constituents",&gen_jets_constituents_); 

  bgoodmatch = 0;
  bnomatch = 0;
  bmultmatch = 0;
  taugoodmatch = 0;
  taunomatch = 0;
  taumultmatch = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GeneratorAna::endJob() 
{
    delete gen_particle_;
	 delete gen_jet_;
	 
	 cout << "\n\nFINAL MATCHING STATUS:\n\n";
	 cout << "B JETS:\n" << "  good: " << bgoodmatch
	 									<< "  failed: " << bnomatch
										<< "  multiple: " << bmultmatch
										<< endl; 
	 cout << "TAU JETS:\n"		<< "  good: " << taugoodmatch
	 									<< "  failed: " << taunomatch
										<< "  multiple: " << taumultmatch
										<< endl; 	 

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


// decides if particle is signal or noise loooking at his "ancestors"
// signal if is son of an Higgs
// in other cases background
bool GeneratorAna::IsSignal (const reco::Candidate* particle)
{
	int pId = particle->pdgId();
	const Candidate* pp = particle;
	bool fIsSignal = false;
	
	while (pp-> numberOfMothers() != 0)
	{
		// not son of the same type of particle
		if (pp->mother(0)->pdgId() != pId)
		{
			// son of an higgs --> signal
			if (pp->mother(0)->pdgId() == 25)
			{
				fIsSignal = true;
				break;
			}
			else
			{break;}
		}
		
		// is the same particle (radiation) --> see the ancestor
		else
		{pp = pp->mother(0);}
	}
	
	return fIsSignal;
}



// decides the if particle is the last of the sequence
// just before hadronization
// by looking at her sons; if it has no sons of the same type
// (es. no b -> b) it means it has hadronized
bool GeneratorAna::IsLast (const reco::Candidate* particle)
{
	int pId = particle->pdgId();
	const Candidate* pp = particle;
	const Candidate* pdau;
	bool fIsLast = true;
	
	// loop on the daughters
	for (unsigned int dau = 0; dau < pp-> numberOfDaughters(); dau++)
	{
		pdau = pp->daughter(dau);
		if (pdau->pdgId() == pId)
		{
			fIsLast = false;
			break;
		}
	}
	
	return fIsLast;
}






// decides if the particle (used for lept and nus) is the son of
// a signal or of a background tau
bool GeneratorAna::IsSonOfSignalTau (const reco::Candidate* particle)
{
	int pId = particle->pdgId();
	const Candidate* pp = particle;
	bool fIsFromSignal;
	
	// get the ancestor (if there is radiation)
	while (pp->mother(0)->pdgId() == pId)
	{pp = pp->mother(0);}
	
	// the mother is a tau?
	if ( abs( pp->mother(0)->pdgId() ) != 15)
	{fIsFromSignal = false;}
	
	else
	{fIsFromSignal = IsSignal (pp->mother(0));}
	
	return fIsFromSignal;
}






// prints info about current particle
void GeneratorAna::PrintParticleInfo (const reco::Candidate* p, int counter)
{
	cout << counter-1
			 << " PdgId=" << p->pdgId() << " status=" << p->status() << " Mass=" << myvector.M()
			 << " p(X,Y,Z)=(" << p->px() << "," << p->py() << "," << p->pz() << ") pt=" <<  p->pt()
			 << " eta=" << p->eta() <<  " phi=" << p->phi() << " Mother=" << p->mother(0)->pdgId()
			 << endl;
}





// print all sons ID + status
void GeneratorAna::PrintDaughters(const reco::Candidate* p)
{
	cout << "  N daughters: " << p->numberOfDaughters() << endl;
	for (unsigned int d = 0; d < p->numberOfDaughters(); d++)
	{
		cout <<"  " << d << " ID: " << p->daughter(d)->pdgId() << "  , Status: "
				 << p->daughter(d)->status() << "  , Pt: " << p->daughter(d)->pt()
				 <<"  , eta: " <<  p->daughter(d)->eta() << "  , phi: " << p->daughter(d)->phi()
				 << endl;
	}
}



// find the jet associated to the particle p; returns a pointer to the jet and
// stores the match result in the flag
// 0: good match
// -1: mult match, found the best
// -2: no match found
const reco::GenJet* GeneratorAna::GetAssociatedJet (
							const reco::Candidate* p,
							int& flag,
							edm::Handle<GenJetCollection>* genJets)
{
	double Rtolerance = 0.35;
	
	std::vector<const GenJet*> MatchingJets;
	
	bool CompFlag = true; // flag on the jets constituents (is there a nu?)
								// used only in case of taus
	
	// TLorentzVector for deltaR
	TLorentzVector vp;
	TLorentzVector vjet;
	
	setMomentum( vp ,p->p4() );
	
	// loop on all jets to find the best match
	// if p is a tau, also require the nu in the jet (algorithm ak5GenJets)
	for(
			GenJetCollection::const_iterator jets = (*genJets)->begin();
			jets != (*genJets)->end();
			++jets
		)
	{
		// if p is a tau, compute CompFlag
		if (abs(p->pdgId()) == 15)
		{
			CompFlag = false;
			for (unsigned int c = 0; c < (jets->getGenConstituents()).size(); c++)
			{
				if ( (jets->getGenConstituent(c))->pdgId() == 16*sgn(p->pdgId()) )
				{
					CompFlag = true;
					break;
				}
			}
		}
	
		setMomentum( vjet, jets->p4() );
		if ( vp.DeltaR(vjet) < Rtolerance && CompFlag == true)
		{MatchingJets.push_back(&(*jets));}
	}	
	
	// no match found --> set the flag and return NULL
	if (MatchingJets.size() == 0)
	{
		flag = -2;
		return NULL;
	}
		
	
	// if only one match, return the element; else find the best match
	if (MatchingJets.size() == 1)
	{
		flag = 0;
		return MatchingJets.at(0);
	}
	
	// multiple match --> find best
	else
	{
		int best = 0;
		setMomentum( vjet, (MatchingJets.at(0))->p4() );
		double Delta = vp.DeltaR(vjet);
		
		for (unsigned int i = 1; i < MatchingJets.size(); i++)
		{
			setMomentum( vjet, (MatchingJets.at(i))->p4() );
			if ( vp.DeltaR(vjet) < Delta )
			{
				best = i;
				Delta = vp.DeltaR(vjet);
			}
		}
	
		flag = -1;
		return MatchingJets.at(best);
	}

}


//define this as a plug-in
DEFINE_FWK_MODULE(GeneratorAna);
