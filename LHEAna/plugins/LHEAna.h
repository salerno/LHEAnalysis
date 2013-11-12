#ifndef LHEAna_H
#define LHEAna_H


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ROOT include
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

//
// class declaration
//

using namespace std;

class LHEAna : public edm::EDAnalyzer {
public:
  explicit LHEAna(const edm::ParameterSet&);
  ~LHEAna();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  //TFile * OutputFile;

  TH1F* h_genH_mass;
  TH1F* h_genH_pt;
  TH1F* h_genH_rap;
  
  TH1F* h_genZ_mass;

  TH1F* h_genZ1_mass;
  TH1F* h_genZ2_mass;
  TH1F* h_genZ1_massnoH;
  TH1F* h_genZ2_massnoH;
  
  TH1F* h_genlep_pt;
  TH1F* h_genlep_eta;
  TH1F* h_genlep_phi;
    
  TH2F * h2_genZ1_mass_vs_genZ2_mass;
  TH2F * h2_genZ1_mass_vs_genZ2_mass_noH; 
  
  TH1F* h_costhetastar;
  TH1F* h_costheta1;
  TH1F* h_costheta2;
  TH1F* h_Phi;
  TH1F* h_Phi1;

  //
  TH1F* h_ptmin;
  TH1F* h_etamax;

  // 
  TString * cutdes; //[40];
  
  int _icut;
  int _icut_MAX;
  TH1F* ICUT;
  
  int _ninitial;
  int _nselected;

  int _n4e, _n4mu, _n2e2mu;

  // ----------member data ---------------------------
};
#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//
