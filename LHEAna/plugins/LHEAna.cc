#include "LHEAna.h"

// CMSSW include
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT include
#include "TLorentzVector.h"


using namespace edm;

// ------------ method called for each event  ------------
// ===============================================================================
void LHEAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ===============================================================================
{

  _ninitial++;

  for(int iloop=0;iloop<1;iloop++) {  // UGLY !!!!!

  edm::Handle <LHEEventProduct> LHEevt;
  iEvent.getByType( LHEevt );
  const lhef::HEPEUP hepeup = LHEevt->hepeup();
  
  int genLHEparton_pdgid[20];//      = 0;
  int genLHEparton_status[20];//     = 0; 
  double genLHEparton_mother1[20];// = 0.;
  double genLHEparton_mother2[20]; // = 0.;
  double genLHEparton_col1[20]; //    = 0.;
  double genLHEparton_col2[20]; //    = 0.;
  double   genLHEparton_px[20]; // = 0;
  double genLHEparton_py[20]; //   = 0;
  double genLHEparton_pz[20]; //   = 0; 
  double genLHEparton_e[20]; //    = 0;
  double genLHEparton_m[20]; //    = 0;

  TLorentzVector gen_H;
  TLorentzVector gen_Z1; 
  TLorentzVector gen_Z2;
  TLorentzVector gen_lep1;
  TLorentzVector gen_lep2;
  TLorentzVector gen_lep3;
  TLorentzVector gen_lep4;
  
  int n_Z   = 0;
  int n_lep = 0;

  int n_ele = 0;
  int n_mu  = 0;

  bool isHiggs = false;

  _icut = 0;

  for(int i = 0; i < hepeup.NUP; i++) {
    genLHEparton_pdgid[i]   = hepeup.IDUP[i];
    genLHEparton_status[i]  = hepeup.ISTUP[i];
    genLHEparton_mother1[i] = hepeup.MOTHUP[i].first;
    genLHEparton_mother2[i] = hepeup.MOTHUP[i].second;
    genLHEparton_col1[i]    = hepeup.ICOLUP[i].first;
    genLHEparton_col2[i]    = hepeup.ICOLUP[i].second;
    genLHEparton_px[i]      = hepeup.PUP[i][0];
    genLHEparton_py[i]      = hepeup.PUP[i][1];
    genLHEparton_pz[i]      = hepeup.PUP[i][2];
    genLHEparton_e[i]       = hepeup.PUP[i][3];
    genLHEparton_m [i]      = hepeup.PUP[i][4];

    //   cout << "pdgid = " << genLHEparton_pdgid[i] 
    // 	 << " mother 1 = " << genLHEparton_mother1[i]
    // 	 << " mother 2 = " << genLHEparton_mother2[i]
    // 	 << " col1     = " << genLHEparton_col1[i]
    // 	 << " col2     = " << genLHEparton_col2[i]
    // 	 << " mass     = " << genLHEparton_m[i]
    // 	 << endl;

    if(fabs(genLHEparton_pdgid[i])==5000000) { 
    //if(fabs(genLHEparton_pdgid[i])==25) {
    //if(fabs(genLHEparton_pdgid[i])==39) {
      gen_H.SetPxPyPzE(genLHEparton_px[i], genLHEparton_py[i], genLHEparton_pz[i], genLHEparton_e[i]);
      isHiggs = true;
    } // if Higgs

    if(fabs(genLHEparton_pdgid[i])==23) { 
      if(n_Z==0) { gen_Z1.SetPxPyPzE(genLHEparton_px[i], genLHEparton_py[i], genLHEparton_pz[i], genLHEparton_e[i]); n_Z++;}
      else { gen_Z2.SetPxPyPzE(genLHEparton_px[i], genLHEparton_py[i], genLHEparton_pz[i], genLHEparton_e[i]); n_Z++;}
    } // if pdgid == Z
    

    // cout << "nlep = " << n_lep <<  " genLHEparton_pdgid = "<< genLHEparton_pdgid[i] << endl;
    if(fabs(genLHEparton_pdgid[i])==11 || fabs(genLHEparton_pdgid[i])==13) {

      if(fabs(genLHEparton_pdgid[i])==11) n_ele++;
      if(fabs(genLHEparton_pdgid[i])==13) n_mu++; 

      if(n_lep==0) { gen_lep1.SetPxPyPzE(genLHEparton_px[i], genLHEparton_py[i], genLHEparton_pz[i], genLHEparton_e[i]);} 
      if(n_lep==1) { gen_lep2.SetPxPyPzE(genLHEparton_px[i], genLHEparton_py[i], genLHEparton_pz[i], genLHEparton_e[i]);} 
      if(n_lep==2) { gen_lep3.SetPxPyPzE(genLHEparton_px[i], genLHEparton_py[i], genLHEparton_pz[i], genLHEparton_e[i]);}  
      if(n_lep==3) { gen_lep4.SetPxPyPzE(genLHEparton_px[i], genLHEparton_py[i], genLHEparton_pz[i], genLHEparton_e[i]);}
      n_lep ++;
    } // if pdgid == e, mu
        
  } // for loop on partons
  
  //cout << " pt = " << gen_lep1.Pt() << " " << gen_lep2.Pt() << " " << gen_lep3.Pt() << " " << gen_lep4.Pt() << endl;

  // ----------------------------------------------------------------------
  //                       Build TLorentzVector
  // ----------------------------------------------------------------------
  
  TLorentzVector myZ1; myZ1 = gen_lep1 + gen_lep2;
  TLorentzVector myZ2; myZ2 = gen_lep3 + gen_lep4;
  
  TLorentzVector myH; myH = myZ1 + myZ2;
  
  // ----------------------------------------------------------------------
  //                       Select Channel
  // ----------------------------------------------------------------------
  //cout << "isHiggs : " << isHiggs  << endl;
  if(isHiggs==false) continue;  cutdes[_icut] = "At least one H"; ICUT->Fill((Float_t)_icut); _icut++;

  //cout << "nZ : " << n_Z  << endl;

  if(n_Z!=2) continue;  cutdes[_icut] = "Two Z's"; ICUT->Fill((Float_t)_icut); _icut++;

  // Channel
  bool is4e    = false;
  bool is4mu   = false;
  bool is2e2mu = false;
  
  if(n_ele>=4)            { is4e    = true; _n4e++;  }
  if(n_mu>=4)             { is4mu   = true; _n4mu++; }
  if(n_ele>=2 && n_mu>=2) { is2e2mu = true; _n2e2mu++; }

  if(is2e2mu==false) continue; cutdes[_icut] = "Channel 2e2mu"; ICUT->Fill((Float_t)_icut); _icut++;
  
  // ----------------------------------------------------------------------
  //                       Compute Angles
  // ----------------------------------------------------------------------

  float costhetastar = 0.;
  float costheta1    = 0.; 
  float costheta2    = 0.; 
  float Phi          = 0.;
  float Phi1         = 0.;

  //// costhetastar
  TVector3 boostX = -(myH.BoostVector());
  TLorentzVector thep4Z1inXFrame( myZ1 );
  TLorentzVector thep4Z2inXFrame( myZ2 );
  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );    
  costhetastar = theZ1X_p3.CosTheta();
  
  //// --------------------------- costheta1
  TVector3 boostV1 = -(myZ1.BoostVector());
  TLorentzVector p4M11_BV1( gen_lep1 );
  TLorentzVector p4M12_BV1( gen_lep2 );
  TLorentzVector p4M21_BV1( gen_lep3 );
  TLorentzVector p4M22_BV1( gen_lep4 );
  p4M11_BV1.Boost( boostV1 );
  p4M12_BV1.Boost( boostV1 );
  p4M21_BV1.Boost( boostV1 );
  p4M22_BV1.Boost( boostV1 );
  
  TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
  //// costheta1
  costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
  
  //// --------------------------- costheta2
  TVector3 boostV2 = -(myZ2.BoostVector());
  if (boostV2.Mag()>=1.) {
    cout << "Warning: Mela::computeAngles: Z2 boost with beta=1, scaling down" << endl;
    boostV2*=0.9999;
  }
  TLorentzVector p4M11_BV2( gen_lep1 );
  TLorentzVector p4M12_BV2( gen_lep2 );
  TLorentzVector p4M21_BV2( gen_lep3 );
  TLorentzVector p4M22_BV2( gen_lep4 );
  p4M11_BV2.Boost( boostV2 );
  p4M12_BV2.Boost( boostV2 );
  p4M21_BV2.Boost( boostV2 );
  p4M22_BV2.Boost( boostV2 );
  
  TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
  //// costheta2
  costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
  
  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  //    TVector3 boostX = -(themyH.BoostVector());
  TLorentzVector p4M11_BX( gen_lep1 );
  TLorentzVector p4M12_BX( gen_lep2 );
  TLorentzVector p4M21_BX( gen_lep3 );
  TLorentzVector p4M22_BX( gen_lep4 );
  
  p4M11_BX.Boost( boostX );
  p4M12_BX.Boost( boostX );
  p4M21_BX.Boost( boostX );
  p4M22_BX.Boost( boostX );
  
  TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
  TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    
  
  TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
  TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 
  
  //// Phi
  TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
  float tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
  float sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
  
  
  //////////////
  
  TVector3 beamAxis(0,0,1);
  TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
  
  TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
  TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
  TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
  
  //// Phi1
  float tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
  float sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
  Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    
  
  // ----------------------------------------------------------------------
  //               Compute Pt, eta min/max of leptons
  // ----------------------------------------------------------------------
  
  double mass_Z1, mass_Z2;
  double rap_Z1, rap_Z2;
  double pt_Z1, pt_Z2;
  if(myZ1.M() >= myZ2.M()) { mass_Z1 = myZ1.M(); mass_Z2 = myZ2.M(); }
  else { mass_Z1 = myZ2.M(); mass_Z2 = myZ1.M(); }

  double * ptlep = new double[4];
  ptlep[0] = gen_lep1.Pt();
  ptlep[1] = gen_lep2.Pt();
  ptlep[2] = gen_lep3.Pt();
  ptlep[3] = gen_lep4.Pt();

  double * etalep = new double[4];
  etalep[0] = gen_lep1.Eta();
  etalep[1] = gen_lep2.Eta();
  etalep[2] = gen_lep3.Eta();
  etalep[3] = gen_lep4.Eta();

  double pt_min  = ptlep[0];
  double eta_max = etalep[0];
  
  for(int ilep=0;ilep<3;ilep++) {
    //cout << "ptmin = " << pt_min << " pT = " << ptlep[ilep+1] << endl;
    if(pt_min>ptlep[ilep+1])   pt_min  = ptlep[ilep+1];
    if(fabs(eta_max) < fabs(etalep[ilep+1])) eta_max = etalep[ilep+1]; 
  } // for loop on leptons
  //double eta1, eta2, eta3, eta4;
  //cout << " pT    = " << ptlep[0] << " " << ptlep[1] << " " << ptlep[2] << " " << ptlep[3] << endl;
  //cout << " pTmin = " << pt_min << endl;

    
  if(isHiggs) {
    
    // Fill Histos
    h_genH_mass->Fill(gen_H.M());
    h_genH_rap->Fill(gen_H.Rapidity());
    h_genH_pt->Fill(gen_H.Pt());

    h_genZ_mass->Fill(gen_Z1.M());
    
    h_genZ1_mass->Fill(mass_Z1);
    h_genZ2_mass->Fill(mass_Z2);

    h2_genZ1_mass_vs_genZ2_mass->Fill(mass_Z1, mass_Z2);

    // Angles
    h_costhetastar->Fill(costhetastar);
    h_costheta1->Fill(costheta1);
    h_costheta2->Fill(costheta2);
    h_Phi->Fill(Phi);
    h_Phi1->Fill(Phi1);

    // Pt, Eta
    h_ptmin->Fill(pt_min);
    h_etamax->Fill(eta_max);

  } // if Higgs was generated
  else {
    
    h_genZ1_massnoH->Fill(mass_Z1);
    h_genZ2_massnoH->Fill(mass_Z2);

    h2_genZ1_mass_vs_genZ2_mass_noH->Fill(mass_Z1, mass_Z2);

  } // else


// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
   
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif



  // -----------------------------
  _icut_MAX = _icut;

  _nselected++;


  } // FOR intial statement...
 
}


// ------------ method called once each job just before starting event loop  ------------
void 
LHEAna::beginJob()
{

  edm::Service<TFileService> fs;
  
  double pi = acos(-1);
   
  //OutputFile = new TFile("outfile.root","recreate");
  
  //  TH1F* h_genH_mass;
  //   TH1F* h_genZ_mass;
  
  //   TH1F* h_genlep_pt;
  //   TH1F* h_genlep_eta;
  //   TH1F* h_genlep_phi;
   
  h_genH_mass  = fs->make<TH1F>("genH_mass", "genH_mass", 200,0,200);
  h_genH_pt    = fs->make<TH1F>("genH_pt", "genH_pt", 200,0,200);
  h_genH_rap   = fs->make<TH1F>("genH_rap", "genH_rap", 200,-5, 5);
 
  h_genZ_mass  = fs->make<TH1F>("genZ_mass", "genZ_mass", 200,0,200);
  h_genZ1_mass = fs->make<TH1F>("genZ1_mass", "genZ1_mass", 200,0,200);
  h_genZ2_mass = fs->make<TH1F>("genZ2_mass", "genZ2_mass", 200,0,200);
  
  h_genZ1_massnoH = fs->make<TH1F>("genZ1_massnoH", "genZ1_massnoH", 200,0,200);
  h_genZ2_massnoH = fs->make<TH1F>("genZ2_massnoH", "genZ2_massnoH", 200,0,200);

  h2_genZ1_mass_vs_genZ2_mass     = fs->make<TH2F>("m1_vs_m2", "m1_vs_m2", 200, 0, 200, 200, 0, 200); 
  h2_genZ1_mass_vs_genZ2_mass_noH = fs->make<TH2F>("m1_vs_m2_noH", "m1_vs_m2_noH", 200, 0, 200, 200, 0, 200); 
    
  h_costhetastar = fs->make<TH1F>("costhetastar", "costhetastar", 100, -1, 1);
  h_costheta1    = fs->make<TH1F>("costheta1", "costheta1", 100, -1, 1);
  h_costheta2    = fs->make<TH1F>("costheta2", "costheta2", 100, -1, 1);
  h_Phi          = fs->make<TH1F>("Phi", "Phi", 100, -pi, pi);
  h_Phi1         = fs->make<TH1F>("Phi1", "Phi1", 100, -pi, pi);

  //
  h_ptmin        = fs->make<TH1F>("ptmin","ptmin",200,0,200);
  h_etamax       = fs->make<TH1F>("etamax","etamax",500,-5, 5);

  // 
  ICUT 		 = fs->make<TH1F>("ICUT","ICUT",50,0.,50.);  
  
  cutdes = new TString[40];
  
  _icut_MAX  = 0;
  _nselected = 0;
  _ninitial  = 0;

  _n4e    = 0;
  _n4mu   = 0;
  _n2e2mu = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
LHEAna::endJob() 
{

  cout << "icut_max = " << _icut_MAX << endl;
  cout << "=============================================================" << endl;
  cout << "Total Number of INITIAL events  = "  << _ninitial                       << endl;
  cout << "=============================================================" << endl;
  for(int i=0;i<_icut_MAX+1;i++) 
    cout << "Cut " << setw(3) << i << " : " << setw(30) << cutdes[i] << " : " << ICUT->GetBinContent(i+1) << endl;
  cout << "=============================================================" << endl;
  cout << "Total Number of events selected = "  << _nselected                       << endl;
  cout << "=============================================================" << endl;
  cout << "Number of 4e    : " <<  _n4e << endl;
  cout << "Number of 4mu   : " <<  _n4mu << endl;
  cout << "Number of 2e2mu : " <<  _n2e2mu << endl;
  
   

 //  OutputFile->Write();
//   OutputFile->Close();

}

// ------------ method called when starting to processes a run  ------------
void 
LHEAna::beginRun(edm::Run const&, edm::EventSetup const&)
{

 

}

// ------------ method called when ending the processing of a run  ------------
void 
LHEAna::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
LHEAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
LHEAna::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LHEAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//
// constructors and destructor
//
LHEAna::LHEAna(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
}


LHEAna::~LHEAna()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//




//define this as a plug-in
DEFINE_FWK_MODULE(LHEAna);
