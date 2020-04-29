#include "jetspectra.hh"

// eicsmear
#include "eicsmear/erhic/EventBase.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/Particle.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/smear/EventS.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/ParticleMCS.h"

// ROOT
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TStyle.h>

// C++
#include <vector>
#include <map>

// fastjet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

// Convenience
using namespace fastjet;
using namespace std;
  
int main(int argc, char* argv[]){
  
  // Parse arguments
  // ---------------
  TString truthname="truth.root";
  TString smearedname = "";
  TString outfilebase = "jets";
  int nevents = -1;

  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay=true;
  for ( auto parg = arguments.begin() ; parg!=arguments.end() ; ++parg){
    string arg=*parg;
    if ( arg == "-o" ){
      if (++parg == arguments.end() ){ argsokay=false; break; }
      outfilebase=*parg;
    } else if ( arg == "-t" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      truthname=*parg;
    } else if ( arg == "-s" ){
      if (++parg ==arguments.end() ){ argsokay=false; break; }
      smearedname=*parg;
    } else if ( arg == "-N" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      nevents=std::stoi(parg->data());
    } else {
      argsokay=false;
      break;
    }
  }
  
  if ( !argsokay ) {
    cerr << "usage: " << argv[0] << endl
      	 << " [-t truthname] (root file)"  << endl
	 << " [-s smearedname] (root file)"  << endl
	 << " [-o OutFileBase] (extension will be added)"  << endl
      	 << " [-N Nevents] (<0 for all)" << endl
	 << endl;
    return -1;
  }

  if (smearedname == "" ){
    smearedname = truthname;
    smearedname.ReplaceAll (".root",".smeared.root" );
  }

  // jet definition - could also be made parameters
  // ----------------------------------------------
  double R=0.5;
  JetDefinition jet_def = JetDefinition(kt_algorithm, R);

  // Note: We're doing everything in the lab frame
  // In a real analysis, you may want to boost to for example
  // the Breit frame before and/or after jetfinding
  double etaMax = 4.0;
  double etaMin = -etaMax;
  double jetPtCut=5.0;
  auto select_jet_eta     = fastjet::SelectorEtaRange ( etaMin + R, etaMax -R  );
  auto select_jet_pt      = SelectorPtRange( jetPtCut, 1000 );
  auto select_jet         = select_jet_eta * select_jet_pt;     


  // some cuts
  // ---------
  double xmin = 1e-4;
  double xmax = 0.99;
  double ymin = 1e-4;
  double ymax = 0.95;
  double Q2min = 1;
  double Q2max = 1000000;

  // --------------
  // Load the trees
  // --------------
  TChain* inTree = new TChain("EICTree");
  inTree->Add(truthname);
  inTree->AddFriend("Smeared",smearedname);
  
  // Setup Input Event Buffer
  erhic::EventMC* inEvent(NULL);
  Smear::Event* inEventS(NULL);
  inTree->SetBranchAddress("event",&inEvent);
  inTree->SetBranchAddress("eventS",&inEventS);

  // Open histo root file and book histograms
  // ----------------------------------------
  TFile * outfile = new TFile ( outfilebase + ".root", "RECREATE");

  TH1::SetDefaultSumw2(true);
  float ptmin = 0;
  float ptmax = 20;
  int ptbins = 80;  
  TH1D* truthpt=new TH1D("truthpt",";p_{T}^{truth} [GeV/c];counts", ptbins, ptmin, ptmax );

  
  // -------------
  // Analysis loop
  // -------------
  // Important: We're running over the truth tree and getting smeared information
  // via the friend mechanism. This is only necessary for truth-smeared comparisons.
  // One can run purely on the smeared tree as well without access to the truth.

  // Containers for jet constituents
  // -------------------------------
  vector<PseudoJet> TruthConstituents;
  vector<PseudoJet> SmearedConstituents;

  for(long iEvent=0; iEvent<inTree->GetEntries(); iEvent++){
    
    //Read next Event
    if(inTree->GetEntry(iEvent) <=0) break;
    if(iEvent%10000 == 0) cout << "Event " << iEvent << endl;
    
    // ---------------
    // event-wise cuts
    // ---------------
    bool truthacceptev = true;
    bool smearacceptev = true;

    // truth level
    if ( inEvent->GetX() < xmin   || inEvent->GetX() > xmax ) truthacceptev = false;
    if ( inEvent->GetY() < ymin   || inEvent->GetY() > ymax ) truthacceptev = false;
    if ( inEvent->GetQ2() < Q2min || inEvent->GetQ2() > Q2max ) truthacceptev = false;
    
    if ( !truthacceptev && !smearacceptev ) continue;

    // -------------------
    // Loop over Particles
    // -------------------
    TruthConstituents.clear();
    SmearedConstituents.clear();

    for(int j=0; j<inEventS->GetNTracks(); j++){
      // Skip beam
      if ( j<3 ) continue;

      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle      
      const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle

      // Particle cuts for truth level
      bool truthacceptp = truthacceptev;
      // Skip non-final particles. 
      if ( inParticle->GetStatus() != 1 ) truthacceptp = false;

      // if ( inParticleS->GetStatus() != 1 ) truthacceptp;
      // bool smearacceptp = smearacceptev;
      // // Particle was not smeared
      // if(inParticleS == NULL) continue; 

      // --------------------
      // Collect constituents
      // --------------------
      if ( truthacceptp ){
	TruthConstituents.push_back ( PseudoJet (inParticle->GetPx(),
						 inParticle->GetPy(),
						 inParticle->GetPz(),
						 inParticle->GetE()) );

      } 
    }
    // --------------------
    // Perform jet analysis
    // --------------------
    ClusterSequence truthcs(TruthConstituents, jet_def);
    vector<PseudoJet> truthjets = sorted_by_pt( select_jet(truthcs.inclusive_jets()) );

    // -------------------
    // Extract observables
    // -------------------
    for ( auto j : truthjets ){
      truthpt->Fill ( j.pt() );
    }

    
  }
  outfile->Write();
  
  return 0;
}

// ---------------------------------------------------------------
// void PlotQA ( const qaparameters& qapars, eventqacollection& eventqa, map<int,pidqacollection>& qabook ){

//   // Stat  position and size
//   // -----------------------
//   gStyle->SetStatX(0.25);
//   gStyle->SetStatW(0.15);
//   gStyle->SetStatY(0.9);
//   gStyle->SetStatH(0.15);

//   // Position of the "Missed: " box
//   float missx = 0.55;
//   float missy = 0.2;
//   float missy2 = 0.8;
//   TText t;
//   t.SetNDC();

//   // prep a pdf collection
//   new TCanvas;
//   gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf[" );

//   // event-wise qa

//   // response-style
//   // NM
//   if ( eventqa.y_NM ) {
//     eventqa.y_NM->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_NM, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.x_NM ) {
//     eventqa.x_NM->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_NM, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.Q2_NM ) {
//     eventqa.Q2_NM->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_NM, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   // DA
//   if ( eventqa.y_DA ) {
//     eventqa.y_DA->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_DA, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.x_DA ) {
//     eventqa.x_DA->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_DA, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.Q2_DA ) {
//     eventqa.Q2_DA->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_DA, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   // JB
//   if ( eventqa.y_JB ) {
//     eventqa.y_JB->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_JB, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.x_JB ) {
//     eventqa.x_JB->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_JB, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.Q2_JB ) {
//     eventqa.Q2_JB->Draw("colz");
//     t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_JB, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }

//   // resolution-style
//   // NM
//   if ( eventqa.dely_NM ) {
//     eventqa.dely_NM->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedy_NM, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.delx_NM ) {
//     eventqa.delx_NM->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedx_NM, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.delQ2_NM ) {
//     eventqa.delQ2_NM->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedQ2_NM, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   // DA
//   if ( eventqa.dely_DA ) {
//     eventqa.dely_DA->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedy_DA, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.delx_DA ) {
//     eventqa.delx_DA->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedx_DA, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.delQ2_DA ) {
//     eventqa.delQ2_DA->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedQ2_DA, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   // JB
//   if ( eventqa.dely_JB ) {
//     eventqa.dely_JB->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedy_JB, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.delx_JB ) {
//     eventqa.delx_JB->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedx_JB, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   if ( eventqa.delQ2_JB ) {
//     eventqa.delQ2_JB->Draw("colz");
//     t.DrawText( missx,missy2, Form("Missed: %ld / %ld",eventqa.missedQ2_JB, qapars.usedevents));
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }

  

//   // particle QA
//   // -----------
//   gStyle->SetStatX(0.55); // reposition stat box
//   for ( auto& pidcoll : qabook ){
//     auto& pid = pidcoll.first;
//     auto& coll = pidcoll.second;

//     // option "s" in Profile shows rms
    
//     coll.DelP_th->Draw("colz");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

//     coll.DelP_eta->Draw("colz");
//     coll.DelP_eta->ProfileX("_px",1,-1,"s")->Draw("same");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

//     coll.DelE_th->Draw("colz");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
    
//     coll.DelE_eta->Draw("colz");
//     coll.DelE_eta->ProfileX("_px",1,-1,"s")->Draw("same");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

//     coll.DelE_E->Draw("colz");
//     coll.DelE_E->ProfileX("_px",1,-1,"s")->Draw("same");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

//     coll.dTh_p->Draw("colz");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
    
//     coll.dEta_p->Draw("colz");
//     coll.dEta_p->ProfileX("_px",1,-1,"s")->Draw("same");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );

//     coll.dPhi_p->Draw("colz");
//     coll.dPhi_p->ProfileX("_px",1,-1,"s")->Draw("same");
//     gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf" );
//   }
//   // close the pdf collection
//   gPad->SaveAs( qapars.outfilebase + qapars.detstring + ".pdf]" );
// }
// // ---------------------------------------------------------------
