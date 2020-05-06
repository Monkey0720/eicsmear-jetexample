

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
#include <cmath>

// fastjet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

// Convenience
using namespace fastjet;
using namespace std;
  
int main(int argc, char* argv[]){
  
  // Parse arguments
  // ---------------
  TString smearedname = "truth.smeared.root";
  TString outfilebase = "simplejet";
  int nevents = -1;

  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay=true;
  for ( auto parg = arguments.begin() ; parg!=arguments.end() ; ++parg){
    string arg=*parg;
    if ( arg == "-o" ){
      if (++parg == arguments.end() ){ argsokay=false; break; }
      outfilebase=*parg;
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
	 << " [-s smearedname] (root file)"  << endl
	 << " [-o OutFileBase] (extension will be added)"  << endl
      	 << " [-N Nevents] (<0 for all)" << endl
	 << endl;
    return -1;
  }

  // jet definition - could also be made parameters
  // ----------------------------------------------
  // Note: We're doing everything in the lab frame
  // In a real analysis, you may want to boost to for example
  // the Breit frame before and/or after jetfinding

  double R=0.5;
  JetDefinition jet_def = JetDefinition(kt_algorithm, R);

  // Only do jetfinding in the region covered by tracking AND calorimetry
  double ConstEtaMax = 3.5;
  double ConstEtaMin = -ConstEtaMax;

  double JetEtaMax = ConstEtaMax - R;
  double JetEtaMin = -JetEtaMax;
  double jetPtCut=4.0;
  auto select_jet_eta     = fastjet::SelectorEtaRange ( JetEtaMin + R, JetEtaMax -R  );
  auto select_jet_pt      = SelectorPtRange( jetPtCut, 1000 );
  auto select_jet         = select_jet_eta * select_jet_pt;     


  // Some cuts
  // ---------
  double xmin = 1e-5;
  double xmax = 0.99;
  double ymin = 0.01;
  double ymax = 0.95;
  double Q2min = 1;
  double Q2max = 1000000;

  // Afterburners etc.
  // -----------------
  // Here, we demonstrate a toy pT acceptance / efficiency combination
  // PURELY a toy, NO connection to reality
  TF1* eff = new TF1("eff","(x>[2]) * [0]*TMath::Erf(x-[1])",0, 3);

  // mostly 99%, dropping toward small pT, sharp cutoff at 0.2
  eff->SetParameters (0.99,-0.8, 0.2);
  
  // --------------
  // Load the trees
  // --------------
  TChain* inTreeS = new TChain("Smeared");
  inTreeS->Add(smearedname);

  // Setup Input Event Buffer
  Smear::Event* inEventS(NULL);
  inTreeS->SetBranchAddress("eventS",&inEventS);

  // Open histo root file and book histograms
  // ----------------------------------------
  TFile * outfile = new TFile ( outfilebase + ".root", "RECREATE");

  TH1::SetDefaultSumw2(true);
  float ptmin = 0;
  float ptmax = 20;
  int ptbins = 40;  
  TH1D* smearedpt=new TH1D("smearedpt",";p_{T}^{smeared} [GeV/c];counts", ptbins, ptmin, ptmax );
  
  // -------------
  // Analysis loop
  // -------------
  // Here we demonstrate using the smeared tree only.
  // Comparisons or afterburners will need access to the truth level
  // which can be accessed via the friend mechanism --> See the extended example.

  // Containers for jet constituents
  // -------------------------------
  vector<PseudoJet> SmearedConstituents;

  for(long iEvent=0; iEvent<inTreeS->GetEntries(); iEvent++){
    
    //Read next Event
    if(inTreeS->GetEntry(iEvent) <=0) break;
    if(iEvent%10000 == 0) cout << "Event " << iEvent << endl;
    
    // ---------------
    // event-wise cuts
    // ---------------
    bool smearacceptev=true;

    // For the high-Q^2 example, the double angle method does pretty well
    auto xS = inEventS->GetXDoubleAngle();
    auto yS = inEventS->GetYDoubleAngle();
    auto Q2S = inEventS->GetQ2DoubleAngle();
    if ( xS < xmin   || xS > xmax ) smearacceptev = false;
    if ( yS < ymin   || yS > ymax ) smearacceptev = false;
    if ( Q2S < Q2min || Q2S > Q2max ) smearacceptev = false;
    
    if ( !smearacceptev ) continue;

    // -------------------
    // Loop over Particles
    // -------------------
    SmearedConstituents.clear();

    for(int j=0; j<inEventS->GetNTracks(); j++){

      // first three particles are always beam and the virtual photon. skip.
      if ( j<3) continue;

      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle
      if ( !inParticleS ) continue; // lost this particle

      // Skip non-final particles. 
      if ( inParticleS->GetStatus() != 1 ) continue;

      // Now we need to work like an experimentalist.
      // For most particles, only partial information is available,
      // and we need to use something "reasonable" to supply the rest
      // This is complicated by the fact that in the default handbook, we don't have
      // PID information. We could pull that from the truth tree or turn on
      // a perfect PID detector, but here we do it the hard way, case by case.
      // 
      // Non-measured is indicated by a 0 in the field
      // (this is unfortunate, will be fixed in the future)

      double epsilon = 1e-7;
      double E = inParticleS->GetE();
      double px = inParticleS->GetPx();
      double py = inParticleS->GetPy();
      double pz = inParticleS->GetPz();
      double P = inParticleS->GetP(); // not used, but will help identification

      bool charged=( fabs(P)>epsilon );
      // Here is the right place for the pT efficiency afterburner.
      // If the particle is lost, we set P=0 and the logic below will treat it right
      if ( charged ){
	auto ran = gRandom->Uniform (0,1);
	auto pt = std::sqrt ( px*px + py*py );  // GetPt() may be better or worse - for now the same
	if ( ran > eff->Eval ( pt ) ) {
	  // if ( pt > 0.9 ) cout << inParticleS->GetPt () << "  " << pt << "  rejected by efficiency, ran = " << ran << " eff = " << eff->Eval ( pt ) << endl;
	  P=px=py=pz=0;
	}
      }

      // the following logic should be safe. a particle should be treated
      // by exactly one of the of cases.
      
      // Nothing measured:
      if ( fabs(P) <= epsilon && fabs(E) <= epsilon ){
	// This can happen when e.g. a low-p particle gets smeared below 0
	// Should be rare but should count as if inParticleS==0 in the first place
	continue;
      }

      // Both measured
      if ( fabs(P) > epsilon && fabs(E) > epsilon ){
	// seen by tracker and a calorimeter
	// either charged hadron or electron
	
	// Here an analyzer would have to use their judgement and decide
	// whether for some phase space discarding one measurement may make sense
	// because the other one has a better resolution.
	// In this example, we use both.

	/* nothing to do */
      }

      // Tracker only
      if ( fabs(P) > epsilon && fabs(E) <= epsilon ){
	// particle seems to be seen only by a tracker (or the calorimeter smeared E to <=0)	
	// For the handbook, this can only happen for the second case, but in principle
	// it indicates a region not covered by calorimetry.
	// What we should do is "cheat" and find out from the truth level if this particle
	// is hadronic or EM. Not a real cheat since we'd know in reality which detector is missing.
	// However, to keep this example free of truth information, we instead
	// assume it's hadronic (because an EIC detector w/o EMCal coverage is silly).
	// Typical assumptions are m=0 or m=m_pi. The latter is more realistic, so do that

	double m =  0.13957;  // GeV
	E = std::sqrt(P*P + m*m);	
      }

      // Calo only
      if ( fabs(P) <= epsilon && fabs(E) > epsilon ){
	// - neutral species (neutron, photon) that's covered only by calorimetry
	// - region that's only covered by calorimetry (_should_ imply EMCAL in a sane detector)
	// - missed track (or failed match)
	// This is the most difficult case and definitely more realistically treated by
	// getting at least the genre from the truth tree.
	// In the absence of such info here, we will reject all of these.
	// That's quite reasonable for HCAL info anyway, since without additional geometry information,
	// any HCAL hit without a pointing track is impossible to treat together with charged constituents.
	// BUT note that any direct-photon analysis is killed. See the extended example instead

	continue;
      }

      // Combine into a Pseudojet now. In general, a more featureful class may
      // be preferable (such as TParticle or better) to track pid etc.
      auto pj = PseudoJet (px, py, pz, E);

      // Particle cuts
      // -------------
      // The previous section made up for missing information.
      // Here, we can add cuts that define or refine our analysis

      // constituent eta cut
      if ( pj.eta() < ConstEtaMin || pj.eta() > ConstEtaMax) continue;


      // --------------------
      // Collect constituents
      // --------------------
      SmearedConstituents.push_back ( pj );
    }
    // --------------------
    // Perform jet analysis
    // --------------------
    ClusterSequence smearcs(SmearedConstituents, jet_def);
    vector<PseudoJet> smearjets = sorted_by_pt( select_jet(smearcs.inclusive_jets()) );

    // -------------------
    // Extract observables
    // -------------------
    for ( auto j : smearjets ){
      smearedpt->Fill ( j.pt() );
    }
    
  }

  new TCanvas;
  gPad->SetLogy();
  smearedpt->Draw();
  gPad->SaveAs( outfilebase + "_pt.png");
  
  outfile->Write();
  
  return 0;
}
