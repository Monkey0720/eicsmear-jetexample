

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

#ifdef GROOMING
// fjcontrib
#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"
#endif


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


  // some cuts
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
  TF1* eff = new TF1("eff","(x>[2]) * [0]*TMath::Erf(x-[1])",0, 100);

  // mostly 99%, dropping toward small pT, sharp cutoff at 0.2
  eff->SetParameters (0.99,-0.8, 0.2);

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
  int ptbins = 40;  
  TH1D* truthpt=new TH1D("truthpt",";p_{T}^{truth} [GeV/c];counts", ptbins, ptmin, ptmax );
  TH1D* smearpt=new TH1D("smearpt",";p_{T}^{smear} [GeV/c];counts", ptbins, ptmin, ptmax );
  TH2D* smearvtruept=new TH2D("smearvtruept",";p_{T}^{truth} [GeV/c];p_{T}^{smear} [GeV/c];counts", ptbins, ptmin, ptmax, ptbins, ptmin, ptmax );

  // brief example for substructure and other fjcontrib tools
#ifdef GROOMING
  int zgbins = 30;
  double zgmin = 0;
  double zgmax = 0.6;
  TH2D* smearvtruezg=new TH2D("smearvtruezg",";z_{g}^{truth};z_{g}^{smear};counts", zgbins, zgmin, zgmax, zgbins, zgmin, zgmax );
  double beta = 0;
  double z_cut = 0.1;
  contrib::SoftDrop sd( beta, z_cut);
#endif

  
  // -------------
  // Analysis loop
  // -------------
  // Important: We're running over the truth tree and getting smeared information
  // via the friend mechanism. This is only necessary for truth-smeared comparisons.
  // One can run purely on the smeared tree as well without access to the truth.

  // Book-keping
  long rejectevent=0;
  long keptevent=0;
  long lostevent=0;
  long fakeevent=0;
  

  // Containers for jet constituents
  // -------------------------------
  vector<PseudoJet> TruthConstituents;
  vector<PseudoJet> SmearedConstituents;
  if ( nevents < 0 )  nevents = inTree->GetEntries();

  for(long iEvent=0; iEvent<nevents; iEvent++){
    
    //Read next Event
    if(inTree->GetEntry(iEvent) <=0) break;
    if(iEvent%10000 == 0) cout << "Event " << iEvent << endl;
    
    // ---------------
    // event-wise cuts
    // ---------------
    bool truthacceptev = true;
    bool smearacceptev = true;

    // truth level
    // ----------- 
    if ( inEvent->GetX() < xmin   || inEvent->GetX() > xmax ) truthacceptev = false;
    if ( inEvent->GetY() < ymin   || inEvent->GetY() > ymax ) truthacceptev = false;
    if ( inEvent->GetQ2() < Q2min || inEvent->GetQ2() > Q2max ) truthacceptev = false;

    // Detector level
    // --------------
    // For the high-Q^2 example, the double angle method does pretty well
    auto xS = inEventS->GetXDoubleAngle();
    auto yS = inEventS->GetYDoubleAngle();
    auto Q2S = inEventS->GetQ2DoubleAngle();
    if ( xS < xmin   || xS > xmax ) smearacceptev = false;
    if ( yS < ymin   || yS > ymax ) smearacceptev = false;
    if ( Q2S < Q2min || Q2S > Q2max ) smearacceptev = false;

    if ( !truthacceptev ) rejectevent++;
     
    if ( !truthacceptev && !smearacceptev ) continue;

    // ----------------
    // Categorize event
    // ----------------
    // We have three classes
    // - Accepted at truth and detector level
    // - Accepted at truth level, LOST via smearing
    // - Rejected at truth level, GAINED (FAKE) via smearing
    // How to handle each depends on details of the analysis,
    // and proper interpretation also requires understanding of the input.
    // Here, all we do is book-keeping
    if (  truthacceptev &&  smearacceptev ) keptevent++;
    if (  truthacceptev && !smearacceptev ) lostevent++;
    if ( !truthacceptev &&  smearacceptev ) fakeevent++;

    // -------------------
    // Loop over Particles
    // -------------------
    TruthConstituents.clear();
    SmearedConstituents.clear();

    for(int j=0; j<inEventS->GetNTracks(); j++){
      // first three particles are always beam and the virtual photon. skip.
      if ( j<3) continue;

      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle      
      const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle

      // Skip non-final particles. 
      if ( inParticle->GetStatus() != 1 ) continue;

      // Particle cuts for truth level
      // -----------------------------
      bool truthacceptp = truthacceptev;
      // only imposing constituent eta cut here
      if ( inParticle->GetEta() < ConstEtaMin || inParticle->GetEta() > ConstEtaMax) truthacceptp=false;

      // Particle cuts for smeared (detector) level
      // ------------------------------------------
      bool smearacceptp = smearacceptev;
      if ( !inParticleS ) smearacceptp=false; // lost this particle 
	    
      // Now we need to work like an experimentalist.
      // For most particles, only partial information is available,
      // and we need to use something "reasonable" to supply the rest
      // This is complicated by the fact that in the default handbook, we don't have
      // PID information. We can pull some from the truth tree,
      // but make sure to not use information we "shouldn't" have
      
      // Non-measured is indicated by a 0 in the field
      // (this is unfortunate, will be fixed in the future)

      PseudoJet pj;
      if ( inParticleS ) {
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
	    P=px=py=pz=0;
	  }
	}	

	// the following logic should be safe. a particle should be treated
	// by exactly one of the of cases. 
	
	// Nothing measured:
	if ( fabs(P) <= epsilon && fabs(E) <= epsilon ){
	  // This can happen when e.g. a low-p particle gets smeared below 0
	  // Should be rare but should count as if inParticleS==0 in the first place
	  smearacceptp=false;
	}
	
	// Both measured
	if ( fabs(P) > epsilon && fabs(E) > epsilon ){
	  // seen by tracker and a calorimeter
	  // either charged hadron or electron
	  
	  // Here an analyzer would have to use their judgement and decide
	  // whether for some phase space discarding one measurement may make sense
	  // because the other one has a better resolution.
	  // In this example, we use both.

	  // Note: In this case, P and E are smeared independently, so in general
	  // the mass-energy relation is not conserved. FastJet does not mind at all.

	  /* nothing to do */
	}
	
	// Tracker only
	if ( fabs(P) > epsilon && fabs(E) <= epsilon ){
	  // particle seems to be seen only by a tracker (or the calorimeter smeared E to <=0)
	  // For the handbook, this can only happen for the second case, but in principle
	  // it indicates a region not covered by calorimetry.

	  // We tap into the truth level information to make up for a framework short-coming:
	  // A real analyser knows whether the (missing!) calo information should be EM or hadronic, but the
	  // ParticleMCS class has lost that information
	  auto abspid=std::abs(inParticle->GetPdgCode());
	  bool IsElectroMagnetic = ( abspid==22 || abspid == 11);	  
	  
	  if (IsElectroMagnetic){
	    // This must be an electron (or in theory a mismatch)
	    double m =  0.000511;  // GeV
	    E = std::sqrt(P*P + m*m);
	  } else {
	    // A charged hadron. A typical assumption is a pion
	    double m =  0.13957;  // GeV
	    E = std::sqrt(P*P + m*m);
	  }
	}
	
	// Calo only
	if ( fabs(P) <= epsilon && fabs(E) > epsilon ){
	  // - neutral species (neutron, photon) that's covered only by calorimetry
	  // - region that's only covered by calorimetry (_should_ imply EMCAL in a sane detector)
	  // - missed track (or failed match)

	  // IMPORTANT: In the following, we assume that we have phi, eta information
	  // (i.e. from the struck detector region).
	  // HOWEVER, this information is unnaturally correct for the case of a charged particle.
	  // The conservative thing to do is to reject all of these particles.
	  smearacceptp = false;

	  // // If you assume you can somehow ascertain whether this particle is charged,
	  // // you could instead proceed along the following lines
	  // // bool charged= OutsideKnowledgeFunction();
	 
	  // // Again tap into the truth level information to make up for a framework short-coming:
	  // // A real analyser knows whether the calo information is EM or hadronic, but the
	  // // ParticleMCS class has lost that information
	  // auto abspid=std::abs(inParticle->GetPdgCode());
	  // bool IsElectroMagnetic = ( abspid==22 || abspid == 11);
	  // double m=0;
	  
	  // if (IsElectroMagnetic){
	  //   // Either a photon or an electron with a lost pointing track
	  //   // Either way, assume P=E and we're not far off
	  //   m=0;
	  // } else {
	  //   // A neutral hadron or a charged one that lost its track.
	  //   // P^2 = E^2 - m^2
	  //   if ( charged ){
	  //     // charged: most likely a pion
	  //     m =  0.13957;  // GeV
	  //   } else {
	  //     // neutral: neutron or K0L;
	  //     // In the absence of further knowledge, either reject this case or assume
	  //     // as an example neutron mass;
	  //     m = 0.93957;  // GeV
	  //   }
	    
	  // }
	  
	  // auto P2 = E*E- m*m;
	  // if ( P2 < 0 ) {
	  //   // Something's not right. reject.
	  //   smearacceptp = false;
	  // } else {
	  //   P = sqrt ( P2);
	  //   auto phi = inParticleS->GetPhi();
	  //   auto theta = inParticleS->GetTheta();
	  //   px = P * sin(theta) * cos(phi);
	  //   py = P * sin(theta) * sin(phi);
	  //   pz = P * cos(theta);
	  // }

	}
	
	// Combine into a Pseudojet now. In general, a more featureful class may
	// be preferable (such as TParticle or better) to track pid etc.
	pj = PseudoJet (px, py, pz, E);
      }
	
      // Particle cuts
      // -------------
      // The previous section made up for missing information.
      // Here, we can add cuts that define or refine our analysis

      // constituent eta cut
      if ( pj.eta() < ConstEtaMin || pj.eta() > ConstEtaMax) smearacceptp=false;
      
      // --------------------
      // Collect constituents
      // --------------------
      if ( truthacceptp ){
	TruthConstituents.push_back ( PseudoJet (inParticle->GetPx(),
						 inParticle->GetPy(),
						 inParticle->GetPz(),
						 inParticle->GetE()) );
      } 
      if ( smearacceptp ){
	SmearedConstituents.push_back ( pj );
      } 
      
    }
    // --------------------
    // Perform jet analysis
    // --------------------
    ClusterSequence truthcs(TruthConstituents, jet_def);
    vector<PseudoJet> truthjets = sorted_by_pt( select_jet(truthcs.inclusive_jets()) );

    ClusterSequence smearcs(SmearedConstituents, jet_def);
    vector<PseudoJet> smearjets = sorted_by_pt( select_jet(smearcs.inclusive_jets()) );

    // -------------------
    // Extract observables
    // -------------------
    for ( auto tj : truthjets ){
      truthpt->Fill ( tj.pt() );
    }

    for ( auto sj : smearjets ){
      smearpt->Fill ( sj.pt() );
    }      

    // -------------------------
    // Simple geometric matching
    // -------------------------
    fastjet::Selector SelectClose = fastjet::SelectorCircle( R );
    for ( auto tj : truthjets ){
      SelectClose.set_reference (tj);
      vector<PseudoJet> matchedjets = sorted_by_pt( SelectClose (smearjets) );
      // Shouldn't have more than 1 match
      if ( matchedjets.size() > 1 ) {
	cout << "Warning: multiple matches. Skipping." << endl;
	continue;
      }
      if ( matchedjets.size() > 0 ){
	auto sj = matchedjets.at(0);
	smearvtruept->Fill ( tj.pt(), sj.pt() );

#ifdef GROOMING
	PseudoJet sd_truth = sd( tj );
	PseudoJet sd_smear = sd( sj );
	double tzg = sd_truth.structure_of<contrib::SoftDrop>().symmetry();
	double szg = sd_smear.structure_of<contrib::SoftDrop>().symmetry();
	smearvtruezg->Fill ( tzg, szg );
#endif // GROOMING
      }
    }

  }

  cout << " Done." << endl;
  cout << "Processed " << nevents << " events. Of these, " << endl;
  cout << " ---  " << keptevent << " were ACCEPTED in both truth and smeared tree," << endl;
  cout << " ---  " << lostevent << " were LOST due to cuts," << endl;
  cout << " ---  " << rejectevent << " were REJECTED at truth level, and " << endl;
  cout << " ---  " << fakeevent << " were recovered as FAKEs, i.e. they passed the cuts only after smearing." << endl;

  new TCanvas;
  gPad->SetLogy();
  truthpt->SetLineColor( kBlack );
  truthpt->Draw();
  smearpt->SetLineColor( kRed );
  smearpt->Draw("same");
  gPad->SaveAs( outfilebase + "_pt.png");
    
  new TCanvas;
  gPad->SetLogz();
  smearvtruept->Draw("colz");
  gPad->SaveAs( outfilebase + "_responsept.png");

#ifdef GROOMING
  new TCanvas;
  smearvtruezg->Draw("colz");
  gPad->SaveAs( outfilebase + "_responsezg.png");
#endif  
    
  outfile->Write();
  
  return 0;
}
