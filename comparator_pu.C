#define comparator_pu_cxx
#include "comparator_pu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void comparator_pu::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L comparator_pu.C
//      Root > comparator_pu t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


   TFile *outfile = new TFile("joined.root","recreate");
   outfile->cd();
   TTree *outtree = new TTree("SCtree","SCtree");

   // Declaration of leaf types
   Int_t           event_nRecoVertex3;
   Float_t         event_rho3;
   Float_t         GenEnergy3;
   Float_t         GenEta3;
   Float_t         GenPhi3;
   Int_t           pho_isInPhiCrack3;
   Int_t           pho_isInEtaCrack3;
   Int_t           pho_isInEBEEGap3;
   Float_t         phoSC_RawEtCetaCorr3;
   Float_t         phoSC_RawEnergyCetaCorr3;
   Float_t         phoSC_GeomEta3;
   Float_t         phoSC_GeomPhi3;
   Float_t         pho_r93;
   Float_t         pho_Brem3;
   Float_t         pho_SCarea3;

   // Declaration of leaf types
   Int_t           event_nRecoVertex4;
   Float_t         event_rho4;
   Float_t         GenEnergy4;
   Float_t         GenEta4;
   Float_t         GenPhi4;
   Int_t           pho_isInPhiCrack4;
   Int_t           pho_isInEtaCrack4;
   Int_t           pho_isInEBEEGap4;
   Float_t         phoSC_RawEtCetaCorr4;
   Float_t         phoSC_RawEnergyCetaCorr4;
   Float_t         phoSC_GeomEta4;
   Float_t         phoSC_GeomPhi4;
   Float_t         pho_r94;
   Float_t         pho_Brem4;
   Float_t         pho_SCarea4;

   outtree->Branch("event_nRecoVertex_0",&event_nRecoVertex3,"event_nRecoVertex_0/I");
   outtree->Branch("event_rho_0",&event_rho3,"event_rho_0/F");
   outtree->Branch("GenEnergy_0",&GenEnergy3,"GenEnergy_0/F");
   outtree->Branch("GenEta_0",&GenEta3,"GenEta_0/F");
   outtree->Branch("GenPhi_0",&GenPhi3,"GenPhi_0/F");
   outtree->Branch("pho_isInPhiCrack_0",&pho_isInPhiCrack3,"pho_isInPhiCrack_0/I");
   outtree->Branch("pho_isInEtaCrack_0",&pho_isInEtaCrack3,"pho_isInEtaCrack_0/I");
   outtree->Branch("pho_isInEBEEGap_0",&pho_isInEBEEGap3,"pho_isInEBEEGap_0/I");
   outtree->Branch("phoSC_RawEtCetaCorr_0",&phoSC_RawEtCetaCorr3,"phoSC_RawEtCetaCorr_0/F");
   outtree->Branch("phoSC_RawEnergyCetaCorr_0",&phoSC_RawEnergyCetaCorr3,"phoSC_RawEnergyCetaCorr_0/F");
   outtree->Branch("phoSC_GeomEta_0",&phoSC_GeomEta3,"phoSC_GeomEta_0/F");
   outtree->Branch("phoSC_GeomPhi_0",&phoSC_GeomPhi3,"phoSC_GeomPhi_0/F");
   outtree->Branch("pho_r9_0",&pho_r93,"pho_r9_0/F");
   outtree->Branch("pho_Brem_0",&pho_Brem3,"pho_Brem_0/F");
   outtree->Branch("pho_SCarea_0",&pho_SCarea3,"pho_SCarea_0/F");

   outtree->Branch("event_nRecoVertex_40",&event_nRecoVertex4,"event_nRecoVertex_40/I");
   outtree->Branch("event_rho_40",&event_rho4,"event_rho_40/F");
   outtree->Branch("GenEnergy_40",&GenEnergy4,"GenEnergy_40/F");
   outtree->Branch("GenEta_40",&GenEta4,"GenEta_40/F");
   outtree->Branch("GenPhi_40",&GenPhi4,"GenPhi_40/F");
   outtree->Branch("pho_isInPhiCrack_40",&pho_isInPhiCrack4,"pho_isInPhiCrack_40/I");
   outtree->Branch("pho_isInEtaCrack_40",&pho_isInEtaCrack4,"pho_isInEtaCrack_40/I");
   outtree->Branch("pho_isInEBEEGap_40",&pho_isInEBEEGap4,"pho_isInEBEEGap_40/I");
   outtree->Branch("phoSC_RawEtCetaCorr_40",&phoSC_RawEtCetaCorr4,"phoSC_RawEtCetaCorr_40/F");
   outtree->Branch("phoSC_RawEnergyCetaCorr_40",&phoSC_RawEnergyCetaCorr4,"phoSC_RawEnergyCetaCorr_40/F");
   outtree->Branch("phoSC_GeomEta_40",&phoSC_GeomEta4,"phoSC_GeomEta_40/F");
   outtree->Branch("phoSC_GeomPhi_40",&phoSC_GeomPhi4,"phoSC_GeomPhi_40/F");
   outtree->Branch("pho_r9_40",&pho_r94,"pho_r9_40/F");
   outtree->Branch("pho_Brem_40",&pho_Brem4,"pho_Brem_40/F");
   outtree->Branch("pho_SCarea_40",&pho_SCarea4,"pho_SCarea_40/F");


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;



      //      cout << "+++ev" << endl;

      bool retry=false;

      if (nb<0) retry=true;
      if (DeltaR(phoSC_GeomPhi,phoSC_GeomPhi2,phoSC_GeomEta,phoSC_GeomEta2)>0.2) retry=true;
      if (retry) nb = DoRetry();

      bool problem=false;
      if (nb<0) problem=true;
      if (DeltaR(phoSC_GeomPhi,phoSC_GeomPhi2,phoSC_GeomEta,phoSC_GeomEta2)>0.2) problem=true;


      if (problem) {
	//	cout << "---" << endl;
	//	cout << DoRetry(1) << endl;
	//	Printer();
	//	cout << DoRetry(1) << endl;
	//	Printer();
	//	cout << "---" << endl;
	continue;
      }

      if (pho_isInPhiCrack || pho_isInEtaCrack || pho_isInEBEEGap) continue;
      if (pho_isInPhiCrack2 || pho_isInEtaCrack2 || pho_isInEBEEGap2) continue;

      event_nRecoVertex3=event_nRecoVertex;
      event_rho3=event_rho;
      GenEnergy3=GenEnergy;
      GenEta3=GenEta;
      GenPhi3=GenPhi;
      pho_isInPhiCrack3=pho_isInPhiCrack;
      pho_isInEtaCrack3=pho_isInEtaCrack;
      pho_isInEBEEGap3=pho_isInEBEEGap;
      phoSC_RawEtCetaCorr3=phoSC_RawEtCetaCorr;
      phoSC_RawEnergyCetaCorr3=phoSC_RawEnergyCetaCorr;
      phoSC_GeomEta3=phoSC_GeomEta;
      phoSC_GeomPhi3=phoSC_GeomPhi;
      pho_r93=pho_r9;
      pho_Brem3=pho_Brem;
      pho_SCarea3=pho_SCarea;

      event_nRecoVertex4=event_nRecoVertex2;
      event_rho4=event_rho2;
      GenEnergy4=GenEnergy2;
      GenEta4=GenEta2;
      GenPhi4=GenPhi2;
      pho_isInPhiCrack4=pho_isInPhiCrack2;
      pho_isInEtaCrack4=pho_isInEtaCrack2;
      pho_isInEBEEGap4=pho_isInEBEEGap2;
      phoSC_RawEtCetaCorr4=phoSC_RawEtCetaCorr2;
      phoSC_RawEnergyCetaCorr4=phoSC_RawEnergyCetaCorr2;
      phoSC_GeomEta4=phoSC_GeomEta2;
      phoSC_GeomPhi4=phoSC_GeomPhi2;
      pho_r94=pho_r92;
      pho_Brem4=pho_Brem2;
      pho_SCarea4=pho_SCarea2;

      outtree->Fill();


   } // end event loop

   outfile->Write();
   outfile->Close();


};

void comparator_pu::Printer(){
  std::cout << phoSC_RawEnergy << " " << phoSC_RawEnergy2 << " " << phoSC_GeomEta << " " << phoSC_GeomEta2 << " " << phoSC_GeomPhi << " " << phoSC_GeomPhi2 << std::endl;
};

float comparator_pu::DeltaR(double phi1, double phi2, double eta1, double eta2){
  
  double dphi=phi2-phi1;
  if (dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
  if (dphi<-TMath::Pi()) dphi=-2*TMath::Pi()-dphi;
  double dR=sqrt(dphi*dphi+(eta2-eta1)*(eta2-eta1));

  return dR;
};

int comparator_pu::DoRetry(bool debug){
  int newindex = event_SCindex2;
  if (event_SCindex2==0) newindex=1; else newindex=0;
  if (debug) {
    cout << "ls,ev " << event_ls << "," << event_number << endl;
    cout << "Getting " << event_ls << "," << event_number*(newindex*2-1) << " in friend" << endl;
  }
  return Ftree->GetEntryWithIndex(event_ls,event_number*(newindex*2-1));
};
