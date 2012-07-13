//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 10 17:52:53 2012 by ROOT version 5.32/01
// from TTree Tree/PileUp info
// found on file: outfile_pu0.root
//////////////////////////////////////////////////////////

#ifndef comparator_pu_h
#define comparator_pu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class comparator_pu {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TTree *Ftree;

   // Declaration of leaf types
   Int_t           event_number;
   Int_t           event_run;
   Int_t           event_ls;
   Int_t           event_SCindex;
   Int_t           event_nRecoVertex;
   Float_t         event_rho;
   Float_t         event_sigma;
   Int_t           pho_isPrompt;
   Float_t         GenEnergy;
   Float_t         GenEta;
   Float_t         GenPhi;
   Int_t           pho_isEB;
   Int_t           pho_isEE;
   Int_t           pho_isInPhiCrack;
   Int_t           pho_isInEtaCrack;
   Int_t           pho_isInEBEEGap;
   Float_t         phoSC_energy;
   Float_t         phoSC_RawEnergy;
   Float_t         phoSC_et;
   Float_t         phoSC_RawEt;
   Float_t         phoSC_RawEtCetaCorr;
   Float_t         phoSC_RawEnergyCetaCorr;
   Int_t           pho_hasConvTracks;
   Float_t         pho_PosEcalX;
   Float_t         pho_PosEcalY;
   Float_t         pho_PosEcalZ;
   Float_t         pho_GeomEta;
   Float_t         pho_GeomPhi;
   Float_t         phoSC_GeomEta;
   Float_t         phoSC_GeomPhi;
   Float_t         pho_vx;
   Float_t         pho_vy;
   Float_t         pho_vz;
   Float_t         pho_r9;
   Float_t         pho_sigmaIetaIeta;
   Float_t         pho_sigmaEtaEta;
   Float_t         pho_e1x5;
   Float_t         pho_e2x5;
   Float_t         pho_e3x3;
   Float_t         pho_e5x5;
   Float_t         pho_maxEnergyXtal;
   Float_t         pho_etaWidth;
   Float_t         pho_phiWidth;
   Float_t         pho_Brem;
   Float_t         pho_hoe;
   Float_t         pho_SCarea;
   Float_t         pho_SCarea_shoelace;
   Int_t           pho_SCnbc;
   Int_t           pho_SCnxtals;
   Float_t         pho_BCenergy[100];   //[pho_SCnbc]
   Int_t           pho_BCnXtals[100];   //[pho_SCnbc]
   Int_t           pho_BCisSeed[100];   //[pho_SCnbc]

   // List of branches
   TBranch        *b_event_number;   //!
   TBranch        *b_event_run;   //!
   TBranch        *b_event_ls;   //!
   TBranch        *b_event_SCindex;   //!
   TBranch        *b_event_nRecoVertex;   //!
   TBranch        *b_event_rho;   //!
   TBranch        *b_event_sigma;   //!
   TBranch        *b_pho_isPrompt;   //!
   TBranch        *b_GenEnergy;   //!
   TBranch        *b_GenEta;   //!
   TBranch        *b_GenPhi;   //!
   TBranch        *b_pho_isEB;   //!
   TBranch        *b_pho_isEE;   //!
   TBranch        *b_pho_isInPhiCrack;   //!
   TBranch        *b_pho_isInEtaCrack;   //!
   TBranch        *b_pho_isInEBEEGap;   //!
   TBranch        *b_phoSC_energy;   //!
   TBranch        *b_phoSC_RawEnergy;   //!
   TBranch        *b_phoSC_et;   //!
   TBranch        *b_phoSC_RawEt;   //!
   TBranch        *b_phoSC_RawEtCetaCorr;   //!
   TBranch        *b_phoSC_RawEnergyCetaCorr;   //!
   TBranch        *b_pho_hasConvTracks;   //!
   TBranch        *b_pho_PosEcalX;   //!
   TBranch        *b_pho_PosEcalY;   //!
   TBranch        *b_pho_PosEcalZ;   //!
   TBranch        *b_pho_GeomEta;   //!
   TBranch        *b_pho_GeomPhi;   //!
   TBranch        *b_phoSC_GeomEta;   //!
   TBranch        *b_phoSC_GeomPhi;   //!
   TBranch        *b_pho_vx;   //!
   TBranch        *b_pho_vy;   //!
   TBranch        *b_pho_vz;   //!
   TBranch        *b_pho_r9;   //!
   TBranch        *b_pho_sigmaIetaIeta;   //!
   TBranch        *b_pho_sigmaEtaEta;   //!
   TBranch        *b_pho_e1x5;   //!
   TBranch        *b_pho_e2x5;   //!
   TBranch        *b_pho_e3x3;   //!
   TBranch        *b_pho_e5x5;   //!
   TBranch        *b_pho_maxEnergyXtal;   //!
   TBranch        *b_pho_etaWidth;   //!
   TBranch        *b_pho_phiWidth;   //!
   TBranch        *b_pho_Brem;   //!
   TBranch        *b_pho_hoe;   //!
   TBranch        *b_pho_SCarea;   //!
   TBranch        *b_pho_SCarea_shoelace;   //!
   TBranch        *b_pho_SCnbc;   //!
   TBranch        *b_pho_SCnxtals;   //!
   TBranch        *b_pho_BCenergy;   //!
   TBranch        *b_pho_BCnXtals;   //!
   TBranch        *b_pho_BCisSeed;   //!

   // Declaration of leaf types
   Int_t           event_number2;
   Int_t           event_run2;
   Int_t           event_ls2;
   Int_t           event_SCindex2;
   Int_t           event_nRecoVertex2;
   Float_t         event_rho2;
   Float_t         event_sigma2;
   Int_t           pho_isPrompt2;
   Float_t         GenEnergy2;
   Float_t         GenEta2;
   Float_t         GenPhi2;
   Int_t           pho_isEB2;
   Int_t           pho_isEE2;
   Int_t           pho_isInPhiCrack2;
   Int_t           pho_isInEtaCrack2;
   Int_t           pho_isInEBEEGap2;
   Float_t         phoSC_energy2;
   Float_t         phoSC_RawEnergy2;
   Float_t         phoSC_et2;
   Float_t         phoSC_RawEt2;
   Float_t         phoSC_RawEtCetaCorr2;
   Float_t         phoSC_RawEnergyCetaCorr2;
   Int_t           pho_hasConvTracks2;
   Float_t         pho_PosEcalX2;
   Float_t         pho_PosEcalY2;
   Float_t         pho_PosEcalZ2;
   Float_t         pho_GeomEta2;
   Float_t         pho_GeomPhi2;
   Float_t         phoSC_GeomEta2;
   Float_t         phoSC_GeomPhi2;
   Float_t         pho_vx2;
   Float_t         pho_vy2;
   Float_t         pho_vz2;
   Float_t         pho_r92;
   Float_t         pho_sigmaIetaIeta2;
   Float_t         pho_sigmaEtaEta2;
   Float_t         pho_e1x52;
   Float_t         pho_e2x52;
   Float_t         pho_e3x32;
   Float_t         pho_e5x52;
   Float_t         pho_maxEnergyXtal2;
   Float_t         pho_etaWidth2;
   Float_t         pho_phiWidth2;
   Float_t         pho_Brem2;
   Float_t         pho_hoe2;
   Float_t         pho_SCarea2;
   Float_t         pho_SCarea_shoelace2;
   Int_t           pho_SCnbc2;
   Int_t           pho_SCnxtals2;
   Float_t         pho_BCenergy2[100];   //[pho_SCnbc]
   Int_t           pho_BCnXtals2[100];   //[pho_SCnbc]
   Int_t           pho_BCisSeed2[100];   //[pho_SCnbc]

   // List of branches
   TBranch        *b_event_number2;   //!
   TBranch        *b_event_run2;   //!
   TBranch        *b_event_ls2;   //!
   TBranch        *b_event_SCindex2;   //!
   TBranch        *b_event_nRecoVertex2;   //!
   TBranch        *b_event_rho2;   //!
   TBranch        *b_event_sigma2;   //!
   TBranch        *b_pho_isPrompt2;   //!
   TBranch        *b_GenEnergy2;   //!
   TBranch        *b_GenEta2;   //!
   TBranch        *b_GenPhi2;   //!
   TBranch        *b_pho_isEB2;   //!
   TBranch        *b_pho_isEE2;   //!
   TBranch        *b_pho_isInPhiCrack2;   //!
   TBranch        *b_pho_isInEtaCrack2;   //!
   TBranch        *b_pho_isInEBEEGap2;   //!
   TBranch        *b_phoSC_energy2;   //!
   TBranch        *b_phoSC_RawEnergy2;   //!
   TBranch        *b_phoSC_et2;   //!
   TBranch        *b_phoSC_RawEt2;   //!
   TBranch        *b_phoSC_RawEtCetaCorr2;   //!
   TBranch        *b_phoSC_RawEnergyCetaCorr2;   //!
   TBranch        *b_pho_hasConvTracks2;   //!
   TBranch        *b_pho_PosEcalX2;   //!
   TBranch        *b_pho_PosEcalY2;   //!
   TBranch        *b_pho_PosEcalZ2;   //!
   TBranch        *b_pho_GeomEta2;   //!
   TBranch        *b_pho_GeomPhi2;   //!
   TBranch        *b_phoSC_GeomEta2;   //!
   TBranch        *b_phoSC_GeomPhi2;   //!
   TBranch        *b_pho_vx2;   //!
   TBranch        *b_pho_vy2;   //!
   TBranch        *b_pho_vz2;   //!
   TBranch        *b_pho_r92;   //!
   TBranch        *b_pho_sigmaIetaIeta2;   //!
   TBranch        *b_pho_sigmaEtaEta2;   //!
   TBranch        *b_pho_e1x52;   //!
   TBranch        *b_pho_e2x52;   //!
   TBranch        *b_pho_e3x32;   //!
   TBranch        *b_pho_e5x52;   //!
   TBranch        *b_pho_maxEnergyXtal2;   //!
   TBranch        *b_pho_etaWidth2;   //!
   TBranch        *b_pho_phiWidth2;   //!
   TBranch        *b_pho_Brem2;   //!
   TBranch        *b_pho_hoe2;   //!
   TBranch        *b_pho_SCarea2;   //!
   TBranch        *b_pho_SCarea_shoelace2;   //!
   TBranch        *b_pho_SCnbc2;   //!
   TBranch        *b_pho_SCnxtals2;   //!
   TBranch        *b_pho_BCenergy2;   //!
   TBranch        *b_pho_BCnXtals2;   //!
   TBranch        *b_pho_BCisSeed2;   //!


   comparator_pu(TTree *tree=0);
   virtual ~comparator_pu();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void Printer();
   float DeltaR(double,double,double,double);
   int DoRetry(bool debug=false);
};

#endif

#ifdef comparator_pu_cxx
comparator_pu::comparator_pu(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outfile_pu0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outfile_pu0.root");
      }
      f->GetObject("Tree",tree);

      TFile *f2 = new TFile("outfile_pu40.root");
      f2->GetObject("Tree",Ftree);

   }
   Init(tree);
}

comparator_pu::~comparator_pu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t comparator_pu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t comparator_pu::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void comparator_pu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   fChain->SetBranchAddress("event_run", &event_run, &b_event_run);
   fChain->SetBranchAddress("event_ls", &event_ls, &b_event_ls);
   fChain->SetBranchAddress("event_SCindex", &event_SCindex, &b_event_SCindex);
   fChain->SetBranchAddress("event_nRecoVertex", &event_nRecoVertex, &b_event_nRecoVertex);
   fChain->SetBranchAddress("event_rho", &event_rho, &b_event_rho);
   fChain->SetBranchAddress("event_sigma", &event_sigma, &b_event_sigma);
   fChain->SetBranchAddress("pho_isPrompt", &pho_isPrompt, &b_pho_isPrompt);
   fChain->SetBranchAddress("GenEnergy", &GenEnergy, &b_GenEnergy);
   fChain->SetBranchAddress("GenEta", &GenEta, &b_GenEta);
   fChain->SetBranchAddress("GenPhi", &GenPhi, &b_GenPhi);
   fChain->SetBranchAddress("pho_isEB", &pho_isEB, &b_pho_isEB);
   fChain->SetBranchAddress("pho_isEE", &pho_isEE, &b_pho_isEE);
   fChain->SetBranchAddress("pho_isInPhiCrack", &pho_isInPhiCrack, &b_pho_isInPhiCrack);
   fChain->SetBranchAddress("pho_isInEtaCrack", &pho_isInEtaCrack, &b_pho_isInEtaCrack);
   fChain->SetBranchAddress("pho_isInEBEEGap", &pho_isInEBEEGap, &b_pho_isInEBEEGap);
   fChain->SetBranchAddress("phoSC_energy", &phoSC_energy, &b_phoSC_energy);
   fChain->SetBranchAddress("phoSC_RawEnergy", &phoSC_RawEnergy, &b_phoSC_RawEnergy);
   fChain->SetBranchAddress("phoSC_et", &phoSC_et, &b_phoSC_et);
   fChain->SetBranchAddress("phoSC_RawEt", &phoSC_RawEt, &b_phoSC_RawEt);
   fChain->SetBranchAddress("phoSC_RawEtCetaCorr", &phoSC_RawEtCetaCorr, &b_phoSC_RawEtCetaCorr);
   fChain->SetBranchAddress("phoSC_RawEnergyCetaCorr", &phoSC_RawEnergyCetaCorr, &b_phoSC_RawEnergyCetaCorr);
   fChain->SetBranchAddress("pho_hasConvTracks", &pho_hasConvTracks, &b_pho_hasConvTracks);
   fChain->SetBranchAddress("pho_PosEcalX", &pho_PosEcalX, &b_pho_PosEcalX);
   fChain->SetBranchAddress("pho_PosEcalY", &pho_PosEcalY, &b_pho_PosEcalY);
   fChain->SetBranchAddress("pho_PosEcalZ", &pho_PosEcalZ, &b_pho_PosEcalZ);
   fChain->SetBranchAddress("pho_GeomEta", &pho_GeomEta, &b_pho_GeomEta);
   fChain->SetBranchAddress("pho_GeomPhi", &pho_GeomPhi, &b_pho_GeomPhi);
   fChain->SetBranchAddress("phoSC_GeomEta", &phoSC_GeomEta, &b_phoSC_GeomEta);
   fChain->SetBranchAddress("phoSC_GeomPhi", &phoSC_GeomPhi, &b_phoSC_GeomPhi);
   fChain->SetBranchAddress("pho_vx", &pho_vx, &b_pho_vx);
   fChain->SetBranchAddress("pho_vy", &pho_vy, &b_pho_vy);
   fChain->SetBranchAddress("pho_vz", &pho_vz, &b_pho_vz);
   fChain->SetBranchAddress("pho_r9", &pho_r9, &b_pho_r9);
   fChain->SetBranchAddress("pho_sigmaIetaIeta", &pho_sigmaIetaIeta, &b_pho_sigmaIetaIeta);
   fChain->SetBranchAddress("pho_sigmaEtaEta", &pho_sigmaEtaEta, &b_pho_sigmaEtaEta);
   fChain->SetBranchAddress("pho_e1x5", &pho_e1x5, &b_pho_e1x5);
   fChain->SetBranchAddress("pho_e2x5", &pho_e2x5, &b_pho_e2x5);
   fChain->SetBranchAddress("pho_e3x3", &pho_e3x3, &b_pho_e3x3);
   fChain->SetBranchAddress("pho_e5x5", &pho_e5x5, &b_pho_e5x5);
   fChain->SetBranchAddress("pho_maxEnergyXtal", &pho_maxEnergyXtal, &b_pho_maxEnergyXtal);
   fChain->SetBranchAddress("pho_etaWidth", &pho_etaWidth, &b_pho_etaWidth);
   fChain->SetBranchAddress("pho_phiWidth", &pho_phiWidth, &b_pho_phiWidth);
   fChain->SetBranchAddress("pho_Brem", &pho_Brem, &b_pho_Brem);
   fChain->SetBranchAddress("pho_hoe", &pho_hoe, &b_pho_hoe);
   fChain->SetBranchAddress("pho_SCarea", &pho_SCarea, &b_pho_SCarea);
   fChain->SetBranchAddress("pho_SCarea_shoelace", &pho_SCarea_shoelace, &b_pho_SCarea_shoelace);
   fChain->SetBranchAddress("pho_SCnbc", &pho_SCnbc, &b_pho_SCnbc);
   fChain->SetBranchAddress("pho_SCnxtals", &pho_SCnxtals, &b_pho_SCnxtals);
   fChain->SetBranchAddress("pho_BCenergy", pho_BCenergy, &b_pho_BCenergy);
   fChain->SetBranchAddress("pho_BCnXtals", pho_BCnXtals, &b_pho_BCnXtals);
   fChain->SetBranchAddress("pho_BCisSeed", pho_BCisSeed, &b_pho_BCisSeed);

   Notify();

   Ftree->BuildIndex("event_ls","event_number*(event_SCindex*2-1)");

   fChain->AddFriend(Ftree,"FTree");

   fChain->SetBranchAddress("FTree.event_number", &event_number2, &b_event_number2);
   fChain->SetBranchAddress("FTree.event_run", &event_run2, &b_event_run2);
   fChain->SetBranchAddress("FTree.event_ls", &event_ls2, &b_event_ls2);
   fChain->SetBranchAddress("FTree.event_SCindex", &event_SCindex2, &b_event_SCindex2);
   fChain->SetBranchAddress("FTree.event_nRecoVertex", &event_nRecoVertex2, &b_event_nRecoVertex2);
   fChain->SetBranchAddress("FTree.event_rho", &event_rho2, &b_event_rho2);
   fChain->SetBranchAddress("FTree.event_sigma", &event_sigma2, &b_event_sigma2);
   fChain->SetBranchAddress("FTree.pho_isPrompt", &pho_isPrompt2, &b_pho_isPrompt2);
   fChain->SetBranchAddress("FTree.GenEnergy", &GenEnergy2, &b_GenEnergy2);
   fChain->SetBranchAddress("FTree.GenEta", &GenEta2, &b_GenEta2);
   fChain->SetBranchAddress("FTree.GenPhi", &GenPhi2, &b_GenPhi2);
   fChain->SetBranchAddress("FTree.pho_isEB", &pho_isEB2, &b_pho_isEB2);
   fChain->SetBranchAddress("FTree.pho_isEE", &pho_isEE2, &b_pho_isEE2);
   fChain->SetBranchAddress("FTree.pho_isInPhiCrack", &pho_isInPhiCrack2, &b_pho_isInPhiCrack2);
   fChain->SetBranchAddress("FTree.pho_isInEtaCrack", &pho_isInEtaCrack2, &b_pho_isInEtaCrack2);
   fChain->SetBranchAddress("FTree.pho_isInEBEEGap", &pho_isInEBEEGap2, &b_pho_isInEBEEGap2);
   fChain->SetBranchAddress("FTree.phoSC_energy", &phoSC_energy2, &b_phoSC_energy2);
   fChain->SetBranchAddress("FTree.phoSC_RawEnergy", &phoSC_RawEnergy2, &b_phoSC_RawEnergy2);
   fChain->SetBranchAddress("FTree.phoSC_et", &phoSC_et2, &b_phoSC_et2);
   fChain->SetBranchAddress("FTree.phoSC_RawEt", &phoSC_RawEt2, &b_phoSC_RawEt2);
   fChain->SetBranchAddress("FTree.phoSC_RawEtCetaCorr", &phoSC_RawEtCetaCorr2, &b_phoSC_RawEtCetaCorr2);
   fChain->SetBranchAddress("FTree.phoSC_RawEnergyCetaCorr", &phoSC_RawEnergyCetaCorr2, &b_phoSC_RawEnergyCetaCorr2);
   fChain->SetBranchAddress("FTree.pho_hasConvTracks", &pho_hasConvTracks2, &b_pho_hasConvTracks2);
   fChain->SetBranchAddress("FTree.pho_PosEcalX", &pho_PosEcalX2, &b_pho_PosEcalX2);
   fChain->SetBranchAddress("FTree.pho_PosEcalY", &pho_PosEcalY2, &b_pho_PosEcalY2);
   fChain->SetBranchAddress("FTree.pho_PosEcalZ", &pho_PosEcalZ2, &b_pho_PosEcalZ2);
   fChain->SetBranchAddress("FTree.pho_GeomEta", &pho_GeomEta2, &b_pho_GeomEta2);
   fChain->SetBranchAddress("FTree.pho_GeomPhi", &pho_GeomPhi2, &b_pho_GeomPhi2);
   fChain->SetBranchAddress("FTree.phoSC_GeomEta", &phoSC_GeomEta2, &b_phoSC_GeomEta2);
   fChain->SetBranchAddress("FTree.phoSC_GeomPhi", &phoSC_GeomPhi2, &b_phoSC_GeomPhi2);
   fChain->SetBranchAddress("FTree.pho_vx", &pho_vx2, &b_pho_vx2);
   fChain->SetBranchAddress("FTree.pho_vy", &pho_vy2, &b_pho_vy2);
   fChain->SetBranchAddress("FTree.pho_vz", &pho_vz2, &b_pho_vz2);
   fChain->SetBranchAddress("FTree.pho_r9", &pho_r92, &b_pho_r92);
   fChain->SetBranchAddress("FTree.pho_sigmaIetaIeta", &pho_sigmaIetaIeta2, &b_pho_sigmaIetaIeta2);
   fChain->SetBranchAddress("FTree.pho_sigmaEtaEta", &pho_sigmaEtaEta2, &b_pho_sigmaEtaEta2);
   fChain->SetBranchAddress("FTree.pho_e1x5", &pho_e1x52, &b_pho_e1x52);
   fChain->SetBranchAddress("FTree.pho_e2x5", &pho_e2x52, &b_pho_e2x52);
   fChain->SetBranchAddress("FTree.pho_e3x3", &pho_e3x32, &b_pho_e3x32);
   fChain->SetBranchAddress("FTree.pho_e5x5", &pho_e5x52, &b_pho_e5x52);
   fChain->SetBranchAddress("FTree.pho_maxEnergyXtal", &pho_maxEnergyXtal2, &b_pho_maxEnergyXtal2);
   fChain->SetBranchAddress("FTree.pho_etaWidth", &pho_etaWidth2, &b_pho_etaWidth2);
   fChain->SetBranchAddress("FTree.pho_phiWidth", &pho_phiWidth2, &b_pho_phiWidth2);
   fChain->SetBranchAddress("FTree.pho_Brem", &pho_Brem2, &b_pho_Brem2);
   fChain->SetBranchAddress("FTree.pho_hoe", &pho_hoe2, &b_pho_hoe2);
   fChain->SetBranchAddress("FTree.pho_SCarea", &pho_SCarea2, &b_pho_SCarea2);
   fChain->SetBranchAddress("FTree.pho_SCarea_shoelace", &pho_SCarea_shoelace2, &b_pho_SCarea_shoelace2);
   fChain->SetBranchAddress("FTree.pho_SCnbc", &pho_SCnbc2, &b_pho_SCnbc2);
   fChain->SetBranchAddress("FTree.pho_SCnxtals", &pho_SCnxtals2, &b_pho_SCnxtals2);
   fChain->SetBranchAddress("FTree.pho_BCenergy", pho_BCenergy2, &b_pho_BCenergy2);
   fChain->SetBranchAddress("FTree.pho_BCnXtals", pho_BCnXtals2, &b_pho_BCnXtals2);
   fChain->SetBranchAddress("FTree.pho_BCisSeed", pho_BCisSeed2, &b_pho_BCisSeed2);





}

Bool_t comparator_pu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void comparator_pu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t comparator_pu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef comparator_pu_cxx
