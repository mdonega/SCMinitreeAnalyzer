#define Looper_cxx
#include "Looper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void Looper::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Looper.C
//      Root > Looper t
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
   Long64_t nbytes = 0, nb = 0;

   Long64_t jentry = 0;

   TH1F *histo = new TH1F("h","h",200,50,150);

   do {

     std::vector<ParticleObject> v;

     {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       if (jentry==nentries) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       {ParticleObject ob; ob=obj_; v.push_back(ob);}
     }

     do {
       jentry++;
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       if (jentry==nentries) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       if (obj_.event_run != v.back().event_run) break;
       if (obj_.event_ls != v.back().event_ls) break;
       if (obj_.event_number != v.back().event_number) break;
       {ParticleObject ob; ob=obj_; v.push_back(ob);}
     }
     while (true);


     for (vector<ParticleObject>::iterator it = v.begin(); it != v.end(); ){
       bool pass=1;
       //       if (!(it->ele_tightID)) pass=0;
       //       if (fabs(it->phoSC_GeomEta)>1.0) pass=0;
       //       if (it->pho_r9<0.94) pass=0;
       if (!pass) it=v.erase(it); else it++;
     }

     //     for (unsigned int i=0; i<v.size(); i++) std::cout << v.at(i).event_SCindex ;
     //     std::cout << std::endl;
       
     TLorentzVector el1;
     TLorentzVector el2;
     if (v.size()<2) continue;
     el1.SetPtEtaPhiM(v.at(0).phoSC_RawEtCetaCorr,v.at(0).phoSC_GeomEta,v.at(0).phoSC_GeomPhi,0.51);
     el2.SetPtEtaPhiM(v.at(1).phoSC_RawEtCetaCorr,v.at(1).phoSC_GeomEta,v.at(1).phoSC_GeomPhi,0.51);
     histo->Fill((el1+el2).M());
     
   } // end event loop
   while (true);

   histo->Draw();

} // end function
