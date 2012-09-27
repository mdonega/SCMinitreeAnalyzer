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

   Long64_t events = 0;

   TH1F *histo_nrecvtx;
   TH1F *histo_nrecvtx_norew;
   TH1F *histo_rho;

   TH1F *histo_r9[2][n_templates_max+1];
   TH1F *histo_r25[2][n_templates_max+1];
   TH1F *histo_brem[2][n_templates_max+1];
   TH1F *histo_sieie[2][n_templates_max+1];
   TH1F *histo_erecegen[2][n_templates_max+1];
   TH1F *histo_mass[2];

   histo_nrecvtx = new TH1F("histo_nrecvtx","histo_nrecvtx",50,0,50);
   histo_nrecvtx_norew = new TH1F("histo_nrecvtx_norew","histo_nrecvtx_norew",50,0,50);
   histo_rho = new TH1F("histo_rho","histo_rho",100,0,20);

   for (int i=0; i<n_templates_max+1; i++){
     for (int j=0; j<2; j++){
       TString reg = (j==0) ? "EB" : "EE";
       histo_r9[j][i] = new TH1F(Form("h_r9_%s_b%d",reg.Data(),i),Form("h_r9_%s_b%d",reg.Data(),i),130*4,0,1.3);
       histo_r25[j][i] = new TH1F(Form("h_r25_%s_b%d",reg.Data(),i),Form("h_r25_%s_b%d",reg.Data(),i),130*4,0,1.3);
       histo_brem[j][i] = new TH1F(Form("h_brem_%s_b%d",reg.Data(),i),Form("h_brem_%s_b%d",reg.Data(),i),100,0,10);
       histo_sieie[j][i] = new TH1F(Form("h_sieie_%s_b%d",reg.Data(),i),Form("h_sieie_%s_b%d",reg.Data(),i),360,0,0.045);
       histo_erecegen[j][i] = new TH1F(Form("h_erecegen_%s_b%d",reg.Data(),i),Form("h_erecegen_%s_b%d",reg.Data(),i),130*4,0,1.3);
     }
   }

   for (int j=0; j<2; j++){
     TString reg = (j==0) ? "EBEB" : "notEBEB";
     histo_mass[j] = new TH1F(Form("h_mass_%s",reg.Data()),Form("h_mass_%s",reg.Data()),100,81.2,101.2);
   }


   do { // start event loop

     events++;
     if (events%100000==0) std::cout << "Processing event " << events << std::endl;

     //     if (events>=100000) break;

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
       if (!(it->ele_tightID)) pass=0;

       float eta = fabs(it->phoSC_GeomEta);
       if (eta>=1.44 && eta<1.56) pass=0;
       if (eta>=2.5) pass=0;
       //       if (fabs(it->phoSC_GeomEta)>1.44) pass=0;
       //       if (it->pho_r9<0.94) pass=0;
       if (!pass) it=v.erase(it); else it++;
     }

     //     for (unsigned int i=0; i<v.size(); i++) std::cout << v.at(i).event_SCindex ;
     //     std::cout << std::endl;
       
     TLorentzVector el1;
     TLorentzVector el2;
     if (v.size()<2) continue;
     v.resize(2);
     el1.SetPtEtaPhiM(v.at(0).phoSC_RawEtCetaCorr,v.at(0).phoSC_GeomEta,v.at(0).phoSC_GeomPhi,0.51);
     el2.SetPtEtaPhiM(v.at(1).phoSC_RawEtCetaCorr,v.at(1).phoSC_GeomEta,v.at(1).phoSC_GeomPhi,0.51);

     float mass = (el1+el2).M();
     if (fabs(mass-91.2)>10) continue;

     float weight=GetWeight(v.at(0).event_nRecoVertex);

     histo_nrecvtx_norew->Fill(v.at(0).event_nRecoVertex);
     histo_nrecvtx->Fill(v.at(0).event_nRecoVertex,weight);
     histo_rho->Fill(v.at(0).event_rho,weight);

     for (int i=0; i<2; i++){
       int reg = (fabs(v.at(i).phoSC_GeomEta)<1.44) ? 0 : 1;
       int binstofill[2] = {Choose_bin_eta(v.at(i).phoSC_GeomEta,reg),n_templates_max};
       for (int j=0; j<2; j++){
       int bin = binstofill[j];
       histo_r9[reg][bin]->Fill(v.at(i).pho_r9,weight);
       histo_r25[reg][bin]->Fill(v.at(i).pho_e5x5/v.at(i).phoSC_RawEnergy,weight);
       if (v.at(i).pho_Brem!=-1) histo_brem[reg][bin]->Fill(v.at(i).pho_Brem,weight);
       histo_sieie[reg][bin]->Fill(v.at(i).pho_sigmaIetaIeta,weight);
       histo_erecegen[reg][bin]->Fill(v.at(i).phoSC_RawEnergyCetaCorr/v.at(i).GenEnergy,weight);
       }
     }
     
     {
       int Zreg = ((fabs(v.at(0).phoSC_GeomEta)<1.44) && (fabs(v.at(1).phoSC_GeomEta)<1.44)) ? 0 : 1;
       histo_mass[Zreg]->Fill(mass,weight);
     }

     
   } // end event loop
   while (true);


   TFile *f = new TFile("output.root","recreate");
   f->cd();
   for (int j=0; j<2; j++){
     histo_mass[j]->Scale(1.0/histo_mass[j]->Integral());
     histo_mass[j]->Write();
   }
   for (int j=0; j<2; j++){
     for (int i=0; i<n_templates_max+1; i++){
       histo_r9[j][i]->Scale(1.0/histo_r9[j][i]->Integral());
       histo_r25[j][i]->Scale(1.0/histo_r25[j][i]->Integral());
       histo_brem[j][i]->Scale(1.0/histo_brem[j][i]->Integral());
       histo_sieie[j][i]->Scale(1.0/histo_sieie[j][i]->Integral());
       histo_erecegen[j][i]->Scale(1.0/histo_erecegen[j][i]->Integral());
     }
     for (int i=0; i<n_templates_max+1; i++) if (histo_r9[j][i]->GetEntries()>0) histo_r9[j][i]->Write();
     for (int i=0; i<n_templates_max+1; i++) if (histo_r25[j][i]->GetEntries()>0) histo_r25[j][i]->Write();
     for (int i=0; i<n_templates_max+1; i++) if (histo_brem[j][i]->GetEntries()>0) histo_brem[j][i]->Write();
     for (int i=0; i<n_templates_max+1; i++) if (histo_sieie[j][i]->GetEntries()>0) histo_sieie[j][i]->Write();
     for (int i=0; i<n_templates_max+1; i++) if (histo_erecegen[j][i]->GetEntries()>0) histo_erecegen[j][i]->Write();
   }

   histo_nrecvtx_norew->Write();
   histo_nrecvtx->Write();
   histo_rho->Write();


   f->Close();



} // end function
