// -*- C++ -*-
//
// Package:    PhotonPUAnalyzer
// Class:      PhotonPUAnalyzer
// 
/**\class PhotonPUAnalyzer PhotonPUAnalyzer.cc Analysis/PhotonPUAnalyzer/src/PhotonPUAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nicolas Pierre Chanon,32 2-C13,+41227674539,
//         Created:  Tue Mar  8 12:23:45 CET 2011
// $Id: PhotonPUAnalyzer.cc,v 1.4 2012/09/21 15:36:10 peruzzi Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "Geometry/CaloTopology/interface/CaloTopology.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
//#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//ParticleFlowReco/interface/PFCluster.h
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"
#include "CommonTools/ParticleFlow/plugins/PFPileUp.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "applyScCorrections_PHOTON.C"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"


//#include "fLocalCorr.h"

#include "CondFormats/EcalObjects/interface/EcalClusterLocalContCorrParameters.h"
#include "CondFormats/DataRecord/interface/EcalClusterLocalContCorrParametersRcd.h"
#include "CondFormats/EcalObjects/interface/EcalClusterCrackCorrParameters.h"
#include "CondFormats/DataRecord/interface/EcalClusterCrackCorrParametersRcd.h"

//#include "RecoEcal/EgammaCoreTools/plugins/EcalClusterEnergyCorrectionObjectSpecificBaseClass.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <Math/VectorUtil.h>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1I.h"
//#include "TBranch.h"

//
// class declaration
//


using namespace ROOT::Math::VectorUtil;
typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;
typedef math::XYZVector Vector;

using namespace std;

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;

class PhotonPUAnalyzer : public edm::EDAnalyzer {
   public:
  explicit PhotonPUAnalyzer(const edm::ParameterSet&);
      ~PhotonPUAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


  bool isInPhiCracks(double phi, double eta);
  bool isInEtaCracks(double eta);

  double f5x5(double iEta);

  typedef std::pair<float,int> OrderPair;
  struct IndexByPt {
    const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
      return j1.first > j2.first;
    }
  };

  IndexByPt indexComparator;

  TH1F *GenParticles_Phi;
  //  TH1F *RecoParticles_Phi;
  //  TH1I *n_RecoParticles;

  edm::InputTag srcRho_;
  edm::InputTag srcSigma_;  

  edm::InputTag               conversionsInputTag_;
  edm::InputTag               beamSpotInputTag_;
  edm::InputTag               rhoIsoInputTag;
  std::vector<edm::InputTag>  isoValInputTags_;
  
  edm::InputTag vertexProducer_;

  edm::InputTag photonsProducer_;

  edm::InputTag SCProducer_;

  edm::InputTag genParticlesProducer_;
  bool doMC_;
  bool doElectrons_;

  edm::InputTag electronsProducer_;
  edm::InputTag reducedBarrelEcalRecHitCollection_;
  edm::InputTag reducedEndcapEcalRecHitCollection_;
  std::string OutputFile_;



  EcalClusterFunctionBaseClass *f;

  TFile* fOutput;
  TTree* myTree_;

  //event info
  int event_number;
  int event_run;
  int event_ls;
  int event_SCindex;
  int event_nRecoVertex;
  int event_genPU;
  int event_genPUtrue;
  
  int coll_size;

  float event_rho;
  float event_sigma;

  float pho_et;
  float pho_energy;
  float pho_px;
  float pho_py;
  float pho_pz;
  float pho_eta;
  float pho_phi;

  float GenEnergy;
  float GenEta;
  float GenPhi;

  int ele_looseID;
  int ele_mediumID;
  int ele_tightID;

  float track_pt;
  float track_eta;
  float track_phi;

  int pixelseed;

  float phoSC_GeomEta;
  float phoSC_GeomPhi;

  float phoSC_et;
  float phoSC_RawEt;
  float phoSC_RawEtCetaCorr;
  float phoSC_RawEnergyCetaCorr;
  float phoSC_correctedenergy;

  //photon tags
  int pho_isEB;
  int pho_isEE;

  //int pho_isInCrack;
  int pho_isInPhiCrack;
  int pho_isInEtaCrack;
  int pho_isInEBEEGap;

  //photon Id
  float phoSC_r9;
  float pho_objectr9;
  float pho_sigmaIetaIeta;
  float pho_e1x5 ;
  float pho_e2x5 ;
  float pho_e3x3 ;
  float pho_e5x5 ;
  float pho_etaWidth;
  float pho_phiWidth;
  float pho_Brem;
  float phoSC_energy;
  float phoSC_RawEnergy;
  float pho_hoe;

  int pho_hasPixelSeed;
  

  EcalClusterLazyTools* lazyTools;  

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PhotonPUAnalyzer::PhotonPUAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  srcRho_ = iConfig.getParameter<edm::InputTag>("srcRho");
  srcSigma_ = iConfig.getParameter<edm::InputTag>("srcSigma");
  vertexProducer_ = iConfig.getParameter<edm::InputTag>("vertexProducer");
  photonsProducer_ = iConfig.getParameter<edm::InputTag>("photonsProducer");
  genParticlesProducer_ = iConfig.getParameter<edm::InputTag>("genParticlesProducer");
  doElectrons_ = iConfig.getParameter<bool>("doElectrons");
  doMC_ = iConfig.getParameter<bool>("doMC");
  electronsProducer_ = iConfig.getParameter<edm::InputTag>("electronsProducer");
  reducedBarrelEcalRecHitCollection_ = iConfig.getParameter<edm::InputTag>("reducedBarrelEcalRecHitCollection");
  reducedEndcapEcalRecHitCollection_ = iConfig.getParameter<edm::InputTag>("reducedEndcapEcalRecHitCollection");
  OutputFile_ = iConfig.getParameter<std::string>("OutputFile");

  conversionsInputTag_    = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
  beamSpotInputTag_       = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
  rhoIsoInputTag          = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
  isoValInputTags_        = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

  GenParticles_Phi = new TH1F("GenParticles_Phi","GenParticles_Phi",128,-3.2,+3.2);
  //  RecoParticles_Phi = new TH1F("RecoParticles_Phi","RecoParticles_Phi",128,-3.2,+3.2);
  //  n_RecoParticles = new TH1I("n_RecoParticles","n_RecoParticles",5,0,5);

  //f = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection", iConfig);
  //  f = EcalClusterFunctionFactory::get()->create("EcalClusterEnergyCorrectionObjectSpecific", iConfig);
  //f = EcalClusterFunctionFactory::get()->create("EcalClusterEnergyCorrection", iConfig);

  cout << "Creating "<<OutputFile_.c_str()<<endl;
  fOutput = new TFile(OutputFile_.c_str(),"RECREATE");
  myTree_ = new TTree("Tree","Minitree");

  myTree_->Branch("event_number",&event_number,"event_number/I");
  myTree_->Branch("event_run",&event_run,"event_run/I");
  myTree_->Branch("event_ls",&event_ls,"event_ls/I");
  myTree_->Branch("event_SCindex",&event_SCindex,"event_SCindex/I");
  myTree_->Branch("event_coll_size",&coll_size,"event_coll_size/I");
  myTree_->Branch("event_nRecoVertex",&event_nRecoVertex,"event_nRecoVertex/I");
  myTree_->Branch("event_rho",&event_rho,"event_rho/F");
  myTree_->Branch("event_sigma",&event_sigma,"event_sigma/F");

  myTree_->Branch("GenEnergy",&GenEnergy,"GenEnergy/F");
  myTree_->Branch("GenEta",&GenEta,"GenEta/F");
  myTree_->Branch("GenPhi",&GenPhi,"GenPhi/F");

  myTree_->Branch("event_genPU",&event_genPU,"event_genPU/I");
  myTree_->Branch("event_genPUtrue",&event_genPUtrue,"event_genPUtrue/I");

  myTree_->Branch("ele_looseID",&ele_looseID,"ele_looseID/I");
  myTree_->Branch("ele_mediumID",&ele_mediumID,"ele_mediumID/I");
  myTree_->Branch("ele_tightID",&ele_tightID,"ele_tightID/I");

  myTree_->Branch("track_pt",&track_pt,"track_pt/F");
  myTree_->Branch("track_eta",&track_eta,"track_eta/F");
  myTree_->Branch("track_phi",&track_phi,"track_phi/F");

  myTree_->Branch("pixelseed",&pixelseed,"pixelseed/I");

  //  myTree_->Branch("pho_isEB",&pho_isEB,"pho_isEB/I");
  //  myTree_->Branch("pho_isEE",&pho_isEE,"pho_isEE/I");

  myTree_->Branch("pho_isInPhiCrack",&pho_isInPhiCrack,"pho_isInPhiCrack/I");
  myTree_->Branch("pho_isInEtaCrack",&pho_isInEtaCrack,"pho_isInEtaCrack/I");
  myTree_->Branch("pho_isInEBEEGap",&pho_isInEBEEGap,"pho_isInEBEEGap/I");

  myTree_->Branch("pho_energy",&pho_energy,"pho_energy/F");
  myTree_->Branch("pho_et",&pho_et,"pho_et/F");

  myTree_->Branch("phoSC_energy",&phoSC_energy,"phoSC_energy/F");
  myTree_->Branch("phoSC_RawEnergy",&phoSC_RawEnergy,"phoSC_RawEnergy/F");
  myTree_->Branch("phoSC_et",&phoSC_et,"phoSC_et/F");
  myTree_->Branch("phoSC_RawEt",&phoSC_RawEt,"phoSC_RawEt/F");
  myTree_->Branch("phoSC_RawEtCetaCorr",&phoSC_RawEtCetaCorr,"phoSC_RawEtCetaCorr/F");
  myTree_->Branch("phoSC_RawEnergyCetaCorr",&phoSC_RawEnergyCetaCorr,"phoSC_RawEnergyCetaCorr/F");
  myTree_->Branch("phoSC_correctedenergy",&phoSC_correctedenergy,"phoSC_correctedenergy/F");

  myTree_->Branch("phoSC_GeomEta", &phoSC_GeomEta, "phoSC_GeomEta/F");
  myTree_->Branch("phoSC_GeomPhi", &phoSC_GeomPhi, "phoSC_GeomPhi/F");

  myTree_->Branch("phoSC_r9",&phoSC_r9,"phoSC_r9/F");
  myTree_->Branch("pho_objectr9",&pho_objectr9,"pho_objectr9/F");
  myTree_->Branch("pho_sigmaIetaIeta",&pho_sigmaIetaIeta,"pho_sigmaIetaIeta/F");
  //  myTree_->Branch("pho_e1x5",&pho_e1x5,"pho_e1x5/F");
  //  myTree_->Branch("pho_e2x5",&pho_e2x5,"pho_e2x5/F");
  myTree_->Branch("pho_e3x3",&pho_e3x3,"pho_e3x3/F");
  myTree_->Branch("pho_e5x5",&pho_e5x5,"pho_e5x5/F");
  myTree_->Branch("pho_etaWidth",&pho_etaWidth,"pho_etaWidth/F");
  myTree_->Branch("pho_phiWidth",&pho_phiWidth,"pho_phiWidth/F");
  myTree_->Branch("pho_Brem",&pho_Brem,"pho_Brem/F");
  //  myTree_->Branch("pho_hoe",&pho_hoe,"pho_hoe/F");


}


PhotonPUAnalyzer::~PhotonPUAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  //  fOutput->WriteObject(myTree_,myTree_->GetName());
  //  fOutput->WriteObject(GenParticles_Phi,GenParticles_Phi->GetName());
  //  fOutput->WriteObject(RecoParticles_Phi,RecoParticles_Phi->GetName());
  //  fOutput->WriteObject(n_RecoParticles,n_RecoParticles->GetName());
  fOutput->Write();
  fOutput->Close();
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PhotonPUAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   
   //Rho, sigma
   edm::Handle<double> rhoHandle;
   iEvent.getByLabel(srcRho_,rhoHandle);
   double rho = *rhoHandle;
   //cout << "rho=" << rho << std::endl;
   edm::Handle<double> sigmaHandle;
   iEvent.getByLabel(srcSigma_,sigmaHandle);
   double sigma = *sigmaHandle;
   //cout << "sigma="<<sigma<<endl;
   edm::Handle<reco::VertexCollection> vertexHandle;
   iEvent.getByLabel(vertexProducer_, vertexHandle);

   // conversions
   edm::Handle<reco::ConversionCollection> conversions_h;
   iEvent.getByLabel(conversionsInputTag_, conversions_h);

   // beam spot
   edm::Handle<reco::BeamSpot> beamspot_h;
   iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
   const reco::BeamSpot &beamSpot = *(beamspot_h.product());

   // iso deposits
   IsoDepositVals isoVals(isoValInputTags_.size());
   for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
     iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
   }
   
   // rho for isolation
   edm::Handle<double> rhoIso_h;
   iEvent.getByLabel(rhoIsoInputTag, rhoIso_h);
   double rhoIso = *(rhoIso_h.product());

   //Photon collection
   edm::Handle<reco::PhotonCollection> photonHandle;
   iEvent.getByLabel(photonsProducer_,photonHandle);
   //cout << "nPhotons=" << photonHandle->size() << std::endl;
   
   //   edm::Handle<reco::SuperClusterCollection> SCHandle;
   //   iEvent.getByLabel(SCProducer_,SCHandle);

   //MC truth
   edm::Handle <reco::GenParticleCollection> genParticles;
   if (doMC_){

     edm::Handle<std::vector<PileupSummaryInfo> > pileupInfo;
     iEvent.getByLabel("addPileupInfo", pileupInfo);
     std::vector<PileupSummaryInfo>::const_iterator PVI;

     for (PVI = pileupInfo->begin(); PVI !=pileupInfo->end(); ++PVI){
       if( PVI->getBunchCrossing() == 0 ){ // in-time PU
	 event_genPU  = PVI->getPU_NumInteractions();
	 event_genPUtrue = PVI->getTrueNumInteractions();
       }
     }
     
     iEvent.getByLabel( genParticlesProducer_, genParticles );
     //cout << "nGenParticles="<<genParticles->size()<<std::endl;
   }
   else {
     event_genPU = -999;
     event_genPUtrue = -999;
   }

   //Electron collection
   edm::Handle<reco::GsfElectronCollection> electronHandle;
   iEvent.getByLabel(electronsProducer_,electronHandle);
   //cout << "nElectrons=" << electronHandle->size() << std::endl;   

   edm::Handle<EcalRecHitCollection> EBRecHits_;
   edm::Handle<EcalRecHitCollection> EERecHits_;
   iEvent.getByLabel(reducedBarrelEcalRecHitCollection_, EBRecHits_);
   iEvent.getByLabel(reducedEndcapEcalRecHitCollection_, EERecHits_);

   lazyTools = new EcalClusterLazyTools( iEvent, iSetup, reducedBarrelEcalRecHitCollection_, reducedEndcapEcalRecHitCollection_ );

   if (doMC_){
     for(unsigned int j=0; j<genParticles->size(); ++j ){
       const reco::GenParticle & p = (*genParticles)[ j ];
       if (p.status()==1 && (p.pdgId()==22 || p.pdgId()==11 || p.pdgId()==-11)) GenParticles_Phi->Fill(p.phi());
     }
   }

   int times_filled=0;

   if (!doElectrons_) coll_size = photonHandle->size(); else coll_size=electronHandle->size();

   std::vector<OrderPair> ordered;
   for (int i=0; i<coll_size; i++){
     reco::SuperClusterRef SCIter;
     if (!doElectrons_) SCIter = (*photonHandle)[i].superCluster();
     else SCIter = (*electronHandle)[i].superCluster();
     if (SCIter.isNull()) continue;
     ordered.push_back(std::make_pair(SCIter->rawEnergy() * sin(SCIter->position().theta()),i));
   }
   std::sort(ordered.begin(),ordered.end(),indexComparator);


   for (unsigned int orderindex=0; orderindex<ordered.size(); orderindex++){

     int index = ordered.at(orderindex).second;

     bool matched=false;

     reco::SuperClusterRef SCIter;
     if (!doElectrons_) SCIter = (*photonHandle)[index].superCluster();
     else SCIter = (*electronHandle)[index].superCluster();
     if (SCIter.isNull()) continue;

     //Event variables
     event_number = iEvent.id().event();
     event_run = iEvent.id().run(); 
     event_ls = iEvent.luminosityBlock();
     
     event_rho=rho;
     event_sigma=sigma;
     event_nRecoVertex=vertexHandle->size();

     if (!doElectrons_) {
       pho_energy = (*photonHandle)[index].energy();
       pho_et = pho_energy / TMath::CosH((*photonHandle)[index].eta());
       pho_objectr9 = (*photonHandle)[index].r9();
     }
     else {
       pho_energy = (*electronHandle)[index].energy();
       pho_et = pho_energy / TMath::CosH((*electronHandle)[index].eta());
       pho_objectr9 = (*electronHandle)[index].r9();
     }

     //SC default
     phoSC_energy = SCIter->energy();
     phoSC_et = SCIter->energy() * sin(SCIter->position().theta());

     // Our corrections (+crack+local): for photons is object energy (SC is electron corrected), for electrons is SC energy (object has track info)
     phoSC_correctedenergy = (!doElectrons_) ? pho_energy : phoSC_energy;

     //SC Raw
     phoSC_RawEnergy = SCIter->rawEnergy();
     phoSC_RawEt = SCIter->rawEnergy() * sin(SCIter->position().theta());

     //SC Raw + C(eta)
     if (fabs(SCIter->eta())<1.5) phoSC_RawEnergyCetaCorr = SCIter->rawEnergy()/f5x5((int)(fabs(SCIter->eta())*(5/0.087)));
     else phoSC_RawEnergyCetaCorr = SCIter->rawEnergy() + SCIter->preshowerEnergy();
     phoSC_RawEtCetaCorr = phoSC_RawEnergyCetaCorr * sin(SCIter->position().theta());

     phoSC_GeomEta = SCIter->eta();
     phoSC_GeomPhi = SCIter->phi();

     pho_isInPhiCrack = isInPhiCracks(SCIter->phi(), SCIter->eta());
     pho_isInEtaCrack = isInEtaCracks(SCIter->eta());

     if (fabs(SCIter->eta())<1.44 || fabs(SCIter->eta())>1.56) pho_isInEBEEGap = false;
     else pho_isInEBEEGap = true;

     pho_sigmaIetaIeta = sqrt((lazyTools->scLocalCovariances(*SCIter)).at(0));

     const reco::CaloClusterPtr phoSeed = SCIter->seed();
     pho_e3x3 = lazyTools->e3x3(*phoSeed);


     phoSC_r9 = pho_e3x3/phoSC_RawEnergy;

     pho_e5x5 = lazyTools->e5x5(*phoSeed);
     pho_etaWidth = SCIter->etaWidth();
     pho_phiWidth = SCIter->phiWidth();
     if (pho_etaWidth!=0) pho_Brem = pho_phiWidth/pho_etaWidth;
     else pho_Brem = -1; 



     if (doMC_){  
       double bestPtdiff=500.0;
       int igpsl=-1;
       
       for (unsigned int j=0; j<genParticles->size(); ++j ) {
	 const reco::GenParticle & p = (*genParticles)[ j ];
	 if (p.status()==1 && reco::deltaR(phoSC_GeomEta, phoSC_GeomPhi, p.eta(), p.phi())<0.2 && fabs(p.pt()-phoSC_RawEtCetaCorr)<bestPtdiff) {
	   bestPtdiff=fabs(p.pt()-phoSC_RawEtCetaCorr);
	   igpsl=j;
	 }
       }
      
       if (igpsl!=-1){
	 const reco::GenParticle & myGenPart = (*genParticles)[ igpsl ];
	 int pdg_to_match;
	 if (!doElectrons_) pdg_to_match=22; else pdg_to_match=11;
	 if(myGenPart.pdgId()==pdg_to_match || myGenPart.pdgId()==-pdg_to_match){
	   matched = true;
	   GenEnergy = myGenPart.energy();
	   GenEta = myGenPart.eta();
	   GenPhi = myGenPart.phi();
	 }
       }
     }


     pixelseed = (!doElectrons_) ? (*photonHandle)[index].hasPixelSeed() : -999;

     track_pt =  (doElectrons_) ? (*electronHandle)[index].gsfTrack()->pt() : -999;
     track_eta = (doElectrons_) ? (*electronHandle)[index].gsfTrack()->eta() : -999;
     track_phi = (doElectrons_) ? (*electronHandle)[index].gsfTrack()->phi() : -999;

     if (!doElectrons_) {
       ele_looseID=-999;
       ele_mediumID=-999;
       ele_tightID=-999;
     }
     else {
       reco::GsfElectronRef ele(electronHandle,index);
       double iso_ch =  (*(isoVals)[0])[ele];
       double iso_em = (*(isoVals)[1])[ele];
       double iso_nh = (*(isoVals)[2])[ele];
       ele_looseID = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele, conversions_h, beamSpot, vertexHandle, iso_ch, iso_em, iso_nh, rhoIso) ? 1 : 0;
       ele_mediumID = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele, conversions_h, beamSpot, vertexHandle, iso_ch, iso_em, iso_nh, rhoIso) ? 1 : 0;
       ele_tightID = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, ele, conversions_h, beamSpot, vertexHandle, iso_ch, iso_em, iso_nh, rhoIso) ? 1 : 0;
     }

     event_SCindex = times_filled;

     if (doMC_ && !matched) continue;

     myTree_->Fill();
     times_filled++;

   } // end photon/electron loop

   delete lazyTools;

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonPUAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonPUAnalyzer::endJob() {
}


bool PhotonPUAnalyzer::isInPhiCracks(double phi, double eta){


  // tranform radiants [-pi,pi] in degrees [0,360]
  phi = (phi+TMath::Pi()) *180/TMath::Pi();

  // each supermodule is 20 degrees wide
  Double_t moduleWidth = 20;

  // the first module is centered at phi=0, so the first cracks are +10 and -10
  Double_t phi0 = 10.;

  // set a fiducial cut around the crack of +2-2 degrees
  Double_t fiducialCut = 2.;

  bool OK = false;
  if (fabs(eta)<1.44){
  for (Int_t i = 0 ; i < 18; ++i){
    if ((phi0 + moduleWidth*i -fiducialCut) <= phi && phi <= (phi0 + moduleWidth*i + fiducialCut)) OK = true;
    //        cout << " PHI " << (phi0 + moduleWidth*i -fiducialCut) << " " << phi << " " <<  (phi0 + moduleWidth*i + fiducialCut)  << " " << OK << endl ;
  }
  }

  //  cout << "is in phi crack ? " << OK << endl;
  return OK;
}

bool PhotonPUAnalyzer::isInEtaCracks(double eta){

     /*
       Yuri Maravin eta cracks def :
     double emAbsEta = fabs(phoSC_GeomEta);
     pho_isInCrack = 0;
     if ( emAbsEta < 0.018 ||
	 (emAbsEta > 0.423 && emAbsEta < 0.461) || 
	 (emAbsEta > 0.770 && emAbsEta < 0.806) || 
	 (emAbsEta > 1.127 && emAbsEta < 1.163) || 
	 (emAbsEta > 1.460 && emAbsEta < 1.558) )
       pho_isInCrack = 1;
     */

  const Int_t nBinsEta = 5;
  Double_t leftEta [nBinsEta]       = {0.00, 0.42, 0.77, 1.13, 1.46};
  Double_t rightEta[nBinsEta]       = {0.02, 0.46, 0.81, 1.16, 9999.};

  bool OK = false;
  if (TMath::Abs(eta)<1.44) {
          for (Int_t i = 0; i< nBinsEta; ++i){
                  if (leftEta[i] < TMath::Abs(eta) && TMath::Abs(eta) < rightEta[i] ) OK = true;
          }
  }
  else if (TMath::Abs(eta)>1.44 && TMath::Abs(eta)<1.56) OK = true;
  else if (TMath::Abs(eta)>1.56) OK = false;
    //    cout << leftEta[i] << " " << TMath::Abs(eta) << " " << rightEta[i] <<  " " << OK << endl;

  //  cout << "IS IN CRACK ? " << OK << endl;
  return OK;
}

double PhotonPUAnalyzer::f5x5( double iEta ) {
  if ( iEta < 40.2198 ) return 1;
  return 1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198);
}





//define this as a plug-in
DEFINE_FWK_MODULE(PhotonPUAnalyzer);
