import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = cms.string('START53_V7D::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/mc/Summer12_DR53X/SingleElectronPt35/AODSIM/PU_Fix10_START53_V7A-v2/00000/D432689E-F7F8-E111-8719-0026189438AA.root')
   fileNames = cms.untracked.vstring('dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/mc/Summer12_DR53X/SingleGammaPt35/AODSIM/PU_Fix20_START53_V7A-v2/00000//92485DB3-F6F8-E111-8A4F-0025905822B6.root')
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.demo = cms.EDAnalyzer("PhotonPUAnalyzer",
                              doElectrons = cms.bool(False),
                              doMC = cms.bool(True),
                              photonsProducer = cms.InputTag("photons"),
                              genParticlesProducer = cms.InputTag("genParticles"),
                              electronsProducer = cms.InputTag("gsfElectrons"),
                              reducedBarrelEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                              reducedEndcapEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                              OutputFile = cms.string("outfile.root"),
                              srcRho = cms.InputTag('kt6PFJets','rho'),
                              srcSigma = cms.InputTag('kt6PFJets','sigma'),
                              vertexProducer = cms.InputTag('offlinePrimaryVerticesWithBS')
                              )

process.p = cms.Path(process.demo)
