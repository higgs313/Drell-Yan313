import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

IsMC=True
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")    

if IsMC:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6' #80X_mcRun2_asymptotic_2016_miniAODv2_v1'
else :
    # process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0' # ICHEP
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4' # Run B-G sept rereco 2016
    # process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v14' # Run H prompt-reco 2016
print process.GlobalTag.globaltag
    
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.GoodVertexFilter = cms.EDFilter("VertexSelector",
                                        src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                        filter = cms.bool(True)   # otherwise it won't filter the events, just produce an empty vertex collection.
                                        )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
		#		    '/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/72446D9C-D89C-E611-9060-002590A3C984.root',
#				    '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root'
				    '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/02456AB6-BBDF-E611-97A4-0CC47AD98F70.root'
                            ),
#                            skipEvents = cms.untracked.uint32(0)                       
                            )



process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root')
                               )



process.demo = cms.EDAnalyzer('DPSwithMuMu',
                              bits          = cms.InputTag("TriggerResults", "", "HLT"),
                              outputFile       = cms.string("patTuple.root"),
                              printDebug        = cms.bool(False),
                              isMC =cms.bool(IsMC),
                              isMPI =cms.bool(True),
			      genjetCollection = cms.InputTag("slimmedGenJets"),
			      src = cms.InputTag("prunedGenParticles"),
			      jetCollection = cms.InputTag("slimmedJets"),
			      beamSpot = cms.InputTag("offlineBeamSpot"),
			      muonCollection = cms.InputTag("slimmedMuons"),
			      metCollection = cms.InputTag("slimmedMETs"),
			      vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
			      triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
			      rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
			      puCollection = cms.InputTag("slimmedAddPileupInfo"),
			      genCollection = cms.InputTag("generator"),
			      lhepCollection = cms.InputTag("externalLHEProducer"),
			      triggerSet = cms.InputTag("selectedPatTrigger"),
                              )

process.p= cms.Path(process.GoodVertexFilter*
                    process.demo)             
                    
process.TFileService = cms.Service("TFileService",fileName = cms.string('patTuple.root'))#for embeded anlzr
#process.outpath = cms.EndPath(process.out)
