import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

#PU = sys.argv[3]
#NAME = sys.argv[2]

process = cms.Process("Match")
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12')
#############   Set the number of events #############

#############   Define the source file ###############
#process.load('RPVSt100tojj_pythia8_13TeV_Asympt25ns_AOD_cfi')
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/D03C6078-BDF6-E611-99CC-FA163EA77046.root',
		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/D66FF887-86F7-E611-B6F6-0090FAA57FC4.root',
		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/E21BE3AC-8CF7-E611-9E89-001E67E6F8C3.root',
		) )
process.TFileService=cms.Service("TFileService", fileName=cms.string('simpleMatching_RPVStop.root'))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 10000 ) )
    
process.selectedPatJetsAK4 = cms.EDFilter("PATJetSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt > 80 && abs(eta) < 2.5") )
		
process.selectedPatJetsAK8 = cms.EDFilter("PATJetSelector",
		src = cms.InputTag("slimmedJetsAK8"),
		cut = cms.string("pt > 200 && abs(eta) < 2.5") )
		


#############   User analyzer (PF jets) ##
process.histos = cms.EDAnalyzer("Matching",
		AK8jets = cms.InputTag( "selectedPatJetsAK8" ),
		AK4jets = cms.InputTag( "selectedPatJetsAK4" ),
		genParticles = cms.InputTag( 'prunedGenParticles' ),
		#### stop
		particle1 = cms.int32( 1000002 ),
		particle2 = cms.int32( 0 ),
		particle3 = cms.int32( 0 ),
		#### squark
		#particle1 = cms.int32( 1000005 ),
		#particle2 = cms.int32( 1000035 ),
		#particle3 = cms.int32( 1000025 ),
		### gluino
		#particle1 = cms.int32( 1000021 ),
		#particle2 = cms.int32( 6 ),
		#particle3 = cms.int32( 24 ),
		#### RSGraviton
		#particle1 = cms.int32( 24 ),
		#particle2 = cms.int32( 0 ),
		#particle3 = cms.int32( 0 ),
		)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
		maxEventsToPrint = cms.untracked.int32(1),
		printVertex = cms.untracked.bool(False),
		src = cms.InputTag("prunedGenParticles")
		)


#############   Path       ###########################
process.p = cms.Path(
	process.selectedPatJetsAK4
	* process.selectedPatJetsAK8
	* process.histos
	* process.printTree
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100

