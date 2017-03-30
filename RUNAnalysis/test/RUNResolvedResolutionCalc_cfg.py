import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

### Example how to run
### cmsRun RUNFullAnalysis_cfg.py PROC=RPVStopStopToJets_UDD312_M-100 jecVersion=supportFiles/Spring16_23Sep2016V2 namePUFile=supportFiles/PileupData2016BCDEFGH_ReReco_69200.root maxEvents=10000 CSVFile=supportFiles/CSVv2_Moriond17_B_H.csv cMVAFile=supportFiles/cMVAv2_Moriond17_B_H.csv

###############################
####### Parameters ############
options = VarParsing ('python')

### General Options
options.register('PROC', 
		'RPVStopStopToJets_UDD312_M-120',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"name"
		)
options.register('local', 
		False,
		VarParsing.multiplicity.singleton,
		VarParsing.varType.bool,
		"Run locally or crab"
		)
options.register('jecVersion', 
		'Fall15_25nsV2',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"Version of the analysis to run. (Full, Resolved, Boosted)"
		)
options.register('namePUFile', 
		'PileupData2016BCDEFGH_ReReco_69200.root',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"namePUFile"
		)


options.parseArguments()

process = cms.Process("RUNAnalysis")
process.load("FWCore.MessageService.MessageLogger_cfi")
NAME = options.PROC

if options.local:
	process.load(NAME+'_RUNA_cfi')
	#process.load('RPVSt100tojj_13TeV_pythia8_RUNtuples_cfi')
else:
	process.source = cms.Source("PoolSource",
		fileNames = cms.untracked.vstring(
			'/store/user/jsomalwa/B2GAnaFW_80X_V2p1/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161018_211709/0000/B2GEDMNtuple_3.root',
			'/store/user/jsomalwa/B2GAnaFW_80X_V2p1/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161018_211709/0000/B2GEDMNtuple_2.root',

	    )
	)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (options.maxEvents) )

##############################
#####   Resolved analysis

ResolvedTriggers = [  'HLT_PFHT800', 'HLT_PFHT900', 'HLT_PFHT750_4Jet', 'HLT_PFHT800_4Jet50', 'HLT_PFJet450' ]

process.ResolvedAnalysisPlots = cms.EDAnalyzer('RUNResolvedResolutionCalc',
		cutAK4jetPt 		= cms.double( 80.0 ),	# default 80.0
		cutAK4HT 		= cms.double( 900.0 ),	# default 800.0
		cutAK4MassAsym		= cms.double( 0.1 ),	# default 0.2
		cutDelta 		= cms.double( 200 ),	# default 180.0
		cutDeltaEtaDijetSyst	= cms.double( 1.0 ),	# default .75
		dataPUFile		= cms.string( options.namePUFile  ),
		jecVersion		= cms.string( options.jecVersion ),
)



process.p = cms.Path( process.ResolvedAnalysisPlots )

process.TFileService=cms.Service("TFileService",fileName=cms.string( ( 'Rootfiles/' if options.local else '' )+'RUNResolvedResolutionCalc_'+NAME+'.root' ) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
