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
options.register('version', 
		'Full',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"Version of the analysis to run. (Full, Resolved, Boosted)"
		)
options.register('systematics', 
		False,		
		VarParsing.multiplicity.singleton,
		VarParsing.varType.bool,
		"Run systematics, default false."
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
options.register('CSVFile', 
		'CSVv2_Moriond17_B_H.csv',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"CSVFile"
		)

options.register('cMVAFile', 
		'cMVAv2_Moriond17_B_H.csv',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"cMVAFile"
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
			#'/store/user/jsomalwa/B2GAnaFW_80X_V2p1/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161018_211211/0000/B2GEDMNtuple_1.root',
			#'/store/user/algomez/RPVStopStopToJets_UDD323_M-120_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/170131_140142/0000/B2GEDMNtuple_5.root',
			#'/store/group/phys_b2g/B2GAnaFW_80X_V2p3/JetHT/Run2016C/JetHT/Run2016C-23Sep2016-v1_B2GAnaFW_80X_V2p3/161216_220503/0000/B2GEDMNtuple_10.root',
			#'/store/user/jsomalwa/B2GAnaFW_80X_V2p1/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161018_211709/0000/B2GEDMNtuple_3.root',
			#'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161222_110143/0000/B2GEDMNtuple_736.root',
			#'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p4/161222_105211/0000/B2GEDMNtuple_146.root',

			#'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3_B2GAnaFW_80X_V2p4/161222_105345/0000/B2GEDMNtuple_20.root',
			#'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3_B2GAnaFW_80X_V2p4/161222_105345/0000/B2GEDMNtuple_21.root',
			'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3_B2GAnaFW_80X_V2p4/161222_105345/0000/B2GEDMNtuple_10.root',
			#'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3_B2GAnaFW_80X_V2p4/161222_105345/0000/B2GEDMNtuple_4.root',

	    ),
	#	lumisToProcess = cms.untracked.VLuminosityBlockRange('1:2083-1:max'),
	)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (options.maxEvents) )

from RUNA.RUNAnalysis.scaleFactors import scaleFactor

bjsample = True if 'bj' in NAME else False
if 'JetHT' in NAME: 
	SF = 1
	isData=True
	options.systematics = False
else:
	isData=False
	SF = scaleFactor(NAME)


##############################
#####   Resolved analysis

ResolvedTriggers = [  'HLT_PFHT800', 'HLT_PFHT900', 'HLT_PFHT750_4Jet', 'HLT_PFHT800_4Jet50', 'HLT_PFJet450' ]

process.ResolvedAnalysisPlots = cms.EDAnalyzer('RUNResolvedAnalysis',
		cutAK4jetPt 		= cms.double( 80.0 ),	# default 80.0
		cutAK4HT 		= cms.double( 900.0 ),	# default 800.0
		cutAK4MassAsym		= cms.double( 0.1 ),	# default 0.2
		cutDelta 		= cms.double( 200 ),	# default 180.0
		cutDeltaEtaDijetSyst	= cms.double( 1.0 ),	# default .75
		triggerPass 		= cms.vstring( ResolvedTriggers ),
		scale 			= cms.double( SF ),
		dataPUFile		= cms.string( options.namePUFile  ),
		jecVersion		= cms.string( options.jecVersion ),
		btagCSVFile		= cms.string( options.CSVFile  ),
		isData			= cms.bool( isData ),
		LHEcont			= cms.bool( True if 'QCD_Pt' in NAME else False ), ## logic is oposite
		pairingMethod		= cms.string( "deltaR" ),
		mkTree			= cms.bool( True ),
)

#process.ResolvedAnalysisPlotsScouting = process.ResolvedAnalysisPlots.clone( cutAK4jetPt = cms.double( 50.0 ), cutAK4HT = cms.double( 450 ), mkTree = cms.bool( True ) )
process.ResolvedAnalysisPlotsMassPairing = process.ResolvedAnalysisPlots.clone( pairingMethod = cms.string( "mass" ), mkTree = cms.bool( False ) )
process.ResolvedAnalysisPlotsChi2Pairing = process.ResolvedAnalysisPlots.clone( pairingMethod = cms.string( "chi2" ), mkTree = cms.bool( False ) )

process.ResolvedAnalysisPlotsJESUp = process.ResolvedAnalysisPlots.clone( systematics = cms.string( 'JESUp' ), mkTree = cms.bool( False ) )
process.ResolvedAnalysisPlotsJESDown = process.ResolvedAnalysisPlots.clone( systematics = cms.string( 'JESDown' ), mkTree = cms.bool( False ) )
process.ResolvedAnalysisPlotsJERUp = process.ResolvedAnalysisPlots.clone( systematics = cms.string( 'JERUp' ), mkTree = cms.bool( False ) )
process.ResolvedAnalysisPlotsJERDown = process.ResolvedAnalysisPlots.clone( systematics = cms.string( 'JERDown' ), mkTree = cms.bool( False ) )

############################################################


##############################
#####   Boosted analysis
BoostedTriggers =  [ 'HLT_PFHT800', 
			'HLT_PFHT900', 
			'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50',
			'HLT_AK8PFJet360_TrimMass30', 'HLT_PFJet450',
			'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
			'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087',
			'HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20' ] 

process.BoostedAnalysisPlots = cms.EDAnalyzer('RUNBoostedAnalysis',
		#cutAK8jetPt 		= cms.double( 150.0 ),	# default 150.0
		#cutAK8HT 		= cms.double( 900.0 ),	# default 900.0
		#cutAK8MassAsym		= cms.double( 0.1 ),	# default 0.1
		#cutTau21 		= cms.double( 0.45 ),	# default 0.45
		#cutDeltaEtaDijet	= cms.double( 1.5 ),	# default 1.5
		triggerPass 		= cms.vstring( BoostedTriggers ),
		dataPUFile		= cms.string( options.namePUFile  ),
		jecVersion		= cms.string( options.jecVersion ),
		isData			= cms.bool( isData ),
		isTTbar			= cms.bool( True if 'TTJets' in NAME else False ),
		btagCSVFile		= cms.string( options.CSVFile  ),
		btagMVAFile		= cms.string( options.cMVAFile  ),
		LHEcont			= cms.bool( True if 'QCD_Pt' in NAME else False ), ## logic is oposite
		scale 			= cms.double( SF ),
		mkTree			= cms.bool( True ),
)

#process.BoostedAnalysisPlotsSortInMass = process.BoostedAnalysisPlots.clone( sortInMass = cms.bool( True ), mkTree = cms.bool( False ) )
#process.BoostedAnalysisPlotsSortInTau21 = process.BoostedAnalysisPlots.clone( sortInTau21 = cms.bool( True ), mkTree = cms.bool( False ) )

process.BoostedAnalysisPlotsJESUp = process.BoostedAnalysisPlots.clone( systematics = cms.string( 'JESUp' ), mkTree = cms.bool( False ) )
process.BoostedAnalysisPlotsJESDown = process.BoostedAnalysisPlots.clone( systematics = cms.string( 'JESDown' ), mkTree = cms.bool( False ) )
process.BoostedAnalysisPlotsJERUp = process.BoostedAnalysisPlots.clone( systematics = cms.string( 'JERUp' ), mkTree = cms.bool( False ) )
process.BoostedAnalysisPlotsJERDown = process.BoostedAnalysisPlots.clone( systematics = cms.string( 'JERDown' ), mkTree = cms.bool( False ) )

process.BoostedAnalysisPlotsPuppi = process.BoostedAnalysisPlots.clone( 
		PUMethod		= cms.string('Puppi'),
		cutDeltaEtaDijet	= cms.double( 1. ),	# default 1.5
		jetPt 			= cms.InputTag('jetsAK8Puppi:jetAK8PuppiPt'),
		jetEta			= cms.InputTag('jetsAK8Puppi:jetAK8PuppiEta'),
		jetPhi 			= cms.InputTag('jetsAK8Puppi:jetAK8PuppiPhi'),
		jetE 			= cms.InputTag('jetsAK8Puppi:jetAK8PuppiE'),
		jetTrimmedMass 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppitrimmedMass'),
		jetPrunedMass 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppisoftDropMass'),
		jetFilteredMass 	= cms.InputTag('jetsAK8Puppi:jetAK8PuppifilteredMass'),
		jetSoftDropMass		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiprunedMass'), 	#### Change name for puppi+softdrop
		jetTau1 		= cms.InputTag('jetsAK8Puppi:jetAK8Puppitau1'),
		jetTau2 		= cms.InputTag('jetsAK8Puppi:jetAK8Puppitau2'),
		jetTau3 		= cms.InputTag('jetsAK8Puppi:jetAK8Puppitau3'),
		jetNSubjets 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppinSubJets'),
		jetSubjetIndex0 	= cms.InputTag('jetsAK8Puppi:jetAK8PuppivSubjetIndex0'),
		jetSubjetIndex1 	= cms.InputTag('jetsAK8Puppi:jetAK8PuppivSubjetIndex1'),
		jetSubjetIndex2 	= cms.InputTag('jetsAK8Puppi:jetAK8PuppivSubjetIndex0'),
		jetSubjetIndex3 	= cms.InputTag('jetsAK8Puppi:jetAK8PuppivSubjetIndex1'),
		jetKeys 		= cms.InputTag('jetKeysAK8Puppi'),
		jetCSVv2 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiCSVv2'),
		jetCMVAv2 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiCMVAv2'),
		jetDoubleB	 	= cms.InputTag('jetsAK8Puppi:jetAK8PuppiDoubleBAK8'),
		jetArea 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppijetArea'),
		jetGenPt 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiGenJetPt'),
		jetGenEta		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiGenJetEta'),
		jetGenPhi 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiGenJetPhi'),
		jetGenE 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiGenJetE'),
		jetHadronFlavour 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppiHadronFlavour'),
		jecFactor 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppijecFactor0'),
		neutralHadronEnergyFrac		= cms.InputTag('jetsAK8Puppi:jetAK8PuppineutralHadronEnergyFrac'),
		neutralEmEnergyFrac 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppineutralEmEnergyFrac'),
		chargedEmEnergyFrac 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppichargedEmEnergyFrac'),
		muonEnergyFrac 			= cms.InputTag('jetsAK8Puppi:jetAK8PuppiMuonEnergy'),
		chargedHadronEnergyFrac		= cms.InputTag('jetsAK8Puppi:jetAK8PuppichargedHadronEnergyFrac'),
		neutralMultiplicity 	= cms.InputTag('jetsAK8Puppi:jetAK8PuppineutralMultiplicity'),
		chargedMultiplicity 		= cms.InputTag('jetsAK8Puppi:jetAK8PuppichargedMultiplicity'),
		#### Subjets
		subjetPt 		= cms.InputTag('subjetsAK8Puppi:subjetAK8PuppiPt'),
		subjetEta 		= cms.InputTag('subjetsAK8Puppi:subjetAK8PuppiEta'),
		subjetPhi 		= cms.InputTag('subjetsAK8Puppi:subjetAK8PuppiPhi'),
		subjetE 		= cms.InputTag('subjetsAK8Puppi:subjetAK8PuppiE'),
		subjetMass 		= cms.InputTag('subjetsAK8Puppi:subjetAK8PuppiMass'),
		subjetCSVv2 		= cms.InputTag('subjetsAK8Puppi:subjetAK8PuppiCSVv2'),
		subjetCMVAv2 		= cms.InputTag('subjetsAK8Puppi:subjetAK8PuppiCMVAv2'),
		mkTree 			= cms.bool( True )
		)

process.BoostedAnalysisPlotsPuppiJESUp = process.BoostedAnalysisPlotsPuppi.clone( systematics = cms.string( 'JESUp' ) )
process.BoostedAnalysisPlotsPuppiJESDown = process.BoostedAnalysisPlotsPuppi.clone( systematics = cms.string( 'JESDown' ) )
process.BoostedAnalysisPlotsPuppiJERUp = process.BoostedAnalysisPlotsPuppi.clone( systematics = cms.string( 'JERUp' ) )
process.BoostedAnalysisPlotsPuppiJERDown = process.BoostedAnalysisPlotsPuppi.clone( systematics = cms.string( 'JERDown' ) )

############################################################


process.p = cms.Path()
if 'Resolved' in options.version:
	outputNAME = 'ResolvedAnalysis_'
	process.p += process.ResolvedAnalysisPlots
	process.p += process.ResolvedAnalysisPlotsMassPairing
	process.p += process.ResolvedAnalysisPlotsChi2Pairing
	if options.systematics:
		process.p += process.ResolvedAnalysisPlotsJESUp
		process.p += process.ResolvedAnalysisPlotsJESDown
		process.p += process.ResolvedAnalysisPlotsJERUp
		process.p += process.ResolvedAnalysisPlotsJERDown

elif 'Boosted' in options.version:
	outputNAME = 'BoostedAnalysis_'
	process.p += process.BoostedAnalysisPlots
	process.p += process.BoostedAnalysisPlotsPuppi
	#process.p += process.BoostedAnalysisPlotsSortInMass
	#process.p += process.BoostedAnalysisPlotsSortInTau21
	if options.systematics:
		process.p += process.BoostedAnalysisPlotsJESUp
		process.p += process.BoostedAnalysisPlotsJESDown
		process.p += process.BoostedAnalysisPlotsJERUp
		process.p += process.BoostedAnalysisPlotsJERDown
		#process.p += process.BoostedAnalysisPlotsPuppiJESUp
		#process.p += process.BoostedAnalysisPlotsPuppiJESDown
else: 
	outputNAME = 'FullAnalysis_'
	process.p += process.ResolvedAnalysisPlots
	#process.p += process.ResolvedAnalysisPlotsScouting
	process.p += process.ResolvedAnalysisPlotsMassPairing
	process.p += process.ResolvedAnalysisPlotsChi2Pairing
	process.p += process.BoostedAnalysisPlots
	#process.p += process.BoostedAnalysisPlotsSortInMass
	#process.p += process.BoostedAnalysisPlotsSortInTau21
	process.p += process.BoostedAnalysisPlotsPuppi

	if options.systematics:
		process.p += process.BoostedAnalysisPlotsJESUp
		process.p += process.BoostedAnalysisPlotsJESDown
		process.p += process.BoostedAnalysisPlotsJERUp
		process.p += process.BoostedAnalysisPlotsJERDown
		process.p += process.ResolvedAnalysisPlotsJESUp
		process.p += process.ResolvedAnalysisPlotsJESDown
		process.p += process.ResolvedAnalysisPlotsJERUp
		process.p += process.ResolvedAnalysisPlotsJERDown

process.TFileService=cms.Service("TFileService",fileName=cms.string( ( 'Rootfiles/' if options.local else '' )+'RUN'+outputNAME+NAME+'.root' ) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
