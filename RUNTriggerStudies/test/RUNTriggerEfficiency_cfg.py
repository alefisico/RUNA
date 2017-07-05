import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('PROC', 
		'JetHT_Run2015C',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"name"
		)
options.register('version', 
		'Full',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"Version of the analysis to run. (Full, Resolved, Boosted)"
		)
options.register('jecVersion', 
		'../../RUNAnalysis/test/supportFiles/Spring16_25nsV7BCD',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"Version of the analysis to run. (Full, Resolved, Boosted)"
		)

options.register('miniAOD', 
		False,
		VarParsing.multiplicity.singleton,
		VarParsing.varType.bool,
		"Run from miniAOD or B2GNtuples"
		)

options.parseArguments()

process = cms.Process("TriggerEfficiency")
process.load("FWCore.MessageService.MessageLogger_cfi")
NAME = options.PROC

#process.load('JetHT_Run2016C_B2GNtuple_cfi')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('275657:1-275657:5')
process.source = cms.Source("PoolSource",
		fileNames = cms.untracked.vstring(
			'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016B-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_133516/0000/B2GEDMNtuple_10.root',
			'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/SingleMuon/Run2016B-23Sep2016-v1_B2GAnaFW_80X_v2p4/161221_133516/0000/B2GEDMNtuple_11.root',
			),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (options.maxEvents) )


if 'MET' in NAME: basedline = 'HLT_PFMET170_HBHECleaned'
elif 'SingleMu' in NAME: basedline =  'HLT_Mu50' #'HLT_IsoMu17_eta2p1_v'
#elif 'SingleElectron' in NAME: basedline =  'HLT_Mu50' #'HLT_IsoMu17_eta2p1_v'
else: basedline = 'HLT_PFHT475'

####### Adding Corrections
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

####################
#### Resolved 
process.ResolvedTriggerEfficiency = cms.EDAnalyzer(('RUNResolvedMiniAODTriggerEfficiency' if options.miniAOD else 'RUNResolvedTriggerEfficiency'),
		#cutAK4jetPt	= cms.double( 50.0 ), 	# default 50.
		#cutAK4jet4Pt	= cms.double( 80.0 ), 	# default 80.
		#cutAK4HT	= cms.double( 800.0 ), 	# default 800.
		baseTrigger		= cms.string(basedline),
		triggerPass		= cms.vstring( [  'HLT_PFHT800', 'HLT_PFHT900', 'HLT_PFHT750_4Jet', 'HLT_PFHT800_4Jet50' ] ),
)

process.ResolvedTriggerEfficiencyPFHT7504Jet = process.ResolvedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT750_4JetPt', 'HLT_PFHT800_4Jet50' ] ) )
process.ResolvedTriggerEfficiencyPFHT800 = process.ResolvedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900' ] ) )
process.ResolvedTriggerEfficiencyPFHT650WideJetMJJ900 = process.ResolvedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5', 'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5' ] ) )
process.ResolvedTriggerEfficiencyPFJet450 = process.ResolvedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [  'HLT_PFHT800', 'HLT_PFHT900', 'HLT_PFHT750_4Jet', 'HLT_PFHT800_4Jet50', 'HLT_PFJet450' ] ) )
process.ResolvedTriggerEfficiencySeveralTriggers = process.ResolvedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [  'HLT_PFHT800', 'HLT_PFHT900', 'HLT_PFHT750_4Jet', 'HLT_PFHT800_4Jet50', 'HLT_PFJet450', 'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5', 'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5' ] ) )
############################################################



####################
#### Boosted
process.BoostedTriggerEfficiency = cms.EDAnalyzer( ('RUNBoostedMiniAODTriggerEfficiency' if options.miniAOD else 'RUNBoostedTriggerEfficiency'),
		cutAK8jetPt	= cms.double( 170.0 ),		# default 150.
		#cutAK8HT	= cms.double( 900.0 ),		# default 900.
		#cutAK8jet1Pt	= cms.double( 500.0 ),		# default 500.
		#cutAK8jet2Pt	= cms.double( 450.0 ),		# default 450.
		#cutAK8jet1Mass	= cms.double( 60.0 ),		# default 60.
		baseTrigger	= cms.string(basedline),
		jecVersion	= cms.string( options.jecVersion ),
		triggerPass	= cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50' ] ),
)

process.BoostedTriggerEfficiencyPFHT800 = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900' ] ) )

process.BoostedTriggerEfficiencyAK8PFHT700TrimMass50 = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50'] ) )

process.BoostedTriggerEfficiencyAK8PFPt360TrimMass30 = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8PFJet360_TrimMass30' ] ) )

process.BoostedTriggerEfficiencyAK8PFHT7504Jet = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT750_4JetPt', 'HLT_PFHT800_4Jet50' ] ) )

process.BoostedTriggerEfficiencyAK8DiPFJet280220TrimMass30Btagp20 = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20', 'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087', 'HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20' ] ) )

process.BoostedTriggerEfficiencyAK8DiPFJet280220TrimMass30 = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8DiPFJet280_200_TrimMass30'] ) )

process.BoostedTriggerEfficiencyPFJet450 = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50', 'HLT_PFJet450' ] ) )

process.BoostedTriggerEfficiencyPFHT650WideJetMJJ900 = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5', 'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5' ] ) )


process.BoostedTriggerEfficiencySeveralTriggers = process.BoostedTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 
			'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50',
			'HLT_AK8PFJet360_TrimMass30', 'HLT_PFJet450',
			'HLT_PFHT750_4JetPt', 'HLT_PFHT800_4Jet50',
			'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
			'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087',
			'HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20' ] ) )

if not options.miniAOD:
	process.BoostedTriggerEfficiencyPuppi = process.BoostedTriggerEfficiency.clone(
			PUMethod = cms.string('Puppi'),
			jetPt = cms.InputTag("jetsAK8Puppi:jetAK8PuppiPt"),
			jetEta = cms.InputTag("jetsAK8Puppi:jetAK8PuppiEta"),
			jetPhi = cms.InputTag("jetsAK8Puppi:jetAK8PuppiPhi"),
			jetE = cms.InputTag("jetsAK8Puppi:jetAK8PuppiE"),
			jetArea = cms.InputTag("jetsAK8Puppi:jetAK8PuppijetArea"),
			jetTrimmedMass = cms.InputTag("jetsAK8Puppi:jetAK8PuppitrimmedMass"),
			jetPrunedMass = cms.InputTag("jetsAK8Puppi:jetAK8PuppiprunedMass"),
			jetSoftDropMass = cms.InputTag("jetsAK8Puppi:jetAK8PuppisoftDropMass"),
			jetCSVv2 = cms.InputTag("jetsAK8Puppi:jetAK8PuppiCSVv2"),
			jecFactor = cms.InputTag("jetsAK8Puppi:jetAK8PuppijecFactor0"),
			neutralHadronEnergyFrac = cms.InputTag("jetsAK8Puppi:jetAK8PuppineutralHadronEnergyFrac"),
			neutralEmEnergyFrac = cms.InputTag("jetsAK8Puppi:jetAK8PuppineutralEmEnergyFrac"),
			chargedEmEnergyFrac = cms.InputTag("jetsAK8Puppi:jetAK8PuppichargedEmEnergyFrac"),
			muonEnergy = cms.InputTag("jetsAK8Puppi:jetAK8PuppiMuonEnergy"),
			chargedHadronEnergyFrac = cms.InputTag("jetsAK8Puppi:jetAK8PuppichargedHadronEnergyFrac"),
			neutralMultiplicity = cms.InputTag("jetsAK8Puppi:jetAK8PuppineutralMultiplicity"),
			chargedMultiplicity = cms.InputTag("jetsAK8Puppi:jetAK8PuppichargedMultiplicity"),
			)

	process.BoostedTriggerEfficiencyPuppiPFHT800 = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900' ] ) )

	process.BoostedTriggerEfficiencyPuppiAK8PFHT700TrimMass50 = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50' ] ) )

	process.BoostedTriggerEfficiencyPuppiAK8PFPt360TrimMass30 = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_AK8PFJet360_TrimMass30' ] ) )

	process.BoostedTriggerEfficiencyPuppiAK8PFHT7504Jet = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_PFHT750_4JetPt'] ) )

	process.BoostedTriggerEfficiencyPuppiAK8DiPFJet280220TrimMass30Btagp20 = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20'] ) )

	process.BoostedTriggerEfficiencyPuppiAK8DiPFJet280220TrimMass30 = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_AK8DiPFJet280_200_TrimMass30'] ) )

	process.BoostedTriggerEfficiencyPuppiPFJet450 = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50', 'HLT_PFJet450' ] ) )

	process.BoostedTriggerEfficiencyPuppiPFHT650WideJetMJJ900 = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5', 'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5' ] ) )

	process.BoostedTriggerEfficiencyPuppiSeveralTriggers = process.BoostedTriggerEfficiencyPuppi.clone( 
			triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 
				'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50',
				'HLT_AK8PFJet360_TrimMass30', 'HLT_PFJet450',
				'HLT_PFHT750_4JetPt', 'HLT_PFHT800_4Jet50',
				'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
				'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087',
				'HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20' ] ) )
############################################################



####################
#### Dijet 
process.DijetTriggerEfficiency = cms.EDAnalyzer( ('RUNDijetMiniAODTriggerEfficiency' if options.miniAOD else 'RUNDijetTriggerEfficiency'),
		#cutAK8jetPt	= cms.double( 150.0 ),		# default 150.
		#cutAK8jet1Pt	= cms.double( 500.0 ),		# default 500.
		#cutAK8jet1Mass	= cms.double( 60.0 ),		# default 60.
		PUMethod 	= cms.string('Puppi'),
		jecVersion	= cms.string( options.jecVersion ),
		baseTrigger	= cms.string(basedline),
		triggerPass	= cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50', 'HLT_AK8PFJet360_TrimMass30' ] ),
)

process.DijetTriggerEfficiencyPFHT800 = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900' ] ) )

process.DijetTriggerEfficiencyAK8PFHT700TrimMass50 = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50', ] ) )

process.DijetTriggerEfficiencyAK8PFPt360TrimMass30 = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8PFJet360_TrimMass30' ] ) )

process.DijetTriggerEfficiencyAK8PFHT7504Jet = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT750_4JetPt', 'HLT_PFHT800_4Jet' ] ) )

process.DijetTriggerEfficiencyAK8DiPFJet280220TrimMass30Btagp20 = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20', 'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087', 'HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20' ] ) )

process.DijetTriggerEfficiencyAK8DiPFJet280220TrimMass30 = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_AK8DiPFJet280_200_TrimMass30'] ) )

process.DijetTriggerEfficiencyPFHT650WideJetMJJ900 = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5', 'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5' ] ) )

process.DijetTriggerEfficiencyPFJet450 = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50', 'HLT_AK8PFJet360_TrimMass30', 'HLT_PFJet450' ] ) )

process.DijetTriggerEfficiencySeveralTriggers = process.DijetTriggerEfficiency.clone( 
		triggerPass = cms.vstring( [ 'HLT_PFHT800', 'HLT_PFHT900', 
			'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50', 'HLT_AK8PFHT750_TrimMass50',
			'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5', 'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5',
			'HLT_AK8PFJet360_TrimMass30', 'HLT_PFJet450',
			'HLT_PFHT750_4JetPt', 'HLT_PFHT800_4Jet50',
			'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
			'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087',
			'HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20' ] ) )
############################################################


process.p = cms.Path()
if 'Resolved' in options.version:
	outputNAME = 'Resolved'
	process.p += process.ResolvedTriggerEfficiency
	process.p += process.ResolvedTriggerEfficiencyPFHT7504Jet
	process.p += process.ResolvedTriggerEfficiencyPFHT800
	process.p += process.ResolvedTriggerEfficiencyPFHT650WideJetMJJ900 
	process.p += process.ResolvedTriggerEfficiencyPFJet450

elif 'Boosted' in options.version:
	outputNAME = 'Boosted'
	process.p += process.BoostedTriggerEfficiency
	process.p += process.BoostedTriggerEfficiencyPFHT800
	process.p += process.BoostedTriggerEfficiencyAK8PFHT700TrimMass50
	process.p += process.BoostedTriggerEfficiencyAK8PFPt360TrimMass30
	process.p += process.BoostedTriggerEfficiencyAK8PFHT7504Jet
	process.p += process.BoostedTriggerEfficiencyAK8DiPFJet280220TrimMass30Btagp20
	process.p += process.BoostedTriggerEfficiencyAK8DiPFJet280220TrimMass30
	process.p += process.BoostedTriggerEfficiencyPFHT650WideJetMJJ900 
	process.p += process.BoostedTriggerEfficiencyPFJet450
	process.p += process.BoostedTriggerEfficiencySeveralTriggers
	if not options.miniAOD:
		process.p += process.BoostedTriggerEfficiencyPuppi
		process.p += process.BoostedTriggerEfficiencyPuppiPFHT800
		process.p += process.BoostedTriggerEfficiencyPuppiAK8PFHT700TrimMass50
		process.p += process.BoostedTriggerEfficiencyPuppiAK8PFPt360TrimMass30
		process.p += process.BoostedTriggerEfficiencyPuppiAK8PFHT7504Jet
		process.p += process.BoostedTriggerEfficiencyPuppiAK8DiPFJet280220TrimMass30Btagp20
		process.p += process.BoostedTriggerEfficiencyPuppiAK8DiPFJet280220TrimMass30
		process.p += process.BoostedTriggerEfficiencyPuppiPFJet450
		process.p += process.BoostedTriggerEfficiencyPuppiSeveralTriggers

elif 'Dijet' in options.version:
	outputNAME = 'Dijet'
	process.p += process.DijetTriggerEfficiency
	process.p += process.DijetTriggerEfficiencyPFHT800
	process.p += process.DijetTriggerEfficiencyAK8PFHT700TrimMass50
	process.p += process.DijetTriggerEfficiencyAK8PFPt360TrimMass30
	process.p += process.DijetTriggerEfficiencyAK8PFHT7504Jet
	process.p += process.DijetTriggerEfficiencyAK8DiPFJet280220TrimMass30Btagp20
	process.p += process.DijetTriggerEfficiencyAK8DiPFJet280220TrimMass30
	process.p += process.DijetTriggerEfficiencyPFJet450
	process.p += process.DijetTriggerEfficiencyPFHT650WideJetMJJ900
	process.p += process.DijetTriggerEfficiencySeveralTriggers
else: 
	outputNAME = 'Full'
	process.p += process.ResolvedTriggerEfficiency
	process.p += process.ResolvedTriggerEfficiencyPFHT7504Jet
	process.p += process.ResolvedTriggerEfficiencyPFHT800
	process.p += process.ResolvedTriggerEfficiencyPFJet450
	process.p += process.ResolvedTriggerEfficiencyPFHT650WideJetMJJ900 
	process.p += process.BoostedTriggerEfficiency
	process.p += process.BoostedTriggerEfficiencyPFHT800
	process.p += process.BoostedTriggerEfficiencyAK8PFHT700TrimMass50
	process.p += process.BoostedTriggerEfficiencyAK8PFPt360TrimMass30
	process.p += process.BoostedTriggerEfficiencyAK8PFHT7504Jet
	process.p += process.BoostedTriggerEfficiencyAK8DiPFJet280220TrimMass30Btagp20
	process.p += process.BoostedTriggerEfficiencyAK8DiPFJet280220TrimMass30
	process.p += process.BoostedTriggerEfficiencyPFJet450
	process.p += process.BoostedTriggerEfficiencyPFHT650WideJetMJJ900 
	process.p += process.BoostedTriggerEfficiencySeveralTriggers
	if not options.miniAOD:
		process.p += process.BoostedTriggerEfficiencyPuppi
		process.p += process.BoostedTriggerEfficiencyPuppiPFHT800
		process.p += process.BoostedTriggerEfficiencyPuppiAK8PFHT700TrimMass50
		process.p += process.BoostedTriggerEfficiencyPuppiAK8PFPt360TrimMass30
		process.p += process.BoostedTriggerEfficiencyPuppiAK8PFHT7504Jet
		process.p += process.BoostedTriggerEfficiencyPuppiAK8DiPFJet280220TrimMass30Btagp20
		process.p += process.BoostedTriggerEfficiencyPuppiAK8DiPFJet280220TrimMass30
		process.p += process.BoostedTriggerEfficiencyPuppiPFJet450
		process.p += process.BoostedTriggerEfficiencyPuppiSeveralTriggers
	process.p += process.DijetTriggerEfficiency
	process.p += process.DijetTriggerEfficiencyPFHT800
	process.p += process.DijetTriggerEfficiencyAK8PFHT700TrimMass50
	process.p += process.DijetTriggerEfficiencyAK8PFPt360TrimMass30
	process.p += process.DijetTriggerEfficiencyAK8PFHT7504Jet
	process.p += process.DijetTriggerEfficiencyAK8DiPFJet280220TrimMass30Btagp20
	process.p += process.DijetTriggerEfficiencyAK8DiPFJet280220TrimMass30
	process.p += process.DijetTriggerEfficiencyPFHT650WideJetMJJ900
	process.p += process.DijetTriggerEfficiencyPFJet450
	process.p += process.DijetTriggerEfficiencySeveralTriggers

if options.miniAOD: outputNAME = outputNAME+'miniAOD'
process.TFileService=cms.Service("TFileService",fileName=cms.string( 'RUN'+outputNAME+'TriggerEfficiency_'+NAME+'.root' ) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
