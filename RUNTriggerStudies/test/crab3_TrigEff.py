##################################################################
########   TO RUN THIS: python crab3_QCD.py
########   DO NOT DO: crab submit crab3_QCD.py
##################################################################

from CRABClient.UserUtilities import config
import argparse, sys
from httplib import HTTPException
from CRABAPI.RawCommand import crabCommand
from multiprocessing import Process
import glob


config = config()

config.General.requestName = ''
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RUNTriggerEfficiency_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = ''
config.Data.outLFNDirBase = '/store/user/algomez/'
config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10

#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_US_Rutgers'
config.Data.outLFNDirBase = '/store/user/algomez/myArea/EOS/TriggerEfficiency/'


def submit(config):
	try:
		crabCommand('submit', config = config)
	except HTTPException, hte:
		print 'Cannot execute commend'
		print hte.headers

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--sample', action='store', default='all', dest='sample', help='Sample to process. Example: QCD, RPV, TTJets.' )
	parser.add_argument('-v', '--version', action='store', default='v01p0', dest='version', help='Version: v01, v02.' )

	try: args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	Samples = {}

	Samples[ 'B2GJetHTB' ] = [ '/JetHT/algomez-Run2016B-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 10, 'Spring16_23Sep2016BCDV2BCD' ]
	Samples[ 'B2GJetHTC' ] = [ '/JetHT/algomez-Run2016C-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 10, 'Spring16_23Sep2016BCDV2BCD' ]
	Samples[ 'B2GJetHTD' ] = [ '/JetHT/algomez-Run2016D-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 10, 'Spring16_23Sep2016BCDV2BCD' ]
	Samples[ 'B2GJetHTE' ] = [ '/JetHT/algomez-Run2016E-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 10, 'Spring16_23Sep2016BCDV2E' ]
	Samples[ 'B2GJetHTF' ] = [ '/JetHT/algomez-Run2016F-PromptReco-v1_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 10, 'Spring16_23Sep2016BCDV2F' ]
	Samples[ 'B2GJetHTG' ] = [ '/JetHT/algomez-Run2016G-PromptReco-v1_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 10, 'Spring16_23Sep2016BCDV2p2' ]
	Samples[ 'B2GJetHTH' ] = [ '/JetHT/algomez-Run2016H-PromptReco-v2_B2GAnaFW_80X_V2p01p1-3e507461fd667ac6961fa4af5b123b09/USER', 10, 'Spring16_23Sep2016BCDV2p2' ]
	Samples[ 'miniAODJetHTB' ] = [ '/JetHT/Run2016B-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODJetHTC' ] = [ '/JetHT/Run2016C-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODJetHTD' ] = [ '/JetHT/Run2016D-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODJetHTE' ] = [ '/JetHT/Run2016E-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODJetHTF' ] = [ '/JetHT/Run2016F-PromptReco-v1/MINIAOD', 10 ]
	Samples[ 'miniAODJetHTG' ] = [ '/JetHT/Run2016G-PromptReco-v1/MINIAOD', 10 ]
	Samples[ 'miniAODJetHTH' ] = [ '/JetHT/Run2016H-PromptReco-v2/MINIAOD', 10 ]

	Samples[ 'B2GSingleMuonB1' ] = [ '/SingleMuon/vorobiev-Run2016B-23Sep2016-v1_B2GAnaFW_80X_v2p4-6f92d8e9b2717da1daa65a0e07f84d5f/USER', 10, 'Spring16_23Sep2016BCDV2', '' ]
	Samples[ 'B2GSingleMuonB3' ] = [ '/SingleMuon/vorobiev-Run2016B-23Sep2016-v3_B2GAnaFW_80X_v2p4-6f92d8e9b2717da1daa65a0e07f84d5f/USER', 10, 'Spring16_23Sep2016BCDV2', '' ]
	Samples[ 'B2GSingleMuonC' ] = [ '/SingleMuon/vorobiev-Run2016C-23Sep2016-v1_B2GAnaFW_80X_v2p4-6f92d8e9b2717da1daa65a0e07f84d5f/USER', 10, 'Spring16_23Sep2016BCDV2', '' ]
	Samples[ 'B2GSingleMuonD' ] = [ '/SingleMuon/vorobiev-Run2016D-23Sep2016-v1_B2GAnaFW_80X_v2p4-6f92d8e9b2717da1daa65a0e07f84d5f/USER', 10, 'Spring16_23Sep2016BCDV2', '' ]
	Samples[ 'B2GSingleMuonE' ] = [ '/SingleMuon/vorobiev-Run2016E-23Sep2016-v1_B2GAnaFW_80X_v2p4-961c7d882d8721e72fac616aaa90ecc1/USER', 10, 'Spring16_23Sep2016EFV2', '' ]
	Samples[ 'B2GSingleMuonF1' ] = [ '/SingleMuon/vorobiev-Run2016F-23Sep2016-v1_B2GAnaFW_80X_v2p4-e26b2444814d18badce3899570108664/USER', 10, 'Spring16_23Sep2016EFV2', '_F1' ]
	Samples[ 'B2GSingleMuonF2' ] = [ '/SingleMuon/vorobiev-Run2016F-23Sep2016-v1_B2GAnaFW_80X_v2p4-e26b2444814d18badce3899570108664/USER', 10, 'Spring16_23Sep2016GV2', '_F2' ]
	Samples[ 'B2GSingleMuonG' ] = [ '/SingleMuon/vorobiev-Run2016G-23Sep2016-v1_B2GAnaFW_80X_v2p4-b45603bfa955d854bdcb9af322c0b037/USER', 10, 'Spring16_23Sep2016GV2', '' ]
	Samples[ 'B2GSingleMuonH1' ] = [ '/SingleMuon/vorobiev-Run2016H-PromptReco-v1_B2GAnaFW_80X_v2p4-376a23645e94877b22a7f32873431514/USER', 10, 'Spring16_23Sep2016HV2', '' ]
	Samples[ 'B2GSingleMuonH2' ] = [ '/SingleMuon/vorobiev-Run2016H-PromptReco-v2_B2GAnaFW_80X_v2p4-376a23645e94877b22a7f32873431514/USER', 10, 'Spring16_23Sep2016HV2', '' ]
	Samples[ 'B2GSingleMuonH3' ] = [ '/SingleMuon/vorobiev-Run2016H-PromptReco-v3_B2GAnaFW_80X_v2p4-376a23645e94877b22a7f32873431514/USER', 10, 'Spring16_23Sep2016HV2', '' ]

	Samples[ 'miniAODSingleMuonB' ] = [ '/SingleMuon/Run2016B-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODSingleMuonC' ] = [ '/SingleMuon/Run2016C-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODSingleMuonD' ] = [ '/SingleMuon/Run2016D-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODSingleMuonE' ] = [ '/SingleMuon/Run2016E-PromptReco-v2/MINIAOD', 10 ]
	Samples[ 'miniAODSingleMuonF' ] = [ '/SingleMuon/Run2016F-PromptReco-v1/MINIAOD', 10 ]
	Samples[ 'miniAODSingleMuonG' ] = [ '/SingleMuon/Run2016G-PromptReco-v1/MINIAOD', 10 ]
	Samples[ 'miniAODSingleMuonH' ] = [ '/SingleMuon/Run2016H-PromptReco-v2/MINIAOD', 10 ]

	processingSamples = {}
	if 'all' in args.sample: 
		for sam in Samples: processingSamples[ sam ] = Samples[ sam ]
	else:
		for sam in Samples: 
			if sam.startswith( args.sample ): processingSamples[ sam ] = Samples[ sam ]

	if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'
		
	for sam in processingSamples:
		dataset = processingSamples[sam][0]
		procName = dataset.split('/')[1]+dataset.split('/')[2].replace( dataset.split('/')[2].split('-')[0], '').split('_')[0]+processingSamples[sam][3]+'_'+args.version
		#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
		config.Data.lumiMask = '/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNAnalysis/test/supportFiles/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON'+processingSamples[sam][3]+'.txt'

		config.Data.inputDataset = dataset
		config.General.requestName = procName
		if 'B2G' in args.sample:
			config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader'
			supportFiles = glob.glob('../../RUNAnalysis/test/supportFiles/'+processingSamples[sam][2]+'*txt')
			config.JobType.inputFiles = supportFiles

		config.JobType.pyCfgParams = [ 'PROC='+procName, ( 'miniAOD=True' if 'miniAOD' in args.sample else 'jecVersion='+processingSamples[sam][2]) ]
		#crabCommand('submit', config = config)
		p = Process(target=submit, args=(config,))
		p.start()
		p.join()
