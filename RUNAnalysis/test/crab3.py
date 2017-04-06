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
config.JobType.psetName = 'RUNFullAnalysis_cfg.py'
#config.JobType.psetName = 'RUNBoostedResolutionCalc_cfg.py'
#config.JobType.psetName = 'RUNResolvedResolutionCalc_cfg.py'
#config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = ''
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
#config.Data.ignoreLocality = True

#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_US_Rutgers'
config.Data.outLFNDirBase = '/store/user/algomez/myArea/EOS/B2GAnaFW_80X_V2p3/'
#config.Data.outLFNDirBase = '/store/user/algomez/myArea/EOS/ResolvedCalc/'

def submit(config):
	try:
		if args.dryrun: crabCommand('submit', '--dryrun', config = config)
		else: crabCommand('submit', config = config)
	except HTTPException, hte:
		print 'Cannot execute commend'
		print hte.headers



#######################################################################################
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--sample', action='store', default='all', dest='sample', help='Sample to process. Example: QCD, RPV, TTJets.' )
	parser.add_argument('-v', '--version', action='store', default='v01p0', dest='version', help='Version: v01, v02.' )
	parser.add_argument('-d', '--dryrun', action='store_true', default=False, help='To run dryrun crab mode.' )

	try: args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	Samples = {}

	Samples[ 'JetHTB' ] = [ '/JetHT/eschmitz-Run2016B-23Sep2016-v3_B2GAnaFW_80X_V2p3-f7a2f235ab93c0f14514c8048ff181ad/USER', 10, 'Spring16_23Sep2016BCDV2', '' ]
	Samples[ 'JetHTC' ] = [ '/JetHT/eschmitz-Run2016C-23Sep2016-v1_B2GAnaFW_80X_V2p3-f7a2f235ab93c0f14514c8048ff181ad/USER', 10, 'Spring16_23Sep2016BCDV2', '' ]
	Samples[ 'JetHTD' ] = [ '/JetHT/eschmitz-Run2016D-23Sep2016-v1_B2GAnaFW_80X_V2p3-f7a2f235ab93c0f14514c8048ff181ad/USER', 10, 'Spring16_23Sep2016BCDV2', '' ]
	Samples[ 'JetHTE' ] = [ '/JetHT/eschmitz-Run2016E-23Sep2016-v1_B2GAnaFW_80X_V2p3-be8c3a8c93082fc7e9e68ba6c89c878a/USER', 10, 'Spring16_23Sep2016EFV2', '' ]
	Samples[ 'JetHTF1' ] = [ '/JetHT/eschmitz-Run2016F-23Sep2016-v1_B2GAnaFW_80X_V2p3-5505a629280c3916514864e51c616eff/USER', 10, 'Spring16_23Sep2016EFV2', '_F1' ]
	Samples[ 'JetHTF2' ] = [ '/JetHT/eschmitz-Run2016F-23Sep2016-v1_B2GAnaFW_80X_V2p3-5505a629280c3916514864e51c616eff/USER', 10, 'Spring16_23Sep2016GV2', '_F2' ]
	Samples[ 'JetHTG' ] = [ '/JetHT/eschmitz-Run2016G-23Sep2016-v1_B2GAnaFW_80X_V2p3-63fc13854a838efc75945a252e802fa8/USER', 10, 'Spring16_23Sep2016GV2', '' ]
	Samples[ 'JetHTH2' ] = [ '/JetHT/eschmitz-Run2016H-PromptReco-v2_B2GAnaFW_80X_V2p3-2b4e6db5c7765185184b3dfd7d730aea/USER', 10, 'Spring16_23Sep2016HV2', '' ]
	Samples[ 'JetHTH3' ] = [ '/JetHT/eschmitz-Run2016H-PromptReco-v3_B2GAnaFW_80X_V2p3-2b4e6db5c7765185184b3dfd7d730aea/USER', 10, 'Spring16_23Sep2016HV2', '' ]


	############### Moriond17 samples
	Samples[ 'QCDHT300' ] = [ '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 40, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDHT500' ] = [ '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 40, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDHT700' ] = [ '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 40, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDHT1000' ] = [ '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 30, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDHT1500' ] = [ '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 30, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDHT2000' ] = [ '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 30, 'Spring16_23Sep2016V2' ]
	
	#Samples[ 'QCDPt120' ] = [ '', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt170' ] = [ '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p04-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt300' ] = [ '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 20, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt470' ] = [ '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p04-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 20, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt600' ] = [ '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 20, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt800' ] = [ '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 20, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt1000' ] = [ '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt1400' ] = [ '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt1800' ] = [ '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt2400' ] = [ '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'QCDPt3200' ] = [ '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 10, 'Spring16_23Sep2016V2' ]

	############### Moriond17 samples
	Samples[ 'QCDHerwig' ] = [ '/QCD_Pt-15to7000_TuneCUETHS1_FlatP6_13TeV_herwigpp/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 10, 'Spring16_23Sep2016V2' ]

	############### Moriond17 samples
	Samples[ 'TTJets' ] = [ '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 100, 'Spring16_23Sep2016V2' ]
	Samples[ 'WJets' ] = [ '/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 10, 'Spring16_23Sep2016V2' ]
	#Samples[ 'ZJets' ] = [ '', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'WZ' ] = [ '/WZ_TuneCUETP8M1_13TeV-pythia8/vorobiev-B2GAnaFW_Spring16MiniAODv2_Moriond17_v80x_ext1_v2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'WWJets' ] = [ '/WWTo4Q_13TeV-powheg/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-9705383e86c3aabba8a30865b5958dd7/USER', 50, 'Spring16_23Sep2016V2' ]
	Samples[ 'ZZJets' ] = [ '/ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 100, 'Spring16_23Sep2016V2' ]
	
#	Samples[ 'RPV80' ] = [ '/RPVStopStopToJets_UDD312_M-80_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
#	Samples[ 'RPV100' ] = [ '/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
#	Samples[ 'RPV120' ] = [ '/RPVStopStopToJets_UDD312_M-120_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
#	Samples[ 'RPV140' ] = [ '/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
#	Samples[ 'RPV180' ] = [ '/RPVStopStopToJets_UDD312_M-180_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV240' ] = [ '/RPVStopStopToJets_UDD312_M-240_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV300' ] = [ '/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV350' ] = [ '/RPVStopStopToJets_UDD312_M-350_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV400' ] = [ '/RPVStopStopToJets_UDD312_M-400_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV450' ] = [ '/RPVStopStopToJets_UDD312_M-450_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV500' ] = [ '/RPVStopStopToJets_UDD312_M-500_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV550' ] = [ '/RPVStopStopToJets_UDD312_M-550_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV600' ] = [ '/RPVStopStopToJets_UDD312_M-600_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV650' ] = [ '/RPVStopStopToJets_UDD312_M-650_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV700' ] = [ '/RPVStopStopToJets_UDD312_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV750' ] = [ '/RPVStopStopToJets_UDD312_M-750_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV800' ] = [ '/RPVStopStopToJets_UDD312_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV850' ] = [ '/RPVStopStopToJets_UDD312_M-850_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2_B2GAnaFW_80X_V2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV900' ] = [ '/RPVStopStopToJets_UDD312_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV950' ] = [ '/RPVStopStopToJets_UDD312_M-950_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2_B2GAnaFW_80X_V2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV1000' ] = [ '/RPVStopStopToJets_UDD312_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2_B2GAnaFW_80X_V2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV1100' ] = [ '/RPVStopStopToJets_UDD312_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV1200' ] = [ '/RPVStopStopToJets_UDD312_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV1300' ] = [ '/RPVStopStopToJets_UDD312_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1-3435c96292394b4f7f0ff1725e305a15/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV1400' ] = [ '/RPVStopStopToJets_UDD312_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2_B2GAnaFW_80X_V2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV1500' ] = [ '/RPVStopStopToJets_UDD312_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/jsomalwa-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2_B2GAnaFW_80X_V2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER', 10, 'Spring16_23Sep2016V2' ]

	Samples[ 'RPV323120' ] = [ '/RPVStopStopToJets_UDD323_M-120_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323180' ] = [ '/RPVStopStopToJets_UDD323_M-180_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323200' ] = [ '/RPVStopStopToJets_UDD323_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323220' ] = [ '/RPVStopStopToJets_UDD323_M-220_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323240' ] = [ '/RPVStopStopToJets_UDD323_M-240_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323280' ] = [ '/RPVStopStopToJets_UDD323_M-280_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323300' ] = [ '/RPVStopStopToJets_UDD323_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323350' ] = [ '/RPVStopStopToJets_UDD323_M-350_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323500' ] = [ '/RPVStopStopToJets_UDD323_M-500_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323600' ] = [ '/RPVStopStopToJets_UDD323_M-600_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323650' ] = [ '/RPVStopStopToJets_UDD323_M-650_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323700' ] = [ '/RPVStopStopToJets_UDD323_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323750' ] = [ '/RPVStopStopToJets_UDD323_M-750_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323800' ] = [ '/RPVStopStopToJets_UDD323_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV323950' ] = [ '/RPVStopStopToJets_UDD323_M-950_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV3231000' ] = [ '/RPVStopStopToJets_UDD323_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV3231100' ] = [ '/RPVStopStopToJets_UDD323_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV3231200' ] = [ '/RPVStopStopToJets_UDD323_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]
	Samples[ 'RPV3231300' ] = [ '/RPVStopStopToJets_UDD323_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/algomez-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4-6b29e1707fe76ab19c1ba543e7f6f24b/USER', 10, 'Spring16_23Sep2016V2' ]



	processingSamples = {}
	if 'all' in args.sample: 
		for sam in Samples: processingSamples[ sam ] = Samples[ sam ]
	else:
		for sam in Samples: 
			if sam.startswith( args.sample ): processingSamples[ sam ] = Samples[ sam ]

	if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'
		
	for sam in processingSamples:
		dataset = processingSamples[sam][0]
		jecVersion = processingSamples[sam][2]
		supportFiles = glob.glob('supportFiles/'+jecVersion+'_*txt') + glob.glob('supportFiles/*csv') + glob.glob('supportFiles/*root')
		config.Data.inputDataset = dataset
		config.Data.unitsPerJob = processingSamples[sam][1]
		if 'JetHT' in dataset: 
			procName = dataset.split('/')[1]+dataset.split('/')[2].replace( dataset.split('/')[2].split('-')[0], '').split('_')[0]+processingSamples[sam][3]+'_'+args.version
			config.Data.lumiMask = '/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNAnalysis/test/supportFiles/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON'+processingSamples[sam][3]+'.txt'
			#config.Data.lumiMask = 'crab_projects/crab_JetHT-Run2016G-23Sep2016-v1_v06p1/results/notFinishedLumis.json'
			config.JobType.pyCfgParams = [ 'PROC='+procName, 'jecVersion='+jecVersion ] 

		else:
			procName = dataset.split('/')[1].split('_TuneCUETP8M1')[0]
			config.JobType.pyCfgParams = ( [ 'PROC='+procName, 'systematics='+('1' if 'RPV' in sam else '0'), 'jecVersion='+jecVersion ] )
			#config.JobType.pyCfgParams = ( [ 'PROC='+procName, 'jecVersion='+jecVersion ] )
			#config.Data.lumiMask = 'crab_projects/crab_QCD_Pt_1800to2400v06/results/notFinishedLumis.json'

		config.JobType.inputFiles =  supportFiles
		config.General.requestName = procName+args.version
		print config
		print '|--- Submmiting sample: ', procName
		p = Process(target=submit, args=(config,))
		p.start()
		p.join()
