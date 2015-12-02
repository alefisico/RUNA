##################################################################
########   TO RUN THIS: python crab3_QCD.py
########   DO NOT DO: crab submit crab3_QCD.py
##################################################################

from CRABClient.UserUtilities import config
from httplib import HTTPException

config = config()

version = 'v01'

config.General.requestName = ''
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RUNFullAnalysis_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = ''
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.ignoreLocality = True
config.JobType.inputFiles = [ 'PileupData2015D_JSON_10-28-2015.root' ]

config.Site.storageSite = 'T3_US_FNALLPC'

def submit(config):
	try:
		crabCommand('submit', config = config)
	except HTTPException, hte:
		print 'Cannot execute commend'
		print hte.headers

if __name__ == '__main__':

	from CRABAPI.RawCommand import crabCommand

	Samples = [ 
	#		'/JetHT/algomez-Run2015D-PromptReco-v4_RUNA_v09-b28b60459a529d7878ced41cd64e24b0/USER',
#			'/JetHT/algomez-Run2015D-05Oct2015-v1_RUNA_v09-0a9621e584219a8f6f3188411a1374b3/USER',
#
#			'/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/jsomalwa-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v12-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/jsomalwa-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v12-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/jsomalwa-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v12-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/algomez-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/mmorris-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-7bcba442d2602d92a19f8ca61d13bc6c/USER',
#			'/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/mmorris-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-7bcba442d2602d92a19f8ca61d13bc6c/USER',
			'/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/mmorris-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-7bcba442d2602d92a19f8ca61d13bc6c/USER',
#
#			'/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/mmorris-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-7bcba442d2602d92a19f8ca61d13bc6c/USER',

#			'/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/dsheffie-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-6cd9a37acb7fba8686d9247b86713620/USER',

#			'/RPVSt100tojj_13TeV_pythia8/algomez-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-6cd9a37acb7fba8686d9247b86713620/USER',
#			'/RPVSt350tojj_13TeV_pythia8/algomez-RunIISpring15MiniAODv2-74X_RUNA_Asympt25ns_v09-6cd9a37acb7fba8686d9247b86713620/USER',


			]

	
	from multiprocessing import Process
	for dataset in Samples:
		config.Data.inputDataset = dataset
		procName = dataset.split('/')[1]+dataset.split('/')[2].replace('algomez-','').replace('RUNA','')+'_'+version
		if 'JetHT' in dataset: 
			#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
			config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt' 
		config.General.requestName = procName
		config.JobType.pyCfgParams = [ 'PROC='+procName, 'local=0' ]
		p = Process(target=submit, args=(config,))
		p.start()
		p.join()
