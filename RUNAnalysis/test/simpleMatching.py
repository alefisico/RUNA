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
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/5C4E54A4-86F7-E611-95B5-549F35AD8BA2.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/00A5C142-C7F6-E611-8F03-002590A370B2.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/E21BE3AC-8CF7-E611-9E89-001E67E6F8C3.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/F4D36532-8DF7-E611-B310-002590E3A212.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/6E9031D5-FCF6-E611-8D0E-02163E019CBE.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/C8D24BCB-09F7-E611-80B7-02163E01A5EF.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/6E62D1F5-19F7-E611-BE19-02163E0145A8.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/AC0CCA16-36F7-E611-BB3E-02163E01231E.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/6ABB8236-28F7-E611-89F5-02163E019B49.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/E2C65F72-95F7-E611-AB40-02163E01A27A.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/EEA683DE-7DF7-E611-A25C-C81F66B7EBF5.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-100_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/3C2574FE-97F6-E611-87CA-FA163E0D94A1.root',

		#### 140
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/8C5D69EF-8CF7-E611-9016-02163E019D67.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/4865AEF4-74F7-E611-8BA1-D067E5F90F2A.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/C6A6B5AE-5EF6-E611-BF04-FA163E07656F.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/AEDDF290-29F7-E611-978B-FA163E77B33C.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/6A5FC54E-33F7-E611-A2AF-FA163EBE61B0.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/D495CC13-76F7-E611-9076-02163E017710.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/92FFB7CF-74F7-E611-B64C-A0369F7FC954.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/16F7776C-A7F6-E611-A71D-24BE05CE2EC1.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/B26818E1-74F7-E611-98AE-A0000420FE80.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/C8C9C8B2-6CF7-E611-B334-0CC47A7C34B0.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/8A82FCD2-75F7-E611-B00F-0025905A60A8.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/B04B1A8F-66F7-E611-B2A4-1866DAEA7A40.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/A4F278D7-74F7-E611-A404-0090FAA58224.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/88AC664A-71F7-E611-BD1F-001E67E6F7BA.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-140_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/B4AAE22B-6CF7-E611-92EB-0025901AC0F8.root',

		#### 200
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/4C4305F4-96F7-E611-8E03-FA163EEFC88A.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/F6AE23A6-A1F7-E611-8624-02163E013A0D.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/F4AEDF09-B1F7-E611-BB94-0CC47A4D76A0.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/06A7E30F-B1F7-E611-966A-0025905B85B6.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/56B3A943-A1F7-E611-AAB1-B083FECFF2BF.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/56A74F55-A0F7-E611-BB6D-549F35AE4F95.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/A8FFD447-A1F7-E611-8210-001E674FCAE9.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/D211C647-A1F7-E611-8AC9-002590D9D896.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-200_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/9C3C781C-49F7-E611-852B-02163E01A627.root',


		#### 300 UDD312
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/42CDB861-AAF7-E611-833B-0025901AC182.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/7A5103A1-CDF7-E611-907B-D4AE526DF2E3.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/406D0E0B-D2F6-E611-96EC-FA163E5F989F.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/6C1BC744-D8F6-E611-9BF6-02163E0148B6.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/A81926A1-3CF7-E611-A727-FA163E60EC40.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/2A13DF2A-CEF7-E611-B7B0-FA163E02238B.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/7806C635-CDF7-E611-95CD-FA163ECF00EE.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/3A8E6A7E-CDF7-E611-BE64-549F35AF44C9.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/EE46B28D-CDF7-E611-86CA-0242AC130002.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/A42BAA3F-34F7-E611-83C1-02163E0143A8.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/8EE6D70C-44F7-E611-A46D-02163E01A793.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/7E3132EB-57F7-E611-9C17-02163E013709.root',
		#### UDD323 300
		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD323_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/18975FE7-45D1-E611-9761-6CC2173BC1A0.root',
		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD323_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4AF776F7-7DD0-E611-9119-047D7BD6DD9A.root',
		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD323_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BCF42E07-46D1-E611-801D-047D7BD6DE5C.root',
		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD323_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5E8D9C57-43D1-E611-A536-848F69FD3D0D.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD323_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/A028ACC7-45D1-E611-AE77-00259029ED1A.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD323_M-300_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/06BB7BB9-45D1-E611-8B62-FA163EF17138.root',

		#### 500
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-500_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/3CC1EB8E-7AF7-E611-A101-FA163EE55D2F.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-500_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/9C2396EE-76F7-E611-BF7C-001E67E6F7BA.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-500_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/6E68BB71-71F7-E611-9067-02163E01A61B.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-500_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/C21D60F4-77F7-E611-B561-02163E01379D.root',
		#### 800
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/D059EB4C-C3F6-E611-9C7C-28924A35056E.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/5000BB60-C3F6-E611-85E4-D8D385AF8ACC.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/00662851-C3F6-E611-851D-848F69FD8933.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/F84A6E23-C4F6-E611-AD78-00266CFB991C.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/0E1FD325-C3F6-E611-9765-6CC2173C9150.root',
#		'/store/mc/RunIISummer16MiniAODv2/RPVStopStopToJets_UDD312_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/130000/B65B5F22-C4F6-E611-8EC3-0025905C54C6.root',
		) )
process.TFileService=cms.Service("TFileService", fileName=cms.string('simpleMatching_RPVStop323300.root'))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 50000 ) )
    
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

