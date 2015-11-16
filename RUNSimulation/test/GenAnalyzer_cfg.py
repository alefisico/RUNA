import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RPVSt100tojj_13TeV_pythia8_GENSIM_TESTHT500_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

'''
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:RPVSt200tojj_13TeV_AODSIM_test.root'
    )
)
'''

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('GenAnalyzer.root')
)

process.demo = cms.EDAnalyzer('GenAnalyzer',
		src = cms.InputTag( 'GenParticles' ),
		stop1Mass = cms.double( 100 ),
		stop2Mass = cms.double( 100 ),
		st1decay = cms.double( 1 ),
)


process.p = cms.Path(process.demo)