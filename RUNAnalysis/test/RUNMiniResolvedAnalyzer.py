#!/usr/bin/env python

'''
File: MyAnalyzer.py --mass 50 (optional) --debug -- final --jetAlgo AK5
Author: Alejandro Gomez Espinosa
Email: gomez@physics.rutgers.edu
Description: My Analyzer 
'''

import sys,os,time
import argparse
from collections import OrderedDict
from multiprocessing import Process
from ROOT import *
from array import array
import numpy as np
from random import randint
try: 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
	from RUNA.RUNAnalysis.commonFunctions import *
except ImportError: 
	sys.path.append('../python') 
	import tdrstyle as tdrstyle
	from commonFunctions import *

gROOT.SetBatch()

######################################
def myAnalyzer( dictSamples, preselection, cuts, signalName, UNC ):


	#outputFileName = 'Rootfiles/RUNMiniScoutingResolvedAnalysis_'+signalName+UNC+'_'+( '' if 'JetHT' in signalName else '80X_')+'V2p1_'+args.version+'p1.root' 
	outputFileName = 'Rootfiles/RUNMiniResolvedAnalysis_'+signalName+UNC+'_'+( '' if 'JetHT' in signalName else '80X_')+'V2p1_'+args.version+'p1.root' 
	outputFile = TFile( outputFileName, 'RECREATE' )


	################################################################################################## Trigger Histos
	nBinsMass	= 200
	maxMass		= 2000
	nBinsHT		= 150
	maxHT		= 1500

	for sam in dictSamples:
		allHistos[ "HT_cutBestPair_"+sam ] = TH1F( "HT_cutBestPair_"+sam, "HT_cutBestPair_"+sam, 5000, 0., 5000 )
		allHistos[ "NPV_cutBestPair_"+sam ] = TH1F( "NPV_cutBestPair_"+sam, "NPV_cutBestPair_"+sam, 100, 0., 100 )
		allHistos[ "jet1Pt_cutBestPair_"+sam ] = TH1F( "jet1Pt_cutBestPair_"+sam, "jet1Pt_cutBestPair_"+sam, 3000, 0., 3000 )
		allHistos[ "jet2Pt_cutBestPair_"+sam ] = TH1F( "jet2Pt_cutBestPair_"+sam, "jet2Pt_cutBestPair_"+sam, 3000, 0., 3000 )
		allHistos[ "jet3Pt_cutBestPair_"+sam ] = TH1F( "jet3Pt_cutBestPair_"+sam, "jet3Pt_cutBestPair_"+sam, 3000, 0., 3000 )
		allHistos[ "jet4Pt_cutBestPair_"+sam ] = TH1F( "jet4Pt_cutBestPair_"+sam, "jet4Pt_cutBestPair_"+sam, 3000, 0., 3000 )
		allHistos[ "massAve_cutBestPair_"+sam ] = TH1F( "massAve_cutBestPair_"+sam, "massAve_cutBestPair_"+sam, 3000, 0., 3000 )
		allHistos[ "massAsym_cutBestPair_"+sam ] = TH1F( "massAsym_cutBestPair_"+sam, "massAsym_cutBestPair_"+sam, 20, 0., 1 )
		allHistos[ "deltaEta_cutBestPair_"+sam ] = TH1F( "deltaEta_cutBestPair_"+sam, "deltaEta_cutBestPair_"+sam, 50, 0., 5 )
		allHistos[ 'deltavsMassAve_cutBestPair_'+sam ] = TH2F( 'deltavsMassAvecut_BestPair_'+sam, 'deltavsMassAvecut_BestPair_'+sam, 1000, 0., 1000, 1000, 0., 1000. )

		allHistos[ "massAve_delta_"+sam ] = TH1F( "massAve_delta_"+sam, "massAve_delta_"+sam, 3000, 0., 3000 )
		allHistos[ "HT_delta_"+sam ] = TH1F( "HT_delta_"+sam, "HT_delta_"+sam, 5000, 0., 5000 )
		allHistos[ "massAve_delta50_"+sam ] = TH1F( "massAve_delta50_"+sam, "massAve_delta50_"+sam, 3000, 0., 3000 )
		allHistos[ "massAve_delta100_"+sam ] = TH1F( "massAve_delta100_"+sam, "massAve_delta100_"+sam, 3000, 0., 3000 )
		allHistos[ "massAve_delta150_"+sam ] = TH1F( "massAve_delta150_"+sam, "massAve_delta150_"+sam, 3000, 0., 3000 )
		allHistos[ "massAve_delta200_"+sam ] = TH1F( "massAve_delta200_"+sam, "massAve_delta200_"+sam, 3000, 0., 3000 )
		allHistos[ "massAve_delta250_"+sam ] = TH1F( "massAve_delta250_"+sam, "massAve_delta250_"+sam, 3000, 0., 3000 )

		allHistos[ "massAsym_n-1_"+sam ] = TH1F( "massAsym_n-1_"+sam, "massAsym_n-1_"+sam, 20, 0., 1 )
		allHistos[ "deltaEta_n-1_"+sam ] = TH1F( "deltaEta_n-1_"+sam, "deltaEta_n-1_"+sam, 50, 0., 5 )
		allHistos[ 'deltavsMassAve_n-1_'+sam ] = TH2F( 'deltavsMassAve_n-1_'+sam, 'deltavsMassAve_n-1_'+sam, 1000, 0., 1000, 1000, 0., 1000. )

		allHistos[ "massAve_delta_1qgl_"+sam ] = TH1F( "massAve_delta_1qgl_"+sam, "massAve_delta_1qgl_"+sam, 3000, 0., 3000 )
		allHistos[ "massAve_delta_2qgl_"+sam ] = TH1F( "massAve_delta_2qgl_"+sam, "massAve_delta_2qgl_"+sam, 3000, 0., 3000 )
		allHistos[ "jet1QGL_delta_"+sam ] = TH1F( "jet1QGL_delta_"+sam, "jet1QGL_delta_"+sam, 20, 0., 1 )
		allHistos[ "jet2QGL_delta_"+sam ] = TH1F( "jet2QGL_delta_"+sam, "jet2QGL_delta_"+sam, 20, 0., 1 )
		allHistos[ "jet3QGL_delta_"+sam ] = TH1F( "jet3QGL_delta_"+sam, "jet3QGL_delta_"+sam, 20, 0., 1 )
		allHistos[ "jet4QGL_delta_"+sam ] = TH1F( "jet4QGL_delta_"+sam, "jet4QGL_delta_"+sam, 20, 0., 1 )

		allHistos[ "massAve_delta_4qgl_"+sam ] = TH1F( "massAve_delta_4qgl_"+sam, "massAve_delta_4qgl_"+sam, 3000, 0., 3000 )

		allHistos[ "massAve_delta_1btag_"+sam ] = TH1F( "massAve_delta_1btag_"+sam, "massAve_delta_1btag_"+sam, 3000, 0., 3000 )
		allHistos[ "massAve_delta_2btag_"+sam ] = TH1F( "massAve_delta_2btag_"+sam, "massAve_delta_2btag_"+sam, 3000, 0., 3000 )
		allHistos[ "jet1Btag_delta_"+sam ] = TH1F( "jet1Btag_delta_"+sam, "jet1Btag_delta_"+sam, 20, 0., 1 )
		allHistos[ "jet2Btag_delta_"+sam ] = TH1F( "jet2Btag_delta_"+sam, "jet2Btag_delta_"+sam, 20, 0., 1 )
		allHistos[ "jet3Btag_delta_"+sam ] = TH1F( "jet3Btag_delta_"+sam, "jet3Btag_delta_"+sam, 20, 0., 1 )
		allHistos[ "jet4Btag_delta_"+sam ] = TH1F( "jet4Btag_delta_"+sam, "jet4Btag_delta_"+sam, 20, 0., 1 )

	for h in allHistos: allHistos[h].Sumw2()

	################################################################################################## Running the Analysis
	SF = 'lumiWeight*puWeight'
	presel = TCut( SF ) *  TCut( preselection )
	fullSel = TCut( SF ) * TCut( preselection + ' && ' + cuts ) 
	print '-'*40
	#print '---- Cuts applied: ', fullSel

	#treeName = 'ResolvedAnalysisPlotsScouting/RUNATree'
	treeName = 'ResolvedAnalysisPlots/RUNATree'

	if 'JetHT' in args.samples: 
		for era in [ 'B', 'C', 'D', 'E', 'F', 'G', 'H' ]:
			dictSamples[ 'JetHT_Run2016'+era ] = dictSamples[ 'JetHT_Run2016' ].replace('2016', '2016'+era)
		dictSamples.pop( 'JetHT_Run2016' )
		fullSel = fullSel + TCut('Entry$ % '+str(randint(0,10))+' == 0')

	for sample in dictSamples:

#
		fileSample = dictSamples[ sample ]
		print '--- Sample ', sample
		if 'JetHT' in sample: sample = 'JetHT_Run2016'

		### All selection
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				fullSel,
				allHistos[ 'massAve_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) ) 

		getHistoFromTree( fileSample, treeName,
				'HT', 
				fullSel,
				allHistos[ 'HT_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )


		#### preselection plots
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				presel,
				allHistos[ 'massAve_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'numPV', 
				presel,
				allHistos[ 'NPV_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'HT', 
				presel,
				allHistos[ 'HT_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsPt[0]', 
				presel,
				allHistos[ 'jet1Pt_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsPt[1]', 
				presel,
				allHistos[ 'jet2Pt_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsPt[2]', 
				presel,
				allHistos[ 'jet3Pt_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsPt[3]', 
				presel,
				allHistos[ 'jet4Pt_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		get2DHistoFromTree( fileSample, treeName,
				'massAve', 'delta1',
				presel,
				allHistos[ 'deltavsMassAve_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )
		get2DHistoFromTree( fileSample, treeName,
				'massAve', 'delta2',
				presel,
				allHistos[ 'deltavsMassAve_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'massAsym',
				presel,
				allHistos[ 'massAsym_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'deltaEta',
				presel,
				allHistos[ 'deltaEta_cutBestPair_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		##### checking diff deltas
		getHistoFromTree( fileSample, treeName,
				'massAve',
				presel + TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>50) && (delta2>50)') ),
				allHistos[ 'massAve_delta50_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'massAve',
				presel + TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>100) && (delta2>100)') ),
				allHistos[ 'massAve_delta100_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'massAve',
				presel + TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>150) && (delta2>150)') ),
				allHistos[ 'massAve_delta150_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'massAve',
				fullSel,
				allHistos[ 'massAve_delta200_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'massAve',
				presel + TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>250) && (delta2>250)') ),
				allHistos[ 'massAve_delta250_'+sample ], 
				( True if 'JetHT' in sample else False ) )



		#### n-1 plots
		get2DHistoFromTree( fileSample, treeName,
				'massAve', 'delta1', 
				presel + TCut( cuts.replace('&& (delta1>200)','').replace('&& (delta2>200)','') ), 
				#presel + TCut( cuts.replace('&& (delta1>100)','').replace('&& (delta2>100)',''), 
				allHistos[ 'deltavsMassAve_n-1_'+sample ], 
				( True if 'JetHT' in sample else False ) )
		get2DHistoFromTree( fileSample, treeName,
				'massAve', 'delta2', 
				#presel + TCut( cuts.replace('&& (delta1>100)','').replace('&& (delta2>100)',''), 
				presel + TCut( cuts.replace('&& (delta1>200)','').replace('&& (delta2>200)','') ), 
				allHistos[ 'deltavsMassAve_n-1_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'massAsym',
				#presel + TCut( cuts.replace('&& (massAsym<0.2)',''), 
				presel + TCut( cuts.replace('&& (massAsym<0.1)','') ), 
				allHistos[ 'massAsym_n-1_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'deltaEta',
				#presel + TCut( cuts.replace('&& (deltaEta<1.2)',''), 
				presel + TCut( cuts.replace('&& (deltaEta<1.)','') ), 
				allHistos[ 'deltaEta_n-1_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		### Test QGL
		oneQGL = TCut( SF ) * TCut( preselection + ' && ' + cuts + ' && ( (jetsQGL[0]>0.5) || (jetsQGL[1]>0.5) || (jetsQGL[2]>0.5) || (jetsQGL[3]>0.5) )')
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				oneQGL, 
				allHistos[ 'massAve_delta_1qgl_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		twoQGL = TCut( SF ) * TCut( preselection + ' && ' + cuts + ' && ( (jetsQGL[0]>0.5) || (jetsQGL[1]>0.5)) && ((jetsQGL[2]>0.5) || (jetsQGL[3]>0.5) )')
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				twoQGL, 
				allHistos[ 'massAve_delta_2qgl_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		fourQGL = TCut( SF ) * TCut( preselection + ' && ' + cuts + ' && (jetsQGL[0]>0.5) && (jetsQGL[1]>0.5) && (jetsQGL[2]>0.5) && (jetsQGL[3]>0.5)')
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				fourQGL, 
				allHistos[ 'massAve_delta_4qgl_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsQGL[0]', 
				fullSel,
				allHistos[ 'jet1QGL_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsQGL[1]', 
				fullSel,
				allHistos[ 'jet2QGL_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsQGL[2]', 
				fullSel,
				allHistos[ 'jet3QGL_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsQGL[3]', 
				fullSel,
				allHistos[ 'jet4QGL_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		### Test Btag
		'''
		oneBtag = TCut('( (jetsBtag[0]>0.5) || (jetsBtag[1]>0.5) || (jetsBtag[2]>0.5) || (jetsBtag[3]>0.5) )')
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				fullSel+oneBtag, 
				allHistos[ 'massAve_delta_1qgl_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		twoBtag = TCut('( (jetsBtag[0]>0.5) || (jetsBtag[1]>0.5)) && ((jetsBtag[2]>0.5) || (jetsBtag[3]>0.5) )')
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				fullSel+twoBtag, 
				allHistos[ 'massAve_delta_2qgl_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		fourBtag = TCut('(jetsBtag[0]>0.5) && (jetsBtag[1]>0.5) && (jetsBtag[2]>0.5) && (jetsBtag[3]>0.5)')
		getHistoFromTree( fileSample, treeName,
				'massAve', 
				fullSel+twoBtag, 
				allHistos[ 'massAve_delta_4qgl_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsBtag[0]', 
				fullSel,
				allHistos[ 'jet1Btag_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsBtag[1]', 
				fullSel,
				allHistos[ 'jet2Btag_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsBtag[2]', 
				fullSel,
				allHistos[ 'jet3Btag_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )

		getHistoFromTree( fileSample, treeName,
				'jetsBtag[3]', 
				fullSel,
				allHistos[ 'jet4Btag_delta_'+sample ], 
				( True if 'JetHT' in sample else False ) )
		'''


	outputFile.Write()

	##### Closing
	print 'Writing output file: '+ outputFileName
	outputFile.Close()


def assignDijets( tmpj1, tmpj2, tmpj3, tmpj4, minVarPairing ):
	"""docstring for assignDijets"""

	Pairing = False
	j1 = TLorentzVector()
	j2 = TLorentzVector()
	j3 = TLorentzVector()
	j4 = TLorentzVector()
	if( '1234' in minVarPairing ):
		if ( tmpj1.DeltaR(tmpj2) > ( tmpj3.DeltaR(tmpj4) ) ):
			j1 = tmpj1
			j2 = tmpj2
			j3 = tmpj3
			j4 = tmpj4
		else:
			j1 = tmpj3
			j2 = tmpj4
			j3 = tmpj1
			j4 = tmpj2
		Pairing = True
	elif( '1324' in minVarPairing ):
		if ( tmpj1.DeltaR(tmpj3) > ( tmpj2.DeltaR(tmpj4) ) ):
			j1 = tmpj1
			j2 = tmpj3
			j3 = tmpj2
			j4 = tmpj4
		else:
			j1 = tmpj2
			j2 = tmpj4
			j3 = tmpj1
			j4 = tmpj3
		Pairing = True
	elif( '1423' in minVarPairing ):
		if ( tmpj1.DeltaR(tmpj4) > ( tmpj2.DeltaR(tmpj3) ) ):
			j1 = tmpj1
			j2 = tmpj4
			j3 = tmpj2
			j4 = tmpj3
		else:
			j1 = tmpj2
			j2 = tmpj3
			j3 = tmpj1
			j4 = tmpj4
		Pairing = True

	return [ Pairing, j1, j2, j3, j4 ]

def DeltaRPairing( j1, j2, j3, j4, offset ):
	"""docstring for DeltaRPairing"""

	mindRDeltaRPairing = {}
	mindRDeltaRPairing[ '1234' ] = ( abs( j1.DeltaR(j2) - offset ) + abs( j3.DeltaR(j4) - offset ) )
	mindRDeltaRPairing[ '1324' ] = ( abs( j1.DeltaR(j3) - offset ) + abs( j2.DeltaR(j4) - offset ) )
	mindRDeltaRPairing[ '1423' ] = ( abs( j1.DeltaR(j4) - offset ) + abs( j2.DeltaR(j3) - offset ) )
	minDeltaRPairing = min(mindRDeltaRPairing, key=mindRDeltaRPairing.get)
	DeltaRPairing = assignDijets( j1, j2, j3, j4, minDeltaRPairing  )

	return [ DeltaRPairing[0], DeltaRPairing[1], DeltaRPairing[2], DeltaRPairing[3], DeltaRPairing[4], mindRDeltaRPairing[ minDeltaRPairing ], minDeltaRPairing ]

def MassAsyming( j1, j2, j3, j4 ):
	"""docstring for MassAsyming"""

	massAsymmetry = {}
	massAsymmetry[ '1234' ] = ( abs( j1.M() - j2.M() ) / abs( j3.M() - j4.M() ) )
	massAsymmetry[ '1324' ] = ( abs( j1.M() - j3.M() ) / abs( j2.M() - j4.M() ) )
	massAsymmetry[ '1423' ] = ( abs( j1.M() - j4.M() ) / abs( j2.M() - j3.M() ) )
	minMassAsyming = min(massAsymmetry, key=massAsymmetry.get)
	MassAsyming = assignDijets( j1, j2, j3, j4, minMassAsyming  )

	return [ MassAsyming[0], MassAsyming[1], MassAsyming[2], MassAsyming[3], MassAsyming[4], massAsymmetry[ minMassAsyming ], minMassAsyming ]

def calcCosThetaStar(j1, j2):
	"""docstring for calcCosThetaStar"""

	tmpCM1 = j1 + j2 
	tmpJ1 = TLorentzVector()
	tmpJ2 = TLorentzVector()
	tmpJ1.SetPtEtaPhiE( j1.Pt(), j1.Eta(), j1.Phi(), j1.E() )
	tmpJ2.SetPtEtaPhiE( j2.Pt(), j2.Eta(), j2.Phi(), j2.E() )
	tmpJ1.Boost( -tmpCM1.BoostVector() )
	tmpJ2.Boost( -tmpCM1.BoostVector() )
	tmpV1 = TVector3( tmpJ1.X(), tmpJ1.Y(), tmpJ1.Z() )
	tmpV2 = TVector3( tmpJ2.X(), tmpJ2.Y(), tmpJ2.Z() )
	#cosThetaStar1 = abs( ( ( pairoff08[1].Px() * pairoff08[2].Px() ) + ( pairoff08[1].Py() * pairoff08[2].Py() ) + ( pairoff08[1].Pz() * pairoff08[2].Pz() ) )  / ( pairoff08[1].E() * pairoff08[2].E() ) )
	cosThetaStar = abs( tmpV1.CosTheta() )

	return cosThetaStar

def dijetVar( listJets ):
	"""docstring for dijetVar"""

	dijet1 = listJets[0] + listJets[1]
	dijet2 = listJets[2] + listJets[3]
	deltaEtaDijet1 = abs( listJets[0].Eta() - listJets[1].Eta() )
	deltaEtaDijet2 = abs( listJets[2].Eta() - listJets[3].Eta() )
	deltaEtaAveDijets = ( abs( listJets[0].Eta() - listJets[1].Eta() ) + abs( listJets[2].Eta() - listJets[3].Eta()  ) ) / 2
	massAsymmetry = abs( dijet1.M() - dijet2.M() ) / ( dijet1.M() + dijet2.M() )
	deltaEtaDijets = abs( dijet1.Eta() - dijet2.Eta() )
	massAve = ( dijet1.M() + dijet2.M() ) / 2
	cosThetaStarDijet1 = calcCosThetaStar( listJets[0], listJets[1] )
	cosThetaStarDijet2 = calcCosThetaStar( listJets[2], listJets[3] )
	deltaDijet1 = ( listJets[0].Pt() + listJets[1].Pt() ) - massAve
	deltaDijet2 = ( listJets[2].Pt() + listJets[3].Pt() ) - massAve
	xi1 = ( ( max( [ listJets[0].M() , listJets[1].M() ] ) / dijet1.M() ) * listJets[0].DeltaR( listJets[1])  )
	xi2 = ( ( max( [ listJets[2].M() , listJets[3].M() ] ) / dijet2.M() ) * listJets[2].DeltaR( listJets[3])  )
	deltaRDijet1 = listJets[0].DeltaR( listJets[1] )
	deltaRDijet2 = listJets[2].DeltaR( listJets[3] )

	return [ massAve, deltaEtaDijet1, deltaEtaDijet2, deltaEtaAveDijets, deltaEtaDijets, massAsymmetry, cosThetaStarDijet1, cosThetaStarDijet2, deltaDijet1, deltaDijet2, xi1, xi2, deltaRDijet1, deltaRDijet2 ]


#################################################################################
if __name__ == '__main__':

	usage = 'usage: %prog [options]'
	
	parser = argparse.ArgumentParser()
	parser.add_argument( '-m', '--mass', action='store', dest='mass', default=100, help='Mass of the Stop' )
	parser.add_argument( '-p', '--process', action='store',  dest='process', default='single', help='Process: all or single.' )
	parser.add_argument( '-d', '--decay', action='store',  dest='decay', default='UDD312', help='Decay: UDD312 or UDD323.' )
	parser.add_argument( '-s', '--sample', action='store',   dest='samples', default='RPV', help='Type of sample' )
	parser.add_argument( '-u', '--unc', action='store',  dest='unc', default='', help='Process: all or single.' )
	parser.add_argument( '-b', '--batchSys', action='store_true',  dest='batchSys', default=False, help='Process: all or single.' )
	parser.add_argument( '-v', '--version', action='store', default='v00p3', dest='version', help='Version of the RUNAnalysis file.' )
	#parser.add_argument( '-l', '--lumi', action='store', type=float, default=1, help='Luminosity, example: 1.' )
	parser.add_argument( '-q', '--qcd', action='store', default='Pt', dest='qcd', help='Type of QCD binning, example: HT.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	if args.batchSys: folder = '/cms/gomez/archiveEOS/Archive/v8020/Analysis/'+args.version+'/'
	else: folder = 'Rootfiles/'

	allSamples = {}
	allSamples[ 'JetHT_Run2016'] = folder+'/RUNAnalysis_JetHT_Run2016C_V2p3_'+args.version+'.root'
	allSamples[ 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) ] = folder+'/RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+args.mass+'_80X_V2p3_'+args.version+'.root'
	allSamples[ 'TTJets' ] = folder+'/RUNAnalysis_TT_80X_V2p3_'+args.version+'.root'
    	allSamples[ 'ZJetsToQQ' ] = folder+'/RUNAnalysis_ZJetsToQQ_80X_V2p3_'+args.version+'.root'
    	allSamples[ 'WJetsToQQ' ] = folder+'/RUNAnalysis_WJetsToQQ_80X_V2p3_'+args.version+'.root'
	allSamples[ 'Dibosons' ] = folder+'/RUNAnalysis_Dibosons_80X_V2p3_'+args.version+'.root'
	allSamples[ 'QCD'+args.qcd+'All' ] = folder+'/RUNAnalysis_QCD'+args.qcd+'All_80X_V2p3_'+args.version+'.root'

	if 'RPV' in args.samples: args.samples = 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)
	if 'single' in args.process: 
		for q in allSamples:
			if q in args.samples:
				dictSamples = { q: allSamples[ q ] } 
				signalSample = q
	else: 
		dictSamples = allSamples
		signalSample = 'RPVStopStopToJets_'+args.decay+'_M-'+args.mass+'_All'

	preselection = '(jetsPt[0]>80) && (jetsPt[1]>80) && (jetsPt[2]>80) && (jetsPt[3]>80) && (HT>900)' 
	cuts = '(delta1>200) && (delta2>200) && (massAsym<0.1) && (deltaEta<1.)'
	#preselection = '(jetsPt[0]>50) * (jetsPt[1]>50) *(jetsPt[2]>50) *(jetsPt[3]>50) * (HT>500)' 
	#cuts = '(delta1>200) * (delta2>200) * (massAsym<0.1) * (deltaEta<1.5)'

	allHistos = {}
	if ('RPV' in args.samples) and args.unc:
		for uncType in [ args.unc+'Up', args.unc+'Down' ]: 
			p = Process( target=myAnalyzer, args=( dictSamples, preselection, cuts, signalSample, uncType ) )
			p.start()
			p.join()
	else:
		p = Process( target=myAnalyzer, args=( dictSamples, preselection, cuts, signalSample, '' ) )
		p.start()
		p.join()


