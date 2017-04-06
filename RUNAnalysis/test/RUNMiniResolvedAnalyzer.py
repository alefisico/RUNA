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
from itertools import permutations 
try: 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
	from RUNA.RUNAnalysis.scaleFactors import * 
	from RUNA.RUNAnalysis.commonFunctions import *
except ImportError: 
	sys.path.append('../python') 
	import tdrstyle as tdrstyle
	from scaleFactors import * 
	from commonFunctions import *

gROOT.SetBatch()

######################################
def myPlotAnalyzer( fileSample, preselection, cuts, sample, UNC ):
	"""docstring for myPlotAnalyzer: creates histograms from tree """


	#outputFileName = 'Rootfiles/RUNMiniScoutingResolvedAnalysis_'+sample+UNC+'_'+( '' if 'JetHT' in sample else '80X_')+'V2p1_'+args.version+'p1.root' 
	outputFileName = 'Rootfiles/RUNMiniResolvedAnalysis_'+sample+UNC+'_'+( '' if 'JetHT' in sample else 'Moriond17_')+'80X_V2p4_'+args.version+'p1.root' 
	outputFile = TFile( outputFileName, 'RECREATE' )


	################################################################################################## Trigger Histos
	nBinsMass	= 200
	maxMass		= 2000
	nBinsHT		= 150
	maxHT		= 1500

	print '--- Sample ', sample
	if 'JetHT' in sample: sample = 'JetHT_Run2016'
	allHistos[ "HT_cutBestPair_"+sample ] = TH1F( "HT_cutBestPair_"+sample, "HT_cutBestPair_"+sample, 5000, 0., 5000 )
	allHistos[ "NPV_cutBestPair_"+sample ] = TH1F( "NPV_cutBestPair_"+sample, "NPV_cutBestPair_"+sample, 100, 0., 100 )
	allHistos[ "jet1Pt_cutBestPair_"+sample ] = TH1F( "jet1Pt_cutBestPair_"+sample, "jet1Pt_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "jet2Pt_cutBestPair_"+sample ] = TH1F( "jet2Pt_cutBestPair_"+sample, "jet2Pt_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "jet3Pt_cutBestPair_"+sample ] = TH1F( "jet3Pt_cutBestPair_"+sample, "jet3Pt_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "jet4Pt_cutBestPair_"+sample ] = TH1F( "jet4Pt_cutBestPair_"+sample, "jet4Pt_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_cutBestPair_"+sample ] = TH1F( "massAve_cutBestPair_"+sample, "massAve_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "massAsym_cutBestPair_"+sample ] = TH1F( "massAsym_cutBestPair_"+sample, "massAsym_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "deltaEta_cutBestPair_"+sample ] = TH1F( "deltaEta_cutBestPair_"+sample, "deltaEta_cutBestPair_"+sample, 50, 0., 5 )
	allHistos[ 'deltavsMassAve_cutBestPair_'+sample ] = TH2F( 'deltavsMassAve_cutBestPair_'+sample, 'deltavsMassAve_cutBestPair_'+sample, 1000, 0., 1000, 1000, 0., 1000. )
	allHistos[ "jet1Btag_cutBestPair_"+sample ] = TH1F( "jet1Btag_cutBestPair_"+sample, "jet1Btag_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "jet2Btag_cutBestPair_"+sample ] = TH1F( "jet2Btag_cutBestPair_"+sample, "jet2Btag_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "jet3Btag_cutBestPair_"+sample ] = TH1F( "jet3Btag_cutBestPair_"+sample, "jet3Btag_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "jet4Btag_cutBestPair_"+sample ] = TH1F( "jet4Btag_cutBestPair_"+sample, "jet4Btag_cutBestPair_"+sample, 20, 0., 1 )

	allHistos[ "massAve_delta_"+sample ] = TH1F( "massAve_delta_"+sample, "massAve_delta_"+sample, 3000, 0., 3000 )
	allHistos[ "HT_delta_"+sample ] = TH1F( "HT_delta_"+sample, "HT_delta_"+sample, 5000, 0., 5000 )
	allHistos[ "massAve_delta50_"+sample ] = TH1F( "massAve_delta50_"+sample, "massAve_delta50_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta100_"+sample ] = TH1F( "massAve_delta100_"+sample, "massAve_delta100_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta150_"+sample ] = TH1F( "massAve_delta150_"+sample, "massAve_delta150_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta200_"+sample ] = TH1F( "massAve_delta200_"+sample, "massAve_delta200_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta250_"+sample ] = TH1F( "massAve_delta250_"+sample, "massAve_delta250_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta300_"+sample ] = TH1F( "massAve_delta300_"+sample, "massAve_delta300_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta400_"+sample ] = TH1F( "massAve_delta400_"+sample, "massAve_delta400_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta500_"+sample ] = TH1F( "massAve_delta500_"+sample, "massAve_delta500_"+sample, 3000, 0., 3000 )

	allHistos[ "massAsym_n-1_"+sample ] = TH1F( "massAsym_n-1_"+sample, "massAsym_n-1_"+sample, 20, 0., 1 )
	allHistos[ "deltaEta_n-1_"+sample ] = TH1F( "deltaEta_n-1_"+sample, "deltaEta_n-1_"+sample, 50, 0., 5 )
	allHistos[ 'deltavsMassAve_n-1_'+sample ] = TH2F( 'deltavsMassAve_n-1_'+sample, 'deltavsMassAve_n-1_'+sample, 1000, 0., 1000, 1000, 0., 1000. )

	allHistos[ "massAve_delta_1qgl_"+sample ] = TH1F( "massAve_delta_1qgl_"+sample, "massAve_delta_1qgl_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta_2qgl_"+sample ] = TH1F( "massAve_delta_2qgl_"+sample, "massAve_delta_2qgl_"+sample, 3000, 0., 3000 )
	allHistos[ "jet1QGL_delta_"+sample ] = TH1F( "jet1QGL_delta_"+sample, "jet1QGL_delta_"+sample, 20, 0., 1 )
	allHistos[ "jet2QGL_delta_"+sample ] = TH1F( "jet2QGL_delta_"+sample, "jet2QGL_delta_"+sample, 20, 0., 1 )
	allHistos[ "jet3QGL_delta_"+sample ] = TH1F( "jet3QGL_delta_"+sample, "jet3QGL_delta_"+sample, 20, 0., 1 )
	allHistos[ "jet4QGL_delta_"+sample ] = TH1F( "jet4QGL_delta_"+sample, "jet4QGL_delta_"+sample, 20, 0., 1 )

	allHistos[ "massAve_delta_4qgl_"+sample ] = TH1F( "massAve_delta_4qgl_"+sample, "massAve_delta_4qgl_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_2qgl_"+sample ] = TH1F( "massAve_2qgl_"+sample, "massAve_2qgl_"+sample, 3000, 0., 3000 )

	allHistos[ "massAve_delta_1CSVv2M_"+sample ] = TH1F( "massAve_delta_1CSVv2M_"+sample, "massAve_delta_1CSVv2M_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta_2CSVv2M_"+sample ] = TH1F( "massAve_delta_2CSVv2M_"+sample, "massAve_delta_2CSVv2M_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta_1CSVv2L_"+sample ] = TH1F( "massAve_delta_1CSVv2L_"+sample, "massAve_delta_1CSVv2L_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta_2CSVv2L_"+sample ] = TH1F( "massAve_delta_2CSVv2L_"+sample, "massAve_delta_2CSVv2L_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta_1CSVv2T_"+sample ] = TH1F( "massAve_delta_1CSVv2T_"+sample, "massAve_delta_1CSVv2T_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_delta_2CSVv2T_"+sample ] = TH1F( "massAve_delta_2CSVv2T_"+sample, "massAve_delta_2CSVv2T_"+sample, 3000, 0., 3000 )

	for h in allHistos: allHistos[h].Sumw2()

	################################################################################################## Running the Analysis
	#SF = 'lumiWeight*puWeight'
	lumiWeight = scaleFactor(sample)
	SF = 'puWeight*'+str(lumiWeight)
	presel = TCut( SF ) *  TCut( preselection )
	fullSel = TCut( SF ) * TCut( preselection + ' && ' + cuts ) 
	print '-'*40

	#treeName = 'ResolvedAnalysisPlotsScouting/RUNATree'
	treeName = 'ResolvedAnalysisPlots/RUNATree'

	### All selection
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fullSel,
			allHistos[ 'massAve_delta_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'HT', 
			fullSel,
			allHistos[ 'HT_delta_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )


	#### preselection plots
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			presel,
			allHistos[ 'massAve_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'numPV', 
			presel,
			allHistos[ 'NPV_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'HT', 
			presel,
			allHistos[ 'HT_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsPt[0]', 
			presel,
			allHistos[ 'jet1Pt_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsPt[1]', 
			presel,
			allHistos[ 'jet2Pt_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsPt[2]', 
			presel,
			allHistos[ 'jet3Pt_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsPt[3]', 
			presel,
			allHistos[ 'jet4Pt_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	get2DHistoFromTree( fileSample, treeName,
			'massAve', 'delta1',
			presel,
			allHistos[ 'deltavsMassAve_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )
	get2DHistoFromTree( fileSample, treeName,
			'massAve', 'delta2',
			presel,
			allHistos[ 'deltavsMassAve_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAsym',
			presel,
			allHistos[ 'massAsym_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'deltaEta',
			presel,
			allHistos[ 'deltaEta_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	##### checking diff deltas
	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>50) && (delta2>50)') ),
			allHistos[ 'massAve_delta50_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>100) && (delta2>100)') ),
			allHistos[ 'massAve_delta100_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>150) && (delta2>150)') ),
			allHistos[ 'massAve_delta150_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAve',
			fullSel,
			#presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>200) && (delta2>200)') ),
			allHistos[ 'massAve_delta200_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>250) && (delta2>250)') ),
			allHistos[ 'massAve_delta250_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>300) && (delta2>300)') ),
			allHistos[ 'massAve_delta300_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>400) && (delta2>400)') ),
			allHistos[ 'massAve_delta400_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('(delta1>200) && (delta2>200)', '(delta1>500) && (delta2>500)') ),
			allHistos[ 'massAve_delta500_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )



	#### n-1 plots
	get2DHistoFromTree( fileSample, treeName,
			'massAve', 'delta1', 
			presel * TCut( cuts.replace('&& (delta1>200)','').replace('&& (delta2>200)','') ), 
			allHistos[ 'deltavsMassAve_n-1_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )
	get2DHistoFromTree( fileSample, treeName,
			'massAve', 'delta2', 
			presel * TCut( cuts.replace('&& (delta1>200)','').replace('&& (delta2>200)','') ), 
			allHistos[ 'deltavsMassAve_n-1_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'massAsym',
			presel * TCut( cuts.replace('&& (massAsym<0.1)','') ), 
			allHistos[ 'massAsym_n-1_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'deltaEta',
			presel * TCut( cuts.replace('&& (deltaEta<1.)','') ), 
			allHistos[ 'deltaEta_n-1_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	### Test QGL
	oneQGL = fullSel * TCut( '( (jetsQGL[0]>0.5) || (jetsQGL[1]>0.5) || (jetsQGL[2]>0.5) || (jetsQGL[3]>0.5) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			oneQGL, 
			allHistos[ 'massAve_delta_1qgl_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	twoQGL = fullSel * TCut(' ( (jetsQGL[0]>0.5) || (jetsQGL[1]>0.5)) && ((jetsQGL[2]>0.5) || (jetsQGL[3]>0.5) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			twoQGL, 
			allHistos[ 'massAve_delta_2qgl_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	fourQGL = fullSel * TCut('(jetsQGL[0]>0.5) && (jetsQGL[1]>0.5) && (jetsQGL[2]>0.5) && (jetsQGL[3]>0.5)')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fourQGL, 
			allHistos[ 'massAve_delta_4qgl_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsQGL[0]', 
			fullSel,
			allHistos[ 'jet1QGL_delta_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsQGL[1]', 
			fullSel,
			allHistos[ 'jet2QGL_delta_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsQGL[2]', 
			fullSel,
			allHistos[ 'jet3QGL_delta_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsQGL[3]', 
			fullSel,
			allHistos[ 'jet4QGL_delta_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	twoQGLwithoutDelta = TCut( SF ) * TCut( preselection + ' && ' + cuts.replace('(delta1>200) && (delta2>200) && ', '') ) * TCut(' ( (jetsQGL[0]>0.5) || (jetsQGL[1]>0.5)) && ((jetsQGL[2]>0.5) || (jetsQGL[3]>0.5) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			twoQGLwithoutDelta, 
			allHistos[ 'massAve_2qgl_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) )
	### Test Btag
	oneBtagM = TCut('( (jetsCSVv2[0]>0.8484) || (jetsCSVv2[1]>0.8484) || (jetsCSVv2[2]>0.8484) || (jetsCSVv2[3]>0.8484) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fullSel*oneBtagM, 
			allHistos[ 'massAve_delta_1CSVv2M_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) )

	twoBtagM = TCut('( (jetsCSVv2[0]>0.8484) || (jetsCSVv2[1]>0.8484)) && ((jetsCSVv2[2]>0.8484) || (jetsCSVv2[3]>0.8484) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fullSel*twoBtagM, 
			allHistos[ 'massAve_delta_2CSVv2M_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) )

	oneBtagL = TCut('( (jetsCSVv2[0]>0.5426) || (jetsCSVv2[1]>0.5426) || (jetsCSVv2[2]>0.5426) || (jetsCSVv2[3]>0.5426) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fullSel*oneBtagL, 
			allHistos[ 'massAve_delta_1CSVv2L_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) )

	twoBtagL = TCut('( (jetsCSVv2[0]>0.5426) || (jetsCSVv2[1]>0.5426)) && ((jetsCSVv2[2]>0.5426) || (jetsCSVv2[3]>0.5426) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fullSel*twoBtagL, 
			allHistos[ 'massAve_delta_2CSVv2L_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) )

	oneBtagT = TCut('( (jetsCSVv2[0]>0.9535) || (jetsCSVv2[1]>0.9535) || (jetsCSVv2[2]>0.9535) || (jetsCSVv2[3]>0.9535) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fullSel*oneBtagT, 
			allHistos[ 'massAve_delta_1CSVv2T_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) )

	twoBtagT = TCut('( (jetsCSVv2[0]>0.9535) || (jetsCSVv2[1]>0.9535)) && ((jetsCSVv2[2]>0.9535) || (jetsCSVv2[3]>0.9535) )')
	getHistoFromTree( fileSample, treeName,
			'massAve', 
			fullSel*twoBtagT, 
			allHistos[ 'massAve_delta_2CSVv2T_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) )


	getHistoFromTree( fileSample, treeName,
			'jetsCSVv2[0]', 
			fullSel,
			allHistos[ 'jet1Btag_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsCSVv2[1]', 
			fullSel,
			allHistos[ 'jet2Btag_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsCSVv2[2]', 
			fullSel,
			allHistos[ 'jet3Btag_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	getHistoFromTree( fileSample, treeName,
			'jetsCSVv2[3]', 
			fullSel,
			allHistos[ 'jet4Btag_cutBestPair_'+sample ], 
			( 0.05 if 'JetHT' in sample else 1 ) )

	outputFile.Write()
	##### Closing
	print 'Writing output file: '+ outputFileName
	outputFile.Close()

def myAnalyzer( fileSample, preselection, cuts, sample, UNC ):
	"""docstring for myAnalyzer: creates new variables from tree """

	outputFileName = 'Rootfiles/RUNMiniResolvedAnalysis_'+sample+UNC+'_'+( '' if 'JetHT' in sample else 'Moriond17_')+'80X_V2p4_'+args.version+'p1.root' 
	outputFile = TFile( outputFileName, 'RECREATE' )

	inputFile, events, numEntries = getTree( fileSample, 'ResolvedAnalysisPlots/RUNATree')
	print '-'*40
	print '------> ', sample
	print '------> Number of events: '+str(numEntries)
	d = 0
	#cutFlowList = OrderedDict()
	#cutFlowList[ 'Process' ] = 0
	#cutFlowList[ 'Preselection' ] = 0
	#for k in listCuts: cutFlowList[ k[0] ] = 0

	for i in xrange(numEntries):
		events.GetEntry(i)

		#---- progress of the reading --------
		fraction = 10.*i/(1.*numEntries)
		if TMath.FloorNint(fraction) > d: print str(10*TMath.FloorNint(fraction))+'%' 
		d = TMath.FloorNint(fraction)
		#if ( i > 100000 ): break

		Run		= events.run
		Lumi    	= events.lumi
		NumEvent	= events.event
		puWeight	= events.puWeight
		lumiWeight	= events.lumiWeight
		HT		= events.HT
		MET		= events.MET
		numJets		= events.numJets

		listOfJets = []
		if len(events.jetsPt) > 4: print len(events.jetsPt)
		for j in range( len(events.jetsPt) ):
			tmpTLV = TLorentzVector()
			tmpTLV.SetPtEtaPhiE( events.jetsPt[j], events.jetsEta[j], events.jetsPhi[j], events.jetsE[j] )
			if (j > 3): print j, tmpTLV.M(), events.jetsPt[j], events.jetsEta[j], events.jetsPhi[j], events.jetsE[j] 
			listOfJets.append( tmpTLV )

		#pairingMinChi( listOfJets, 30 )
		#for i in listOfJets: print i.Pt(), len(listOfJets)


def pairingMinChi( listOfJets, sigma ):
	"""docstring for pairingMinChi"""
	
	minChi2 = 9999999999999
	if ( len(listOfJets) > 3 ):
		comb = range(0, len(listOfJets) )
		listOfCombinations = list(permutations(comb, 4))
		print comb, listOfCombinations

		for a in range( 0, len(listOfJets) ):
			for b in range( 0, len(listOfJets) ):
				for c in range( 0, len(listOfJets) ):
					for d in range( 0, len(listOfJets) ):
						#if( (a!=b) and (a!=c) and (a!=d) and (b!=c) and (b!=d) and (c!=d) and (a<b) and (b<c) and (a<d) and (c<d) ):
						if( (a<b) and (a<c) and (a<d) and (b>c) and (b>d) and (c>d) ):
							print a, b, c, d
							dijet1Mass = ( listOfJets[a] + listOfJets[b] ).M()	
							dijet2Mass = ( listOfJets[c] + listOfJets[d] ).M()	
							tmpChi2 = TMath.Power( (dijet1Mass - dijet2Mass ), 2 ) / TMath.Power( sigma, 2 )
							#print a, b, c, d, dijet1Mass, dijet2Mass, tmpChi2
							if( tmpChi2 < minChi2 ):
								minChi2 = tmpChi2
								bestInd = [ a, b, c, d ]
								#print bestInd
	if ( len(bestInd) == 4 ): reorderedJETS = [ listOfJets[ bestInd[0] ], listOfJets[ bestInd[1] ], listOfJets[ bestInd[2] ], listOfJets[ bestInd[3] ] ]


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
	parser.add_argument( '-p', '--process', action='store',  dest='process', default='Plots', help='Process: all or single.' )
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
	allSamples[ 'JetHT_Run2016B'] = folder+'/RUNAnalysis_JetHT_Run2016B_80X_V2p4_'+args.version+'.root'
	allSamples[ 'JetHT_Run2016C'] = folder+'/RUNAnalysis_JetHT_Run2016C_80X_V2p4_'+args.version+'.root'
	allSamples[ 'JetHT_Run2016D'] = folder+'/RUNAnalysis_JetHT_Run2016D_80X_V2p4_'+args.version+'.root'
	allSamples[ 'JetHT_Run2016E'] = folder+'/RUNAnalysis_JetHT_Run2016E_80X_V2p4_'+args.version+'.root'
	allSamples[ 'JetHT_Run2016F'] = folder+'/RUNAnalysis_JetHT_Run2016F_80X_V2p4_'+args.version+'.root'
	allSamples[ 'JetHT_Run2016G'] = folder+'/RUNAnalysis_JetHT_Run2016G_80X_V2p4_'+args.version+'.root'
	allSamples[ 'JetHT_Run2016H'] = folder+'/RUNAnalysis_JetHT_Run2016H_80X_V2p4_'+args.version+'.root'
	allSamples[ 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) ] = folder+'/RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+args.mass+'_80X_V2p4_'+args.version+'.root'
	allSamples[ 'TTJets' ] = folder+'/RUNAnalysis_TT_80X_V2p4_'+args.version+'.root'
    	allSamples[ 'ZJetsToQQ' ] = folder+'/RUNAnalysis_ZJetsToQQ_80X_V2p4_'+args.version+'.root'
    	allSamples[ 'WJetsToQQ' ] = folder+'/RUNAnalysis_WJetsToQQ_80X_V2p4_'+args.version+'.root'
	allSamples[ 'Dibosons' ] = folder+'/RUNAnalysis_Dibosons_80X_V2p4_'+args.version+'.root'
	allSamples[ 'QCD'+args.qcd+'All' ] = folder+'/RUNAnalysis_QCD'+args.qcd+'All_80X_V2p4_'+args.version+'.root'

	if 'RPV' in args.samples: args.samples = 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)
	dictSamples = OrderedDict()
	for sam in allSamples: 
		if sam.startswith( args.samples ): dictSamples[ sam ] = allSamples[ sam ]

	preselection = '(jetsPt[0]>80) && (jetsPt[1]>80) && (jetsPt[2]>80) && (jetsPt[3]>80) && (HT>900)' 
	cuts = '(delta1>200) && (delta2>200) && (massAsym<0.1) && (deltaEta<1.)'
	#preselection = '(jetsPt[0]>50) * (jetsPt[1]>50) *(jetsPt[2]>50) *(jetsPt[3]>50) * (HT>500)' 
	#cuts = '(delta1>200) * (delta2>200) * (massAsym<0.1) * (deltaEta<1.5)'

	allHistos = {}
	for sample in dictSamples:
		if 'Plots' in args.process: 
			if ('RPV' in args.samples) and args.unc:
				for uncType in [ args.unc+'Up', args.unc+'Down' ]: 
					p = Process( target=myPlotAnalyzer, args=( dictSamples[sample], preselection, cuts, sample, uncType ) )
			else: p = Process( target=myPlotAnalyzer, args=( dictSamples[sample], preselection, cuts, sample, '' ) )
		else:
			#p = Process( target=myAnalyzer, args=( dictSamples[sample], preselection, cuts, sample, '' ) )
			p = Process( target=myAnalyzer, args=( 'RUNResolvedAnalysis_QCD_Pt_3200toInf.root', preselection, cuts, sample, '' ) )
	p.start()
	p.join()


