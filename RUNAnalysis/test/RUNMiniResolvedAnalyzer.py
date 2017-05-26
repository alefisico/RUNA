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
from itertools import combinations
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


	outputFileName = 'Rootfiles/RUNMiniResolvedAnalysis_'+sample+UNC+'_'+( '' if 'JetHT' in sample else 'Moriond17_')+'80X_V2p4_'+args.version+'p0.root' 
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
	allHistos[ 'deltavsMassAve_cutBestPair_'+sample ] = TH2F( 'deltavsMassAve_cutBestPair_'+sample, 'deltavsMassAve_cutBestPair_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "jet1Btag_cutBestPair_"+sample ] = TH1F( "jet1Btag_cutBestPair_"+sample, "jet1Btag_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "jet2Btag_cutBestPair_"+sample ] = TH1F( "jet2Btag_cutBestPair_"+sample, "jet2Btag_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "jet3Btag_cutBestPair_"+sample ] = TH1F( "jet3Btag_cutBestPair_"+sample, "jet3Btag_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "jet4Btag_cutBestPair_"+sample ] = TH1F( "jet4Btag_cutBestPair_"+sample, "jet4Btag_cutBestPair_"+sample, 20, 0., 1 )

	allHistos[ "massAve_delta_"+sample ] = TH1F( "massAve_delta_"+sample, "massAve_delta_"+sample, 3000, 0., 3000 )
	allHistos[ "HT_delta_"+sample ] = TH1F( "HT_delta_"+sample, "HT_delta_"+sample, 5000, 0., 5000 )
	allHistos[ "massAve_woMassAsym_"+sample ] = TH1F( "massAve_woMassAsym_"+sample, "massAve_woMassAsym_"+sample, 3000, 0., 3000 )
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
	SF = 'lumiWeight*puWeight'
	#lumiWeight = scaleFactor(sample)
	#SF = 'puWeight*'+str(lumiWeight)
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
	getHistoFromTree( fileSample, treeName,
			'massAve',
			presel * TCut( cuts.replace('&& (massAsym<0.1)', '') ),
			allHistos[ 'massAve_woMassAsym_'+sample ], 
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

	outputFileName = 'Rootfiles/RUNMiniResolvedAnalysis_'+sample+UNC+'_'+( '' if 'JetHT' in sample else 'Moriond17_')+'80X_V2p4_'+args.version+'p2.root' 
	#outputFileName = 'Rootfiles/RUNMiniResolvedAnalysis_'+sample+UNC+'_'+( '' if 'JetHT' in sample else 'Moriond17_')+'80X_V2p4_'+args.version+'p2_4jets.root' 
	outputFile = TFile( outputFileName, 'RECREATE' )

	print '--- Sample ', sample
	if 'JetHT' in sample: sample = 'JetHT_Run2016'
	allHistos[ "HT_cutBestPair_"+sample ] = TH1F( "HT_cutBestPair_"+sample, "HT_cutBestPair_"+sample, 5000, 0., 5000 )
	allHistos[ "massAve_cutBestPair_"+sample ] = TH1F( "massAve_cutBestPair_"+sample, "massAve_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "massAsym_cutBestPair_"+sample ] = TH1F( "massAsym_cutBestPair_"+sample, "massAsym_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ "deltaEta_cutBestPair_"+sample ] = TH1F( "deltaEta_cutBestPair_"+sample, "deltaEta_cutBestPair_"+sample, 50, 0., 5 )
	allHistos[ 'deltavsMassAve_cutBestPair_'+sample ] = TH2F( 'deltavsMassAve_cutBestPair_'+sample, 'deltavsMassAve_cutBestPair_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "minChi2_"+sample ] = TH1F( "minChi2_"+sample, "minChi2_"+sample, 500, 0., 50 )
	allHistos[ "secondMinChi2_"+sample ] = TH1F( "secondMinChi2_"+sample, "secondMinChi2_"+sample, 500, 0., 50 )
	allHistos[ "minVsSecondMinChi2_"+sample ] = TH2F( "minVsSecondMinChi2_"+sample, "minVsSecondMinChi2_"+sample, 500, 0., 50, 500, 0., 50 )
	allHistos[ "massAveVsMinChi2_"+sample ] = TH2F( "massAveVsMinChi2_"+sample, "massAveVsMinChi2_"+sample, 3000, 0., 3000, 500, 0., 50 )
	allHistos[ "massAve_minChi_cutBestPair_"+sample ] = TH1F( "massAve_minChi_cutBestPair_"+sample, "massAve_minChi_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minChi_cutBestPair_"+sample ] = TH1F( "deltaEta_minChi_cutBestPair_"+sample, "deltaEta_minChi_cutBestPair_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minChi_cutBestPair_"+sample ] = TH1F( "massAsym_minChi_cutBestPair_"+sample, "massAsym_minChi_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minChi_cutBestPair_'+sample ] = TH2F( 'deltavsMassAve_minChi_cutBestPair_'+sample, 'deltavsMassAve_minChi_cutBestPair_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "cosThetaStar_minChi_cutBestPair_"+sample ] = TH1F( "cosThetaStar_minChi_cutBestPair_"+sample, "cosThetaStar_minChi_cutBestPair_"+sample, 20, 0., 2)
	allHistos[ 'massAveVscosThetaStar_minChi_cutBestPair_'+sample ] = TH2F( 'massAveVscosThetaStar_minChi_cutBestPair_'+sample, 'massAveVscosThetaStar_minChi_cutBestPair_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "xi_minChi_cutBestPair_"+sample ] = TH1F( "xi_minChi_cutBestPair_"+sample, "xi_minChi_cutBestPair_"+sample, 20, 0., 2)
	allHistos[ 'massAveVsXi_minChi_cutBestPair_'+sample ] = TH2F( 'massAveVsXi_minChi_cutBestPair_'+sample, 'massAveVsXi_minChi_cutBestPair_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "massAve_minChi_cutDeltaEta_"+sample ] = TH1F( "massAve_minChi_cutDeltaEta_"+sample, "massAve_minChi_cutDeltaEta_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minChi_cutDeltaEta_"+sample ] = TH1F( "deltaEta_minChi_cutDeltaEta_"+sample, "deltaEta_minChi_cutDeltaEta_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minChi_cutDeltaEta_"+sample ] = TH1F( "massAsym_minChi_cutDeltaEta_"+sample, "massAsym_minChi_cutDeltaEta_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minChi_cutDeltaEta_'+sample ] = TH2F( 'deltavsMassAve_minChi_cutDeltaEta_'+sample, 'deltavsMassAve_minChi_cutDeltaEta_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "massAve_minChi_cutDelta_"+sample ] = TH1F( "massAve_minChi_cutDelta_"+sample, "massAve_minChi_cutDelta_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minChi_cutDelta_"+sample ] = TH1F( "deltaEta_minChi_cutDelta_"+sample, "deltaEta_minChi_cutDelta_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minChi_cutDelta_"+sample ] = TH1F( "massAsym_minChi_cutDelta_"+sample, "massAsym_minChi_cutDelta_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minChi_cutDelta_'+sample ] = TH2F( 'deltavsMassAve_minChi_cutDelta_'+sample, 'deltavsMassAve_minChi_cutDelta_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "cosThetaStar_minChi_cutDelta_"+sample ] = TH1F( "cosThetaStar_minChi_cutDelta_"+sample, "cosThetaStar_minChi_cutDelta_"+sample, 20, 0., 2)
	allHistos[ 'massAveVscosThetaStar_minChi_cutDelta_'+sample ] = TH2F( 'massAveVscosThetaStar_minChi_cutDelta_'+sample, 'massAveVscosThetaStar_minChi_cutDelta_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "xi_minChi_cutDelta_"+sample ] = TH1F( "xi_minChi_cutDelta_"+sample, "xi_minChi_cutDelta_"+sample, 20, 0., 2)
	allHistos[ 'massAveVsXi_minChi_cutDelta_'+sample ] = TH2F( 'massAveVsXi_minChi_cutDelta_'+sample, 'massAveVsXi_minChi_cutDelta_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "massAve_minChi_cutXi_"+sample ] = TH1F( "massAve_minChi_cutXi_"+sample, "massAve_minChi_cutXi_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minChi_cutDelta75_"+sample ] = TH1F( "massAve_minChi_cutDelta75_"+sample, "massAve_minChi_cutDelta75_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minChi_cutDelta100_"+sample ] = TH1F( "massAve_minChi_cutDelta100_"+sample, "massAve_minChi_cutDelta100_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minChi_cutDelta125_"+sample ] = TH1F( "massAve_minChi_cutDelta125_"+sample, "massAve_minChi_cutDelta125_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minChi_cutDelta150_"+sample ] = TH1F( "massAve_minChi_cutDelta150_"+sample, "massAve_minChi_cutDelta150_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minChi_cutDelta175_"+sample ] = TH1F( "massAve_minChi_cutDelta175_"+sample, "massAve_minChi_cutDelta175_"+sample, 3000, 0., 3000 )

	allHistos[ "minDeltaR_"+sample ] = TH1F( "minDeltaR_"+sample, "minDeltaR_"+sample, 50, 0., 5 )
	allHistos[ "secondMinDeltaR_"+sample ] = TH1F( "secondMinDeltaR_"+sample, "secondMinDeltaR_"+sample, 50, 0., 5 )
	allHistos[ "minVsSecondMinDeltaR_"+sample ] = TH2F( "minVsSecondMinDeltaR_"+sample, "minVsSecondMinDeltaR_"+sample, 50, 0., 5, 50, 0., 5 )
	allHistos[ "massAveVsMinDeltaR_"+sample ] = TH2F( "massAveVsMinDeltaR_"+sample, "massAveVsMinDeltaR_"+sample, 3000, 0., 3000, 50, 0., 5 )
	allHistos[ "massAve_minDeltaR_cutBestPair_"+sample ] = TH1F( "massAve_minDeltaR_cutBestPair_"+sample, "massAve_minDeltaR_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minDeltaR_cutBestPair_"+sample ] = TH1F( "deltaEta_minDeltaR_cutBestPair_"+sample, "deltaEta_minDeltaR_cutBestPair_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minDeltaR_cutBestPair_"+sample ] = TH1F( "massAsym_minDeltaR_cutBestPair_"+sample, "massAsym_minDeltaR_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minDeltaR_cutBestPair_'+sample ] = TH2F( 'deltavsMassAve_minDeltaR_cutBestPair_'+sample, 'deltavsMassAve_minDeltaR_cutBestPair_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "cosThetaStar_minDeltaR_cutBestPair_"+sample ] = TH1F( "cosThetaStar_minDeltaR_cutBestPair_"+sample, "cosThetaStar_minDeltaR_cutBestPair_"+sample, 20, 0., 2)
	allHistos[ 'massAveVscosThetaStar_minDeltaR_cutBestPair_'+sample ] = TH2F( 'massAveVscosThetaStar_minDeltaR_cutBestPair_'+sample, 'massAveVscosThetaStar_minDeltaR_cutBestPair_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "xi_minDeltaR_cutBestPair_"+sample ] = TH1F( "xi_minDeltaR_cutBestPair_"+sample, "xi_minDeltaR_cutBestPair_"+sample, 20, 0., 2)
	allHistos[ 'massAveVsXi_minDeltaR_cutBestPair_'+sample ] = TH2F( 'massAveVsXi_minDeltaR_cutBestPair_'+sample, 'massAveVsXi_minDeltaR_cutBestPair_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "massAve_minDeltaR_cutDeltaEta_"+sample ] = TH1F( "massAve_minDeltaR_cutDeltaEta_"+sample, "massAve_minDeltaR_cutDeltaEta_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minDeltaR_cutDeltaEta_"+sample ] = TH1F( "deltaEta_minDeltaR_cutDeltaEta_"+sample, "deltaEta_minDeltaR_cutDeltaEta_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minDeltaR_cutDeltaEta_"+sample ] = TH1F( "massAsym_minDeltaR_cutDeltaEta_"+sample, "massAsym_minDeltaR_cutDeltaEta_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minDeltaR_cutDeltaEta_'+sample ] = TH2F( 'deltavsMassAve_minDeltaR_cutDeltaEta_'+sample, 'deltavsMassAve_minDeltaR_cutDeltaEta_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "massAve_minDeltaR_cutDelta_"+sample ] = TH1F( "massAve_minDeltaR_cutDelta_"+sample, "massAve_minDeltaR_cutDelta_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minDeltaR_cutDelta_"+sample ] = TH1F( "deltaEta_minDeltaR_cutDelta_"+sample, "deltaEta_minDeltaR_cutDelta_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minDeltaR_cutDelta_"+sample ] = TH1F( "massAsym_minDeltaR_cutDelta_"+sample, "massAsym_minDeltaR_cutDelta_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minDeltaR_cutDelta_'+sample ] = TH2F( 'deltavsMassAve_minDeltaR_cutDelta_'+sample, 'deltavsMassAve_minDeltaR_cutDelta_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "cosThetaStar_minDeltaR_cutDelta_"+sample ] = TH1F( "cosThetaStar_minDeltaR_cutDelta_"+sample, "cosThetaStar_minDeltaR_cutDelta_"+sample, 20, 0., 2)
	allHistos[ 'massAveVscosThetaStar_minDeltaR_cutDelta_'+sample ] = TH2F( 'massAveVscosThetaStar_minDeltaR_cutDelta_'+sample, 'massAveVscosThetaStar_minDeltaR_cutDelta_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "xi_minDeltaR_cutDelta_"+sample ] = TH1F( "xi_minDeltaR_cutDelta_"+sample, "xi_minDeltaR_cutDelta_"+sample, 20, 0., 2)
	allHistos[ 'massAveVsXi_minDeltaR_cutDelta_'+sample ] = TH2F( 'massAveVsXi_minDeltaR_cutDelta_'+sample, 'massAveVsXi_minDeltaR_cutDelta_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "massAve_minDeltaR_cutXi_"+sample ] = TH1F( "massAve_minDeltaR_cutXi_"+sample, "massAve_minDeltaR_cutXi_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minDeltaR_cutMassAsym_"+sample ] = TH1F( "massAve_minDeltaR_cutMassAsym_"+sample, "massAve_minDeltaR_cutMassAsym_"+sample, 3000, 0., 3000 )

	allHistos[ "minMassAsym_"+sample ] = TH1F( "minMassAsym_"+sample, "minMassAsym_"+sample, 20, 0., 1 )
	allHistos[ "secondMinMassAsym_"+sample ] = TH1F( "secondMinMassAsym_"+sample, "secondMinMassAsym_"+sample, 20, 0., 1 )
	allHistos[ "minVsSecondMinMassAsym_"+sample ] = TH2F( "minVsSecondMinMassAsym_"+sample, "minVsSecondMinMassAsym_"+sample, 20, 0., 1, 20, 0., 1 )
	allHistos[ "massAveVsMinMassAsym_"+sample ] = TH2F( "massAveVsMinMassAsym_"+sample, "massAveVsMinMassAsym_"+sample, 3000, 0., 3000, 20, 0., 1 )
	allHistos[ "massAve_minMass_cutBestPair_"+sample ] = TH1F( "massAve_minMass_cutBestPair_"+sample, "massAve_minMass_cutBestPair_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minMass_cutBestPair_"+sample ] = TH1F( "deltaEta_minMass_cutBestPair_"+sample, "deltaEta_minMass_cutBestPair_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minMass_cutBestPair_"+sample ] = TH1F( "massAsym_minMass_cutBestPair_"+sample, "massAsym_minMass_cutBestPair_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minMass_cutBestPair_'+sample ] = TH2F( 'deltavsMassAve_minMass_cutBestPair_'+sample, 'deltavsMassAve_minMass_cutBestPair_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "cosThetaStar_minMass_cutBestPair_"+sample ] = TH1F( "cosThetaStar_minMass_cutBestPair_"+sample, "cosThetaStar_minMass_cutBestPair_"+sample, 20, 0., 2)
	allHistos[ 'massAveVscosThetaStar_minMass_cutBestPair_'+sample ] = TH2F( 'massAveVscosThetaStar_minMass_cutBestPair_'+sample, 'massAveVscosThetaStar_minMass_cutBestPair_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "xi_minMass_cutBestPair_"+sample ] = TH1F( "xi_minMass_cutBestPair_"+sample, "xi_minMass_cutBestPair_"+sample, 20, 0., 2)
	allHistos[ 'massAveVsXi_minMass_cutBestPair_'+sample ] = TH2F( 'massAveVsXi_minMass_cutBestPair_'+sample, 'massAveVsXi_minMass_cutBestPair_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "massAve_minMass_cutDeltaEta_"+sample ] = TH1F( "massAve_minMass_cutDeltaEta_"+sample, "massAve_minMass_cutDeltaEta_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minMass_cutDeltaEta_"+sample ] = TH1F( "deltaEta_minMass_cutDeltaEta_"+sample, "deltaEta_minMass_cutDeltaEta_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minMass_cutDeltaEta_"+sample ] = TH1F( "massAsym_minMass_cutDeltaEta_"+sample, "massAsym_minMass_cutDeltaEta_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minMass_cutDeltaEta_'+sample ] = TH2F( 'deltavsMassAve_minMass_cutDeltaEta_'+sample, 'deltavsMassAve_minMass_cutDeltaEta_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "massAve_minMass_cutDelta_"+sample ] = TH1F( "massAve_minMass_cutDelta_"+sample, "massAve_minMass_cutDelta_"+sample, 3000, 0., 3000 )
	allHistos[ "deltaEta_minMass_cutDelta_"+sample ] = TH1F( "deltaEta_minMass_cutDelta_"+sample, "deltaEta_minMass_cutDelta_"+sample, 50, 0., 5)
	allHistos[ "massAsym_minMass_cutDelta_"+sample ] = TH1F( "massAsym_minMass_cutDelta_"+sample, "massAsym_minMass_cutDelta_"+sample, 20, 0., 1 )
	allHistos[ 'deltavsMassAve_minMass_cutDelta_'+sample ] = TH2F( 'deltavsMassAve_minMass_cutDelta_'+sample, 'deltavsMassAve_minMass_cutDelta_'+sample, 1000, 0., 1000, 2000, -1000., 1000. )
	allHistos[ "cosThetaStar_minMass_cutDelta_"+sample ] = TH1F( "cosThetaStar_minMass_cutDelta_"+sample, "cosThetaStar_minMass_cutDelta_"+sample, 20, 0., 2)
	allHistos[ 'massAveVscosThetaStar_minMass_cutDelta_'+sample ] = TH2F( 'massAveVscosThetaStar_minMass_cutDelta_'+sample, 'massAveVscosThetaStar_minMass_cutDelta_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "xi_minMass_cutDelta_"+sample ] = TH1F( "xi_minMass_cutDelta_"+sample, "xi_minMass_cutDelta_"+sample, 20, 0., 2)
	allHistos[ 'massAveVsXi_minMass_cutDelta_'+sample ] = TH2F( 'massAveVsXi_minMass_cutDelta_'+sample, 'massAveVsXi_minMass_cutDelta_'+sample, 1000, 0., 1000, 20, 0., 2. )
	allHistos[ "massAve_minMass_cutXi_"+sample ] = TH1F( "massAve_minMass_cutXi_"+sample, "massAve_minMass_cutXi_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minMass_cutDelta75_"+sample ] = TH1F( "massAve_minMass_cutDelta75_"+sample, "massAve_minMass_cutDelta75_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minMass_cutDelta100_"+sample ] = TH1F( "massAve_minMass_cutDelta100_"+sample, "massAve_minMass_cutDelta100_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minMass_cutDelta125_"+sample ] = TH1F( "massAve_minMass_cutDelta125_"+sample, "massAve_minMass_cutDelta125_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minMass_cutDelta150_"+sample ] = TH1F( "massAve_minMass_cutDelta150_"+sample, "massAve_minMass_cutDelta150_"+sample, 3000, 0., 3000 )
	allHistos[ "massAve_minMass_cutDelta175_"+sample ] = TH1F( "massAve_minMass_cutDelta175_"+sample, "massAve_minMass_cutDelta175_"+sample, 3000, 0., 3000 )

	allHistos[ "minChi2vsminDeltaR_"+sample ] = TH2F( "minChi2vsminDeltaR_"+sample, "minChi2vsminDeltaR_"+sample, 500, 0., 50, 50, 0., 5 )
	allHistos[ "massAveminChi2vsminDeltaR_"+sample ] = TH2F( "massAveminChi2vsminDeltaR_"+sample, "massAveminChi2vsminDeltaR_"+sample, 3000, 0., 3000, 3000, 0., 3000 )

	for h in allHistos: allHistos[h].Sumw2()

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
		#print Run, Lumi, NumEvent
		SF = puWeight * lumiWeight

		listOfJets = []
		for j in range( len(events.jetsPt) ):
			tmpTLV = TLorentzVector()
			tmpTLV.SetPtEtaPhiE( events.jetsPt[j], events.jetsEta[j], events.jetsPhi[j], events.jetsE[j] )
			listOfJets.append( tmpTLV )

		if len(listOfJets) > 3:
		#if len(listOfJets) == 4:
			comb = range(0, len(listOfJets) )
			possibleCombinations = []
			if len(comb) == 4 : possibleCombinations = [ [ 0, 1, 2, 3 ], [ 0, 2, 1, 3 ], [ 0, 3, 1, 2 ] ]
			else:
				for item in combinations(comb, 4): 
					possibleCombinations.append( sorted(item) )

			listMinChi, listOfPairsMinChi = bestPairing( listOfJets, possibleCombinations, method='minChi2' )
			listOldMinChi, listOfPairsOldMinChi = bestPairing( listOfJets, possibleCombinations, method='minChi2', offset=10 ) ## offset is dummy
			listDeltaR, listOfPairsDeltaR = bestPairing( listOfJets, possibleCombinations, method='deltaR', offset=0.8 )
			listMassAsym, listOfPairsMassAsym = bestPairing( listOfJets, possibleCombinations, method='mass' )

			############ Min Chi2
			varMinChi = dijetVar( listOfPairsMinChi ) 
			allHistos[ 'minChi2_'+sample ].Fill( listMinChi[0], SF )
			allHistos[ 'massAveVsMinChi2_'+sample ].Fill( varMinChi[0], listMinChi[0], SF )
			allHistos[ "massAve_minChi_cutBestPair_"+sample ].Fill( varMinChi[0], SF )
			allHistos[ "deltaEta_minChi_cutBestPair_"+sample ].Fill( varMinChi[4], SF )
			allHistos[ "massAsym_minChi_cutBestPair_"+sample ].Fill( varMinChi[5], SF )
			allHistos[ 'deltavsMassAve_minChi_cutBestPair_'+sample ].Fill( varMinChi[0], varMinChi[8], SF )
			allHistos[ 'deltavsMassAve_minChi_cutBestPair_'+sample ].Fill( varMinChi[0], varMinChi[9], SF )
			allHistos[ 'cosThetaStar_minChi_cutBestPair_'+sample ].Fill( varMinChi[10], SF )
			allHistos[ 'massAveVscosThetaStar_minChi_cutBestPair_'+sample ].Fill( varMinChi[0], varMinChi[10], SF )
			allHistos[ 'xi_minChi_cutBestPair_'+sample ].Fill( varMinChi[11], SF )
			allHistos[ 'massAveVsXi_minChi_cutBestPair_'+sample ].Fill( varMinChi[0], varMinChi[11], SF )
			if len(listMinChi)>1: 
				allHistos[ 'secondMinChi2_'+sample ].Fill( listMinChi[1], SF )
				allHistos[ 'minVsSecondMinChi2_'+sample ].Fill( listMinChi[0], listMinChi[1], SF )
			if ( varMinChi[4] < 1.0 ):
				allHistos[ "massAve_minChi_cutDeltaEta_"+sample ].Fill( varMinChi[0], SF )
				allHistos[ "deltaEta_minChi_cutDeltaEta_"+sample ].Fill( varMinChi[4], SF )
				allHistos[ "massAsym_minChi_cutDeltaEta_"+sample ].Fill( varMinChi[5], SF )
				allHistos[ 'deltavsMassAve_minChi_cutDeltaEta_'+sample ].Fill( varMinChi[0], varMinChi[8], SF )
				allHistos[ 'deltavsMassAve_minChi_cutDeltaEta_'+sample ].Fill( varMinChi[0], varMinChi[9], SF )
				if (( varMinChi[8] > 200 ) and ( varMinChi[9] > 200 )):
					allHistos[ "massAve_minChi_cutDelta_"+sample ].Fill( varMinChi[0], SF )
					allHistos[ "deltaEta_minChi_cutDelta_"+sample ].Fill( varMinChi[4], SF )
					allHistos[ "massAsym_minChi_cutDelta_"+sample ].Fill( varMinChi[5], SF )
					allHistos[ 'deltavsMassAve_minChi_cutDelta_'+sample ].Fill( varMinChi[0], varMinChi[8], SF )
					allHistos[ "massAsym_minChi_cutDelta_"+sample ].Fill( varMinChi[5], SF )
					allHistos[ 'cosThetaStar_minChi_cutDelta_'+sample ].Fill( varMinChi[10], SF )
					allHistos[ 'massAveVscosThetaStar_minChi_cutDelta_'+sample ].Fill( varMinChi[0], varMinChi[10], SF )
					allHistos[ 'xi_minChi_cutDelta_'+sample ].Fill( varMinChi[11], SF )
					allHistos[ 'massAveVsXi_minChi_cutDelta_'+sample ].Fill( varMinChi[0], varMinChi[11], SF )
					if ( varMinChi[11] > .5 ):
						allHistos[ "massAve_minChi_cutXi_"+sample ].Fill( varMinChi[0], SF )
				if (( varMinChi[8] > 75 ) and ( varMinChi[9] > 75 )):
					allHistos[ "massAve_minChi_cutDelta75_"+sample ].Fill( varMinChi[0], SF )
				if (( varMinChi[8] > 100 ) and ( varMinChi[9] > 100 )):
					allHistos[ "massAve_minChi_cutDelta100_"+sample ].Fill( varMinChi[0], SF )
				if (( varMinChi[8] > 125 ) and ( varMinChi[9] > 125 )):
					allHistos[ "massAve_minChi_cutDelta125_"+sample ].Fill( varMinChi[0], SF )
				if (( varMinChi[8] > 150 ) and ( varMinChi[9] > 150 )):
					allHistos[ "massAve_minChi_cutDelta150_"+sample ].Fill( varMinChi[0], SF )
				if (( varMinChi[8] > 175 ) and ( varMinChi[9] > 175 )):
					allHistos[ "massAve_minChi_cutDelta175_"+sample ].Fill( varMinChi[0], SF )



			############ Delta R
			varDeltaR = dijetVar( listOfPairsDeltaR ) 
			allHistos[ 'minDeltaR_'+sample ].Fill( listDeltaR[0], SF )
			allHistos[ 'massAveVsMinDeltaR_'+sample ].Fill( varDeltaR[0], listDeltaR[0], SF )
			allHistos[ "massAve_minDeltaR_cutBestPair_"+sample ].Fill( varDeltaR[0], SF )
			allHistos[ "deltaEta_minDeltaR_cutBestPair_"+sample ].Fill( varDeltaR[4], SF )
			allHistos[ "massAsym_minDeltaR_cutBestPair_"+sample ].Fill( varDeltaR[5], SF )
			allHistos[ 'deltavsMassAve_minDeltaR_cutBestPair_'+sample ].Fill( varDeltaR[0], varDeltaR[8], SF )
			allHistos[ 'deltavsMassAve_minDeltaR_cutBestPair_'+sample ].Fill( varDeltaR[0], varDeltaR[9], SF )
			allHistos[ 'cosThetaStar_minDeltaR_cutBestPair_'+sample ].Fill( varDeltaR[10], SF )
			allHistos[ 'massAveVscosThetaStar_minDeltaR_cutBestPair_'+sample ].Fill( varDeltaR[0], varDeltaR[10], SF )
			allHistos[ 'xi_minDeltaR_cutBestPair_'+sample ].Fill( varDeltaR[11], SF )
			allHistos[ 'massAveVsXi_minDeltaR_cutBestPair_'+sample ].Fill( varDeltaR[0], varDeltaR[11], SF )
			if len(listDeltaR)>1: 
				allHistos[ 'secondMinDeltaR_'+sample ].Fill( listDeltaR[1], SF )
				allHistos[ 'minVsSecondMinDeltaR_'+sample ].Fill( listDeltaR[0], listDeltaR[1], SF )
			if ( varDeltaR[4] < 1.0 ):
				allHistos[ "massAve_minDeltaR_cutDeltaEta_"+sample ].Fill( varDeltaR[0], SF )
				allHistos[ "deltaEta_minDeltaR_cutDeltaEta_"+sample ].Fill( varDeltaR[4], SF )
				allHistos[ "massAsym_minDeltaR_cutDeltaEta_"+sample ].Fill( varDeltaR[5], SF )
				allHistos[ 'deltavsMassAve_minDeltaR_cutDeltaEta_'+sample ].Fill( varDeltaR[0], varDeltaR[8], SF )
				allHistos[ 'deltavsMassAve_minDeltaR_cutDeltaEta_'+sample ].Fill( varDeltaR[0], varDeltaR[9], SF )
				if (( varDeltaR[8] > 200 ) and ( varDeltaR[9] > 200 )):
					allHistos[ "massAve_minDeltaR_cutDelta_"+sample ].Fill( varDeltaR[0], SF )
					allHistos[ "deltaEta_minDeltaR_cutDelta_"+sample ].Fill( varDeltaR[4], SF )
					allHistos[ "massAsym_minDeltaR_cutDelta_"+sample ].Fill( varDeltaR[5], SF )
					allHistos[ 'deltavsMassAve_minDeltaR_cutDelta_'+sample ].Fill( varDeltaR[0], varDeltaR[8], SF )
					allHistos[ 'deltavsMassAve_minDeltaR_cutDelta_'+sample ].Fill( varDeltaR[0], varDeltaR[9], SF )
					allHistos[ 'cosThetaStar_minDeltaR_cutDelta_'+sample ].Fill( varDeltaR[10], SF )
					allHistos[ 'massAveVscosThetaStar_minDeltaR_cutDelta_'+sample ].Fill( varDeltaR[0], varDeltaR[10], SF )
					allHistos[ 'xi_minDeltaR_cutDelta_'+sample ].Fill( varDeltaR[11], SF )
					allHistos[ 'massAveVsXi_minDeltaR_cutDelta_'+sample ].Fill( varDeltaR[0], varDeltaR[11], SF )
					if ( varDeltaR[11] > .5 ):
						allHistos[ "massAve_minDeltaR_cutXi_"+sample ].Fill( varDeltaR[0], SF )
					if ( varDeltaR[5] < .1 ):
						allHistos[ "massAve_minDeltaR_cutMassAsym_"+sample ].Fill( varDeltaR[0], SF )

			############ Mass Asym
			varMassAsym = dijetVar( listOfPairsMassAsym ) 
			allHistos[ 'minMassAsym_'+sample ].Fill( listMassAsym[0], SF )
			allHistos[ 'massAveVsMinMassAsym_'+sample ].Fill( varMassAsym[0], listMassAsym[0], SF )
			allHistos[ "massAve_minMass_cutBestPair_"+sample ].Fill( varMassAsym[0], SF )
			allHistos[ "deltaEta_minMass_cutBestPair_"+sample ].Fill( varMassAsym[4], SF )
			allHistos[ "massAsym_minMass_cutBestPair_"+sample ].Fill( varMassAsym[5], SF )
			allHistos[ 'deltavsMassAve_minMass_cutBestPair_'+sample ].Fill( varMassAsym[0], varMassAsym[8], SF )
			allHistos[ 'deltavsMassAve_minMass_cutBestPair_'+sample ].Fill( varMassAsym[0], varMassAsym[9], SF )
			allHistos[ 'cosThetaStar_minMass_cutBestPair_'+sample ].Fill( varMassAsym[10], SF )
			allHistos[ 'massAveVscosThetaStar_minMass_cutBestPair_'+sample ].Fill( varMassAsym[0], varMassAsym[10], SF )
			allHistos[ 'xi_minMass_cutBestPair_'+sample ].Fill( varMassAsym[11], SF )
			allHistos[ 'massAveVsXi_minMass_cutBestPair_'+sample ].Fill( varMassAsym[0], varMassAsym[11], SF )
			if len(listMassAsym)>1: 
				allHistos[ 'secondMinMassAsym_'+sample ].Fill( listMassAsym[1], SF )
				allHistos[ 'minVsSecondMinMassAsym_'+sample ].Fill( listMassAsym[0], listMassAsym[1], SF )
			if ( varMassAsym[4] < 1.0 ):
				allHistos[ "massAve_minMass_cutDeltaEta_"+sample ].Fill( varMassAsym[0], SF )
				allHistos[ "deltaEta_minMass_cutDeltaEta_"+sample ].Fill( varMassAsym[4], SF )
				allHistos[ "massAsym_minMass_cutDeltaEta_"+sample ].Fill( varMassAsym[5], SF )
				allHistos[ 'deltavsMassAve_minMass_cutDeltaEta_'+sample ].Fill( varMassAsym[0], varMassAsym[8], SF )
				allHistos[ 'deltavsMassAve_minMass_cutDeltaEta_'+sample ].Fill( varMassAsym[0], varMassAsym[9], SF )
				if (( varMassAsym[8] > 200 ) and ( varMassAsym[9] > 200 )):
					allHistos[ "massAve_minMass_cutDelta_"+sample ].Fill( varMassAsym[0], SF )
					allHistos[ "deltaEta_minMass_cutDelta_"+sample ].Fill( varMassAsym[4], SF )
					allHistos[ "massAsym_minMass_cutDelta_"+sample ].Fill( varMassAsym[5], SF )
					allHistos[ 'deltavsMassAve_minMass_cutDelta_'+sample ].Fill( varMassAsym[0], varMassAsym[8], SF )
					allHistos[ 'deltavsMassAve_minMass_cutDelta_'+sample ].Fill( varMassAsym[0], varMassAsym[9], SF )
					allHistos[ 'cosThetaStar_minMass_cutDelta_'+sample ].Fill( varMassAsym[10], SF )
					allHistos[ 'massAveVscosThetaStar_minMass_cutDelta_'+sample ].Fill( varMassAsym[0], varMassAsym[10], SF )
					allHistos[ 'xi_minMass_cutDelta_'+sample ].Fill( varMassAsym[11], SF )
					allHistos[ 'massAveVsXi_minMass_cutDelta_'+sample ].Fill( varMassAsym[0], varMassAsym[11], SF )
					if ( varMassAsym[11] > .5 ):
						allHistos[ "massAve_minMass_cutXi_"+sample ].Fill( varMassAsym[0], SF )
				if (( varMassAsym[8] > 75 ) and ( varMassAsym[9] > 75 )):
					allHistos[ "massAve_minMass_cutDelta75_"+sample ].Fill( varMassAsym[0], SF )
				if (( varMassAsym[8] > 100 ) and ( varMassAsym[9] > 100 )):
					allHistos[ "massAve_minMass_cutDelta100_"+sample ].Fill( varMassAsym[0], SF )
				if (( varMassAsym[8] > 125 ) and ( varMassAsym[9] > 125 )):
					allHistos[ "massAve_minMass_cutDelta125_"+sample ].Fill( varMassAsym[0], SF )
				if (( varMassAsym[8] > 150 ) and ( varMassAsym[9] > 150 )):
					allHistos[ "massAve_minMass_cutDelta150_"+sample ].Fill( varMassAsym[0], SF )
				if (( varMassAsym[8] > 175 ) and ( varMassAsym[9] > 175 )):
					allHistos[ "massAve_minMass_cutDelta175_"+sample ].Fill( varMassAsym[0], SF )


			allHistos[ "minChi2vsminDeltaR_"+sample ].Fill( listMinChi[0], listDeltaR[0], SF )
			allHistos[ "massAveminChi2vsminDeltaR_"+sample ].Fill( varMinChi[0], varDeltaR[0], SF )

	outputFile.Write()
	##### Closing
	print 'Writing output file: '+ outputFileName
	outputFile.Close()



def bestPairing( listOfJets, possibleCombinations, method='', offset=0.8 ):
	"""docstring for pairing"""
	
	bestMin = 9999999999999
	if ( len(listOfJets) > 3 ):
		bestInd = ''
		listMin = []
		for a in possibleCombinations:
			if 'minChi2' in method:
				dijet1Mass = ( listOfJets[a[0]] + listOfJets[a[1]] ).M()	
				dijet2Mass = ( listOfJets[a[2]] + listOfJets[a[3]] ).M()	
				aveMass = (dijet1Mass + dijet2Mass )/2 ;
				sigma = ( (0.0724281 - (4.14372e-05 * aveMass) ) if (offset < 1) else 30 )
				tmpMin = TMath.Power( (dijet1Mass - dijet2Mass )/(aveMass if (offset<1) else 1 ), 2 ) / TMath.Power( sigma, 2 ) 
				#print tmpMin, TMath.Power( (dijet1Mass - dijet2Mass ), 2 ) , TMath.Power( sigma, 2 )
				#print dijet1Mass, dijet2Mass, tmpChi2
			elif 'deltaR' in method:
				tmpMin = abs( listOfJets[a[0]].DeltaR( listOfJets[a[1]] ) - offset ) + abs( listOfJets[a[2]].DeltaR( listOfJets[a[3]] ) - offset )
			elif 'mass' in method:
				dijet1Mass = ( listOfJets[a[0]] + listOfJets[a[1]] ).M()	
				dijet2Mass = ( listOfJets[a[2]] + listOfJets[a[3]] ).M()	
				tmpMin = abs( dijet1Mass - dijet2Mass ) / ( dijet1Mass + dijet2Mass )
			else: 
				tmpMin = 9999999
				print 'Error in pairing.'
				sys.exit(0)

			if( tmpMin < bestMin ):
				#print 'TMPMIN', tmpMin
				bestMin = tmpMin
				listMin.append( tmpMin )
				bestInd = a 
		reorderedJETS = [ listOfJets[ bestInd[0] ], listOfJets[ bestInd[1] ], listOfJets[ bestInd[2] ], listOfJets[ bestInd[3] ] ]
		listMin.sort()
		#print bestInd, listMin[0]
		return listMin, reorderedJETS

	else: return None



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
	dijetCosTheta = cosThetaStarDijet1 + cosThetaStarDijet2 
	deltaDijet1 = ( listJets[0].Pt() + listJets[1].Pt() ) - massAve
	deltaDijet2 = ( listJets[2].Pt() + listJets[3].Pt() ) - massAve
	xi1 = ( ( max( [ listJets[0].M() , listJets[1].M() ] ) / dijet1.M() ) * listJets[0].DeltaR( listJets[1])  )
	xi2 = ( ( max( [ listJets[2].M() , listJets[3].M() ] ) / dijet2.M() ) * listJets[2].DeltaR( listJets[3])  )
	dijetxi = ( xi1 + xi2 )
	#deltaRDijet1 = listJets[0].DeltaR( listJets[1] )
	#deltaRDijet2 = listJets[2].DeltaR( listJets[3] )

	return [ massAve, deltaEtaDijet1, deltaEtaDijet2, deltaEtaAveDijets, deltaEtaDijets, massAsymmetry, cosThetaStarDijet1, cosThetaStarDijet2, deltaDijet1, deltaDijet2, dijetCosTheta, dijetxi ]


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
			p = Process( target=myAnalyzer, args=( dictSamples[sample], preselection, cuts, sample, '' ) )
	p.start()
	p.join()


