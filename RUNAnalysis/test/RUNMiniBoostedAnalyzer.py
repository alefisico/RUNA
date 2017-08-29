#!/usr/bin/env python

'''
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
from random import randint
try: 
	from RUNA.RUNAnalysis.commonFunctions import *
	from RUNA.RUNAnalysis.cuts import selection
	from RUNA.RUNAnalysis.scaleFactors import scaleFactor
except ImportError: 
	sys.path.append('../python') 
	from commonFunctions import *
	from cuts import selection
	from scaleFactors import scaleFactor

gROOT.SetBatch()
######################################
def myPlotAnalyzer( fileSample, listCuts, sample, UNC ):

	outputFileName = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_'+sample+UNC+'_'+( '' if 'JetHT' in sample else 'Moriond17_')+'80X_V2p4_'+args.version+'p2.root' 
	outputFile = TFile( outputFileName, 'RECREATE' )


	################################################################################################## Histos
	massBins = 500
	massXmin = 0.
	massXmax = 500.
	listOfOptions = [ [ j,k] for j in range(len(listCuts)-1) for k in range(1, len(listCuts) ) if k > j ]

	print '--- Sample ', sample
	sf = scaleFactor(sample)
	if 'JetHT' in sample: sample = 'JetHT_Run2016'
	#elif 'QCD_HT' in sample: sample = 'QCDHTAll'
	#elif 'QCD_Pt' in sample: sample = 'QCDPtAll'
	
	allHistos[ "massAve_preSel_"+sample ] = TH1F( "massAve_preSel_"+sample, "massAve_preSel_"+sample, massBins, massXmin, massXmax )
	allHistos[ "HT_preSel_"+sample ] = TH1F( "HT_preSel_"+sample, "HT_preSel_"+sample, 5000, 0, 5000 )
	allHistos[ "numJets_preSel_"+sample ] = TH1F( "numJets_preSel_"+sample, "numJets_preSel_"+sample, 10, 0, 10)
	allHistos[ "deltaEtaDijet_preSel_"+sample ] = TH1F( "deltaEtaDijet_preSel_"+sample, "deltaEtaDijet_preSel_"+sample, 50, 0., 5 )
	allHistos[ "prunedMassAsym_preSel_"+sample ] = TH1F( "prunedMassAsym_preSel_"+sample, "prunedMassAsym_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau21_preSel_"+sample ] = TH1F( "jet1Tau21_preSel_"+sample, "jet1Tau21_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau21_preSel_"+sample ] = TH1F( "jet2Tau21_preSel_"+sample, "jet2Tau21_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau32_preSel_"+sample ] = TH1F( "jet1Tau32_preSel_"+sample, "jet1Tau32_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau32_preSel_"+sample ] = TH1F( "jet2Tau32_preSel_"+sample, "jet2Tau32_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet1btagCSVv2_preSel_"+sample ] = TH1F( "jet1btagCSVv2_preSel_"+sample, "jet1btagCSVv2_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet2btagCSVv2_preSel_"+sample ] = TH1F( "jet2btagCSVv2_preSel_"+sample, "jet2btagCSVv2_preSel_"+sample, 20, 0., 1 )

	allHistos[ "deltaEtaDijet_n-1_"+sample ] = TH1F( "deltaEtaDijet_n-1_"+sample, "deltaEtaDijet_n-1_"+sample, 50, 0., 5 )
	allHistos[ "prunedMassAsym_n-1_"+sample ] = TH1F( "prunedMassAsym_n-1_"+sample, "prunedMassAsym_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau21_n-1_"+sample ] = TH1F( "jet1Tau21_n-1_"+sample, "jet1Tau21_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau21_n-1_"+sample ] = TH1F( "jet2Tau21_n-1_"+sample, "jet2Tau21_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau32_n-1_"+sample ] = TH1F( "jet1Tau32_n-1_"+sample, "jet1Tau32_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau32_n-1_"+sample ] = TH1F( "jet2Tau32_n-1_"+sample, "jet2Tau32_n-1_"+sample, 20, 0., 1 )
	allHistos[ "deltaEtaDijet_n-1_2btag_"+sample ] = TH1F( "deltaEtaDijet_n-1_2btag_"+sample, "deltaEtaDijet_n-1_2btag_"+sample, 50, 0., 5 )
	allHistos[ "prunedMassAsym_n-1_2btag_"+sample ] = TH1F( "prunedMassAsym_n-1_2btag_"+sample, "prunedMassAsym_n-1_2btag_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau21_n-1_2btag_"+sample ] = TH1F( "jet1Tau21_n-1_2btag_"+sample, "jet1Tau21_n-1_2btag_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau21_n-1_2btag_"+sample ] = TH1F( "jet2Tau21_n-1_2btag_"+sample, "jet2Tau21_n-1_2btag_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau32_n-1_2btag_"+sample ] = TH1F( "jet1Tau32_n-1_2btag_"+sample, "jet1Tau32_n-1_2btag_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau32_n-1_2btag_"+sample ] = TH1F( "jet2Tau32_n-1_2btag_"+sample, "jet2Tau32_n-1_2btag_"+sample, 20, 0., 1 )
	allHistos[ "jet1btagCSVv2_n-1_2btag_"+sample ] = TH1F( "jet1btagCSVv2_n-1_2btag_"+sample, "jet1btagCSVv2_n-1_2btag_"+sample, 20, 0., 1 )
	allHistos[ "jet2btagCSVv2_n-1_2btag_"+sample ] = TH1F( "jet2btagCSVv2_n-1_2btag_"+sample, "jet2btagCSVv2_n-1_2btag_"+sample, 20, 0., 1 )
	listCuts.append( [ '1btag' ] )
	listCuts.append( [ '2btag' ] )

	for var in listCuts:
		if 'deltaEta' in var[0]: 
			allHistos[ var[0]+'_'+sample ] = TH1F( var[0]+'_'+sample, var[0]+'_'+sample, 50, 0., 5. )
			for var1 in listCuts: allHistos[ var[0]+'_'+var1[0]+"_"+sample ] = TH1F( var[0]+'_'+var1[0]+"_"+sample, var[0]+'_'+var1[0]+"_"+sample, 50, 0., 5. )
		else: 
			allHistos[ var[0]+'_'+sample ] = TH1F( var[0]+'_'+sample, var[0]+'_'+sample, 20, 0., 1. )
			for var1 in listCuts: allHistos[ var[0]+'_'+var1[0]+"_"+sample ] = TH1F( var[0]+'_'+var1[0]+"_"+sample, var[0]+'_'+var1[0]+"_"+sample, 20, 0., 1. )
		allHistos[ "massAve_"+var[0]+'_'+sample ] = TH1F( "massAve_"+var[0]+'_'+sample, "massAve_"+var[0]+'_'+sample, massBins, massXmin, massXmax )

		#allHistos[ "HT_"+var[0]+"_"+sample ] = TH1F( "HT_"+var[0]+"_"+sample, "HT_"+var[0]+"_"+sample, 5000, 0., 5000 )
		#allHistos[ "MET_"+var[0]+"_"+sample ] = TH1F( "MET_"+var[0]+"_"+sample, "MET_"+var[0]+"_"+sample, 500, 0., 500 )
		#allHistos[ "numJets_"+var[0]+"_"+sample ] = TH1F( "numJets_"+var[0]+"_"+sample, "numJets_"+var[0]+"_"+sample, 20, 0., 20 )
		allHistos[ "jet1Pt_"+var[0]+"_"+sample ] = TH1F( "jet1Pt_"+var[0]+"_"+sample, "jet1Pt_"+var[0]+"_"+sample, 2000, 0., 2000 )
		allHistos[ "jet2Pt_"+var[0]+"_"+sample ] = TH1F( "jet2Pt_"+var[0]+"_"+sample, "jet2Pt_"+var[0]+"_"+sample, 2000, 0., 2000 )
	allHistos[ 'massAve_jet2Tau32WOTau21_'+sample ] = TH1F( 'massAve_jet2Tau32WOTau21_'+sample, 'massAve_jet2Tau32WOTau21_'+sample, massBins, massXmin, massXmax )
	listCuts.remove( ['1btag'] )
	listCuts.remove( ['2btag'] )

	for ind in listOfOptions:
		tmpName = listCuts[ind[0]][0]+'Vs'+listCuts[ind[1]][0]+'_'+sample
		allHistos[ tmpName ] = TH2F( tmpName, tmpName, 
				(50 if 'deltaEta' in listCuts[ind[0]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[0]][0] else 1. ),
				(50 if 'deltaEta' in listCuts[ind[1]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[1]][0] else 1. ) 
				)

	tmpNameSam = listCuts[-2][0]+'Vs'+listCuts[-1][0]+'_'+sample
	for k in [ 'A', 'B', 'C', 'D' ]:
		allHistos[ "massAve_"+tmpNameSam+'_'+k ] = TH1F( "massAve_"+tmpNameSam+'_'+k, "massAve_"+tmpNameSam+'_'+k, massBins, massXmin, massXmax )
		allHistos[ "massAve_"+tmpNameSam+'_1btag_'+k ] = TH1F( "massAve_"+tmpNameSam+'_1btag_'+k, "massAve_"+tmpNameSam+'_1btag_'+k, massBins, massXmin, massXmax )
		allHistos[ "massAve_"+tmpNameSam+'_2btag_'+k ] = TH1F( "massAve_"+tmpNameSam+'_2btag_'+k, "massAve_"+tmpNameSam+'_2btag_'+k, massBins, massXmin, massXmax )
		allHistos[ tmpNameSam+'_'+k ] = TH2F( tmpNameSam+'_'+k, tmpNameSam+'_'+k, 
				(50 if 'deltaEta' in listCuts[-2][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-2][0] else 1. ),
				(50 if 'deltaEta' in listCuts[-1][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-1][0] else 1. ) 
				)
		allHistos[ tmpNameSam+'_1btag_'+k ] = TH2F( tmpNameSam+'_1btag_'+k, tmpNameSam+'_1btag_'+k, 
				(50 if 'deltaEta' in listCuts[-2][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-2][0] else 1. ),
				(50 if 'deltaEta' in listCuts[-1][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-1][0] else 1. ) 
				)

		allHistos[ tmpNameSam+'_2btag_'+k ] = TH2F( tmpNameSam+'_2btag_'+k, tmpNameSam+'_2btag_'+k, 
				(50 if 'deltaEta' in listCuts[-2][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-2][0] else 1. ),
				(50 if 'deltaEta' in listCuts[-1][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-1][0] else 1. ) 
				)

	for h in allHistos: allHistos[h].Sumw2()

	################################################################################################## Running the Analysis
	print '-'*40
	#lumiWeight = ( 1 if 'JetHT' in sample else sf  )
	SF = TCut( '1' if 'JetHT' in sample else 'lumiWeight * puWeight') 
	preselection = TCut('HT>900') + TCut("numJets==2") 
	stringSel = '' 
	for var in listCuts: stringSel = stringSel+'('+var[0]+('>' if '32' in var[0] else '<')+str(var[1])+')'
	stringSel = stringSel.replace(')(',') && (')
	print '---- Cuts applied: ', stringSel, preselection + TCut( stringSel )

	ABCDRegions = {}
	ABCDRegions[ '_A' ] = preselection + TCut( stringSel )
	ABCDRegions[ '_B' ] = preselection + TCut( stringSel.replace('deltaEtaDijet<', 'deltaEtaDijet>') )
	ABCDRegions[ '_C' ] = preselection + TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>') )
	ABCDRegions[ '_D' ] = preselection + TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>').replace('deltaEtaDijet<', 'deltaEtaDijet>') )

	#btag1Selection = TCut('(jet1btagCSVv2 > 0.8484) || (jet2btagCSVv2 > 0.8484)')
	btag1Selection = TCut('(jet1btagCSVv2 > 0.5426) || (jet2btagCSVv2 > 0.5426)')
	ABCDRegions1Btag = {}
	ABCDRegions1Btag[ '_A' ] = preselection + TCut( stringSel ) + TCut( btag1Selection )
	ABCDRegions1Btag[ '_B' ] = preselection + TCut( stringSel.replace('deltaEtaDijet<', 'deltaEtaDijet>') ) + TCut( btag1Selection )
	ABCDRegions1Btag[ '_C' ] = preselection + TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>') ) + TCut( btag1Selection )
	ABCDRegions1Btag[ '_D' ] = preselection + TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>').replace('deltaEtaDijet<', 'deltaEtaDijet>') ) + TCut( btag1Selection )

	#btag2Selection = '(jet1btagCSVv2 > 0.8484) && (jet2btagCSVv2 > 0.8484)'
	btag2Selection = '(jet1btagCSVv2 > 0.5426) && (jet2btagCSVv2 > 0.5426)'
	ABCDRegions2Btag = {}
	ABCDRegions2Btag[ '_A' ] = preselection + TCut( stringSel ) + TCut( btag2Selection )
	ABCDRegions2Btag[ '_B' ] = preselection + TCut( stringSel.replace('deltaEtaDijet<', 'deltaEtaDijet>') ) + TCut( btag2Selection )
	ABCDRegions2Btag[ '_C' ] = preselection + TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>') ) + TCut( btag2Selection )
	ABCDRegions2Btag[ '_D' ] = preselection + TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>').replace('deltaEtaDijet<', 'deltaEtaDijet>') ) + TCut( btag2Selection )

	sel = preselection + TCut( stringSel )
	btag2Sel = sel + TCut( btag2Selection ) + TCut( 'btagWeight' )
	btag1Sel = sel + TCut( btag1Selection ) + TCut( 'btagWeight' )

	treeName = 'BoostedAnalysisPlots'+('Puppi' if 'Puppi' in args.grooming else '')+'/RUNATree'


	### Preselection
	getHistoFromTree( fileSample, treeName,
			'HT', 
			SF,
			preselection,
			allHistos[ 'HT_preSel_'+sample ], 
			1 )

	getHistoFromTree( fileSample, treeName,
			'numJets', 
			SF,
			preselection,
			allHistos[ 'numJets_preSel_'+sample ], 
			1 )
	
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			preselection, 
			allHistos[ 'massAve_preSel_'+sample ], 
			1 )

	getHistoFromTree( fileSample, treeName,
			'prunedMassAsym', 
			SF,
			preselection,
			allHistos[ 'prunedMassAsym_preSel_'+sample ], 
			1 )

	getHistoFromTree( fileSample, treeName,
			'deltaEtaDijet', 
			SF,
			preselection,
			allHistos[ 'deltaEtaDijet_preSel_'+sample ], 
			1 )

	getHistoFromTree( fileSample, treeName,
			'jet1Tau21', 
			SF,
			preselection,
			allHistos[ 'jet1Tau21_preSel_'+sample ], 
			1 )

	getHistoFromTree( fileSample, treeName,
			'jet2Tau21', 
			SF,
			preselection,
			allHistos[ 'jet2Tau21_preSel_'+sample ], 
			1 )

	getHistoFromTree( fileSample, treeName,
			'jet1Tau32', 
			SF,
			preselection,
			allHistos[ 'jet1Tau32_preSel_'+sample ], 
			1 )

	getHistoFromTree( fileSample, treeName,
			'jet2Tau32', 
			SF,
			preselection,
			allHistos[ 'jet2Tau32_preSel_'+sample ], 
			1 )
	getHistoFromTree( fileSample, treeName,
			'jet1btagCSVv2', 
			SF,
			preselection,
			allHistos[ 'jet1btagCSVv2_preSel_'+sample ], 
			1 )
	getHistoFromTree( fileSample, treeName,
			'jet2btagCSVv2', 
			SF,
			preselection,
			allHistos[ 'jet2btagCSVv2_preSel_'+sample ], 
			1 )

	### All selection
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			sel, 
			allHistos[ 'massAve_deltaEtaDijet_'+sample ], 
			1 ) #( 0.10 if 'JetHT' in sample else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'jet1Pt', 
			SF,
			sel, 
			allHistos[ 'jet1Pt_deltaEtaDijet_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'jet2Pt', 
			SF,
			sel, 
			allHistos[ 'jet2Pt_deltaEtaDijet_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) ) 

	### Partial selection
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			preselection + TCut('(jet1Tau21<0.45) && (jet2Tau21<0.45)'), 
			allHistos[ 'massAve_jet2Tau21_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			preselection * TCut('(jet1Tau21<0.45) && (jet2Tau21<0.45) && (jet1Tau32>0.57) && (jet2Tau32>0.57)'), 
			allHistos[ 'massAve_jet1Tau21_'+sample ], #### just the label
			( 0.10 if 'JetHT' in sample else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			preselection + TCut('(jet1Tau21<0.45) && (jet2Tau21<0.45) && (jet1Tau32>0.57) && (jet2Tau32>0.57) && (prunedMassAsym<0.1)'), 
			allHistos[ 'massAve_prunedMassAsym_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) ) 


	### ttbar selection inclusive
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			preselection + TCut( stringSel.replace('(jet1Tau32>0.57)','(jet1Tau32<0.57)').replace('(jet2Tau32>0.57)','(jet2Tau32<0.57)') ), 
			#preselection + TCut( stringSel ) + TCut( '(jet1Tau32<0.57) && (jet2Tau32<0.57)') , 
			allHistos[ 'massAve_jet2Tau32_'+sample ], 
			1 ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			preselection + TCut( stringSel.replace('(jet1Tau32>0.57)','(jet1Tau32<0.57)').replace('(jet2Tau32>0.57)','(jet2Tau32<0.57)').replace('(jet1Tau21<0.45)','(jet1Tau21>0.45)').replace('(jet2Tau21<0.45)','(jet2Tau21>0.45)') ), 
			#preselection + TCut( stringSel.replace('(jet1Tau21<0.45)','(jet1Tau21>0.45)').replace('(jet2Tau21<0.45)','(jet2Tau21>0.45)') ) + TCut( '(jet1Tau32<0.57) && (jet2Tau32<0.57)'), 
			allHistos[ 'massAve_jet1Tau32_'+sample ], 
			1 ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			preselection + TCut( stringSel.replace('(jet1Tau32>0.57)','(jet1Tau32<0.57)').replace('(jet2Tau32>0.57)','(jet2Tau32<0.57)').replace('&& (jet1Tau21<0.45)','').replace('&& (jet2Tau21<0.45)','') ), 
			#preselection + TCut( stringSel.replace('&& (jet1Tau21<0.45)','').replace('&& (jet2Tau21<0.45)','') ) + TCut( '(jet1Tau32<0.57) && (jet2Tau32<0.57)'), 
			allHistos[ 'massAve_jet2Tau32WOTau21_'+sample ], 
			1 ) 

	### n-1 selection
	getHistoFromTree( fileSample, treeName,
			'prunedMassAsym', 
			SF,
			preselection + TCut( stringSel.replace('&& (prunedMassAsym<0.1)','') ), 
			allHistos[ 'prunedMassAsym_n-1_'+sample ], 
			1 ) 

	getHistoFromTree( fileSample, treeName,
			'deltaEtaDijet', 
			SF,
			preselection + TCut( stringSel.replace(('&& (deltaEtaDijet<1.5)' if 'pruned' in args.grooming else '&& (deltaEtaDijet<1.0)' ),'') ), 
			allHistos[ 'deltaEtaDijet_n-1_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet1Tau21', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau21<0.45)','').replace('&& (jet1Tau21<0.45)','') ), 
			allHistos[ 'jet1Tau21_n-1_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet2Tau21', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau21<0.45)','').replace('&& (jet1Tau21<0.45)','') ), 
			allHistos[ 'jet2Tau21_n-1_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet1Tau32', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau32<0.57)','').replace('&& (jet1Tau32<0.57)','') ), 
			allHistos[ 'jet1Tau32_n-1_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet2Tau32', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau32<0.57)','').replace('&& (jet1Tau32<0.57)','') ), 
			allHistos[ 'jet2Tau32_n-1_'+sample ], 
			1 ) 

	### ABCD plots
	for region, selABCD in ABCDRegions.items():
		getHistoFromTree( fileSample, treeName,
				'prunedMassAve', 
				SF,
				selABCD, 
				allHistos[ 'massAve_prunedMassAsymVsdeltaEtaDijet_'+sample+region ], 
				( 0.10 if 'JetHT' in sample else 1 ) ) 
	
		get2DHistoFromTree( fileSample, treeName,
				'prunedMassAsym', 'deltaEtaDijet',
				SF,
				selABCD, 
				allHistos[ 'prunedMassAsymVsdeltaEtaDijet_'+sample+region ],
				1 ) 

	## Btagging
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			btag1Sel, 
			allHistos[ 'massAve_1btag_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) ) 

	for region, selABCD in ABCDRegions1Btag.items():
		getHistoFromTree( fileSample, treeName,
				'prunedMassAve', 
				selABCD, 
				SF,
				allHistos[ 'massAve_prunedMassAsymVsdeltaEtaDijet_'+sample+'_1btag'+region ], 
				( 0.10 if 'JetHT' in sample else 1 ) ) 
	
		get2DHistoFromTree( fileSample, treeName,
				'prunedMassAsym', 'deltaEtaDijet',
				SF,
				selABCD, 
				allHistos[ 'prunedMassAsymVsdeltaEtaDijet_'+sample+'_1btag'+region ],
				1 ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			SF,
			btag2Sel, 
			allHistos[ 'massAve_2btag_'+sample ], 
			( 0.10 if 'JetHT' in sample else 1 ) ) 

	for region, selABCD in ABCDRegions2Btag.items():
		getHistoFromTree( fileSample, treeName,
				'prunedMassAve', 
				SF,
				selABCD, 
				allHistos[ 'massAve_prunedMassAsymVsdeltaEtaDijet_'+sample+'_2btag'+region ], 
				( 0.10 if 'JetHT' in sample else 1 ) ) 
	
		get2DHistoFromTree( fileSample, treeName,
				'prunedMassAsym', 'deltaEtaDijet',
				SF,
				selABCD, 
				allHistos[ 'prunedMassAsymVsdeltaEtaDijet_'+sample+'_2btag'+region ],
				1 ) 

	### n-1 selection btagged
	getHistoFromTree( fileSample, treeName,
			'prunedMassAsym', 
			SF,
			preselection + TCut( stringSel.replace('&& (prunedMassAsym<0.1)','') ) + TCut( btag2Selection ), 
			allHistos[ 'prunedMassAsym_n-1_2btag_'+sample ], 
			1 ) 

	getHistoFromTree( fileSample, treeName,
			'deltaEtaDijet', 
			SF,
			preselection + TCut( stringSel.replace(('&& (deltaEtaDijet<1.5)' if 'pruned' in args.grooming else '&& (deltaEtaDijet<1.0)' ),'') ) + TCut( btag2Selection ), 
			allHistos[ 'deltaEtaDijet_n-1_2btag_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet1Tau21', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau21<0.45)','').replace('&& (jet1Tau21<0.45)','') ) + TCut( btag2Selection ), 
			allHistos[ 'jet1Tau21_n-1_2btag_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet2Tau21', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau21<0.45)','').replace('&& (jet1Tau21<0.45)','') ) + TCut( btag2Selection ), 
			allHistos[ 'jet2Tau21_n-1_2btag_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet1Tau32', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau32<0.57)','').replace('&& (jet1Tau32<0.57)','') ) + TCut( btag2Selection ), 
			allHistos[ 'jet1Tau32_n-1_2btag_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet2Tau32', 
			SF,
			preselection + TCut( stringSel.replace('&& (jet2Tau32<0.57)','').replace('&& (jet1Tau32<0.57)','') ) + TCut( btag2Selection ), 
			allHistos[ 'jet2Tau32_n-1_2btag_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet1btagCSVv2', 
			SF,
			sel,
			allHistos[ 'jet1btagCSVv2_n-1_2btag_'+sample ], 
			1 ) 
	getHistoFromTree( fileSample, treeName,
			'jet2btagCSVv2', 
			SF,
			sel,
			allHistos[ 'jet2btagCSVv2_n-1_2btag_'+sample ], 
			1 ) 


	outputFile.Write()
	##### Closing
	print 'Writing output file: '+ outputFileName
	outputFile.Close()


######################################
def myAnalyzer( fileSample, listCuts, sample, UNC ):

	outputFileName = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_'+sample+UNC+'_'+( '' if 'JetHT' in sample else 'Moriond17_')+'80X_V2p4_'+args.version+'p3.root' 
	outputFile = TFile( outputFileName, 'RECREATE' )


	################################################################################################## Histos
	massBins = 500
	massXmin = 0.
	massXmax = 500.
	listOfOptions = [ [ j,k] for j in range(len(listCuts)-1) for k in range(1, len(listCuts) ) if k > j ]

	#print '--- Sample ', sample
	#sf = scaleFactor(sample)
	if 'JetHT' in sample: sample = 'JetHT_Run2016'
	elif 'QCD_HT' in sample: sample = 'QCDHTAll'
	elif 'QCD_Pt' in sample: sample = 'QCDPtAll'
	allHistos[ "massAve_preSel_"+sample ] = TH1F( "massAve_preSel_"+sample, "massAve_preSel_"+sample, massBins, massXmin, massXmax )
	allHistos[ "HT_preSel_"+sample ] = TH1F( "HT_preSel_"+sample, "HT_preSel_"+sample, 5000, 0, 5000 )
	allHistos[ "deltaEtaDijet_preSel_"+sample ] = TH1F( "deltaEtaDijet_preSel_"+sample, "deltaEtaDijet_preSel_"+sample, 50, 0., 5 )
	allHistos[ "prunedMassAsym_preSel_"+sample ] = TH1F( "prunedMassAsym_preSel_"+sample, "prunedMassAsym_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau21_preSel_"+sample ] = TH1F( "jet1Tau21_preSel_"+sample, "jet1Tau21_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau21_preSel_"+sample ] = TH1F( "jet2Tau21_preSel_"+sample, "jet2Tau21_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau31_preSel_"+sample ] = TH1F( "jet1Tau31_preSel_"+sample, "jet1Tau31_preSel_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau31_preSel_"+sample ] = TH1F( "jet2Tau31_preSel_"+sample, "jet2Tau31_preSel_"+sample, 20, 0., 1 )

	allHistos[ "deltaEtaDijet_n-1_"+sample ] = TH1F( "deltaEtaDijet_n-1_"+sample, "deltaEtaDijet_n-1_"+sample, 50, 0., 5 )
	allHistos[ "prunedMassAsym_n-1_"+sample ] = TH1F( "prunedMassAsym_n-1_"+sample, "prunedMassAsym_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau21_n-1_"+sample ] = TH1F( "jet1Tau21_n-1_"+sample, "jet1Tau21_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau21_n-1_"+sample ] = TH1F( "jet2Tau21_n-1_"+sample, "jet2Tau21_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet1Tau31_n-1_"+sample ] = TH1F( "jet1Tau31_n-1_"+sample, "jet1Tau31_n-1_"+sample, 20, 0., 1 )
	allHistos[ "jet2Tau31_n-1_"+sample ] = TH1F( "jet2Tau31_n-1_"+sample, "jet2Tau31_n-1_"+sample, 20, 0., 1 )
	listCuts.append( [ 'btag' ] )
	for var in listCuts:
		if 'deltaEta' in var[0]: 
			allHistos[ var[0]+'_'+sample ] = TH1F( var[0]+'_'+sample, var[0]+'_'+sample, 50, 0., 5. )
			for var1 in listCuts: allHistos[ var[0]+'_'+var1[0]+"_"+sample ] = TH1F( var[0]+'_'+var1[0]+"_"+sample, var[0]+'_'+var1[0]+"_"+sample, 50, 0., 5. )
		else: 
			allHistos[ var[0]+'_'+sample ] = TH1F( var[0]+'_'+sample, var[0]+'_'+sample, 20, 0., 1. )
			for var1 in listCuts: allHistos[ var[0]+'_'+var1[0]+"_"+sample ] = TH1F( var[0]+'_'+var1[0]+"_"+sample, var[0]+'_'+var1[0]+"_"+sample, 20, 0., 1. )
		allHistos[ "massAve_"+var[0]+'_'+sample ] = TH1F( "massAve_"+var[0]+'_'+sample, "massAve_"+var[0]+'_'+sample, massBins, massXmin, massXmax )

		allHistos[ "HT_"+var[0]+"_"+sample ] = TH1F( "HT_"+var[0]+"_"+sample, "HT_"+var[0]+"_"+sample, 5000, 0., 5000 )
		allHistos[ "MET_"+var[0]+"_"+sample ] = TH1F( "MET_"+var[0]+"_"+sample, "MET_"+var[0]+"_"+sample, 500, 0., 500 )
		allHistos[ "numJets_"+var[0]+"_"+sample ] = TH1F( "numJets_"+var[0]+"_"+sample, "numJets_"+var[0]+"_"+sample, 20, 0., 20 )
		allHistos[ "jet1Pt_"+var[0]+"_"+sample ] = TH1F( "jet1Pt_"+var[0]+"_"+sample, "jet1Pt_"+var[0]+"_"+sample, 2000, 0., 2000 )
		allHistos[ "jet2Pt_"+var[0]+"_"+sample ] = TH1F( "jet2Pt_"+var[0]+"_"+sample, "jet2Pt_"+var[0]+"_"+sample, 2000, 0., 2000 )
	listCuts.remove( ['btag'] )

	for ind in listOfOptions:
		tmpName = listCuts[ind[0]][0]+'Vs'+listCuts[ind[1]][0]+'_'+sample
		allHistos[ tmpName ] = TH2F( tmpName, tmpName, 
				(50 if 'deltaEta' in listCuts[ind[0]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[0]][0] else 1. ),
				(50 if 'deltaEta' in listCuts[ind[1]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[1]][0] else 1. ) 
				)

	tmpNameSam = listCuts[-2][0]+'Vs'+listCuts[-1][0]+'_'+sample
	#tmpNameSam = tmpName #listCuts[-2][0]+'Vs'+listCuts[-1][0]+'_'+sample
	for k in [ 'A', 'B', 'C', 'D' ]:
		#allHistos[ "massAve_"+tmpNameSam+'_'+k ] = TH1F( "massAve_"+tmpNameSam+'_'+k, "massAve_"+tmpNameSam+'_'+k,  len(boostedMassAveBins)-1, boostedMassAveBins )
		#allHistos[ "massAve_"+tmpNameSam+'_btag_'+k ] = TH1F( "massAve_"+tmpNameSam+'_btag_'+k, "massAve_"+tmpNameSam+'_btag_'+k,  len(boostedMassAveBins)-1, boostedMassAveBins )
		allHistos[ "massAve_"+tmpNameSam+'_'+k ] = TH1F( "massAve_"+tmpNameSam+'_'+k, "massAve_"+tmpNameSam+'_'+k, massBins, massXmin, massXmax )
		allHistos[ "massAve_"+tmpNameSam+'_btag_'+k ] = TH1F( "massAve_"+tmpNameSam+'_btag_'+k, "massAve_"+tmpNameSam+'_btag_'+k, massBins, massXmin, massXmax )
		allHistos[ tmpNameSam+'_'+k ] = TH2F( tmpNameSam+'_'+k, tmpNameSam+'_'+k, 
				(50 if 'deltaEta' in listCuts[-2][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-2][0] else 1. ),
				(50 if 'deltaEta' in listCuts[-1][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-1][0] else 1. ) 
				#(50 if 'deltaEta' in listCuts[ind[0]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[0]][0] else 1. ),
				#(50 if 'deltaEta' in listCuts[ind[1]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[1]][0] else 1. ) 
				)
		allHistos[ tmpNameSam+'_btag_'+k ] = TH2F( tmpNameSam+'_btag_'+k, tmpNameSam+'_btag_'+k, 
				(50 if 'deltaEta' in listCuts[-2][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-2][0] else 1. ),
				(50 if 'deltaEta' in listCuts[-1][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[-1][0] else 1. ) 
				#(50 if 'deltaEta' in listCuts[ind[0]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[0]][0] else 1. ),
				#(50 if 'deltaEta' in listCuts[ind[1]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[1]][0] else 1. ) 
				)

	for h in allHistos: allHistos[h].Sumw2()

	####### Get GenTree 
	inputFile, events, numEntries = getTree( fileSample, ('BoostedAnalysisPlotsPuppi'+( '' if 'PDF' in UNC else UNC)+'/RUNATree' if 'Puppi' in args.grooming else 'BoostedAnalysisPlots'+( '' if 'PDF' in UNC else UNC)+'/RUNATree' ) )
	print '-'*40
	print '------> ', sample
	print '------> Number of events: '+str(numEntries)
	d = 0
	cutFlowList = OrderedDict()
	cutFlowScaledList = OrderedDict()
	cutFlowScaledListWeights = OrderedDict()
	cutFlowList[ 'Process' ] = 0
	cutFlowList[ 'Preselection' ] = 0
	cutFlowScaledList[ 'Process' ] = 0
	cutFlowScaledList[ 'Preselection' ] = 0
	cutFlowScaledListWeights[ 'Process' ] = 0
	cutFlowScaledListWeights[ 'Preselection' ] = 0
	for k in listCuts: 
		cutFlowList[ k[0] ] = 0
		cutFlowScaledList[ k[0] ] = 0
		cutFlowScaledListWeights[ k[0] ] = 0
	cutFlowList[ 'btag' ] = 0
	cutFlowScaledList[ 'btag' ] = 0
	cutFlowScaledListWeights[ 'btag' ] = 0

	for i in xrange(numEntries):
		events.GetEntry(i)

		#---- progress of the reading --------
		fraction = 10.*i/(1.*numEntries)
		if TMath.FloorNint(fraction) > d: print str(10*TMath.FloorNint(fraction))+'%' 
		d = TMath.FloorNint(fraction)
		#if ( i > 100000 ): break

		Run      = events.run
		Lumi     = events.lumi
		NumEvent = events.event
		puWeight	= events.puWeight
		pdfWeight	= events.pdfWeight
		lumiWeight	= events.lumiWeight
		HT		= events.HT
		MET		= events.MET
		numJets		= events.numJets
		massAve		= getattr( events, (args.grooming+"MassAve").replace('Puppi','') )
		jet1Mass	= getattr( events, 'jet1'+(args.grooming+"Mass").replace('pruned','Pruned').replace('soft','Soft').replace('Puppi',''))
		jet2Mass	= getattr( events, 'jet2'+(args.grooming+"Mass").replace('pruned','Pruned').replace('soft','Soft').replace('Puppi',''))
		jet1Pt          = events.jet1Pt
		jet2Pt          = events.jet2Pt
		jet1Eta          = events.jet1Eta
		jet2Eta          = events.jet2Eta
		#jet1CosThetaStar	= events.jet1CosThetaStar
		#jet2CosThetaStar	= events.jet2CosThetaStar
		jet1BtagCSV		= ( events.jet1btagCSVv2 > 0.5426 )
		jet2BtagCSV		= ( events.jet2btagCSVv2 > 0.5426 )
		
		#print 'Entry ', Run, ':', Lumi, ':', NumEvent

		if 'JetHT' in sample: scale = 1
		#elif 'RPV' in sample: scale = 2606 * puWeight * SF
		else: scale = args.lumi * puWeight * lumiWeight

		if 'PDF' in UNC:
			if 'Up' in UNC: scale = scale*(1+pdfWeight)
			else: scale = scale*(1-pdfWeight)

		cutFlowList[ 'Process' ] += 1
		cutFlowScaledList[ 'Process' ] += scale
		cutFlowScaledList[ 'Process' ] += (puWeight*puWeight)

		########## DDT
		##jet1RhoDDT = TMath.Log( jet1Mass*jet1Mass/jet1Pt )
		##jet2RhoDDT = TMath.Log( jet2Mass*jet2Mass/jet2Pt )
		##jet1Tau21DDT = events.jet1Tau21 + 0.063 * jet1RhoDDT 
		##jet2Tau21DDT = events.jet2Tau21 + 0.063 * jet2RhoDDT 
		
		#### Pre-selection
		HTCut = ( HT > 900 )
		dijetCut =  ( numJets == 2 )
		#jetPtCut =  ( jet1Pt > 500 ) and ( jet2Pt > 450 )
		jetPtCut =  ( jet1Pt > 150 ) and ( jet2Pt > 150 )
		
		#if HTCut and dijetCut and jetPtCut:
		if HTCut and dijetCut :
			if ( events.jet1Tau32 > 0.57 ) and  ( events.jet2Tau32 > 0.57 ) and ( events.jet1Tau21 < 0.45 ) and  ( events.jet2Tau21 < 0.45 ) and ( events.prunedMassAsym < 0.10 ) and  ( events.deltaEtaDijet < 1.5 ):
					allHistos[ 'massAve_deltaEtaDijet_'+sample ].Fill( massAve, scale )
					if ( jet1BtagCSV and jet2BtagCSV ): 
						print str(Run)+':'+str(Lumi)+':'+str(NumEvent)

#				cutFlowList[ 'Preselection' ] += 1
#				cutFlowScaledList[ 'Preselection' ] += scale
#				cutFlowScaledList[ 'Preselection' ] += (puWeight*puWeight)
#				sigCutsList = []
#				allHistos[ "HT_"+sample ].Fill( HT, scale )
#				allHistos[ "MET_"+sample ].Fill( MET, scale )
#				allHistos[ "massAve_"+sample ].Fill( massAve, scale )
#				allHistos[ "numJets_"+sample ].Fill( numJets, scale )
#				allHistos[ "jet1Pt_"+sample ].Fill( jet1Pt, scale )
#				allHistos[ "jet2Pt_"+sample ].Fill( jet2Pt, scale )
#				allHistos[ "jet1RhoDDT_"+sample ].Fill( jet1RhoDDT, scale )
#				allHistos[ "jet2RhoDDT_"+sample ].Fill( jet2RhoDDT, scale )
#				allHistos[ "jet1Tau21VsRhoDDT_"+sample ].Fill( events.jet1Tau21, jet1RhoDDT, scale )
#				allHistos[ "jet2Tau21VsRhoDDT_"+sample ].Fill( events.jet2Tau21, jet2RhoDDT, scale )
#				allHistos[ "jet1Tau21DDT_"+sample ].Fill( jet1Tau21DDT, scale )
#				allHistos[ "jet2Tau21DDT_"+sample ].Fill( jet2Tau21DDT, scale )
#				allHistos[ "jet1Tau21DDTVsRhoDDT_"+sample ].Fill( jet1Tau21DDT, jet1RhoDDT, scale )
#				allHistos[ "jet2Tau21DDTVsRhoDDT_"+sample ].Fill( jet2Tau21DDT, jet2RhoDDT, scale )
#				allHistos[ "prunedMassAsym_"+sample ].Fill( events.prunedMassAsym, scale )
#				allHistos[ "deltaEtaDijet_"+sample ].Fill( events.deltaEtaDijet, scale )
#				allHistos[ "jet1CosThetaStar_"+sample ].Fill( jet1CosThetaStar, scale )
#				allHistos[ "jet2CosThetaStar_"+sample ].Fill( jet2CosThetaStar, scale )
#				allHistos[ "jet1Tau21_"+sample ].Fill( events.jet1Tau21, scale )
#				allHistos[ "jet2Tau21_"+sample ].Fill( events.jet2Tau21, scale )
#				allHistos[ "jet1Tau31_"+sample ].Fill( events.jet1Tau31, scale )
#				allHistos[ "jet2Tau31_"+sample ].Fill( events.jet2Tau31, scale )
#				allHistos[ "jet1Tau32_"+sample ].Fill( events.jet1Tau32, scale )
#				allHistos[ "jet2Tau32_"+sample ].Fill( events.jet2Tau32, scale )
#				allHistos[ "jet1SubjetPtRatio_"+sample ].Fill( events.jet1SubjetPtRatio, scale )
#				allHistos[ "jet2SubjetPtRatio_"+sample ].Fill( events.jet2SubjetPtRatio, scale )
#				allHistos[ "jet1BtagCSV_"+sample ].Fill( 1 if jet1BtagCSV else 0 )
#				allHistos[ "jet2BtagCSV_"+sample ].Fill( 1 if jet1BtagCSV else 0 )
#
#				bothBtag = ( jet1BtagCSV and jet2BtagCSV )
#				oneBtag = ( jet1BtagCSV or jet2BtagCSV )
#				if bothBtag: allHistos[ "jetsBtagCSV_"+sample ].Fill( 2 )
#				elif oneBtag: allHistos[ "jetsBtagCSV_"+sample ].Fill( 1 )
#				else: allHistos[ "jetsBtagCSV_"+sample ].Fill( 0 )
#
#				for var in listCuts:
#					#allHistos[ var[0]+'_'+sample ].Fill( getattr( events, var[0] ), scale )
#					nextCut = False
#					if ( getattr( events, var[0] ) < var[1] ): nextCut = True 
#					else: nextCut = False
#					sigCutsList.append( nextCut )
#
#					if all(sigCutsList): 
#						allHistos[ 'massAve_'+var[0]+'_'+sample ].Fill( massAve, scale )  ### adding two prong scale factor
#						allHistos[ 'jet1Tau21_'+var[0]+'_'+sample ].Fill( events.jet1Tau21, scale )
#						allHistos[ 'jet2Tau21_'+var[0]+'_'+sample ].Fill( events.jet2Tau21, scale )
#						allHistos[ 'prunedMassAsym_'+var[0]+'_'+sample ].Fill( events.prunedMassAsym, scale )
#						allHistos[ 'deltaEtaDijet_'+var[0]+'_'+sample ].Fill( events.deltaEtaDijet, scale )
#						allHistos[ "HT_"+var[0]+"_"+sample ].Fill( HT, scale )
#						allHistos[ "MET_"+var[0]+"_"+sample ].Fill( MET, scale )
#						allHistos[ "numJets_"+var[0]+"_"+sample ].Fill( numJets, scale )
#						allHistos[ "jet1Pt_"+var[0]+"_"+sample ].Fill( jet1Pt, scale )
#						allHistos[ "jet2Pt_"+var[0]+"_"+sample ].Fill( jet2Pt, scale )
#						cutFlowList[ var[0] ] += 1
#						cutFlowScaledList[ var[0] ] += scale
#						cutFlowScaledList[ var[0] ] += (puWeight*puWeight)
#
#				if oneBtag and all(sigCutsList): 
#					allHistos[ 'massAve_btag_'+sample ].Fill( massAve, scale )  ### adding two prong scale factor
#					allHistos[ 'jet1Tau21_btag_'+sample ].Fill( events.jet1Tau21, scale )
#					allHistos[ 'jet2Tau21_btag_'+sample ].Fill( events.jet2Tau21, scale )
#					allHistos[ 'prunedMassAsym_btag_'+sample ].Fill( events.prunedMassAsym, scale )
#					allHistos[ 'deltaEtaDijet_btag_'+sample ].Fill( events.deltaEtaDijet, scale )
#					allHistos[ "HT_btag_"+sample ].Fill( HT, scale )
#					allHistos[ "MET_btag_"+sample ].Fill( MET, scale )
#					allHistos[ "numJets_btag_"+sample ].Fill( numJets, scale )
#					allHistos[ "jet1Pt_btag_"+sample ].Fill( jet1Pt, scale )
#					allHistos[ "jet2Pt_btag_"+sample ].Fill( jet2Pt, scale )
#					cutFlowList[ 'btag' ] += 1
#					cutFlowScaledList[ 'btag' ] += scale
#					cutFlowScaledList[ 'btag' ] += (puWeight*puWeight)
#				#### n-1 plots
#				'''
#				if ( 'low' in args.RANGE ):
#					if ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ) and (  getattr( events, listCuts[1][0] ) < listCuts[1][1] ) and ( getattr( events, listCuts[2][0] ) < listCuts[2][1] ) and ( getattr( events, listCuts[3][0] ) < listCuts[3][1] ) and ( getattr( events, listCuts[4][0] ) < listCuts[4][1] ): allHistos[ 'deltaEtaDijet_n-1_'+sample ].Fill( events.deltaEtaDijet, scale )
#					if ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ) and (  getattr( events, listCuts[1][0] ) < listCuts[1][1] ) and ( getattr( events, listCuts[2][0] ) < listCuts[2][1] ) and ( getattr( events, listCuts[3][0] ) < listCuts[3][1] ) and ( getattr( events, listCuts[5][0] ) < listCuts[5][1] ): allHistos[ 'prunedMassAsym_n-1_'+sample ].Fill( events.prunedMassAsym, scale )
#					if ( getattr( events, listCuts[2][0] ) < listCuts[2][1] ) and ( getattr( events, listCuts[3][0] ) < listCuts[3][1] ) and ( getattr( events, listCuts[4][0] ) < listCuts[4][1] ) and ( getattr( events, listCuts[5][0] ) < listCuts[5][1] ): 
#						allHistos[ 'jet1Tau21_n-1_'+sample ].Fill( events.jet1Tau21, scale )
#						allHistos[ 'jet2Tau21_n-1_'+sample ].Fill( events.jet2Tau21, scale )
#					if ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ) and ( getattr( events, listCuts[1][0] ) < listCuts[1][1] ) and ( getattr( events, listCuts[4][0] ) < listCuts[4][1] ) and ( getattr( events, listCuts[5][0] ) < listCuts[5][1] ): 
#						allHistos[ 'jet1Tau31_n-1_'+sample ].Fill( events.jet1Tau31, scale )
#						allHistos[ 'jet2Tau31_n-1_'+sample ].Fill( events.jet2Tau31, scale )
#				else:
#				'''
#				if ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ) and (  getattr( events, listCuts[1][0] ) < listCuts[1][1] ) and ( getattr( events, listCuts[3][0] ) < listCuts[3][1] ): allHistos[ 'prunedMassAsym_n-1_'+sample ].Fill( events.prunedMassAsym, scale )
#				if ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ) and (  getattr( events, listCuts[1][0] ) < listCuts[1][1] ) and ( getattr( events, listCuts[2][0] ) < listCuts[2][1] ): allHistos[ 'deltaEtaDijet_n-1_'+sample ].Fill( events.deltaEtaDijet, scale )
#				if ( getattr( events, listCuts[2][0] ) < listCuts[2][1] ) and ( getattr( events, listCuts[3][0] ) < listCuts[3][1] ): 
#					allHistos[ 'jet1Tau21_n-1_'+sample ].Fill( events.jet1Tau21, scale )
#					allHistos[ 'jet2Tau21_n-1_'+sample ].Fill( events.jet2Tau21, scale )
#
#				##########
#
#				for Ind in listOfOptions:
#					allHistos[ listCuts[Ind[0]][0]+'Vs'+listCuts[Ind[1]][0]+'_'+sample ].Fill( getattr( events, listCuts[Ind[0]][0] ), getattr( events, listCuts[Ind[1]][0] ), scale )
#					tmpSigCutsList = [ x for i,x in enumerate(sigCutsList) if i not in Ind ]
#					
#				##### Bkg estimation/ABCD method
#				if ( all(sigCutsList[:-2]) ): # and ( getattr( events, listCuts[5][0] ) > (listCuts[5][1]*2) )): 
#					allHistos[ listCuts[-2][0]+'Vs'+listCuts[-1][0]+'_'+sample+'_Bkg' ].Fill( getattr( events, listCuts[0][0] ), getattr( events, listCuts[1][0] ), scale )
#					plotABCD( [ ( getattr( events, listCuts[-2][0] ) < listCuts[-2][1] ), ( getattr( events, listCuts[-1][0] ) < listCuts[-1][1] ) ], [ listCuts[-2][0], listCuts[-1][0] ], events, massAve, scale, sample )
#
#				####### bkg estimation alternatives
#				if sigCutsList[2]: 
#					allHistos[ 'jet1Tau21VsdeltaEtaDijet_'+sample+'_Bkg' ].Fill( getattr( events, 'jet1Tau21' ), getattr( events, 'deltaEtaDijet' ), scale )
#					allHistos[ 'jet2Tau21VsdeltaEtaDijet_'+sample+'_Bkg' ].Fill( getattr( events, 'jet2Tau21' ), getattr( events, 'deltaEtaDijet' ), scale )
#					plotABCDv2( [ ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ), ( getattr( events, listCuts[1][0] ) < listCuts[1][1] ), ( getattr( events, listCuts[-1][0] ) < listCuts[-1][1] ) ], [ listCuts[0][0], listCuts[-1][0] ], events, massAve, scale, sample )
#					plotABCDv2( [ ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ), ( getattr( events, listCuts[1][0] ) < listCuts[1][1] ), ( getattr( events, listCuts[-1][0] ) < listCuts[-1][1] ) ], [ listCuts[1][0], listCuts[-1][0] ], events, massAve, scale, sample )
#
#				if sigCutsList[-1]: 
#					allHistos[ 'jet1Tau21VsprunedMassAsym_'+sample+'_Bkg' ].Fill( getattr( events, 'jet1Tau21' ), getattr( events, 'prunedMassAsym' ), scale )
#					allHistos[ 'jet2Tau21VsprunedMassAsym_'+sample+'_Bkg' ].Fill( getattr( events, 'jet2Tau21' ), getattr( events, 'prunedMassAsym' ), scale )
#					plotABCDv2( [ ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ), ( getattr( events, listCuts[1][0] ) < listCuts[1][1] ), ( getattr( events, listCuts[-2][0] ) < listCuts[-2][1] ) ], [ listCuts[0][0], listCuts[-2][0] ], events, massAve, scale, sample )
#					plotABCDv2( [ ( getattr( events, listCuts[0][0] ) < listCuts[0][1] ), ( getattr( events, listCuts[1][0] ) < listCuts[1][1] ), ( getattr( events, listCuts[-2][0] ) < listCuts[-2][1] ) ], [ listCuts[1][0], listCuts[-2][0] ], events, massAve, scale, sample )
#
#						
#
#		dummy = 1
#		for q in cutFlowList: 
#			allHistos[ 'cutFlow_'+sample ].SetBinContent( dummy, cutFlowList[q] )
#			allHistos[ 'cutFlow_'+sample ].GetXaxis().SetBinLabel( dummy, q )
#			allHistos[ 'cutFlow_Scaled_'+sample ].SetBinContent( dummy, cutFlowScaledList[q] )
#			allHistos[ 'cutFlow_Scaled_'+sample ].GetXaxis().SetBinLabel( dummy, q )
#			allHistos[ 'cutFlow_Scaled_Weights_'+sample ].SetBinContent( dummy, cutFlowScaledListWeights[q] )
#			allHistos[ 'cutFlow_Scaled_Weights_'+sample ].GetXaxis().SetBinLabel( dummy, q )
#			dummy+=1
#
#	for sample in dictSamples:
#		nameABCD = listCuts[-2][0]+'Vs'+listCuts[-1][0]+'_'+sample
#		allHistos[ 'massAve_'+nameABCD+'_BC' ].Multiply( allHistos[ 'massAve_'+nameABCD+'_B' ], allHistos[ 'massAve_'+nameABCD+'_C' ], 1, 1, '')
#		allHistos[ 'massAve_'+nameABCD+'_ABCDProj' ].Divide( allHistos[ 'massAve_'+nameABCD+'_BC' ], allHistos[ 'massAve_'+nameABCD+'_D' ], 1, 1, '')
#		'''
#		### The two lines above are doing exactly the following:
#		for ibin in range( 0, allHistos[ 'massAve_'+nameABCD+'_B' ].GetNbinsX() ):
#			Bcont = allHistos[ 'massAve_'+nameABCD+'_B' ].GetBinContent( ibin )
#			Berr = allHistos[ 'massAve_'+nameABCD+'_B' ].GetBinError( ibin )
#			Ccont = allHistos[ 'massAve_'+nameABCD+'_C' ].GetBinContent( ibin )
#			Cerr = allHistos[ 'massAve_'+nameABCD+'_C' ].GetBinError( ibin )
#			Dcont = allHistos[ 'massAve_'+nameABCD+'_D' ].GetBinContent( ibin )
#			Derr = allHistos[ 'massAve_'+nameABCD+'_D' ].GetBinError( ibin )
#
#			try: Nbkg = ( Bcont * Ccont ) / Dcont
#			except ZeroDivisionError: Nbkg = 0
#			allHistos[ "massAve_"+nameABCD+'_ABCDProj' ].SetBinContent( ibin, Nbkg )
#			#try: NbkgErr = Nbkg * TMath.Sqrt( TMath.Power( Berr / Bcont, 2 ) + TMath.Power( Cerr / Ccont, 2 ) + TMath.Power( Derr / Dcont, 2 ) )
#			try: NbkgErr = Nbkg * TMath.Sqrt( TMath.Power( TMath.Sqrt(Bcont) / Bcont, 2 ) + TMath.Power( TMath.Sqrt(Ccont) / Ccont, 2 ) + TMath.Power( TMath.Sqrt(Dcont) / Dcont, 2 ) )
#			except ZeroDivisionError: NbkgErr = 0
#			allHistos[ "massAve_"+nameABCD+'_ABCDProj' ].SetBinError( ibin, NbkgErr )
#		'''


	outputFile.Write()
	##### Closing
	print 'Writing output file: '+ outputFileName
	outputFile.Close()


def plotABCD( listSel, var, fromTree, massAve, scale, sample ):
	"""docstring for plotABCD"""

	nameABCD = var[0]+'Vs'+var[1]+'_'+sample
	if listSel[0] and listSel[1]: 
		allHistos[ 'massAve_'+nameABCD+'_A' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_A' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )
	elif listSel[0] and not listSel[1]: 
		allHistos[ 'massAve_'+nameABCD+'_B' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_B' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )
	elif not listSel[0] and listSel[1]: 
		allHistos[ 'massAve_'+nameABCD+'_C' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_C' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )
	else:
		allHistos[ 'massAve_'+nameABCD+'_D' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_D' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )

def plotABCDv2( listSel, var, fromTree, massAve, scale, sample ):
	"""docstring for plotABCD"""

	nameABCD = var[0]+'Vs'+var[1]+'_'+sample
	if (listSel[0] and listSel[1]) and listSel[2]: 
		allHistos[ 'massAve_'+nameABCD+'_A' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_A' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )
	elif (listSel[0] and listSel[1]) and not listSel[2]: 
		allHistos[ 'massAve_'+nameABCD+'_B' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_B' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )
	elif not listSel[0] and not listSel[1] and listSel[2]: 
		allHistos[ 'massAve_'+nameABCD+'_C' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_C' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )
	elif not listSel[0] and not listSel[1] and not listSel[2]: 
		allHistos[ 'massAve_'+nameABCD+'_D' ].Fill( massAve, scale )
		allHistos[ nameABCD+'_D' ].Fill( getattr( fromTree, var[0] ), getattr( fromTree, var[1] ), scale )




#################################################################################
if __name__ == '__main__':

	usage = 'usage: %prog [options]'
	
	parser = argparse.ArgumentParser()
	parser.add_argument( '-m', '--mass', action='store', dest='mass', default=100, help='Mass of the Stop' )
	parser.add_argument( '-g', '--grooming', action='store',  dest='grooming', default='pruned', help='Jet Algorithm' )
	parser.add_argument( '-p', '--process', action='store',  dest='process', default='Plots', help='Process: Plots or run.' )
	parser.add_argument( '-d', '--decay', action='store',  dest='decay', default='UDD312', help='Decay: UDD312 or UDD323.' )
	parser.add_argument( '-s', '--sample', action='store',   dest='samples', default='RPV', help='Type of sample' )
	parser.add_argument( '-u', '--unc', action='store',  dest='unc', default='', help='Process: all or single.' )
	parser.add_argument( '-q', '--qcd', action='store', default='Pt', dest='qcd', help='Type of QCD binning, example: HT.' )
	parser.add_argument( '-l', '--lumi', action='store', type=float, default=1787, help='Luminosity, example: 1.' )
	parser.add_argument( '-b', '--batchSys', action='store_true',  dest='batchSys', default=False, help='Process: all or single.' )
	parser.add_argument( '-v', '--version', action='store', default='v05', dest='version', help='Version of the RUNAnalysis file.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	if 'Pt' in args.qcd: QCDSF = ( 0.86 if 'Puppi' in args.grooming else 0.89 ) 
	else: QCDSF = 1

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
	allSamples[ 'TT' ] = folder+'/RUNAnalysis_TT_80X_V2p4_'+args.version+'.root'
    	allSamples[ 'ZJetsToQQ' ] = folder+'/RUNAnalysis_ZJetsToQQ_80X_V2p4_'+args.version+'.root'
    	allSamples[ 'WJetsToQQ' ] = folder+'/RUNAnalysis_WJetsToQQ_80X_V2p4_'+args.version+'.root'
	allSamples[ 'Dibosons' ] = folder+'/RUNAnalysis_Dibosons_80X_V2p4_'+args.version+'.root'
	#allSamples[ 'WWTo4Q' ] = folder+'/RUNAnalysis_WWTo4Q_80X_V2p4_'+args.version+'.root'
	#allSamples[ 'ZZTo4Q' ] = folder+'/RUNAnalysis_ZZTo4Q_80X_V2p4_'+args.version+'.root'
	#allSamples[ 'WZ' ] = folder+'/RUNAnalysis_WZ_80X_V2p4_'+args.version+'.root'
	allSamples[ 'QCD'+args.qcd+'All' ] = folder+'/RUNAnalysis_QCD'+args.qcd+'All_80X_V2p4_'+args.version+'.root'

	if 'pruned' in args.grooming: 
		cuts = [ 
				[ 'jet1Tau32', 0.57 ], [ 'jet2Tau32', 0.57 ],
				#[ 'jet1Tau32', 0 ], [ 'jet2Tau32', 0 ],
				[ 'jet1Tau21', 0.45 ], [ 'jet2Tau21', 0.45 ], 
				[ 'prunedMassAsym', 0.10 ], 
				[ 'deltaEtaDijet', 1.5 ]
				]
	else:
		cuts = [ 
				[ 'jet1Tau21', 0.45 ], [ 'jet2Tau21', 0.45 ], [ 'prunedMassAsym', 0.10 ], [ 'deltaEtaDijet', 1.5 ] 
				]
		
	if 'RPV' in args.samples: args.samples = 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)
	dictSamples = OrderedDict()
	for sam in allSamples: 
		if sam.startswith( args.samples ): dictSamples[ sam ] = allSamples[ sam ]

	allHistos = {}

	for sample in dictSamples:
		if 'Plots' in args.process: 
			if ('RPV' in args.samples) and args.unc:
				for uncType in [ args.unc+'Up', args.unc+'Down' ]: 
					p = Process( target=myPlotAnalyzer, args=( dictSamples[sample], cuts, sample, uncType ) )
			else: 
				p = Process( target=myPlotAnalyzer, args=( dictSamples[sample], cuts, sample, '' ) )
		else:
			p = Process( target=myAnalyzer, args=( dictSamples[sample], cuts, sample, '' ) )
	p.start()
	p.join()
