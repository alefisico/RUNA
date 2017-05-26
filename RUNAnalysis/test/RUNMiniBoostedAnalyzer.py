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
def myAnalyzer( fileSample, listCuts, signalName, UNC ):

	outputFileName = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_'+signalName+UNC+'_'+( '' if 'JetHT' in signalName else 'Moriond17_')+'80X_V2p4_'+args.version+'p1.root' 
	outputFile = TFile( outputFileName, 'RECREATE' )


	################################################################################################## Histos
	massBins = 500
	massXmin = 0.
	massXmax = 500.
	listOfOptions = [ [ j,k] for j in range(len(listCuts)-1) for k in range(1, len(listCuts) ) if k > j ]

	print '--- Sample ', signalName
	sf = scaleFactor(signalName)
	if 'JetHT' in signalName: signalName = 'JetHT_Run2016'
	elif 'QCD_HT' in signalName: signalName = 'QCDHTAll'
	elif 'QCD_Pt' in signalName: signalName = 'QCDPtAll'
	'''
	allHistos[ "cutFlow_"+signalName ] = TH1F( "cutflow_"+signalName, "cutflow_"+signalName, len(listCuts), 0., len(listCuts) )
	allHistos[ "cutFlow_Scaled_"+signalName ] = TH1F( "cutflow_scaled_"+signalName, "cutflow_scaled_"+signalName, len(listCuts), 0., len(listCuts) )
	allHistos[ "cutFlow_Scaled_Weights_"+signalName ] = TH1F( "cutflow_scaled_weights_"+signalName, "cutflow_scaled_weights_"+signalName, len(listCuts), 0., len(listCuts) )
	allHistos[ "HT_"+signalName ] = TH1F( "HT_"+signalName, "HT_"+signalName, 5000, 0., 5000 )
	allHistos[ "MET_"+signalName ] = TH1F( "MET_"+signalName, "MET_"+signalName, 500, 0., 500 )
	allHistos[ "massAve_"+signalName ] = TH1F( "massAve_"+signalName, "massAve_"+signalName, 500, 0., 500 )
	allHistos[ "numJets_"+signalName ] = TH1F( "numJets_"+signalName, "numJets_"+signalName, 20, 0., 20 )
	allHistos[ "jet1Pt_"+signalName ] = TH1F( "jet1Pt_"+signalName, "jet1Pt_"+signalName, 2000, 0., 2000 )
	allHistos[ "jet2Pt_"+signalName ] = TH1F( "jet2Pt_"+signalName, "jet2Pt_"+signalName, 2000, 0., 2000 )
	allHistos[ "jet1CosThetaStar_"+signalName ] = TH1F( "jet1CosThetaStar_"+signalName, "jet1CosThetaStar_"+signalName, 20, 0., 1 )
	allHistos[ "jet2CosThetaStar_"+signalName ] = TH1F( "jet2CosThetaStar_"+signalName, "jet2CosThetaStar_"+signalName, 20, 0., 1 )
	allHistos[ "jet1Tau32_"+signalName ] = TH1F( "jet1Tau32_"+signalName, "jet1Tau32_"+signalName, 20, 0., 1 )
	allHistos[ "jet2Tau32_"+signalName ] = TH1F( "jet2Tau32_"+signalName, "jet2Tau32_"+signalName, 20, 0., 1 )
	allHistos[ "jet1RhoDDT_"+signalName ] = TH1F( "jet1RhoDDT_"+signalName, "jet1RhoDDT_"+signalName, 200, -10, 10 )
	allHistos[ "jet2RhoDDT_"+signalName ] = TH1F( "jet2RhoDDT_"+signalName, "jet2RhoDDT_"+signalName, 200, -10, 10 )
	allHistos[ "jet1Tau21VsRhoDDT_"+signalName ] = TH2F( "jet1Tau21VsRhoDDT_"+signalName, "jet1Tau21VsRhoDDT_"+signalName, 20, 0., 1., 200, -10, 10 )
	allHistos[ "jet2Tau21VsRhoDDT_"+signalName ] = TH2F( "jet2Tau21VsRhoDDT_"+signalName, "jet2Tau21VsRhoDDT_"+signalName, 20, 0., 1., 200, -10, 10 )
	allHistos[ "jet1Tau21DDT_"+signalName ] = TH1F( "jet1Tau21DDT_"+signalName, "jet1Tau21DDT_"+signalName, 30, -1, 2 )
	allHistos[ "jet2Tau21DDT_"+signalName ] = TH1F( "jet2Tau21DDT_"+signalName, "jet2Tau21DDT_"+signalName, 30, -1, 2 )
	allHistos[ "jet1Tau21DDTVsRhoDDT_"+signalName ] = TH2F( "jet1Tau21DDTVsRhoDDT_"+signalName, "jet1Tau21DDTVsRhoDDT_"+signalName, 20, 0., 1., 200, -10, 10 )
	allHistos[ "jet2Tau21DDTVsRhoDDT_"+signalName ] = TH2F( "jet2Tau21DDTVsRhoDDT_"+signalName, "jet2Tau21DDTVsRhoDDT_"+signalName, 20, 0., 1., 200, -10, 10 )
	allHistos[ "jet1Tau31_"+signalName ] = TH1F( "jet1Tau31_"+signalName, "jet1Tau31_"+signalName, 20, 0., 1 )
	allHistos[ "jet2Tau31_"+signalName ] = TH1F( "jet2Tau31_"+signalName, "jet2Tau31_"+signalName, 20, 0., 1 )
	allHistos[ "jet1SubjetPtRatio_"+signalName ] = TH1F( "jet1SubjetPtRatio_"+signalName, "jet1SubjetPtRatio_"+signalName, 20, 0., 1 )
	allHistos[ "jet2SubjetPtRatio_"+signalName ] = TH1F( "jet2SubjetPtRatio_"+signalName, "jet2SubjetPtRatio_"+signalName, 20, 0., 1 )
	allHistos[ "jet1BtagCSV_"+signalName ] = TH1F( "jet1BtagCSV_"+signalName, "jet1BtagCSV_"+signalName, 5, 0., 5 )
	allHistos[ "jet2BtagCSV_"+signalName ] = TH1F( "jet2BtagCSV_"+signalName, "jet2BtagCSV_"+signalName, 5, 0., 5 )
	allHistos[ "jetsBtagCSV_"+signalName ] = TH1F( "jetsBtagCSV_"+signalName, "jetsBtagCSV_"+signalName, 5, 0., 5 )
	'''
	allHistos[ "massAve_preSel_"+signalName ] = TH1F( "massAve_preSel_"+signalName, "massAve_preSel_"+signalName, massBins, massXmin, massXmax )
	allHistos[ "massAve_preSel_lowNPV_"+signalName ] = TH1F( "massAve_preSel_lowNPV_"+signalName, "massAve_preSel_lowNPV_"+signalName, massBins, massXmin, massXmax )
	allHistos[ "NPV_preSel_lowNPV_"+signalName ] = TH1F( "NPV_preSel_lowNPV_"+signalName, "NPV_preSel_lowNPV_"+signalName, 50, 0, 50 )
	allHistos[ "massAve_preSel_medNPV_"+signalName ] = TH1F( "massAve_preSel_medNPV_"+signalName, "massAve_preSel_medNPV_"+signalName, massBins, massXmin, massXmax )
	allHistos[ "NPV_preSel_medNPV_"+signalName ] = TH1F( "NPV_preSel_medNPV_"+signalName, "NPV_preSel_medNPV_"+signalName, 50, 0, 50 )
	allHistos[ "massAve_preSel_highNPV_"+signalName ] = TH1F( "massAve_preSel_highNPV_"+signalName, "massAve_preSel_highNPV_"+signalName, massBins, massXmin, massXmax )
	allHistos[ "NPV_preSel_highNPV_"+signalName ] = TH1F( "NPV_preSel_highNPV_"+signalName, "NPV_preSel_highNPV_"+signalName, 50, 0, 50 )
	allHistos[ "HT_preSel_"+signalName ] = TH1F( "HT_preSel_"+signalName, "HT_preSel_"+signalName, 5000, 0, 5000 )
	allHistos[ "deltaEtaDijet_preSel_"+signalName ] = TH1F( "deltaEtaDijet_preSel_"+signalName, "deltaEtaDijet_preSel_"+signalName, 50, 0., 5 )
	allHistos[ "prunedMassAsym_preSel_"+signalName ] = TH1F( "prunedMassAsym_preSel_"+signalName, "prunedMassAsym_preSel_"+signalName, 20, 0., 1 )
	allHistos[ "jet1Tau21_preSel_"+signalName ] = TH1F( "jet1Tau21_preSel_"+signalName, "jet1Tau21_preSel_"+signalName, 20, 0., 1 )
	allHistos[ "jet2Tau21_preSel_"+signalName ] = TH1F( "jet2Tau21_preSel_"+signalName, "jet2Tau21_preSel_"+signalName, 20, 0., 1 )
	allHistos[ "jet1Tau31_preSel_"+signalName ] = TH1F( "jet1Tau31_preSel_"+signalName, "jet1Tau31_preSel_"+signalName, 20, 0., 1 )
	allHistos[ "jet2Tau31_preSel_"+signalName ] = TH1F( "jet2Tau31_preSel_"+signalName, "jet2Tau31_preSel_"+signalName, 20, 0., 1 )

	allHistos[ "deltaEtaDijet_n-1_"+signalName ] = TH1F( "deltaEtaDijet_n-1_"+signalName, "deltaEtaDijet_n-1_"+signalName, 50, 0., 5 )
	allHistos[ "prunedMassAsym_n-1_"+signalName ] = TH1F( "prunedMassAsym_n-1_"+signalName, "prunedMassAsym_n-1_"+signalName, 20, 0., 1 )
	allHistos[ "jet1Tau21_n-1_"+signalName ] = TH1F( "jet1Tau21_n-1_"+signalName, "jet1Tau21_n-1_"+signalName, 20, 0., 1 )
	allHistos[ "jet2Tau21_n-1_"+signalName ] = TH1F( "jet2Tau21_n-1_"+signalName, "jet2Tau21_n-1_"+signalName, 20, 0., 1 )
	allHistos[ "jet1Tau31_n-1_"+signalName ] = TH1F( "jet1Tau31_n-1_"+signalName, "jet1Tau31_n-1_"+signalName, 20, 0., 1 )
	allHistos[ "jet2Tau31_n-1_"+signalName ] = TH1F( "jet2Tau31_n-1_"+signalName, "jet2Tau31_n-1_"+signalName, 20, 0., 1 )
	listCuts.append( [ 'btag' ] )
	for var in listCuts:
		if 'deltaEta' in var[0]: 
			allHistos[ var[0]+'_'+signalName ] = TH1F( var[0]+'_'+signalName, var[0]+'_'+signalName, 50, 0., 5. )
			for var1 in listCuts: allHistos[ var[0]+'_'+var1[0]+"_"+signalName ] = TH1F( var[0]+'_'+var1[0]+"_"+signalName, var[0]+'_'+var1[0]+"_"+signalName, 50, 0., 5. )
		else: 
			allHistos[ var[0]+'_'+signalName ] = TH1F( var[0]+'_'+signalName, var[0]+'_'+signalName, 20, 0., 1. )
			for var1 in listCuts: allHistos[ var[0]+'_'+var1[0]+"_"+signalName ] = TH1F( var[0]+'_'+var1[0]+"_"+signalName, var[0]+'_'+var1[0]+"_"+signalName, 20, 0., 1. )
		allHistos[ "massAve_"+var[0]+'_'+signalName ] = TH1F( "massAve_"+var[0]+'_'+signalName, "massAve_"+var[0]+'_'+signalName, massBins, massXmin, massXmax )

		allHistos[ "HT_"+var[0]+"_"+signalName ] = TH1F( "HT_"+var[0]+"_"+signalName, "HT_"+var[0]+"_"+signalName, 5000, 0., 5000 )
		allHistos[ "MET_"+var[0]+"_"+signalName ] = TH1F( "MET_"+var[0]+"_"+signalName, "MET_"+var[0]+"_"+signalName, 500, 0., 500 )
		allHistos[ "numJets_"+var[0]+"_"+signalName ] = TH1F( "numJets_"+var[0]+"_"+signalName, "numJets_"+var[0]+"_"+signalName, 20, 0., 20 )
		allHistos[ "jet1Pt_"+var[0]+"_"+signalName ] = TH1F( "jet1Pt_"+var[0]+"_"+signalName, "jet1Pt_"+var[0]+"_"+signalName, 2000, 0., 2000 )
		allHistos[ "jet2Pt_"+var[0]+"_"+signalName ] = TH1F( "jet2Pt_"+var[0]+"_"+signalName, "jet2Pt_"+var[0]+"_"+signalName, 2000, 0., 2000 )
	listCuts.remove( ['btag'] )
	allHistos[ "massAve_deltaEtaDijet_lowNPV_"+signalName ] = TH1F( "massAve_deltaEtaDijet_lowNPV_"+signalName, "massAve_deltaEtaDijet_lowNPV_"+signalName, massBins, massXmin, massXmax )
	allHistos[ "NPV_deltaEtaDijet_lowNPV_"+signalName ] = TH1F( "NPV_deltaEtaDijet_lowNPV_"+signalName, "NPV_deltaEtaDijet_lowNPV_"+signalName, 50, 0, 50 )
	allHistos[ "massAve_deltaEtaDijet_medNPV_"+signalName ] = TH1F( "massAve_deltaEtaDijet_medNPV_"+signalName, "massAve_deltaEtaDijet_medNPV_"+signalName, massBins, massXmin, massXmax )
	allHistos[ "NPV_deltaEtaDijet_medNPV_"+signalName ] = TH1F( "NPV_deltaEtaDijet_medNPV_"+signalName, "NPV_deltaEtaDijet_medNPV_"+signalName, 50, 0, 50 )
	allHistos[ "massAve_deltaEtaDijet_highNPV_"+signalName ] = TH1F( "massAve_deltaEtaDijet_highNPV_"+signalName, "massAve_deltaEtaDijet_highNPV_"+signalName, massBins, massXmin, massXmax )
	allHistos[ "NPV_deltaEtaDijet_highNPV_"+signalName ] = TH1F( "NPV_deltaEtaDijet_highNPV_"+signalName, "NPV_deltaEtaDijet_highNPV_"+signalName, 50, 0, 50 )

	for ind in listOfOptions:
		tmpName = listCuts[ind[0]][0]+'Vs'+listCuts[ind[1]][0]+'_'+signalName
		allHistos[ tmpName ] = TH2F( tmpName, tmpName, 
				(50 if 'deltaEta' in listCuts[ind[0]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[0]][0] else 1. ),
				(50 if 'deltaEta' in listCuts[ind[1]][0] else 20 ), 0., (5. if 'deltaEta' in listCuts[ind[1]][0] else 1. ) 
				)

	tmpNameSam = listCuts[-2][0]+'Vs'+listCuts[-1][0]+'_'+signalName
	#tmpNameSam = tmpName #listCuts[-2][0]+'Vs'+listCuts[-1][0]+'_'+signalName
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

	################################################################################################## Running the Analysis
	print '-'*40
	#lumiWeight = ( 1 if 'JetHT' in signalName else sf  )
	#SF = 'puWeight*'+str(lumiWeight)
	SF = TCut('lumiWeight * puWeight') # * '+str(args.lumi))
	preselection = TCut('HT>900') * TCut( SF )
	stringSel = '' 
	for var in listCuts: stringSel = stringSel+'('+var[0]+'<'+str(var[1])+')'
	stringSel = stringSel.replace(')(',') * (')
	print '---- Cuts applied: ', stringSel, preselection * TCut( stringSel )

	ABCDRegions = {}
	ABCDRegions[ '_A' ] = TCut( stringSel )
	ABCDRegions[ '_B' ] = TCut( stringSel.replace('deltaEtaDijet<', 'deltaEtaDijet>') )
	ABCDRegions[ '_C' ] = TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>') )
	ABCDRegions[ '_D' ] = TCut( stringSel.replace('prunedMassAsym<', 'prunedMassAsym>').replace('deltaEtaDijet<', 'deltaEtaDijet>') )

	btagSelection = stringSel + '&& (jet1btagCSVv2 > 0.8) && (jet2btagCSVv2 > 0.8)'
	ABCDRegionsBtag = {}
	ABCDRegionsBtag[ '_A' ] = TCut( btagSelection )
	ABCDRegionsBtag[ '_B' ] = TCut( btagSelection.replace('deltaEtaDijet<', 'deltaEtaDijet>') )
	ABCDRegionsBtag[ '_C' ] = TCut( btagSelection.replace('prunedMassAsym<', 'prunedMassAsym>') )
	ABCDRegionsBtag[ '_D' ] = TCut( btagSelection.replace('prunedMassAsym<', 'prunedMassAsym>').replace('deltaEtaDijet<', 'deltaEtaDijet>') )

	sel = preselection * TCut( stringSel )
	btagSel = sel * TCut( btagSelection )

	treeName = 'BoostedAnalysisPlots'+('Puppi' if 'Puppi' in args.grooming else '')+'/RUNATree'

	if 'JetHT' in signalName: SF = 1

	### Preselection
	getHistoFromTree( fileSample, treeName,
			'HT', 
			preselection,
			allHistos[ 'HT_preSel_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 
	
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			preselection, 
			allHistos[ 'massAve_preSel_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			preselection+TCut( '(numPV<12)' ), 
			allHistos[ 'massAve_preSel_lowNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'numPV', 
			preselection+TCut( '(numPV<12)' ), 
			allHistos[ 'NPV_preSel_lowNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			preselection+TCut( '(numPV>12) && (numPV<24)' ), 
			allHistos[ 'massAve_preSel_medNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'numPV', 
			preselection+TCut( '(numPV>12) && (numPV<24)' ), 
			allHistos[ 'NPV_preSel_medNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			preselection+TCut( '(numPV>24)' ), 
			allHistos[ 'massAve_preSel_highNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'numPV', 
			preselection+TCut( '(numPV>24)' ), 
			allHistos[ 'NPV_preSel_highNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAsym', 
			preselection,
			allHistos[ 'prunedMassAsym_preSel_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'deltaEtaDijet', 
			preselection,
			allHistos[ 'deltaEtaDijet_preSel_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 
	getHistoFromTree( fileSample, treeName,
			'jet1Tau21', 
			preselection,
			allHistos[ 'jet1Tau21_preSel_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 
	getHistoFromTree( fileSample, treeName,
			'jet2Tau21', 
			preselection,
			allHistos[ 'jet2Tau21_preSel_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 



	### All selection
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			sel, 
			allHistos[ 'massAve_deltaEtaDijet_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			sel+TCut('(numPV<12)'), 
			allHistos[ 'massAve_deltaEtaDijet_lowNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'numPV', 
			sel+TCut('(numPV<12)'), 
			allHistos[ 'NPV_deltaEtaDijet_lowNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			sel+TCut('(numPV>12) && (numPV<24)'), 
			allHistos[ 'massAve_deltaEtaDijet_medNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'numPV', 
			sel+TCut('(numPV>12) && (numPV<24)'), 
			allHistos[ 'NPV_deltaEtaDijet_medNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			sel+TCut('(numPV>24)'), 
			allHistos[ 'massAve_deltaEtaDijet_highNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'numPV', 
			sel+TCut('(numPV>24)'), 
			allHistos[ 'NPV_deltaEtaDijet_highNPV_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 


	### n-1 selection
	getHistoFromTree( fileSample, treeName,
			'prunedMassAsym', 
			preselection + TCut( stringSel.replace('&& (prunedMassAsym<0.1)','') ), 
			allHistos[ 'prunedMassAsym_n-1_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	getHistoFromTree( fileSample, treeName,
			'deltaEtaDijet', 
			preselection + TCut( stringSel.replace(('&& (deltaEtaDijet<1.5)' if 'pruned' in args.grooming else '&& (deltaEtaDijet<1.0)' ),'') ), 
			allHistos[ 'deltaEtaDijet_n-1_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 
	getHistoFromTree( fileSample, treeName,
			'jet1Tau21', 
			preselection + TCut( stringSel.replace('&& (jet1Tau21<0.45)','') ), 
			allHistos[ 'jet1Tau21_n-1_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 
	getHistoFromTree( fileSample, treeName,
			'jet2Tau21', 
			preselection + TCut( stringSel.replace('&& (jet2Tau21<0.45)','') ), 
			allHistos[ 'jet2Tau21_n-1_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	### ABCD plots
	for region, selABCD in ABCDRegions.items():
		getHistoFromTree( fileSample, treeName,
				'prunedMassAve', 
				selABCD, 
				allHistos[ 'massAve_prunedMassAsymVsdeltaEtaDijet_'+signalName+region ], 
				( 0.05 if 'JetHT' in signalName else 1 ) ) 
	
		get2DHistoFromTree( fileSample, treeName,
				'prunedMassAsym', 'deltaEtaDijet',
				selABCD, 
				allHistos[ 'prunedMassAsymVsdeltaEtaDijet_'+signalName+region ],
				( 0.05 if 'JetHT' in signalName else 1 ) ) 

	## Btagging
	getHistoFromTree( fileSample, treeName,
			'prunedMassAve', 
			btagSel, 
			allHistos[ 'massAve_btag_'+signalName ], 
			( 0.05 if 'JetHT' in signalName else 1 ) ) 

	for region, selABCD in ABCDRegionsBtag.items():
		getHistoFromTree( fileSample, treeName,
				'prunedMassAve', 
				selABCD, 
				allHistos[ 'massAve_prunedMassAsymVsdeltaEtaDijet_'+signalName+'_btag'+region ], 
				( 0.05 if 'JetHT' in signalName else 1 ) ) 
	
		get2DHistoFromTree( fileSample, treeName,
				'prunedMassAsym', 'deltaEtaDijet',
				selABCD, 
				allHistos[ 'prunedMassAsymVsdeltaEtaDijet_'+signalName+'_btag'+region ],
				( 0.05 if 'JetHT' in signalName else 1 ) ) 


#	for sample in dictSamples:
#
#		####### Get GenTree 
#		inputFile, events, numEntries = getTree( dictSamples[ sample ], ('BoostedAnalysisPlotsPuppi'+( '' if 'PDF' in UNC else UNC)+'/RUNATree' if 'Puppi' in args.grooming else 'BoostedAnalysisPlots'+( '' if 'PDF' in UNC else UNC)+'/RUNATree' ) )
#		print '-'*40
#		print '------> ', sample
#		print '------> Number of events: '+str(numEntries)
#		d = 0
#		cutFlowList = OrderedDict()
#		cutFlowScaledList = OrderedDict()
#		cutFlowScaledListWeights = OrderedDict()
#		cutFlowList[ 'Process' ] = 0
#		cutFlowList[ 'Preselection' ] = 0
#		cutFlowScaledList[ 'Process' ] = 0
#		cutFlowScaledList[ 'Preselection' ] = 0
#		cutFlowScaledListWeights[ 'Process' ] = 0
#		cutFlowScaledListWeights[ 'Preselection' ] = 0
#		for k in listCuts: 
#			cutFlowList[ k[0] ] = 0
#			cutFlowScaledList[ k[0] ] = 0
#			cutFlowScaledListWeights[ k[0] ] = 0
#		cutFlowList[ 'btag' ] = 0
#		cutFlowScaledList[ 'btag' ] = 0
#		cutFlowScaledListWeights[ 'btag' ] = 0
#
#		for i in xrange(numEntries):
#			events.GetEntry(i)
#
#			#---- progress of the reading --------
#			fraction = 10.*i/(1.*numEntries)
#			if TMath.FloorNint(fraction) > d: print str(10*TMath.FloorNint(fraction))+'%' 
#			d = TMath.FloorNint(fraction)
#			#if ( i > 100000 ): break
#
#			Run      = events.run
#			Lumi     = events.lumi
#			NumEvent = events.event
#			puWeight	= events.puWeight
#			if 'v06' in args.version: pdfWeight	= events.pdfWeight
#			lumiWeight	= events.lumiWeight
#			HT		= events.HT
#			MET		= events.MET
#			numJets		= events.numJets
#			massAve		= getattr( events, (args.grooming+"MassAve").replace('Puppi','') )
#			jet1Mass	= getattr( events, 'jet1'+(args.grooming+"Mass").replace('pruned','Pruned').replace('soft','Soft').replace('Puppi',''))
#			jet2Mass	= getattr( events, 'jet2'+(args.grooming+"Mass").replace('pruned','Pruned').replace('soft','Soft').replace('Puppi',''))
#			jet1Pt          = events.jet1Pt
#			jet2Pt          = events.jet2Pt
#			jet1Eta          = events.jet1Eta
#			jet2Eta          = events.jet2Eta
#			jet1CosThetaStar	= events.jet1CosThetaStar
#			jet2CosThetaStar	= events.jet2CosThetaStar
#			jet1BtagCSV		= ( events.jet1btagCSVv2 > 0.800 )
#			jet2BtagCSV		= ( events.jet2btagCSVv2 > 0.800 )
#			
#			#print 'Entry ', Run, ':', Lumi, ':', NumEvent
#
#			if 'DATA' in sample: scale = 1
#			#elif 'RPV' in sample: scale = 2606 * puWeight * SF
#			else: scale = 2666 * puWeight * lumiWeight
#
#			if 'PDF' in UNC:
#				if 'Up' in UNC: scale = scale*(1+pdfWeight)
#				else: scale = scale*(1-pdfWeight)
#
#			cutFlowList[ 'Process' ] += 1
#			cutFlowScaledList[ 'Process' ] += scale
#			cutFlowScaledList[ 'Process' ] += (puWeight*puWeight)
#
#			########## DDT
#			jet1RhoDDT = TMath.Log( jet1Mass*jet1Mass/jet1Pt )
#			jet2RhoDDT = TMath.Log( jet2Mass*jet2Mass/jet2Pt )
#			jet1Tau21DDT = events.jet1Tau21 + 0.063 * jet1RhoDDT 
#			jet2Tau21DDT = events.jet2Tau21 + 0.063 * jet2RhoDDT 
#			
#			#### Pre-selection
#			HTCut = ( HT > 900 )
#			dijetCut =  ( numJets > 1 )
#			#jetPtCut =  ( jet1Pt > 500 ) and ( jet2Pt > 450 )
#			jetPtCut =  ( jet1Pt > 150 ) and ( jet2Pt > 150 )
#			
#			#if HTCut and dijetCut and jetPtCut:
#			if HTCut and dijetCut :
#				cutFlowList[ 'Preselection' ] += 1
#				cutFlowScaledList[ 'Preselection' ] += scale
#				cutFlowScaledList[ 'Preselection' ] += (puWeight*puWeight)
#				sigCutsList = []
#				allHistos[ "HT_"+sam ].Fill( HT, scale )
#				allHistos[ "MET_"+sam ].Fill( MET, scale )
#				allHistos[ "massAve_"+sam ].Fill( massAve, scale )
#				allHistos[ "numJets_"+sam ].Fill( numJets, scale )
#				allHistos[ "jet1Pt_"+sam ].Fill( jet1Pt, scale )
#				allHistos[ "jet2Pt_"+sam ].Fill( jet2Pt, scale )
#				allHistos[ "jet1RhoDDT_"+sam ].Fill( jet1RhoDDT, scale )
#				allHistos[ "jet2RhoDDT_"+sam ].Fill( jet2RhoDDT, scale )
#				allHistos[ "jet1Tau21VsRhoDDT_"+sam ].Fill( events.jet1Tau21, jet1RhoDDT, scale )
#				allHistos[ "jet2Tau21VsRhoDDT_"+sam ].Fill( events.jet2Tau21, jet2RhoDDT, scale )
#				allHistos[ "jet1Tau21DDT_"+sam ].Fill( jet1Tau21DDT, scale )
#				allHistos[ "jet2Tau21DDT_"+sam ].Fill( jet2Tau21DDT, scale )
#				allHistos[ "jet1Tau21DDTVsRhoDDT_"+sam ].Fill( jet1Tau21DDT, jet1RhoDDT, scale )
#				allHistos[ "jet2Tau21DDTVsRhoDDT_"+sam ].Fill( jet2Tau21DDT, jet2RhoDDT, scale )
#				allHistos[ "prunedMassAsym_"+sam ].Fill( events.prunedMassAsym, scale )
#				allHistos[ "deltaEtaDijet_"+sam ].Fill( events.deltaEtaDijet, scale )
#				allHistos[ "jet1CosThetaStar_"+sam ].Fill( jet1CosThetaStar, scale )
#				allHistos[ "jet2CosThetaStar_"+sam ].Fill( jet2CosThetaStar, scale )
#				allHistos[ "jet1Tau21_"+sam ].Fill( events.jet1Tau21, scale )
#				allHistos[ "jet2Tau21_"+sam ].Fill( events.jet2Tau21, scale )
#				allHistos[ "jet1Tau31_"+sam ].Fill( events.jet1Tau31, scale )
#				allHistos[ "jet2Tau31_"+sam ].Fill( events.jet2Tau31, scale )
#				allHistos[ "jet1Tau32_"+sam ].Fill( events.jet1Tau32, scale )
#				allHistos[ "jet2Tau32_"+sam ].Fill( events.jet2Tau32, scale )
#				allHistos[ "jet1SubjetPtRatio_"+sam ].Fill( events.jet1SubjetPtRatio, scale )
#				allHistos[ "jet2SubjetPtRatio_"+sam ].Fill( events.jet2SubjetPtRatio, scale )
#				allHistos[ "jet1BtagCSV_"+sam ].Fill( 1 if jet1BtagCSV else 0 )
#				allHistos[ "jet2BtagCSV_"+sam ].Fill( 1 if jet1BtagCSV else 0 )
#
#				bothBtag = ( jet1BtagCSV and jet2BtagCSV )
#				oneBtag = ( jet1BtagCSV or jet2BtagCSV )
#				if bothBtag: allHistos[ "jetsBtagCSV_"+sam ].Fill( 2 )
#				elif oneBtag: allHistos[ "jetsBtagCSV_"+sam ].Fill( 1 )
#				else: allHistos[ "jetsBtagCSV_"+sam ].Fill( 0 )
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
#						allHistos[ "HT_"+var[0]+"_"+sam ].Fill( HT, scale )
#						allHistos[ "MET_"+var[0]+"_"+sam ].Fill( MET, scale )
#						allHistos[ "numJets_"+var[0]+"_"+sam ].Fill( numJets, scale )
#						allHistos[ "jet1Pt_"+var[0]+"_"+sam ].Fill( jet1Pt, scale )
#						allHistos[ "jet2Pt_"+var[0]+"_"+sam ].Fill( jet2Pt, scale )
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
#					allHistos[ "HT_btag_"+sam ].Fill( HT, scale )
#					allHistos[ "MET_btag_"+sam ].Fill( MET, scale )
#					allHistos[ "numJets_btag_"+sam ].Fill( numJets, scale )
#					allHistos[ "jet1Pt_btag_"+sam ].Fill( jet1Pt, scale )
#					allHistos[ "jet2Pt_btag_"+sam ].Fill( jet2Pt, scale )
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
	parser.add_argument( '-p', '--process', action='store',  dest='process', default='single', help='Process: all or single.' )
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

	if 'pruned' in args.grooming: cuts = selection['CHSpruned'] 
	else: cuts = selection['PUPPIsoftDrop']
		
	if 'RPV' in args.samples: args.samples = 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)
	dictSamples = OrderedDict()
	for sam in allSamples: 
		if sam.startswith( args.samples ): dictSamples[ sam ] = allSamples[ sam ]

	allHistos = {}

	for sample in dictSamples:
		if ('RPV' in args.samples) and args.unc:
			for uncType in [ args.unc+'Up', args.unc+'Down' ]: 
				p = Process( target=myAnalyzer, args=( dictSamples[sample], cuts, sample, uncType ) )
		else: 
			p = Process( target=myAnalyzer, args=( dictSamples[sample], cuts, sample, '' ) )
		p.start()
		p.join()
