#!/usr/bin/env python

'''
File: MyAnalyzer.py --mass 50 (optional) --debug -- final --jetAlgo AK5
Author: Alejandro Gomez Espinosa
Email: gomez@physics.rutgers.edu
Description: My Analyzer 
'''

import sys,os,time, re
from math import *
from string import *
from array import array
import argparse
from ROOT import * 
from collections import OrderedDict
from multiprocessing import Process
import numpy as np
try: 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
	from RUNA.RUNAnalysis.commonFunctions import *
except ImportError: 
	sys.path.append('../python') 
	import tdrstyle as tdrstyle
	from commonFunctions import *


TMVA.Tools.Instance()
gROOT.Reset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)


#----------------------------------------------------------------------
### Main Optimization
def calcROCs( BkgSamples, SigSamples, treename, varList, window, cutsList ):
	"""docstring for calcROCs"""
	
	outputFile = TFile( 'test.root', "RECREATE" )
	
	allHistos = {}
	allSumw2 = {}
	for var in varList: 
		allHistos[ var[0]+"_Sig" ] = TH1F( var[0]+"_Sig", var[0]+"_Sig", var[1], var[2], var[3] )
		allHistos[ var[0]+"_Sig" ].Sumw2()
		allSumw2[ var[0]+"_Sig" ] = 0
		for bkgSample in BkgSamples: 
			allHistos[ var[0]+'_'+bkgSample ] = TH1F( var[0]+"_"+bkgSample, var[0]+"_"+bkgSample, var[1], var[2], var[3] )
			allHistos[ var[0]+'_'+bkgSample ].Sumw2()
			allSumw2[ var[0]+"_"+bkgSample  ] = 0
			allHistos[ var[0]+"_"+bkgSample+"_BkgROC"] = TH1F( var[0]+"_"+bkgSample+"_BkgROC", var[0]+"_"+bkgSample+"_BkgROC; "+var[0], var[1], var[2], var[3] )
			allHistos[ var[0]+"_"+bkgSample+"_SigROC"] = TH1F( var[0]+"_"+bkgSample+"_SigROC", var[0]+"_"+bkgSample+"_SigROC; "+var[0], var[1], var[2], var[3] )

	massAve = 'prunedMassAve' if 'Boosted' in args.boosted else 'massAve'
	SF = TCut( 'lumiWeight * puWeight * '+str(args.lumi) )
	if 'Resolved' in args.boosted: preselection = TCut( '(jetsPt[0]>80) && (jetsPt[1]>80) && (jetsPt[2]>80) && (jetsPt[3]>80) && (HT>900)')
	else: preselection = TCut('')
	cuts = SF * preselection * TCut('('+massAve+'>'+str(args.mass-window)+') && ('+massAve+'<'+str(args.mass+window)+') ' )
	if len(cutsList) > 0:
		#for cutVar in cutsList: cuts += ' && '+cutVar[0]+('<' if cutVar[4] else '>' )+str(cutVar[5])
		for cutVar in cutsList: cuts = cuts * TCut( cutVar[0]+('<' if cutVar[4] else '>' )+str(cutVar[5]) )
	print '-'*40
	print '---- Cuts applied: ', cuts

	#### Signal
	print '-'*40
	print '---- RPV mass', args.mass
	for sigVar in varList: 
		allHistos[ sigVar[0]+"_Sig" ] = getHistoFromTree( SigSamples[ args.mass ], treename, sigVar[0], cuts, allHistos[ sigVar[0]+"_Sig" ], 1 ) 
		print '---- processing... ', sigVar[0]


		#allSumw2[ sigVar[0]+"_Sig" ] += sigSF*sigSF


	allROCs = {}
	for bkgSample in BkgSamples: 
		print '-'*40
		print '---- ', bkgSample
		for bkgVar in varList: 
			print '---- processing... ', bkgVar[0]
			allHistos[ bkgVar[0]+"_"+bkgSample ] = getHistoFromTree( BkgSamples[ bkgSample ][0], treename, bkgVar[0], cuts, allHistos[ bkgVar[0]+"_"+bkgSample ], 1 ) 

		for var in varList:
			BkgHisto = allHistos[ var[0]+"_"+bkgSample ]
			BkgTotal = BkgHisto.Integral() 

			SigHisto = allHistos[ var[0]+"_Sig" ]
			SigTotal = SigHisto.Integral() 

			ibinCut = 0
			#if ( len( var ) > 2 ): ibinCut = SigHisto.GetXaxis().FindBin( var[2] )

			BkgROCValues = []	
			BkgNumValues = []	
			BkgErrValues = []	
			SigROCValues = []	
			SigNumValues = []	
			SigErrValues = []	
			SigROCBins = []	
			SigROCLowEdge = []	
			firstBin =  SigHisto.FindFirstBinAbove( 0, 1 )
			lastBin =  SigHisto.FindLastBinAbove( 0, 1 )+1 
			#print 'Bfore', allSumw2[ bkgVar[0]+"_"+bkgSample ], TMath.Sqrt( allSumw2[ bkgVar[0]+"_"+bkgSample ] )
			#errorVal = Double(0)
			#ingral = BkgHisto.IntegralAndError( 0, 1, errorVal )
			#print '0'*10, ingral, errorVal
			for ibin in range( firstBin, lastBin ):
				if var[4]:
					#print ibin, SigHisto.GetXaxis().GetBinLowEdge(ibin), SigHisto.Integral( 0, ibin ), SigHisto.Integral( firstBin, ibin ) #, SigTotal, SigHisto.Integral( 0, ibin+1 ) / SigTotal,  BkgHisto.GetXaxis().GetBinLowEdge(ibin), BkgHisto.Integral( 0, ibin+1 ) , BkgTotal, BkgHisto.Integral( 0, ibin+1 ) / BkgTotal
					try: effBkg = 1 - BkgHisto.Integral( firstBin, ibin ) / BkgTotal 
					except ZeroDivisionError: effBkg = 1
					allHistos[ var[0]+"_"+bkgSample+"_BkgROC" ].SetBinContent( ibin, effBkg )
					BkgROCValues.append( effBkg )
					bkgErrorInt = Double(0)
					bkgIntegral = BkgHisto.IntegralAndError( firstBin, ibin, bkgErrorInt )
					BkgNumValues.append( bkgIntegral )
					BkgErrValues.append( TMath.Sqrt( bkgErrorInt ) )
					#print '1'*10, ingral, errorVal
					#print BkgHisto.GetSumOfWeights()
					#for i in BkgHisto.GetSumw2(): print i
					#print allSumw2[ var[0]+"_"+bkgSample ], TMath.Sqrt( allSumw2[ var[0]+"_"+bkgSample ] )

					try: effSig = SigHisto.Integral( firstBin, ibin ) / SigTotal 
					except ZeroDivisionError: effSig = 0
					allHistos[ var[0]+"_"+bkgSample+"_SigROC" ].SetBinContent( ibin, effSig )
					SigROCValues.append( effSig )
					SigROCBins.append( ibin )
					SigROCLowEdge.append( SigHisto.GetXaxis().GetBinLowEdge(ibin+1) )
					#SigNumValues.append( SigHisto.Integral( firstBin, ibin ) )
					sigErrorInt = Double(0)
					sigIntegral = SigHisto.IntegralAndError( firstBin, ibin, sigErrorInt )
					SigNumValues.append( sigIntegral )
					SigErrValues.append( TMath.Sqrt( sigErrorInt ) )

				else:
					effBkg = 1 - BkgHisto.Integral( ibin, lastBin ) / BkgTotal 
					allHistos[ var[0]+"_"+bkgSample+"_BkgROC" ].SetBinContent( ibin, effBkg )
					BkgROCValues.append( effBkg )
					#BkgNumValues.append( BkgHisto.Integral( ibin, lastBin ) )
					bkgErrorInt = Double(0)
					bkgIntegral = BkgHisto.IntegralAndError( ibin, lastBin, bkgErrorInt )
					BkgNumValues.append( bkgIntegral )
					BkgErrValues.append( TMath.Sqrt( bkgErrorInt ) )

					effSig = SigHisto.Integral( ibin, lastBin ) / SigTotal 
					allHistos[ var[0]+"_"+bkgSample+"_SigROC" ].SetBinContent( ibin, effSig )
					SigROCValues.append( effSig )
					SigROCBins.append( ibin )
					SigROCLowEdge.append( SigHisto.GetXaxis().GetBinLowEdge(ibin+1) )
					#SigNumValues.append( SigHisto.Integral( ibin, lastBin ) )
					sigErrorInt = Double(0)
					sigIntegral = SigHisto.IntegralAndError( ibin, lastBin, sigErrorInt )
					SigNumValues.append( sigIntegral )
					SigErrValues.append( TMath.Sqrt( sigErrorInt) )

			allROCs[ var[0]+"_"+bkgSample+"_ROC" ] = [ BkgROCValues, SigROCValues, SigROCBins, SigROCLowEdge, BkgNumValues, SigNumValues, BkgErrValues, SigErrValues ]

	outputTextFile = 'ROCfiles/ROC'+args.boosted+'Values_QCD'+args.qcd+'_'+args.decay+str(args.mass)+'_cut'+str(len(cutsList))+'_'+args.version+args.grooming+'.txt'
	print '--- Creating ', outputTextFile 
	print >> open(outputTextFile, 'w+'), allROCs
	outputFile.Write()
	outputFile.Close()

def makeROCs( textFile, variables, bkgSamples, perVariable, cutsList, printVar, returnSOB='' ):
	"""docstring for makeROCs"""

	f = open( textFile, 'r' )
	tmpDictVar = eval( f.read() )

	dictVar = OrderedDict()
	for q in tmpDictVar:
		for tmpVar in variables:
			if tmpVar[0] in q: 
				dictVar[ q ] = tmpDictVar[ q ]

	dictROC = OrderedDict()
	dictNumEvents = {}
	for varROC in dictVar:
		BkgROCArray = np.asarray( dictVar[ varROC ][0] )
		SigROCArray = np.asarray( dictVar[ varROC ][1] )
		SigROCBinsArray = np.asarray( dictVar[ varROC ][2] )
		SigROCLowEdgeArray = np.asarray( dictVar[ varROC ][3] )
		BkgNumArray = np.asarray( dictVar[ varROC ][4] )
		SigNumArray = np.asarray( dictVar[ varROC ][5] )
		BkgErrArray = np.asarray( dictVar[ varROC ][6] )
		SigErrArray = np.asarray( dictVar[ varROC ][7] )
		ROC = TGraph( len( BkgROCArray ), SigROCArray, BkgROCArray )		
		dictROC[ varROC ] = ROC
		dictNumEvents[ varROC ] = [ BkgNumArray, SigNumArray, SigROCLowEdgeArray, BkgErrArray, SigErrArray ]
		if (args.printValue and ( printVar in varROC )):
			'''
			if var[0] in varROC:
				x = find_nearest(SigROCLowEdgeArray, var[6] )
				print 'For cut', var[6], 'in', var[0], 'effSig: ', SigROCArray[x], ', effBkg: ', BkgROCArray[x]
			'''
			for Cut in cutsList:
				if printVar in Cut[0]:
					y = find_nearest( SigROCArray, Cut[6] )
					print 'Optimal value for ', varROC, SigROCLowEdgeArray[y], 'is', SigROCArray[y], 'bkgRej = ', BkgROCArray[y]

	if not args.printValue:
		if perVariable:
			dictNumVar = OrderedDict()
			for var in variables:
				dictVariables = { rocs: dictROC[ rocs ] for rocs in dictROC if var[0] in rocs }
				#plotROC( var[0], var[0], dictVariables, len(cutsList), True, False )
				dictVariablesNum = { rocs: dictNumEvents[ rocs ] for rocs in dictNumEvents if var[0] in rocs }

				totalSig = dictVariablesNum[ var[0]+'_QCD'+args.qcd+'All_ROC' ][1] 
				#totalSigErr = np.sqrt( np.multiply( totalSig, dictVariablesNum[ var[0]+'_QCD'+args.qcd+'All_ROC' ][4] ) ) 
				totalSigErr = np.power( dictVariablesNum[ var[0]+'_QCD'+args.qcd+'All_ROC' ][4], 2 ) ### it is power because I sqrt the error

				totalBkg = np.zeros( len( dictVariablesNum[ var[0]+'_QCD'+args.qcd+'All_ROC' ][1] ) )
				totalBkgErr = np.zeros( len( dictVariablesNum[ var[0]+'_QCD'+args.qcd+'All_ROC' ][1] ) )
				for bkg in dictVariablesNum: 
					totalBkg += dictVariablesNum[ bkg ][0] 
					totalBkgErr += np.power( dictVariablesNum[ bkg ][3], 4 )   ### it is power because I sqrt the error
				totalBkgErr = np.sqrt( totalBkgErr ) 

				sqrtSigBkg = np.sqrt( np.add(totalSig, totalBkg ) )
				#sigOverSqrtSigBkg = np.divide( totalSig, sqrtSigBkg )

				sigBkgErr = np.sqrt( np.add( np.power( totalSigErr, 2 ), np.power( totalBkgErr, 2 ) ) )
				sqrtSigBkgErr = np.divide( sigBkgErr, 2*sqrtSigBkg )
				totalErrSigOverSqrtSigBkg = np.sqrt( np.add( np.power( np.divide( totalSigErr, totalSig ), 2), np.power( np.divide( sqrtSigBkgErr, sqrtSigBkg ), 2) ) )
				#print args.mass, sigOverSqrtSigBkg
				#print totalErrSigOverSqrtSigBkg
				#sigOverSqrtSigBkg = np.divide( totalSig, np.sqrt(totalBkg) ) 
				with np.errstate(divide='ignore', invalid='ignore'):
					sigOverSqrtSigBkg = np.true_divide( totalSig, np.sqrt(totalBkg) ) 
					sigOverSqrtSigBkg[ sigOverSqrtSigBkg == np.inf ] = 0
					sigOverSqrtSigBkg = np.nan_to_num( sigOverSqrtSigBkg )
				totalErrSigOverSqrtBkg =  np.add( np.divide( totalBkgErr, np.sqrt( totalBkg ) ), -np.divide( (totalSig*totalSigErr), (2*np.power(totalBkg, 1.5)) ) )
				#sigOverSqrtSigBkg = np.sqrt( 2*( (totalSig+totalBkg)*np.log( 1 + np.divide( totalSig, totalBkg ) ) - totalSig ) )

				if var[0] == 'deltaEta': dummyBin = 5
				elif 'delta' in var[0]: dummyBin = 300.
				else: dummyBin = 1
				binsHisto = (dictVariablesNum[ var[0]+'_QCD'+args.qcd+'All_ROC' ][2])/dummyBin
				'''
				SOB = TGraphErrors( len(totalSig), 
						binsHisto, 
						sigOverSqrtSigBkg, 
						array( 'd', [0]*len(sigOverSqrtSigBkg)), 
						#totalErrSigOverSqrtBkg ) 
						totalErrSigOverSqrtSigBkg ) 
				'''

				SOB = TGraph( len(totalSig), binsHisto, sigOverSqrtSigBkg ) 
				dictNumVar[ var[0] ] = SOB 
				#if returnSOB in var[0]: return SOB
				#else: continue

			plotROC( 'SOB', '', dictNumVar, len(cutsList), False, True )
		else: 
			for bkg in bkgSamples: 
				dictBkg = { rocs: dictROC[ rocs ] for rocs in dictROC if bkg in rocs }
				plotROC( bkg, bkg, dictBkg, len(cutsList), False, False )

def plotROC( name, sample, dictROC, numCuts, varOrbkg, SOB, diffMasses=False ):
	"""docstring for plotROC"""

	dictSOB = OrderedDict()
	dictSOBError = OrderedDict()
	dictSOBLowEdge = OrderedDict()

	if diffMasses:
		dictROC = OrderedDict()
		tmpName = name.replace( 'SOB_', '' )

		#for masses in [80, 90, 100, 110, 120, 130, 140, 150, 170, 180, 190, 210, 220, 230, 240, 300 ]: 
		for masses in [ 240, 300, 350, 450, 500, 550, 600 ]: 
			dictROC[ 'M_{#tilde{t}} = '+str(masses)+' GeV' ] = makeROCs( 'ROCfiles/ROC'+args.boosted+'Values_QCD'+args.qcd+'_'+args.decay+str(masses)+'_'+args.version+args.grooming+'.txt', [ [tmpName] ], '', True, '', '', returnSOB=tmpName)
	
	f1 = TF1("f1","x",0,1)
	f1.SetLineColor(1)
	f1.SetLineStyle(3)
	can = TCanvas('c1', 'c1',  10, 10, 1000, 750 )
	can.SetGrid()
	#can.SetLogy()
	
	PT = TText(0.6, 0.05, sample )
	multiGraph = TMultiGraph()
	#legend=TLegend(0.65,0.15,0.90,0.45)
	legend=TLegend(0.70,0.70,0.90,0.87)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)
	dummy=1
	for h in dictROC:
		#tmpName = h.replace( name+'_', '' ).replace( '_ROC', '' )
		tmpName = str(h)
		legend.AddEntry( dictROC[ h ], tmpName, 'l' )
		#for var in variables:
		#	if var[0] in tmpName: varColor = var[7]
		#dictROC[h].SetLineColor( bkgSamples[ tmpName ][1] if varOrbkg else varColor )
		dictROC[h].SetLineColor( dummy ) 
		dictROC[h].SetLineWidth( 2 )
		dictROC[h].SetMarkerStyle(4)
		multiGraph.Add( dictROC[h] )
		dummy += 1
		if dummy == 5: dummy = 6
		if dummy == 9: dummy = 40
	multiGraph.Draw("ALP")
	PT.Draw()
	if SOB:
		if diffMasses: multiGraph.GetXaxis().SetTitle( 'Leading jet #tau_{21}' )
		else: multiGraph.GetXaxis().SetTitle( 'A.U.' )
		#multiGraph.GetYaxis().SetTitle('S/#sqrt{S+B}')
		multiGraph.GetYaxis().SetTitle('S/#sqrt{B}')
		#multiGraph.GetYaxis().SetTitle('#sqrt{2*((s+b)ln(s/b)-s)}')
	else:
		multiGraph.GetXaxis().SetRangeUser(-0.05,1.05)
		multiGraph.GetYaxis().SetRangeUser(-0.05,1.05)
		multiGraph.GetXaxis().SetTitle('Signal efficiency')
		multiGraph.GetYaxis().SetTitle('Bkg rejection')
		f1.Draw("same")
	multiGraph.GetYaxis().SetTitleOffset(0.95)
	legend.Draw()
	if diffMasses: can.SaveAs('Plots/'+name+'_'+args.boosted+'_QCD'+args.qcd+'_ROC_cut'+str(numCuts)+'_'+args.version+args.grooming+'.'+args.ext)
	else: can.SaveAs('Plots/'+name+'_'+args.boosted+args.decay+'_M'+str(args.mass)+'_QCD'+args.qcd+'_ROC_cut'+str(numCuts)+'_'+args.version+args.grooming+'.'+args.ext)
	del can

def plotSigOverBkg( name, sample, dictNum, numCuts ):
	"""docstring for plotROC"""

	can = TCanvas('c1', 'c1',  10, 10, 1000, 750 )
	can.SetGrid()

	PT = TText(0.7, 1.10, sample )
	#legend=TLegend(0.15,0.15,0.45,0.45)
	#legend.SetFillStyle(0)
	#legend.SetTextSize(0.03)
	SOB.SetLineWidth( 2 )
	SOB.SetMarkerStyle(4)
	SOB.Draw("ALP")
	PT.Draw()
	SOB.GetYaxis().SetTitleOffset(0.95)
	#legend.Draw()
	can.SaveAs('Plots/'+name+'_'+args.boosted+str(args.mass)+'_QCD'+args.qcd+'_SOB_cut'+str(numCuts)+'_'+args.version+args.grooming+'.'+args.ext)
	del can

def calcOptimizeDelta( BkgSamples, SigSamples, nameInRoot, massList ):
	"""docstring for calcOptimizeDelta"""

	QCDHisto = BkgSamples.Get( nameInRoot+'_QCD'+args.qcd+'All' )
	QCDHisto.Scale( QCDSF* args.lumi )

	signalHistos = OrderedDict()
	SOsqrtSB = []
	SOsqrtSBError = []
	for xmass in massList: 
		NAME = 'RPVStopStopToJets_'+args.decay+'_M-'+str( xmass )
		signalInputFile = TFile( SigSamples.replace( str(args.mass), str(xmass) ) )
		signalHistos[ xmass ] = signalInputFile.Get( nameInRoot+'_'+NAME )
		signalHistos[ xmass ].Scale( args.lumi )
		minMass = signalHistos[ xmass ].GetXaxis().FindBin( xmass-50 )
		maxMass = signalHistos[ xmass ].GetXaxis().FindBin( xmass+50 )
		bkgError = Double(0)
		bkgIntegral = QCDHisto.IntegralAndError( minMass, maxMass, bkgError )
		sigError = Double(0)
		sigIntegral = signalHistos[ xmass ].IntegralAndError( minMass, maxMass, sigError )

		sqrtSigBkg = sigIntegral / TMath.Sqrt( bkgIntegral )
		#sigBkgErr = TMath.Sqrt( (bkgError*bkgError) + (sigError*sigError) ) / (2*TMath.Sqrt( bkgIntegral + sigIntegral )) 
		#totalErrSOsqrtSB = TMath.Sqrt( TMath.Power( (sigError/sigIntegral), 2 ) + TMath.Power( (sigBkgErr/TMath.Sqrt( bkgIntegral + sigIntegral )), 2) )
		totalErrSOsqrtSB = sqrtSigBkg * TMath.Sqrt( TMath.Power( sigError/sigIntegral, 2 ) + (0.25 * TMath.Power( bkgError,2 )/bkgIntegral ) )  ### Error in s/sqrt(B)
		#print minMass, maxMass, bkgIntegral, bkgError, sigIntegral, sigError, totalErrSOsqrtSB

		SOsqrtSB.append( sqrtSigBkg )
		SOsqrtSBError.append( totalErrSOsqrtSB )

	#graphSOsqrtSB = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', SOsqrtSB ), array('d', [0]*len(massList) ), array('d', SOsqrtSBError) )
	graphSOsqrtSB = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', SOsqrtSB ), array('d', [0]*len(massList) ), array('d', [0]*len(SOsqrtSB) ) )
	return graphSOsqrtSB
	
def optimizeDelta( BkgSamples, SigSamples, massList ):
	"""docstring for optimizeDelta"""
	
	outputName = 'Plots/massAve_Btag_DeltaOptimization'+args.version+'.'+args.ext
	graphs = OrderedDict()
	'''
	graphs[ 'with delta' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta', massList )
	graphs[ 'without delta' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_deltaEta', massList )
	'''
	'''
	graphs[ 'delta150' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta150', massList )
	graphs[ 'delta200' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta200', massList )
	graphs[ 'delta250' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta250', massList )
	graphs[ 'delta300' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta300', massList )
	graphs[ 'delta350' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta350', massList )
	graphs[ 'delta400' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta400', massList )
	graphs[ 'delta450' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta450', massList )
	graphs[ 'delta500' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta500', massList )
	'''
	graphs[ 'Nominal' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta', massList )
	graphs[ '+ 1CSVL' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_1CSVv2L', massList )
	graphs[ '+ 2CSVL' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_2CSVv2L', massList )
	graphs[ '+ 1CSVM' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_1CSVv2M', massList )
	graphs[ '+ 2CSVM' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_2CSVv2M', massList )
	graphs[ '+ 1CSVT' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_1CSVv2T', massList )
	graphs[ '+ 2CSVT' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_2CSVv2T', massList )
	'''
	graphs[ 'delta150' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta150', massList )
	graphs[ 'delta250' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta250', massList )
	graphs[ 'delta300' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta300', massList )
	graphs[ '+ 1QGL' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_1qgl', massList )
	graphs[ '+ 2QGL' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_2qgl', massList )
	graphs[ '+ 4QGL' ] = calcOptimizeDelta( BkgSamples, SigSamples, 'massAve_delta_4qgl', massList )
	'''

	can = TCanvas('c1', 'c1',  10, 10, 1000, 750 )
	can.SetGrid()

	legend=TLegend(0.70,0.45,0.95,0.90)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)
	dummy=1
	multiGraph = TMultiGraph()
	for g in graphs:
		graphs[g].SetLineWidth( 2 )
		graphs[g].SetMarkerStyle(4)
		graphs[g].SetLineColor(dummy)
		legend.AddEntry( graphs[g], g, 'pl' )
		multiGraph.Add( graphs[g] )
		dummy+=1

	multiGraph.Draw("ALP")
	multiGraph.GetYaxis().SetTitleOffset(0.95)
	multiGraph.GetXaxis().SetTitle('Average dijet mass [GeV]')
	multiGraph.GetYaxis().SetTitle('S/#sqrt{B}')
	legend.Draw()
	can.SaveAs(outputName)
	del can

#----------------------------------------------------------------------
def RUNTMVATraining( BkgSample, SigSample, treename, outputFileName, variables ):

	bkgName = outputFileName.split('_')[1].replace(".root",'')

	print "\n", "="*80
	print ' ---> Running Training for', bkgName
	print "="*80

	# get signal and background data for training/testing
	sigFile, sigTree = getTree( SigSample, treename )
	bkgFile, bkgTree = getTree( BkgSample, treename )

	# everything is done via a TMVA factory
	outputFile = TFile(outputFileName, "recreate")
	#factory = TMVA.Factory("sigbkg", outputFile, "!V:Transformations=I;N;D")
	factory = TMVA.Factory( "TMVATraining"+args.boosted+bkgName+"_RPVSt"+args.mass, outputFile, "!V:!Silent:Transformations=I;D;P;G,D:AnalysisType=Classification")

	# define input var
	factory.AddTarget( 'massAve', 'F' )
	for k in variables: factory.AddVariable( k, "F" )

	#factory.SetWeightExpression( "weight" if ( 'Mini' in inputFiles ) else 'puWeight*lumiWeight' )
	factory.SetWeightExpression( 'puWeight*lumiWeight' )

	# define from which trees data are to be taken
	factory.AddSignalTree( sigTree ) 
	factory.AddBackgroundTree( bkgTree ) 

	# remove problematic events and specify how
	# many events are to be used
	# for training and testing
	if ( ('QCD' in bkgName) or ('TTJets' in bkgName)): bkgTrain = 100000
	else: bkgTrain = 0
	if ('QCD' in bkgName) : bkgTest = 100000
	else: bkgTest = 0
	sigCut = TCut("massAve>"+str(int(args.mass)-30)+" && massAve < "+str(int(args.mass)+30))
	bkgCut = TCut("massAve>"+str(int(args.mass)-30)+" && massAve < "+str(int(args.mass)+30))
	factory.PrepareTrainingAndTestTree(sigCut, bkgCut,
			#"nTrain_Signal=0:nTrain_Background=100000:SplitMode=Random:NormMode=NumEvents:!V")
			"nTrain_Background="+str(bkgTrain)+":nTest_Background="+str(bkgTest)+":SplitMode=Random:NormMode=NumEvents:!V")

	# define multivariate methods to be run
	'''
	factory.BookMethod( TMVA.Types.kCuts, 
			"Cuts", 
			"!H:!V:"\
			"FitMethod=MC:"\
			"EffSel:"\
			"SampleSize=200000:"\
			"VarProp=FSmart" )
	factory.BookMethod( TMVA.Types.kMLP,
			"MLP",
			"!H:!V:"\
			"VarTransform=N:"\
			"HiddenLayers=10:"\
			"TrainingMethod=BFGS")

	factory.BookMethod( TMVA.Types.kBDT,
			"BDT",
			"!V:"\
			"BoostType=AdaBoost:"\
			"NTrees=400:"\
			"nEventsMin=100:"\
			"nCuts=50")
	factory.BookMethod( TMVA.Types.kBDT, 
			"BDTG",
			"!H:!V:NTrees=100:"\
			"BoostType=Grad:"\
			"Shrinkage=0.30:"\
			"UseBaggedGrad:"\
			"GradBaggingFraction=0.5:"\
			"SeparationType=GiniIndex:"\
			"nCuts=20:"\
			"NNodesMax=5" )
	'''
	factory.BookMethod( TMVA.Types.kCuts, 
			"CutsGA",
			"!H:!V:"\
			"FitMethod=GA:"\
			#"CutRangeMin[0]=-1:"\
			#"CutRangeMin[1]=0:"\
			#"CutRangeMax[0]=10:"\
			#"VarProp[1]=FMax:"\  ## varProp for gaussian shape variables
			#"EffSel:Steps=30:"\  ## Default 40 
			"Cycles=3:"\
			#"PopSize=400:"\ ## Default 300
			#"SC_steps=10:"\  ## Default 10 
			#"SC_rate=5:"\  ## default 5
			#"SC_factor=0.95"  ## default
			)

	factory.TrainAllMethods()  
	factory.TestAllMethods()
	factory.EvaluateAllMethods()

	outputFile.Close()
	del factory


def getOptimizeValues( variables, sample ):
	"""docstring for getOptimizeValues"""

	dictVariables = {}
	for k in variables: dictVariables[ k ] = array('f', [0.] )
	
	reader = TMVA.Reader( "!Color:Silent" )
	for k in variables: reader.AddVariable( k, dictVariables[k] )
	reader.BookMVA( "CutsGA method", "weights/TMVATraining"+args.boosted+sample+"_RPVSt"+args.mass+"_CutsGA.weights.xml" )
	#passed = reader.EvaluateMVA( "CutsGA method", 0.9 )
	mcuts = reader.FindMVA( "CutsGA method" )
	cutsMin = array( 'd', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] )        
        cutsMax = array( 'd', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] )
	finalValues = []
	for k in np.arange(0,1,0.01):
		mcuts.GetCuts( k, cutsMin, cutsMax )
		for ivar in range( len(variables) ):
			#print "... %.2g: %.5g < %s <= %.5g" % (k, cutsMin[ivar], variables[ivar], cutsMax[ivar])
			finalValues.append( [ k, variables[ivar], cutsMin[ivar], cutsMax[ivar] ] )
				
	return finalValues


def ApplicationCreateCombinedTree( variables, outputFileName, bkgSamples, SigSample, treename ):
	"""docstring for ApplicationCreateCombinedTree"""

	# Create a new root output file.
	#outfileName = TString( "tmva_example_multiple_backgrounds__applied.root" )
	outputFile = TFile( outputFileName, "RECREATE" )
	outputTree = TTree("multiBkg","multiple backgrounds tree")

	dictVariables = {}
	for k in variables:
		dictVariables[ k ] = array('f', [0.] )
		outputTree.Branch( k, dictVariables[ k ], k+'/F' )
	
	MassAve = array( 'f', [ 0.] )
	outputTree.Branch( 'massAve', MassAve, 'massAve/F' ) 
	Weight = array( 'f', [ 0.] )
	outputTree.Branch( 'weight', Weight, 'weight/F' ) 
	ClassID = array( 'f', [ 0.] )
	outputTree.Branch( 'classID', ClassID, 'classID/F' ) 
	for bkg in bkgSamples:
		dictVariables[ 'cls'+bkg ] = array('f', [0.] )
		outputTree.Branch( 'cls'+bkg, dictVariables[ 'cls'+bkg ], 'cls'+bkg+'/F' )

	# ===== create three readers for the three different signal/background classifications, .. one for each background
	readers = {}
	for bkg in bkgSamples:
		readers[ 'reader'+bkg ] = TMVA.Reader( "!Color:!Silent" )
		for k in variables: readers[ 'reader'+bkg ].AddVariable( k, dictVariables[k] )
		readers[ 'reader'+bkg ].BookMVA( "BDT method", "weights/TMVATraining"+bkg+"_BDTG.weights.xml" )

	sigFile, sigTree, sigNumEntries = getTree( SigSample, treename )
	for k in variables: sigTree.SetBranchAddress( k, dictVariables[ k ] )
	for i in xrange(sigNumEntries):
		sigTree.GetEntry(i)
		ClassID[0] = 0
		MassAve[0] = sigTree.massAve
		Weight[0] = sigTree.lumiWeight * sigTree.puWeight
		for bkg in bkgSamples:
			dictVariables[ 'cls'+bkg ][0] = readers[ 'reader'+bkg ].EvaluateMVA( 'BDT method' )
		outputTree.Fill()
	print ' ---> End of Signal Tree'


	dummyCls = 0
	for sample in bkgSamples: 
		dummyCls += 1
		bkgFile, bkgTree, bkgNumEntries  = getTree( bkgSamples[ sample ], treename )
		for k in variables: bkgTree.SetBranchAddress( k, dictVariables[ k ] )
		for i in xrange(bkgNumEntries):
			bkgTree.GetEntry(i)
			ClassID[0] = dummyCls 
			MassAve[0] = bkgTree.massAve
			Weight[0] = bkgTree.lumiWeight * bkgTree.puWeight
			for bkg in bkgSamples:
				dictVariables[ 'cls'+bkg ][0] = readers[ 'reader'+bkg ].EvaluateMVA( 'BDT method' )
			outputTree.Fill()
		print ' ---> End of', sample, 'Tree'
		
	#  write output tree
	outputFile.Write()
	outputFile.Close()

	del readers
	print "--- Created root file: ",  outputFileName, " containing the MVA output histograms."
	print "==> Application of readers is done! combined tree created"
#----------------------------------------------------------------------

#################################################################################
if __name__ == '__main__':

	usage = 'usage: %prog [options]'
	
	parser = argparse.ArgumentParser()
	parser.add_argument( '-m', '--mass', action='store', type=int, dest='mass', default=100, help='Mass of the Stop' )
	parser.add_argument( '-t', '--typeROC', action='store',  dest='typeROC', default='var', help='Type of ROC: variables only (var) or bkgs only (bkg)')
	parser.add_argument( '-g', '--grom', action='store', default='', dest='grooming', help='Grooming Algorithm, example: Pruned, Puppi.' )
	parser.add_argument( '-P', '--printValue', action='store',  dest='printValue', type=bool, default=False, help='Print values of ROCs')
	parser.add_argument( '-Q', '--quantity', action='store',  dest='quantity', default='massAsym', help='Variable to print ROC.')
	parser.add_argument( '-p', '--process', action='store',  dest='process', default='Simple', help='Process: simple or TMVA' )
	parser.add_argument( '-b', '--boosted', action='store',  dest='boosted', default='Boosted', help='Boosted or Resolved version' )
	parser.add_argument( '-e', '--eff', action='store', dest='effS', type=int, default=0, help='Mass of the Stop' )
	parser.add_argument( '-l', '--lumi', action='store', dest='lumi', type=int, default=35870, help='Mass of the Stop' )
	parser.add_argument( '-E', '--extension', action='store', dest='ext', default='png', help='Extension of plots.' )
	parser.add_argument( '-B', '--batchSys', action='store',  dest='batchSys', type=bool, default=False, help='BatchSys: False (lxplus), True hexfarm.' )
	parser.add_argument( '-v', '--version', action='store', default='v01', help='Version: v01, v02.' )
	parser.add_argument( '-q', '--qcd', action='store', default='HT', dest='qcd', help='Type of QCD binning, example: HT.' )
	parser.add_argument( '-d', '--decay', action='store', default='UDD312', dest='decay', help='Decay, example: UDD312, UDD323.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	if args.batchSys: folder = '/cms/gomez/archiveEOS/Archive/v8020/Analysis/'+args.version+'/'
	else: folder = 'Rootfiles/'
	QCDSF = ( 0.24 if 'Resolved' in args.boosted else ( 1 if 'Puppi' in args.grooming else 0.366 ) )  ### v09

	bkgSamples = OrderedDict()
	bkgSamples[ 'QCD'+args.qcd+'All' ] = [ folder+'/RUNAnalysis_QCD'+args.qcd+'All_80X_V2p4_'+args.version+'.root', kBlue-4 ]
	if 'Boosted' in args.boosted:
		bkgSamples[ 'TTJets' ] = [ folder+'/RUNAnalysis_TTJets_80X_V2p4_'+args.version+'.root', kGreen+2 ]
		bkgSamples[ 'WJetsToQQ' ] = [ folder+'/RUNAnalysis_WJetsToQQ_80X_V2p4_'+args.version+'.root', 38 ]
		bkgSamples[ 'ZJetsToQQ' ] = [ folder+'/RUNAnalysis_ZJetsToQQ_80X_V2p4_'+args.version+'.root', kOrange ]
		bkgSamples[ 'Dibosons' ] = [ folder+'/RUNAnalysis_Dibosons_80X_V2p4_'+args.version+'.root', kMagenta+2 ]
		#bkgSamples[ 'WWTo4Q' ] = [ folder+'/RUNAnalysis_WWTo4Q_80X_V2p4_'+args.version+'.root', kMagenta+2 ]
		#bkgSamples[ 'ZZTo4Q' ] = [ folder+'/RUNAnalysis_ZZTo4Q_80X_V2p4_'+args.version+'.root', kOrange+2 ]
		#bkgSamples[ 'WZ' ] = [ folder+'/RUNAnalysis_WZ_80X_V2p4_'+args.version+'.root', kCyan ]

	sigSamples = {}
	sigSamples[ args.mass ] = folder+'/RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_80X_V2p4_'+args.version+'.root'

	SigSample = folder+'/RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'-madgraph_RunIISpring15MiniAODv2-74X_Asympt25ns_v09_v03.root'
	treename = "ResolvedAnalysisPlots/RUNATree" if ( 'Resolved' in args.boosted ) else 'BoostedAnalysisPlots'+args.grooming+'/RUNATree'

	var = [
		#[ 'Resolved', 'mindR', 50, 0., 5., True, 0., 0.8 ],
		#[ 'Resolved', 'jetsQGL[0]', 50, 0., 1., False, 0., 0., 5, 1 ],
		#[ 'Resolved', 'jetsQGL[1]', 50, 0., 1., False, 0., 0., 5, 1 ],
		#[ 'Resolved', 'jetsQGL[2]', 50, 0., 1., False, 0., 0., 5, 1 ],
		#[ 'Resolved', 'jetsQGL[3]', 50, 0., 1., False, 0., 0., 5, 1 ],
		[ 'Resolved', 'deltaEta', 50, 0., 5., True, 0, 0., 5, 1 ],
		[ 'Resolved', 'massAsym', 20, 0., 1., True, 0., 0., 1, 2 ],
		#[ 'Resolved', 'cosThetaStar1', 20, 0., 1., True, 0., 0., 1, 3 ],
		#[ 'Resolved', 'cosThetaStar2', 20, 0., 1., True, 0., 0., 1, 4 ],
		[ 'Resolved', 'delta1', 300, 0, 300,  False, 0, 0., 300, 5 ],
		[ 'Resolved', 'delta2', 300, 0, 300, False, 0, 0., 300, 6 ],
		#[ 'Resolved', 'xi1', 20, 0., 1., True, 0, 0.6 ],
		#[ 'Resolved', 'xi2', 20, 0., 1., True , 0, 0.6],
		[ 'Boosted', "prunedMassAsym", 20, 0., 1., True, 0., 0.2, 1, 1 ],
		[ 'Boosted', "jet1CosThetaStar", 20, 0., 1, True, 0., 0.8, 1, 2 ],
#		[ 'Boosted', "jet2CosThetaStar", 20, 0., 1, True, 0., 0.8, 1, 3 ],
		[ 'Boosted', "jet1Tau21", 20, 0., 1., True, 0., 0., 1, 4 ],
		[ 'Boosted', "jet2Tau21", 20, 0., 1., True, 0., 0., 1, 5 ],
		[ 'Boosted', "jet1Tau31", 20, 0., 1., True, 0., 0., 1, 6 ],
		[ 'Boosted', "jet2Tau31", 20, 0., 1., True, 0., 0., 1, 7 ],
		[ 'Boosted', "jet1Tau32", 20, 0., 1., True, 0., 0 , 1, 8 ],
		[ 'Boosted', "jet2Tau32", 20, 0., 1., True, 0., 0 , 1, 9 ],
		[ 'Boosted', "deltaEtaDijet", 50, 0., 5., True, 0.0, 0.6, 5, 11 ],
		[ 'Boosted', "jet1SubjetPtRatio", 20, 0., 1., True, 0., 0 , 1, 12 ],
		[ 'Boosted', "jet2SubjetPtRatio", 20, 0., 1., True, 0., 0 , 1, 13 ],
		]
	

	if 'calcROC' in args.process: 
		variables = [ x[1:] for x in var if ( args.boosted in x[0] ) ]
		cuts = [ x[1:] for x in var if ( ( args.boosted in x[0] ) and ( x[6]!=x[3] ) ) ]
		p0 = Process( target=calcROCs, args=( bkgSamples, sigSamples, treename, variables, ( 50 if 'Resolved' in args.boosted else 20 ), cuts ) )
		p0.start()
		p0.join()
	
	elif 'plotROC' in args.process: 
		variables = [ x[1:] for x in var if ( ( args.boosted in x[0] ) and ( x[6]==x[3] ) ) ]
		cuts = [ x[1:] for x in var if ( ( args.boosted in x[0] ) and ( x[6]!=x[3] ) ) ]
		if args.printValue: 
			makeROCs( 'ROCfiles/ROC'+args.boosted+'Values_QCD'+args.qcd+'_'+args.decay+str(args.mass)+'_cut'+str(len(cuts)-1)+'_'+args.version+args.grooming+'.txt', cuts, bkgSamples, True if 'var' in args.typeROC else False, cuts, args.quantity )
		else: 
			makeROCs( 'ROCfiles/ROC'+args.boosted+'Values_QCD'+args.qcd+'_'+args.decay+str(args.mass)+'_cut'+str(len(cuts))+'_'+args.version+args.grooming+'.txt', variables, bkgSamples, True if 'var' in args.typeROC else False, cuts, args.quantity )

	elif 'tmp' in args.process: 	
		#plotROC( 'SOB_jet1Tau21', '', '', '', True, True, diffMasses=True )
		plotROC( 'SOB_massAsym', '', '', '', True, True, diffMasses=True )

	elif 'delta' in args.process:
		optimizeDelta( TFile( folder+'/RUNMiniResolvedAnalysis_QCD'+args.qcd+'All_Moriond17_80X_V2p4_'+args.version+'.root' ),
			folder+'/RUNMiniResolvedAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_Moriond17_80X_V2p4_'+args.version+'.root',
			[ 200, 240, 300, 400, 500, 600, 700, 800, 900] ) 

	elif 'TMVA' in args.process:
		variables = [ x[1] for x in var if ( args.boosted in x[0] ) ]
		for sample in bkgSamples: 
			p0 = Process( target=RUNTMVATraining, args=( bkgSamples[ sample ], SigSample, treename, 'Rootfiles/RUNTMVATraining_'+bkgSamples[ sample ].split('_')[1]+'_RPVSt'+str(args.mass)+'.root', variables ) )
			p0.start()
			p0.join()
		#outputFileName = 'Rootfiles/RUN'+args.boosted+'OptimizationStudiesTMP.root'
		#p1 = Process( target=ApplicationCreateCombinedTree, args=( variables, outputFileName, bkgSamples, SigSample, treename ) )
		#p1.start()
		#p1.join()
	elif 'Print' in args.process:
		results = {}
		variables = [ x[1] for x in var if ( args.boosted in x[0] ) ]
		for sample in bkgSamples: 
			optValues = getOptimizeValues( variables, sample )
			results[ sample ] = optValues

		listSamples = []
		finalValues = {}
		print '\t', '\t'.join(str(p) for p in variables) 
		for k in results:
			listSamples.append( k )
			tmpMax = []
			for q in results[ k ]:
				if ( q[0] == args.effS/100. ): 
					for v in variables:
						if (q[1] == v ): tmpMax.append( q[3] )
			print k, '\t', '\t'.join(str(round(p, 2)) for p in tmpMax) 
				

