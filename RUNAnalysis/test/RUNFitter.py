#!/usr/bin/env python

###################
### Make Fitting
###################

#from ROOT import RooRealVar, RooDataHist, RooArgList, RooArgSet, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf, RooMultiPdf
from ROOT import *
from array import array
import argparse
import glob,sys, os
import warnings
import random
import numpy as np
from collections import OrderedDict
from multiprocessing import Process
from itertools import combinations
try: 
	import RUNA.RUNAnalysis.CMS_lumi as CMS_lumi 
	from RUNA.RUNAnalysis.scaleFactors import * #scaleFactor as SF
	from RUNA.RUNAnalysis.histoLabels import labels, labelAxis 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
except ImportError:
	sys.path.append('../python') 
	import CMS_lumi as CMS_lumi 
	from scaleFactors import * #scaleFactor as SF
	from histoLabels import labels, labelAxis 
	import tdrstyle as tdrstyle

gSystem.SetIncludePath('-I$ROOFITSYS/include')
if os.access('RooPowerFunction.cxx', os.R_OK): ROOT.gROOT.ProcessLine('.L RooPowerFunction.cxx+')

gROOT.Reset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 1
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls( 5000000 )
gRandom.SetSeed(0)

gStyle.SetOptFit()
gStyle.SetStatY(0.91)
gStyle.SetStatX(0.95)
gStyle.SetStatW(0.15)
gStyle.SetStatH(0.15) 
gStyle.SetTextSize(0.5)

xline = array('d', [0,2000])
yline = array('d', [0,0])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)
#massBins = range( 0, 600, 20 ) + range( 600, 840, 30 ) + range( 840, 1080, 40 ) + range( 1080, 1280, 50 ) + range( 1280, 2000, 60 )
massBins = [ 0, 10, 20, 30, 41, 52, 63, 74, 86, 98, 111, 124, 137, 151, 165, 180, 195, 210, 226, 242, 259, 276, 294, 312, 331, 350, 370, 390, 412, 433, 455, 478, 502, 526, 551, 577, 603, 631, 659, 688, 717, 748, 779, 812, 845, 879, 914, 950, 988, 1026, 1066, 1106, 1148, 1191, 1235, 1281, 1328, 1376, 1426, 1477, 1529, 1583, 1639, 1696, 1755, 1816, 1878, 1942, 2008 ]

def createPseudoExperiment( rawFunction, parRawFunction, numEvents, minX, maxX, plot ):
	"""docstring for createPseudoExperiment"""

	randomNumEventsQCD = numEvents #random.randint( int(numEvents-round(TMath.Sqrt(numEvents))), int(numEvents+round(TMath.Sqrt(numEvents))) )
	print "randomNumber Of PseudoExperiment events", randomNumEventsQCD, numEvents, int(numEvents-round(TMath.Sqrt(numEvents))), int(numEvents+round(TMath.Sqrt(numEvents)))
	hMainPSE = TH1D("hbkgPSE", "hbkgPSE", int( (maxX-minX) ) , minX, maxX)
	hMainPSE.SetBinErrorOption(TH1.kPoisson)
	if isinstance( rawFunction, TH1 ):
		hMainPSE.FillRandom( rawFunction, int(randomNumEventsQCD) )
	else:
		fitFunction = rawFunction[0][0].Clone() 
		fitFunction.SetRange( minX, maxX )
		fitFunction.SetName(rawFunction[0][0].GetName()+"PSE")
		#### giving initial values to fit
		for k in range( len(parRawFunction) ): fitFunction.SetParameter(k, parRawFunction[rawFunction[0][0].GetName()][k] )
		hMainPSE.FillRandom( rawFunction[0][0].GetName()+"PSE", int(randomNumEventsQCD) )

	c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	c1.SetLogy()
	gStyle.SetOptFit()
	gStyle.SetStatY(0.94)
	gStyle.SetStatX(0.9)
	gStyle.SetStatW(0.15)
	gStyle.SetStatH(0.15) 
	hMainPSE.GetXaxis().SetTitle( histYaxis )
	hMainPSE.GetYaxis().SetTitle("Events / "+str(rebinX)+" GeV" ) # dN/dM_{bbjj} [GeV^{-1}]")
	#hMainPSE.GetYaxis().SetTitleOffset(1.2);
	hMainPSE.SetTitle("QCD PseudoExperiments")
	hMainPSE.Sumw2()
	hMainPSE.Draw()
	c1.SaveAs(outputDir+hist+"_PseudoExperiment_Fit_ResolvedAnalysis_"+args.version+"."+args.extension)
	del c1

	return hMainPSE

def rootFitter( inFile, hist, scale, fitFunctions, minX, maxX, rebinX, plot, log=True, MCPlot='' ):
	"""Simple rootFitter"""

	#### Initial histo
	if isinstance( inFile, TFile ):
		unbinnedRawHisto = inFile.Get( hist )
		unbinnedRawHisto.SetBinErrorOption(TH1.kPoisson)
	else: unbinnedRawHisto = inFile

	if isinstance( rebinX, int ):
		minBinX = int(minX/rebinX)+1
		maxBinX = int(maxX/rebinX)+1
		rawHisto = unbinnedRawHisto.Rebin( rebinX )
		if MCPlot: 
			MCHisto = MCPlot.Rebin( rebinX )
			print '7'*50, MCHisto.Integral()	
	else: 
		minBinX = rebinX[2].index(minX)+1
		maxBinX = rebinX[2].index(maxX)+1
		rawHisto = unbinnedRawHisto.Rebin( rebinX[0], rebinX[1], array( 'd', rebinX[2] ) )
		if MCPlot: 
			MCHisto = MCPlot.Rebin( rebinX[0], rebinX[1]+"QCD", array( 'd', rebinX[2] ) )
	rawHisto.Scale( scale, ('width' if not isinstance( rebinX, int ) else '' ) ) 		##### dividing each bin by the width
	if MCPlot: MCHisto.Scale( 1, ('width' if not isinstance( rebinX, int ) else '' ) ) 		##### dividing each bin by the width
	print '|----> Raw number of events: ', rawHisto.Integral()#, ', bin size:', rawHisto.GetBinWidth(1)

	#### extracting bin contents 
	print '|----> Creating histograms from bin', minBinX, '(', rawHisto.GetBinLowEdge(minBinX),') to ', maxBinX, '(', rawHisto.GetBinLowEdge(maxBinX),'), initial integral: ', rawHisto.Integral()
	tmpBinContent = []
	tmpBinError = []
	tmpMCBinContent = []
	tmpMCBinError = []
	for ibin in range( minBinX, maxBinX ):
		#print rawHisto.GetBinContent(ibin), rawHisto.GetBinCenter(ibin) #rawHisto.GetXaxis().GetBinWidth(ibin), rawHisto.GetBinContent(ibin) / rawHisto.GetXaxis().GetBinWidth(ibin)
		tmpBinContent.append( rawHisto.GetBinContent(ibin) ) 
		#if ( (ibin==11) or (ibin==12) or (ibin==13) ): tmpBinError.append( 100000 ) 
		#else: tmpBinError.append( rawHisto.GetBinError(ibin) ) 
		#tmpBinError.append( rawHisto.GetBinError(ibin) ) 
		tmpBinError.append( ( rawHisto.GetBinError(ibin) if ( rawHisto.GetBinContent(ibin) > 0 ) else 1.8/( 1 if rawHisto.GetBinWidth(ibin)==1 else rawHisto.GetBinWidth(ibin) ) ) ) 
		#print rawHisto.GetBinLowEdge( ibin ), rawHisto.GetBinContent(ibin) , rebinX, rawHisto.GetBinError(ibin)
		if MCPlot:
			tmpMCBinContent.append( MCHisto.GetBinContent(ibin) ) 
			tmpMCBinError.append( MCHisto.GetBinError(ibin) ) 

	binContents = np.array(tmpBinContent)
	binError = np.array(tmpBinError)
	binfinalError = []
	MCbinContents = np.array(tmpMCBinContent)
	MCbinError = np.array(tmpMCBinError)
	print '|----> bins :', binContents
	print '|----> bins error :', binError
	
	#### creating histo for fitting
	numBins = maxBinX - minBinX
	print '|----> Histogram min: ', minX, '(', minBinX, '), max: ', maxX, '(', maxBinX, '), numBins: ', numBins
	if isinstance( rebinX, int ): 
		finalHisto = TH1F("finalHisto"+hist, "finalHisto"+hist, numBins, minX, maxX)
		finalErrHisto = TH1F("finalErrHisto"+hist, "finalErrHisto"+hist, numBins, minX, maxX)
	else: 
		finalHisto = TH1F("finalHisto"+hist, "finalHisto"+hist, numBins-1, array('d',  rebinX[2][minBinX-1:maxBinX] ) )  
		finalErrHisto = TH1F("finalErrHisto"+hist, "finalErrHisto"+hist, numBins-1, array('d',  rebinX[2][minBinX-1:maxBinX] ) )  
		print rebinX[2][minBinX:maxBinX]
	finalHisto.Sumw2(False)
	for ibin in range( 0, numBins ):
		#print ibin, binContents[ibin], finalHisto.GetBinLowEdge(ibin+1), binError[ibin]
		finalHisto.SetBinContent( ibin+1, binContents[ibin] )
		#finalBinError = binError[ibin] 
		finalBinError = ( TMath.Sqrt( MCbinError[ibin]*MCbinError[ibin] + binError[ibin]*binError[ibin] ) if MCPlot else binError[ibin] )
		#try: print ibin, binContents[ibin], MCbinContents[ibin], finalBinError, binError[ibin], MCbinError[ibin]
		#except IndexError: continue
		finalHisto.SetBinError( ibin+1, finalBinError )
		binfinalError.append( finalBinError )
		finalErrHisto.SetBinContent( ibin+1, binError[ibin] )
	finalHisto.SetBinErrorOption(TH1.kPoisson)
	print '|----> Number of events: ', finalHisto.Integral()

	fitParameters = OrderedDict()
	fitParErrors = OrderedDict()
	for fitFunc in fitFunctions:

		#if ( 'gaus' in fitFunc[0].GetName() ): fitFunc[0].SetName( 'gaus'+hist )
		#### giving initial values to fit
		if( len(fitFunc[1])>0 ):
			for k in range( len(fitFunc[1]) ): fitFunc[0].SetParameter(k, fitFunc[1][k])

		### Fitting
		tmpFit = 1000000000000000000
		keepFitting = True
		while keepFitting:
			#finalHisto.Fit(fitFunc[0],"MIRS","",minX,maxX)
			fitResults = TFitResultPtr( finalHisto.Fit( fitFunc[0],"ESR"+("" if 'QCDPtAll' in hist else "LL"),"",minX,maxX) )
			#### this is just a trick to keep fitting...
			tmpFitFunc = finalHisto.GetFunction( fitFunc[0].GetName() )
			chi2Ndf = tmpFitFunc.GetChisquare() / tmpFitFunc.GetNDF()
			try: tmpStatus = ( tmpFit - fitFunc[0].GetParameter(1) ) / fitFunc[0].GetParameter(1)
			except ZeroDivisionError: tmpStatus = 0
			if ((tmpStatus < 0.001) or (chi2Ndf < 1.5 ) ): 
				keepFitting = False
				print "|----> Final fit for ", fitFunc[0].GetName(), 'with chi2/ndf', tmpFitFunc.GetChisquare() , tmpFitFunc.GetNDF()
				print "|----> Likelihood ratio: ", 2*fitResults.MinFcnValue()

			else: 
				tmpFit = fitFunc[0].GetParameter(1)
				print "|----> Keep fitting...."

		numParam = fitFunc[0].GetNpar()
		fitParameters[ fitFunc[0].GetName() ] = [ fitFunc[0].GetParameter(k) for k in range( numParam ) ]
		fitParErrors[ fitFunc[0].GetName() ] =  [ fitFunc[0].GetParError(k) for k in range( numParam ) ]
		print "|----> Fitter parameters for", fitFunc[0].GetName(), fitParameters[ fitFunc[0].GetName() ], fitParErrors[ fitFunc[0].GetName() ]
		finalFitFunc = fitFunc[0] 

		listBinCenter = []
		listFitValues = []
		listFitErrors = []
		for ibin in range( 0, numBins ):
			binCenter = finalHisto.GetBinCenter( ibin+1 )
			fitVal = fitFunc[0].Eval( binCenter )
			err = array( 'd', [0] )   ### error in fit
			fitResults.GetConfidenceIntervals( 1, 1, 1, array('d',[binCenter]), err, 0.683, False ) 
			listFitValues.append( fitVal )
			listBinCenter.append( binCenter )
			listFitErrors.append( err[0] )

		finalFitUp = TGraph( len(listFitValues), array( 'd', listBinCenter ), np.add( listFitValues, listFitErrors ) ) 
		finalFitDown = TGraph( len(listFitValues), array( 'd', listBinCenter ), np.subtract( listFitValues, listFitErrors ) ) 

		######### Plotting Histograms
		if plot:
			c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
			if log: c1.SetLogy()
			finalHisto.GetXaxis().SetTitle( histYaxis )
			finalHisto.GetYaxis().SetTitle( ( "Events / "+ str(finalHisto.GetBinWidth(1))  +"GeV" if isinstance( rebinX, int ) else "< Events / GeV >" )) 
			finalHisto.GetYaxis().SetTitleOffset(0.9);
			finalHisto.GetXaxis().SetRangeUser( minX-50, maxX+50 )
			finalHisto.SetTitle("")
			finalHisto.Draw()
			#rawHisto.Draw("histe same")
			c1.SaveAs(outputDir+hist.replace('ResolvedAnalysisPlots/','')+"_"+args.process+"_"+('doubleGaus' if args.doubleGaus else fitFunc[0].GetName())+"Fit_ResolvedAnalysis_"+args.version+"."+args.extension)
			#print finalHisto.GetBinWidth(1)
			del c1

			c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
			if log: c1.SetLogy()
			finalErrHisto.GetXaxis().SetTitle( histYaxis )
			finalErrHisto.GetYaxis().SetTitle( ( "Events / "+ str(finalHisto.GetBinWidth(1))  +"GeV" if isinstance( rebinX, int ) else "< Events / GeV >" )) 
			finalErrHisto.GetYaxis().SetTitleOffset(0.9);
			finalErrHisto.GetXaxis().SetRangeUser( minX-50, maxX+50 )
			finalErrHisto.SetTitle("")
			finalErrHisto.Draw()
			#rawHisto.Draw("histe same")
			c1.SaveAs(outputDir+hist.replace('ResolvedAnalysisPlots/','')+"_"+args.process+"_"+('doubleGaus' if args.doubleGaus else fitFunc[0].GetName())+"Fit_Errors_ResolvedAnalysis_"+args.version+"."+args.extension)
			#print finalErrHisto.GetBinWidth(1)
			del c1
	
	return [ fitParameters, fitParErrors, binContents, binfinalError, chi2Ndf, finalHisto, finalFitFunc, finalFitUp, finalFitDown, listFitErrors ]
##############################################################


def histoFunctionFit( nameHisto, initFitFunction, parameters, parErrors, massBin, massBinErr, minX, maxX, rebinX ):
	"""docstring for histoFunctionFit"""

	fitFunction = initFitFunction.Clone()
	for i in range(0, initFitFunction.GetNpar() ): 
		fitFunction.SetParameter( i, parameters[fitFunction.GetName()][i] )
		fitFunction.SetParError( i, parErrors[fitFunction.GetName()][i] )

	if isinstance( rebinX, int ): histoFit = TH1D( nameHisto, nameHisto, len(massBin) , minX, maxX)
	else: 
		minBinX = rebinX[2].index(minX)+1
		maxBinX = rebinX[2].index(maxX)+1
		numBins = maxBinX - minBinX
		histoFit = TH1D( nameHisto, nameHisto, numBins-1, array('d',  rebinX[2][minBinX-1:maxBinX] ) )  
	histoFit.Sumw2(False)

	for ibin in range( 0, len(massBin)):
		#print ibin, massBin[ibin], histoFit.GetBinLowEdge(ibin+1)
		histoFit.SetBinContent( ibin+1, massBin[ibin] )
		#histoFit.SetBinError( ibin, ( massBinErr[ibin] if ( massBin[ibin] > 0 ) else 1.8 ) )
		histoFit.SetBinError( ibin+1, massBinErr[ibin] )

	histoFit.SetBinErrorOption(TH1.kPoisson)
			
	fitResult = histoFit.Fit( fitFunction, "ELLSR", "", minX, maxX )
	print '%'*30, 2*fitResult.MinFcnValue()

	return histoFit, fitFunction


def fitToHisto( nameHisto, initFitFunction, minX, maxX, rebinX ):
	"""docstring for fitToHisto"""

	print '|----> For signal+bkg, parameters:' 
	for k in range( initFitFunction.GetNpar() ) :
		print k, initFitFunction.GetParameter(k) 

	if isinstance( rebinX, int ): newHisto = TH1D( nameHisto, nameHisto, (maxX-minX)/rebinX , minX, maxX)
	else: 
		minBinX = rebinX[2].index(minX)+1
		maxBinX = rebinX[2].index(maxX)+1
		numBins = maxBinX - minBinX
		newHisto = TH1D( nameHisto, nameHisto, numBins-1, array('d',  rebinX[2][minBinX-1:maxBinX] ) )  
	newHisto.Sumw2(False)

	points = []
	errPoints = []
	for ibin in range( 1, newHisto.GetNbinsX()+1 ):
		valIntegral = initFitFunction.Eval( newHisto.GetBinCenter(ibin) ) 
		newHisto.SetBinContent( ibin, valIntegral )
		newHisto.SetBinError( ibin, TMath.Sqrt(valIntegral) )
		points.append( valIntegral )
		errPoints.append( TMath.Sqrt(valIntegral) )

	newHisto.SetBinErrorOption(TH1.kPoisson)
	#c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	#newHisto.Draw()
	#c1.SaveAs(outputDir+hist+"_SignalPlusBkg_ResolvedAnalysis_"+args.version+"."+args.extension)
	#del c1
			
	return newHisto, points, errPoints


def residualAndPulls(dataPoints, dataErrPoints, function, histo, minX, maxX, rebinX, extraName='', altPull='' ):
	"""docstring for residualAndPulls"""

	if isinstance( rebinX, int ): 
		hPull = TH1D("hpull_"+extraName+function.GetName(), "hpull_"+extraName+function.GetName(), len(dataPoints) , minX, maxX)
		if altPull: hAltPull = TH1D("haltpull_"+extraName+function.GetName(), "haltpull_"+extraName+function.GetName(), len(dataPoints) , minX, maxX)
		hResidual = TH1D("hresidual_"+extraName+function.GetName(), "hresidual_"+extraName+function.GetName(), len(dataPoints) , minX, maxX)
	else: 
		minBinX = rebinX[2].index(minX)+1
		maxBinX = rebinX[2].index(maxX)+1
		numBins = maxBinX - minBinX
		hPull = TH1D("hpull_"+extraName+function.GetName(), "hpull_"+extraName+function.GetName(), numBins-1, array('d',  rebinX[2][minBinX-1:maxBinX] ) ) 
		if altPull: hAltPull = TH1D("haltpull_"+extraName+function.GetName(), "haltpull_"+extraName+function.GetName(), numBins-1, array('d',  rebinX[2][minBinX-1:maxBinX] ) ) 
		hResidual = TH1D("hresidual_"+extraName+function.GetName(), "hresidual_"+extraName+function.GetName(), numBins-1, array('d',  rebinX[2][minBinX-1:maxBinX] ) ) 
	hPull.Sumw2()
	if altPull: hAltPull.Sumw2()
	hResidual.Sumw2()

	######## Calculating Pull and Residual
	chi2 = 0 
	RSS = 0
	nof = 0
	for ibin in range(0, len(dataPoints) ):
	
		binCont = dataPoints[ibin]
		binErr = dataErrPoints[ibin]
		histo.GetBinCenter(ibin)
		valIntegral = function.Eval( histo.GetBinCenter(ibin+1) ) 
		diff = (binCont - valIntegral)/ valIntegral
		#errDiff = diff * TMath.Sqrt( TMath.Power( P4Gaus.GetParError(0) / P4Gaus.GetParameter(0),2 ) + TMath.Power( P4Gaus.GetParError(1)/ P4Gaus.GetParameter(1), 2 )  + TMath.Power( P4Gaus.GetParError(2)/ P4Gaus.GetParameter(2), 2 )  + TMath.Power( P4Gaus.GetParError(3)/ P4Gaus.GetParameter(3), 2 ) )
		#errDiff = diff * TMath.Sqrt( TMath.Power( function.GetParError(0) / function.GetParameter(0),2 ) + TMath.Power( function.GetParError(1)/ function.GetParameter(1), 2 )  + TMath.Power( function.GetParError(2)/ function.GetParameter(2), 2 )  + TMath.Power( function.GetParError(3)/ function.GetParameter(3), 2 ) )
		#print binCont, binErr, valIntegral, diff

		if (binCont != 0):
			pull = (binCont - valIntegral)/ binErr
			#print binCont, binErr, valIntegral, diff, pull, histo.GetBinCenter(ibin+1), hPull.GetBinLowEdge(ibin+1)
			chi2 += TMath.Power(pull,2)
			RSS += TMath.Power((binCont - valIntegral),2)
			nof += 1
			
			hPull.SetBinContent(ibin+1, pull)
			hPull.SetBinError(ibin+1, 1.0)
	
			hResidual.SetBinContent(ibin+1, diff)
			hResidual.SetBinError(ibin+1, binErr/valIntegral )

			#if altPull: 
			#altpull = (binCont - valIntegral)/abs(altPull[ibin])
			#hPull.SetBinContent(ibin+1, altpull)
			#hPull.SetBinError(ibin+1, 1.0)

		#print '|---> Significance of high mass bins: ', binCont, valIntegral, binErr, pull
		#print binCont, valIntegral

	NDoF = nof - function.GetNpar() - 1
	print '|----> ############### chi2 and nof: ', chi2, nof

	return hPull, hResidual, RSS, NDoF



def FitterCombination( inFileData, inFileBkg, inFileSignalName, hist, scale, bkgFunction, minX, maxX, rebinX ):
	"""docstring for FitterCombination"""

	### Fit QCD MC
	print "|----> Fitting MC QCD"
	hMCQCD = inFileBkg.Get( hist+('QCD'+args.qcd+'All' if args.miniTree else '') )
	hMCQCD.Scale(scale*QCDSF)
	htmpMCQCD = hMCQCD.Clone()
	BkgParameters = rootFitter( htmpMCQCD, 
					hist+('QCD'+args.qcd+'All' if args.miniTree else ''), 
					1,
					bkgFunction, 
					minX, 
					maxX, 
					rebinX, 
					True ) #( False if args.bkgAsData else True ) )

	bkgParameters = BkgParameters[0]
	bkgParErrors = BkgParameters[1]
	bkgpoints = BkgParameters[2]
	bkgpointsErr = BkgParameters[3]

	if args.final: legend=TLegend(0.60,0.63,0.95,0.88)
	else: legend=TLegend(0.18,0.15,0.50,0.35)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)
	
	if (not args.bkgAsData) or (args.pseudoExperiment):
		print "|----> Fitting Data"
		if args.pseudoExperiment:
			#inFileData = createPseudoExperiment( bkgFunction, bkgParameters, sum(bkgpoints), minX, maxX, True )
			rawHistoForPSE = hMCQCD.Clone()
			numBkgEvents = rawHistoForPSE.Integral()
			RawHistoForPSE = inFileData.Get( hist+('JetHT_Run2016' if args.miniTree else '') )
			numBkgEvents = RawHistoForPSE.Integral()
			#newFile = TFile( "pseudoExperiment"+args.decay+".root", 'recreate' )
			#inFileData = createPseudoExperiment( rawHistoForPSE, bkgParameters, numBkgEvents, 0, 3000, True )
			#newFile.Write()
			newFile = TFile( "pseudoExperiment"+args.decay+".root", 'open' )
			inFileData = newFile.Get( "hbkgPSE" )

		DataParameters = rootFitter( inFileData, 
						hist+('JetHT_Run2016' if args.miniTree else ''), 
						1,
						bkgFunction, 
						minX, maxX, 
						rebinX, 
						True, 
						MCPlot=(rawHistoForPSE if args.pseudoExperiment else "" ))

		points = DataParameters[2]
		pointsErr = DataParameters[3]
		hMain = DataParameters[5]
		mainP4 = DataParameters[6]
		mainP4Up = DataParameters[7]
		mainP4Down = DataParameters[8]
		functionErr = DataParameters[9]

		legend.AddEntry( hMain, ('PseudoExperiment' if args.pseudoExperiment else 'Data - '+('btagged' if 'UDD323' in args.decay else 'inclusive')+' selection'), 'ep' )
		legend.AddEntry( mainP4, ( args.func if args.comparison else 'Background Fit'), 'l' )
		legend.AddEntry( mainP4Up, ( args.func+" Fit Error" if args.comparison else 'Fit Error'), 'l' )

		funcDict = {}
		if args.comparison:
			for iF in range( 1, len(bkgFunction) ):
				tmpDummy, funcDict[ bkgFunction[iF][0].GetName() ] = histoFunctionFit( 'Data'+bkgFunction[iF][0].GetName(), 
													bkgFunction[iF][0], 
													bkgParameters, bkgParErrors, 
													bkgpoints, bkgpointsErr, 
													minX, maxX, rebinX )
				legend.AddEntry( funcDict[ bkgFunction[iF][0].GetName() ], bkgFunction[iF][0].GetName(), 'l' )
		
		#hMCQCD, qcdMCP4 = histoFunctionFit( 'QCD'+args.qcd+'All', 
		#					bkgFunction[0][0], bkgParameters, 
		#					bkgParErrors, bkgpoints, bkgpointsErr, minX, maxX )
		#legend.AddEntry( qcdMCP4, 'Fit to MC QCD pythia', 'l' )

		if args.final:
			signalFuncs = OrderedDict()
			hSBPulls = OrderedDict()
			hSBResiduals = OrderedDict()
			for imass in [ args.mass, 600 ]:
				signalMassWidth = int(2*( 10.48 + ( 0.04426 * imass) ))
				lowEdgeWindow = massBins[min(range(len(massBins)), key=lambda x:abs(massBins[x]-int(imass-2*signalMassWidth)))]
				highEdgeWindow = massBins[min(range(len(massBins)), key=lambda x:abs(massBins[x]-int(imass+2*signalMassWidth)))]
				signalFileName = inFileSignalName.replace(str(args.mass), str(imass) )
				tmpHistName = hist+('RPVStopStopToJets_'+args.decay+'_M-'+str(imass) if args.miniTree else '')
				SignalParameters = rootFitter( TFile( signalFileName ), 
						tmpHistName, 
						scale, 
						[ [ fitFunctions['gaus'][0].Clone(), [ 1, imass, 50 ] ] ], 
						lowEdgeWindow, 
						highEdgeWindow,
						rebinX, 
						True,
						False )
				signalFunction = SignalParameters[6] 
				signalFunction.SetTitle( str(imass) )
				signalFuncs[ imass ] = signalFunction
				legend.AddEntry( signalFuncs[ imass ], 'M_{#tilde{t}} = '+str(imass)+' GeV', 'l' )
		
				##### for plotting purposes, signal + bkg
				signalPlusBkg = TF1( 'signalPlusBkg'+str(imass), 'gaus+'+bkgFunction[0][0].GetName() )
				signalPlusBkg.SetParameter( 0, signalFuncs[ imass ].GetParameter( 0 ) ) 
				signalPlusBkg.SetParError( 0, signalFuncs[ imass ].GetParError( 0 ) ) 
				signalPlusBkg.SetParameter( 1, signalFuncs[ imass ].GetParameter( 1 ) ) 
				signalPlusBkg.SetParError( 1, signalFuncs[ imass ].GetParError( 1 ) ) 
				signalPlusBkg.SetParameter( 2, signalFuncs[ imass ].GetParameter( 2 ) ) 
				signalPlusBkg.SetParError( 2, signalFuncs[ imass ].GetParError( 2 ) ) 
				for i in range(0, len(DataParameters[0][args.func]) ): 
					signalPlusBkg.SetParameter( i+3, mainP4.GetParameter(i) )# DataParameters[1][args.func][i] )
					signalPlusBkg.SetParError( i+3, mainP4.GetParError(i) )# DataParameters[1][args.func][i] )
				hSignalPlusBkg, signalPlusBkgPoints, signalPlusBkgPointsErrors = fitToHisto( 'hSignalPlusBkg'+str(imass), 
														signalPlusBkg, 
														minX, maxX, rebinX )

				hSBPulls[imass], hSBResiduals[imass], sbchi2, sbNDF = residualAndPulls(signalPlusBkgPoints, 
													signalPlusBkgPointsErrors, 
													mainP4, hMain, 
													minX, maxX, rebinX, 
													extraName="SB"+str(imass) )

	else:
		hMain, mainP4 = histoFunctionFit( 'QCD'+args.qcd+'All', 
							bkgFunction[0][0], 
							bkgParameters, bkgParErrors, 
							bkgpoints, bkgpointsErr, 
							minX, maxX, rebinX )

		points = bkgpoints
		pointsErr = bkgpointsErr
		legend.AddEntry( hMain, 'QCD multijet MC', 'ep' )
		legend.AddEntry( mainP4, args.func, 'l' )
		funcDict = {}
		if args.comparison:
			for iF in range( 1, len(bkgFunction) ):
				tmpDummy, funcDict[ bkgFunction[iF][0].GetName() ] = histoFunctionFit( 'QCD'+args.qcd+'All'+bkgFunction[iF][0].GetName(), 
													bkgFunction[iF][0], 
													bkgParameters, bkgParErrors, 
													bkgpoints, bkgpointsErr, 
													minX, maxX, rebinX )
				legend.AddEntry( funcDict[ bkgFunction[iF][0].GetName() ], bkgFunction[iF][0].GetName(), 'l' )

	print '|----> DATA Plotted:', points
	print '|----> DATA Err:', pointsErr

	hPull, hResidual, chi2, NDF = residualAndPulls(points, pointsErr, mainP4, hMain, minX, maxX, rebinX, altPull=functionErr )

	######### Plotting Histograms
	maxXPlot = maxX+500
	tdrStyle.SetPadRightMargin(0.05)
	if args.final or args.comparison: 
		gStyle.SetOptFit(0)
		#sigBkgPoints = bkgpoints+sigpoints
		#hSigPull, hSigResidual, sigchi2, sigNDF = residualAndPulls(sigBkgPoints, pointsErr, mainP4, hMain, minX, maxX )
	else: 
		gStyle.SetOptFit()
		gStyle.SetStatY(0.91)
		gStyle.SetStatX(0.95)
		gStyle.SetStatW(0.15)
		gStyle.SetStatH(0.30) 

	c3 = TCanvas('c1', 'c1',  10, 10, (800 if args.comparison else 1250), 500 )
	if not args.comparison:
  		tdrStyle.SetPadLeftMargin(0.15)
		pad1 = TPad("pad1", "Fit",0,0.03,0.50,1.00,-1)
		pad2 = TPad("pad2", "Pull",0.50,0.45,1.00,0.95,-1);
		pad3 = TPad("pad3", "Residual",0.50,0,1.00,0.507,-1);
		pad1.Draw()
		pad2.Draw()
		pad3.Draw()

		pad1.cd()
		pad1.SetLogy()
		#pad1.SetLogx()
	else: c3.SetLogy()

	hMain.SetMarkerStyle(8)
	hMain.GetYaxis().SetTitle( ( "Events / "+ str(hMain.GetBinWidth(1))  +"GeV" if isinstance( rebinX, int ) else "< Events / GeV >" ) ) 
	hMain.GetXaxis().SetTitle( histYaxis )
	hMain.GetXaxis().SetTitleSize(0.055)
	if not args.comparison: hMain.GetYaxis().SetTitleOffset(1.15);
	hMain.SetTitle("")
	hMain.SetMarkerStyle(20)
	#hMain.SetMaximum( 1.5 * hMain.GetMaximum() )
	#hMain.SetMaximum( 200 )
	#hMain.SetMinimum( 10 )
	hMain.Draw('E0')
	hMain.GetXaxis().SetRangeUser( minX, maxXPlot  )
	if args.final:
		color = 4
		for imass in signalFuncs:
			signalFuncs[imass].SetLineStyle(5)
			signalFuncs[imass].SetLineColor(kRed-color)
			signalFuncs[imass].Draw("same")
			color+=1
	mainP4.SetLineColor(kBlue-4)
	mainP4.Draw("same")
	mainP4Up.SetLineColor(kBlue-4)
	mainP4Up.SetLineStyle(7)
	mainP4Up.Draw("same")
	mainP4Down.SetLineColor(kBlue-4)
	mainP4Down.SetLineStyle(7)
	mainP4Down.Draw("same")
	if args.comparison:
		color = 2
		for iF in funcDict:
			funcDict[iF].SetLineColor(color)
			funcDict[iF].Draw("same")
			color+=1
			if (color==4): color=6
	'''
	#if (not args.bkgAsData) or (args.pseudoExperiment) or (args.comparison):
		#qcdMCP4.SetLineColor( kMagenta )
		#qcdMCP4.Draw("same")	
		#if isinstance( inFileSignal, TFile):
		#	qcdHTMCP4.SetLineColor( kViolet )
		#	qcdHTMCP4.Draw("same")	
	'''
	legend.Draw("same")
	CMS_lumi.relPosX = 0.13
	CMS_lumi.cmsTextSize = 0.60
	CMS_lumi.lumiTextSize = 0.50
	CMS_lumi.CMS_lumi((c3 if args.comparison else pad1), 4, 0)
	#labels( hist, '', '', 0.20, 0.45, 'left' )

	if not args.comparison:
		pad2.cd()
		pad2.SetGrid()
		gStyle.SetOptStat(0)
		hPull.GetYaxis().SetTitle("#frac{(Data - Fit)}{Unc.}")
		hPull.GetXaxis().SetTitleSize(0.06)
		hPull.GetYaxis().SetLabelSize(0.07)
		hPull.GetYaxis().SetTitleSize(0.08)
		hPull.GetYaxis().SetTitleOffset(0.80)
		hPull.GetYaxis().CenterTitle()
		hPull.SetMarkerStyle(7)
		hPull.SetMaximum(4)
		hPull.SetMinimum(-3)
		hPull.GetXaxis().SetRangeUser( minX, maxXPlot )
		hPull.Draw("e")
		if args.final:
			color=4
			for imass in hSBPulls:
				hSBPulls[imass].SetLineStyle(5)
				hSBPulls[imass].SetLineColor(kRed-color)
				hSBPulls[imass].SetFillStyle(3004)
				hSBPulls[imass].SetFillColor(kRed-color)
				hSBPulls[imass].GetXaxis().SetRangeUser( imass-50, imass+50)
				hSBPulls[imass].Draw("same hist")
				color+=1
		line.Draw("same")
		
		pad3.cd()
		pad3.SetGrid()
		pad3.SetTopMargin(0)
		pad3.SetBottomMargin(0.3)
		gStyle.SetOptStat(0)
		hResidual.GetXaxis().SetTitle( histYaxis )
		hResidual.GetYaxis().SetTitle("#frac{(Data - Fit)}{Fit}")
		hResidual.GetXaxis().SetTitleSize(0.10)
		hResidual.GetXaxis().SetLabelSize(0.07)
		hResidual.GetYaxis().SetLabelSize(0.07)
		hResidual.GetYaxis().SetTitleSize(0.08)
		hResidual.GetYaxis().SetTitleOffset(0.80)
		hResidual.GetYaxis().CenterTitle()
		hResidual.SetMarkerStyle(7)
		hResidual.SetMaximum(4)
		hResidual.SetMinimum(-1)
		hResidual.GetXaxis().SetRangeUser( minX, maxXPlot )
		#hResidual.Sumw2()
		hResidual.Draw("e")
		if args.final:
			color=4
			for imass in hSBResiduals:
				hSBResiduals[imass].SetLineStyle(5)
				hSBResiduals[imass].SetLineColor(kRed-color)
				hSBResiduals[imass].SetFillStyle(3004)
				hSBResiduals[imass].SetFillColor(kRed-color)
				hSBResiduals[imass].GetXaxis().SetRangeUser( imass-50, imass+50)
				hSBResiduals[imass].Draw("same hist")
				color+=1
		line.Draw("same")

	c3.SaveAs("Plots/"+hist.replace('ResolvedAnalysisPlots/','')+"_"+('PseudoExperiment' if args.pseudoExperiment else ('MC' if args.bkgAsData else args.process))+"Fit"+( 'diff' if args.comparison else bkgFunction[0][0].GetName())+"_"+( 'final_' if args.final else '' )+"_"+( 'resoBasedBin_' if not isinstance( rebinX, int ) else '' )+"ResolvedAnalysis_"+args.version+"."+args.extension)
	#c3.SaveAs("Plots/"+hist.replace('ResolvedAnalysisPlots/','')+"_"+('PseudoExperiment' if args.pseudoExperiment else ('MC' if args.bkgAsData else args.process))+"Fit"+( 'diff' if args.comparison else bkgFunction[0][0].GetName())+"_"+( 'final_' if args.final else '' )+"_"+( 'resoBasedBin_' if not isinstance( rebinX, int ) else '' )+"ResolvedAnalysis_"+args.version+"_"+str(args.mass)+"."+args.extension)
	del c3



def createCards( dataFile, bkgFile, inFileSignal, listMass, hist, scale, bkgFunctions, minX, maxX, rebinX ):
	"""function to run Roofit and save workspace for RooStats"""
	
	warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='.*class stack<RooAbsArg\*,deque<RooAbsArg\*> >' )

	########## Fitting background and taking parameters
	if args.pseudoExperiment:
		newFile = TFile( "pseudoExperiment"+args.decay+".root", 'open' )
		dataFile = newFile.Get( "hbkgPSE" )

	print '|----> Background'
	bkgFuncParameters = rootFitter( dataFile, 
					hist+( 'QCD'+args.qcd+'All' if args.bkgAsData else ('JetHT_Run2016' if args.miniTree else '')) , 
					( scale*QCDSF if args.bkgAsData else 1 ), 
					bkgFunctions, 
					minX, 
					maxX, 
					rebinX, 
					True ) #( False if args.bkgAsData else True ))

	########## Storing parameters in RooRealVar
	dictPar = OrderedDict()
	for func in bkgFuncParameters[0]:
		for par in range( 1, len( bkgFuncParameters[0][func] ) ):
			tmpNamePar = func+'p'+str(par)
			dictPar[ tmpNamePar ] = RooRealVar( tmpNamePar, 
								tmpNamePar, 
								bkgFuncParameters[0][func][par],
								#-abs(bkgFuncParameters[1][func][par])*100, 
								#abs(bkgFuncParameters[1][func][par])*100)
								(bkgFuncParameters[0][func][par] - (2*bkgFuncParameters[1][func][par])), 
								(bkgFuncParameters[0][func][par] + (2*bkgFuncParameters[1][func][par]))
								)
			#print bkgFuncParameters[0][func][par], (bkgFuncParameters[0][func][par] - (2*bkgFuncParameters[1][func][par])), (bkgFuncParameters[0][func][par] + (2*bkgFuncParameters[1][func][par]))

	############ Creating data plot from input in bkgFuncParameters
	hData = TH1D("hData", "hData", len(bkgFuncParameters[2]) , minX, maxX)
	hData.Sumw2()
	for ibin in range(0, len( bkgFuncParameters[2] ) ):
		#print ibin+1, bkgFuncParameters[2][ibin]
		hData.SetBinContent( ibin+1, bkgFuncParameters[2][ibin] )
		hData.SetBinError( ibin+1, bkgFuncParameters[3][ibin] )
	print '|---> Events in data: ', hData.Integral()

	########### Starting signal fits
	listAcceptance = []
	listAccepError = []
	listChi2 = []
	listMean = []
	listMeanError = []
	listSigma = []
	listSigmaError = []
	listMass.sort()

	for imass in listMass:

		######## Fitting signal and extracting parameters
		TMPResolution = ( 10.48 + ( 0.04426 * imass) )
		print '|----> Signal'
		if (TMPResolution/rebinX) < 1 : tmpResolution = rebinX
		else: tmpResolution = rebinX* round(TMPResolution/rebinX)
		print '|----> tmp resolution', TMPResolution, tmpResolution
		SignalParameters = rootFitter( TFile.Open( inFileSignal.replace( str(args.mass), str(imass) ) ), 
						hist+('RPVStopStopToJets_'+args.decay+'_M-'+str(imass) if args.miniTree else ''), 
						scale, 
						[ [ fitFunctions['gaus'][0], ( [ 1, 300, 100, 1, imass, tmpResolution ] if args.doubleGaus else [ 1, imass, tmpResolution ] ) ] ], 
						( minX if args.doubleGaus else uncDict[imass][(9 if 'UDD312' in args.decay else 10)][0]),  
						( maxX if args.doubleGaus else uncDict[imass][(9 if 'UDD312' in args.decay else 10)][1]), 
						uncDict[imass][(9 if 'UDD312' in args.decay else 10)][2], ##rebinX, 
						True,
						False )

		jmass = int(SignalParameters[0]['gaus'][1])
		massWindow = int(SignalParameters[0]['gaus'][2])*2		### size of search
		#print TMPResolution, jmass, massWindow
		lowerLimitSearch = ( jmass-massWindow if args.window else minX )
		upperLimitSearch = ( jmass+massWindow if args.window else maxX )

		print '|----> mass window', massWindow, ', lowerLimitSearch', lowerLimitSearch, ', upperLimitSearch', upperLimitSearch, ', sigma: ',  SignalParameters[0]['gaus'][2]
		mass = RooRealVar( 'mass', 'mass', lowerLimitSearch, upperLimitSearch )
		#mass.setBins((upperLimitSearch-lowerLimitSearch))
		mypdfs = RooArgList()

		### creating RooStats objects for data
		rooDataHist = RooDataHist('data_obs','data_obs',RooArgList(mass),hData)
		rooDataHist2 = RooDataHist('data_obs','data_obs',RooArgList(mass),hData)
		rooDataHist.Print()

		### creating RooStats objects for signal
		sigMean = RooRealVar('sigMean', 'sigMean', 
				SignalParameters[0]['gaus'][1], 
				(SignalParameters[0]['gaus'][1] - (2*SignalParameters[1]['gaus'][1])),
				(SignalParameters[0]['gaus'][1] + (2*SignalParameters[1]['gaus'][1])),
				)
		sigSigma = RooRealVar('sigSigma', 'sigSigma', 
				SignalParameters[0]['gaus'][2], 
				(SignalParameters[0]['gaus'][2] - (2*SignalParameters[1]['gaus'][2])),
				(SignalParameters[0]['gaus'][2] + (2*SignalParameters[1]['gaus'][2])),
				)
		if args.doubleGaus:
			tmpPer = SignalParameters[0]['gaus'][0] / ( SignalParameters[0]['gaus'][0] + SignalParameters[0]['gaus'][3]) 
			sigConst1 = RooRealVar('sigConst1', 'sigConst1', tmpPer ) 
			signal1 = RooGaussian( 'signal1', 'signal1', mass, sigMean, sigSigma )
			sigConst2 = RooRealVar('sigConst2', 'sigConst2', 1-tmpPer ) 
			sigMean2 = RooRealVar('sigMean2', 'sigMean2', SignalParameters[0]['gaus'][4]) 
			sigSigma2 = RooRealVar('sigSigma2', 'sigSigma2', SignalParameters[0]['gaus'][5]) 
			signal2 = RooGaussian( 'signal2', 'signal2', mass, sigMean2, sigSigma2 )
			signal = RooAddPdf( 'signal', "signal1+signal2", RooArgList(signal1, signal2), RooArgList(sigConst1, sigConst2) )
		else:
			signal = RooGaussian( 'signal', 'signal', mass, sigMean, sigSigma )
		signal.Print()

		if args.doubleGaus: signalGaus = TF1( "RPVStop"+str(imass), "gaus(0)+gaus(3)", minX, maxX )
		else: signalGaus = TF1( "RPVStop"+str(imass), "gaus", lowerLimitSearch, upperLimitSearch )
		signalGaus.SetParameter( 0, SignalParameters[0]['gaus'][0] ) 
		signalGaus.SetParameter( 1, SignalParameters[0]['gaus'][1] ) 
		signalGaus.SetParameter( 2, SignalParameters[0]['gaus'][2] ) 
		if args.doubleGaus:
			signalGaus.SetParameter( 3, SignalParameters[0]['gaus'][3] ) 
			signalGaus.SetParameter( 4, SignalParameters[0]['gaus'][4] ) 
			signalGaus.SetParameter( 5, SignalParameters[0]['gaus'][5] ) 
			sigAcc = round( signalGaus.Integral( minX, maxX )/uncDict[imass][9][2], 2 )
		else: 
			sigAcc = round( signalGaus.Integral( lowerLimitSearch, upperLimitSearch)/uncDict[imass][(9 if 'UDD312' in args.decay else 10)][2], 2 ) 
		print '|----> Signal events:', sigAcc, ', total number of events:', signalGaus.Integral(0, 2000)

		#hBkg = inFileBkg.Get(hist)
		#hData = hBkg.Clone()
		#hData.Add( hSignal )


		### creating RooStats objects for background
		dictRooPdf = {}
		for bkgFunc in bkgFunctions:

			tmpBkgFunct = bkgFunc[0].Clone()
			for p in range( 0, tmpBkgFunct.GetNpar()): 
				tmpBkgFunct.SetParameter( p, bkgFuncParameters[0][ tmpBkgFunct.GetName() ][p] )
			if ( tmpBkgFunct.GetName() in args.func ):
				bkgAcc = round( tmpBkgFunct.Integral( lowerLimitSearch, upperLimitSearch ), 2 )
				print '|----> Background events in mass (', imass, '):', bkgAcc 
			#print round( tmpBkgFunct.Integral( 940, 1300), 2 )
			#sys.exit(0)
			#c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
			#c1.SetLogy()
			#tmpBkgFunct.Draw()
			#hData.Draw("hist same")
			#c1.SaveAs("test.png")
		
			tmpList = RooArgList() 
			tmpList.add( mass ) 
			for par in dictPar: 
				if tmpBkgFunct.GetName() in par: tmpList.add( dictPar[par] )
			
			if 'Bias' in args.process:
				print '|----> Adding ', bkgFunc[0].GetName()+'_pdf', ' to bias test'
				if 'Landau' in bkgFunc[0].GetName(): 
					LandauSigMean = RooRealVar('LandauSigMean', 'LandauSigMean', bkgFuncParameters[0]['Landau'][1] ) #, -abs(bkgFuncParameters[0]['Landau'][1])*100, abs(bkgFuncParameters[0]['Landau'][1])*100)
					LandauSigSigma = RooRealVar('LandauSigSigma', 'LandauSigSigma', bkgFuncParameters[0]['Landau'][2] ) #, -abs(bkgFuncParameters[0]['Landau'][2])*100, abs(bkgFuncParameters[0]['Landau'][2])*100)
					dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ] = RooLandau( bkgFunc[0].GetName()+'_pdf' , bkgFunc[2], mass, LandauSigMean, LandauSigSigma ) 
				else:
					dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ] = RooGenericPdf( bkgFunc[0].GetName()+'_pdf' , bkgFunc[2], tmpList ) 
				dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ].fitTo( rooDataHist, RooFit.SumW2Error(kTRUE) )
				dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ].Print()
				mypdfs.add( dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ] )

			else:
				background = RooGenericPdf('background', bkgFunc[2], tmpList ) 
				background.Print()

		#### S+B model this is to fit in RooFit
		#signal_norm = RooRealVar( 'signal_norm', 'signal_norm', sigAcc , (sigAcc-TMath.Sqrt(sigAcc)), (sigAcc+TMath.Sqrt(sigAcc)) )
		#signal_norm = RooRealVar( 'signal_norm', 'signal_norm', sigAcc, 0, 10000000000 ) #+(20.*TMath.Sqrt(sigAcc)))
		#signal_norm.setConstant()
		#signal_norm.Print()
		#background_norm = RooRealVar('background_norm','background_norm', bkgAcc, (bkgAcc-TMath.Sqrt(bkgAcc)), (bkgAcc+TMath.Sqrt(bkgAcc))) #, 0., 20.*TMath.Sqrt(bkgAcc))
		#background_norm.Print()
		#model = RooAddPdf("model","s+b",RooArgList(background,signal),RooArgList(background_norm,signal_norm))
		#res = model.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Strategy(0), RooFit.SumW2Error(kFALSE))
		#res.Print()
		###################################################

		##### creating workspace
		outputRootFile = '/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/Rootfiles/workspace_RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'_Resolved_'+args.cut+'_'+('BiasTest_' if 'Bias' in args.process else '')+('QCD_' if args.bkgAsData else ( 'PSE_' if args.pseudoExperiment else '' ))+('massWindow_' if args.window else '' )+('doubleGaus_' if args.doubleGaus else '' )+args.version+'.root' 
		if 'Bias' in args.process:
			cat = RooCategory( "pdf_index", "Index of Pdf which is active" )
			multipdf = RooMultiPdf( "roomultipdf", "All Pdfs", cat, mypdfs )
			norm = RooRealVar( "roomultipdf_norm", "Number of background events", 0, 20.*TMath.Sqrt(bkgAcc) )
			
			#### test of pdfs
			'''
			c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
			c1.SetLogy()
			plot = mass.frame()
			dummy=1
			for ifit in dictRooPdf: 
				dictRooPdf[ ifit ].plotOn( plot, RooFit.LineColor(dummy))
				dummy+=1
			plot.Draw()
			c1.SaveAs('Plots/BiasTest_pdfs.'+args.extension)
			'''

			bkgWS = RooWorkspace("bkgWS")
			getattr(bkgWS,'import')(cat)
			getattr(bkgWS,'import')(norm)
			getattr(bkgWS,'import')(multipdf, RooCmdArg())
			bkgWS.Print()
			bkgWS.writeToFile(outputRootFile, True)
		else: 
			myWS = RooWorkspace("myWS")
			getattr(myWS,'import')(background, RooCmdArg())
			getattr(myWS,'import')(signal)
			getattr(myWS,'import')(rooDataHist2) 
			##### IF included _norm, then the rate should be 1
			#getattr(myWS,'import')(background_norm)
			#getattr(myWS,'import')(signal_norm)
			myWS.Print()
			myWS.writeToFile(outputRootFile, True)

		print '|----> Workspace created:', outputRootFile

		##### write a datacard
		datacard = open('/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/Datacards/datacard_RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'_Resolved_'+args.cut+'_'+('BiasTest_' if 'Bias' in args.process else '')+('QCD_' if args.bkgAsData else ( 'PSE_' if args.pseudoExperiment else '' ) )+('massWindow_' if args.window else '' )+('doubleGaus_' if args.doubleGaus else '' )+args.version+'.txt','w')
		datacard.write('imax 1\n')
		datacard.write('jmax 1\n')
		datacard.write('kmax *\n')
		datacard.write('---------------\n')
		if 'Bias' in args.process:
			datacard.write("shapes background * "+outputRootFile+" bkgWS:roomultipdf \n")
			datacard.write("shapes data_obs * "+outputRootFile.replace('BiasTest_', '')+" myWS:data_obs \n")
			datacard.write("shapes signal * "+outputRootFile.replace('BiasTest_', '')+" myWS:signal\n")
		else: datacard.write("shapes * * "+outputRootFile+" myWS:$PROCESS \n")
		datacard.write('---------------\n')
		datacard.write('bin RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'\n')
		datacard.write('observation -1\n')
		datacard.write('------------------------------\n')
		datacard.write('bin          RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'\tRPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'\n')
		datacard.write('process      signal     background\n')
		datacard.write('process      0          1\n')
		datacard.write('rate         '+str(sigAcc)+'      '+str(bkgAcc)+'\n')
		#datacard.write('rate         1		1\n')    ##### rate is 1 is _norm are included
		datacard.write('------------------------------\n')
		if 'Bias' in args.process: datacard.write('pdf_index\tdiscrete\n')
		datacard.write("lumi\tlnN\t1.025\t-\n")
		datacard.write("triggerResolved\tlnN\t1.05\t-\n")
		datacard.write("sigMean\tparam\t"+str(round(SignalParameters[0]['gaus'][1],3))+'\t'+str(round((SignalParameters[0]['gaus'][1])*0.02,3))+"\t### jes shape unc.\n")
		datacard.write("sigSigma\tparam\t"+str(round(SignalParameters[0]['gaus'][2],3))+'\t'+str(round((SignalParameters[0]['gaus'][2])*0.1,3))+"\t### jer shape unc.\n")
		datacard.write("sigMean\tparam\t"+str(round(SignalParameters[0]['gaus'][1],3))+'\t'+str(round(SignalParameters[1]['gaus'][1],3))+"\t### statistical unc.\n")
		datacard.write("sigSigma\tparam\t"+str(round(SignalParameters[0]['gaus'][2],3))+'\t'+str(round(SignalParameters[1]['gaus'][2],3))+"\t### statistical unc.\n")
		datacard.write("SigNormJES\tlnN\t"+str(1+uncDict[imass][(1 if 'UDD312' in args.decay else 5)])+"\t-\n")
		datacard.write("SigNormJER\tlnN\t"+str(1+uncDict[imass][(2 if 'UDD312' in args.decay else 6)])+"\t-\n")
		datacard.write("SigNormPU\tlnN\t"+str(1+uncDict[imass][(3 if 'UDD312' in args.decay else 7)])+"\t-\n")
		datacard.write("SigNormPDF\tlnN\t"+str(1+uncDict[imass][(4 if 'UDD312' in args.decay else 8)])+"\t-\n")
		if 'UDD323' in args.decay: datacard.write('SigNormBtag\tlnN\t'+str(1+uncDict[imass][0])+"\t-\n") 
		datacard.close()
		print '|----> Datacard created:', datacard
	
		if args.acceptance:
			##### For acceptance
			eventsInAccep = sigAcc / ( scaleFactor( 'RPVStopStopToJets_'+args.decay+'_M-'+str( imass ) )* scale ) ### because sigAcc is scaled
			totalEvents = search( dictEvents, 'RPVStopStopToJets_'+args.decay+'_M-'+str( imass ) )[0]
			failedEvents = totalEvents - eventsInAccep
			accXeff = eventsInAccep / totalEvents 
			accXeffErr = TMath.Sqrt( (1/failedEvents) + (1/eventsInAccep) ) * failedEvents * eventsInAccep / TMath.Power( ( totalEvents ), 2 )
			listAcceptance.append( accXeff ) 
			listAccepError.append( accXeffErr )

			listMean.append( SignalParameters[0]['gaus'][1] )
			listMeanError.append( SignalParameters[1]['gaus'][1] )
			listSigma.append( SignalParameters[0]['gaus'][2] )
			listSigmaError.append( SignalParameters[1]['gaus'][2] )
			listChi2.append( SignalParameters[4] )


	if args.acceptance:

		print '|----> Creating acceptance plots'
		print listMass, listAcceptance, listAccepError
		accXeffGraph = TGraphErrors(len(listMass), array('d',listMass), array('d',listAcceptance), array('d',[0]*len(listMass)), array('d',listAccepError) )
		canAccep = TCanvas('canAccep', 'canAccep',  10, 10, 750, 500 )
		canAccep.SetLogy()

		accXeffGraph.SetLineColor(kRed)
		accXeffGraph.SetLineWidth(2)
		accXeffGraph.SetMarkerStyle(8)
		#legend.AddEntry( accXeffGraph, '#tau_{21} < 0.45', 'l' )

		accXeffGraph.Draw("ap")
		accXeffGraph.GetYaxis().SetTitle( 'Acceptance #times efficiency' )
		accXeffGraph.GetXaxis().SetTitle( "Resonance mass [GeV]" )
		accXeffGraph.GetYaxis().SetTitleOffset( 0.8 )
		accXeffGraph.GetYaxis().SetRangeUser( 0.0001, 0.1  )

		#legend.Draw()
		CMS_lumi.extraText = "Simulation Preliminary"
		CMS_lumi.lumi_13TeV = ''
		CMS_lumi.relPosX = 0.12
		CMS_lumi.CMS_lumi(canAccep, 4, 0)
		canAccep.SaveAs( 'Plots/signalAcceptance_massAve_'+args.cut+'_'+args.decay+'_ResolvedAnalysis_'+args.version+'.'+args.extension )
		del canAccep

		meanXeffGraph = TGraphErrors(len(listMass), array('d',listMass), array('d',listMean), array('d',[0]*len(listMass)), array('d',listMeanError) )
		canMean = TCanvas('canMean', 'canMean',  10, 10, 750, 500 )
		meanXeffGraph.SetLineColor(kRed)
		meanXeffGraph.SetLineWidth(2)
		meanXeffGraph.SetMarkerStyle(8)
		meanXeffGraph.Draw("ap")
		meanXeffGraph.GetYaxis().SetTitle( 'Mean' )
		meanXeffGraph.GetXaxis().SetTitle( "Resonance mass [GeV]" )
		meanXeffGraph.GetYaxis().SetTitleOffset( 0.8 )
		#meanXeffGraph.GetYaxis().SetRangeUser( 0.0001, 0.1  )
		CMS_lumi.extraText = "Simulation Preliminary"
		CMS_lumi.lumi_13TeV = ''
		CMS_lumi.relPosX = 0.12
		CMS_lumi.CMS_lumi(canMean, 4, 0)
		canMean.SaveAs( 'Plots/signalMean_massAve_'+args.cut+'_'+args.decay+'_ResolvedAnalysis_'+args.version+'.'+args.extension )
		del canMean

		sigmaXeffGraph = TGraphErrors(len(listMass), array('d',listMass), array('d',listSigma), array('d',[0]*len(listMass)), array('d',listSigmaError) )
		#sigmaXeffGraph = TGraph(len(listMass), array('d',listMass), array('d',listSigma) ) 
		#pol1Sigma = TF1( "pol1Sigma", "TMath::Log(x)", 200, 1400 )
		#print '|----> Fitting sigma'
		#sigmaXeffGraph.Fit(pol1Sigma)
		canSigma = TCanvas('canSigma', 'canSigma',  10, 10, 750, 500 )
		sigmaXeffGraph.SetLineColor(kRed)
		sigmaXeffGraph.SetLineWidth(2)
		sigmaXeffGraph.SetMarkerStyle(8)
		sigmaXeffGraph.Draw("ap")
		sigmaXeffGraph.GetYaxis().SetTitle( 'Sigma' )
		sigmaXeffGraph.GetXaxis().SetTitle( "Resonance mass [GeV]" )
		sigmaXeffGraph.GetYaxis().SetTitleOffset( 0.8 )
		#sigmaXeffGraph.GetYaxis().SetRangeUser( 0.0001, 0.1  )
		CMS_lumi.extraText = "Simulation Preliminary"
		CMS_lumi.lumi_13TeV = ''
		CMS_lumi.relPosX = 0.12
		CMS_lumi.CMS_lumi(canSigma, 4, 0)
		canSigma.SaveAs( 'Plots/signalSigma_massAve_'+args.cut+'_'+args.decay+'_ResolvedAnalysis_'+args.version+'.'+args.extension )
		del canSigma

		chi2XeffGraph = TGraph(len(listMass), array('d',listMass), array('d',listChi2) )
		canchi2 = TCanvas('canchi2', 'canchi2',  10, 10, 750, 500 )
		canchi2.SetLogy()
		chi2XeffGraph.SetLineColor(kRed)
		chi2XeffGraph.SetLineWidth(2)
		chi2XeffGraph.SetMarkerStyle(8)
		chi2XeffGraph.Draw("ap")
		chi2XeffGraph.GetYaxis().SetTitle( 'chi2' )
		chi2XeffGraph.GetXaxis().SetTitle( "Resonance mass [GeV]" )
		chi2XeffGraph.GetYaxis().SetTitleOffset( 0.8 )
		#chi2XeffGraph.GetYaxis().SetRangeUser( 0.0001, 0.1  )
		CMS_lumi.extraText = "Simulation Preliminary"
		CMS_lumi.lumi_13TeV = ''
		CMS_lumi.relPosX = 0.12
		CMS_lumi.CMS_lumi(canchi2, 4, 0)
		canchi2.SaveAs( 'Plots/signalchi2_massAve_'+args.cut+'_'+args.decay+'_ResolvedAnalysis_'+args.version+'.'+args.extension )
		del canchi2


def doftest( RSS1, RSS2, nPar1, nPar2, nBinsFit):
	"""docstring for doftest"""

	Fvalue = abs(( ( RSS1 - RSS2 ) / ( nPar2 - nPar1 ) ) / ( RSS2 / float( nBinsFit - nPar2 - 1 ) ))

	Fdist = TF1("Fdistr","TMath::Sqrt( (TMath::Power([0]*x,[0]) * TMath::Power([1],[1])) / (TMath::Power([0]*x + [1],[0]+[1])) / (x*TMath::Beta([0]/2,[1]/2)) )",0,100)
	Fdist.SetParameter( 0, (nPar2-nPar1) )
	Fdist.SetParameter( 1, (nBinsFit-nPar2) )
	CL = 1 - Fdistr.Integral( 0.00000001, Fvalue )
	altCL =  1. - TMath.FDistI( Fvalue, nPar2-nPar1, nBinsFit-nPar2 )

	return [ RSS1, RSS2, Fvalue, CL, altCL ]

	

def FisherTest( dataFile, hist, bkgFunctions, minX, maxX, rebinX ):
	"""docstring for FisherTest"""

	if args.bkgAsData:
		fitParameters = rootFitter( dataFile, 
				hist+( 'QCD'+args.qcd+'All' if args.miniTree else ''), 
				scale, 
				bkgFunctions,
				minX, 
				maxX, 
				rebinX, 
				False )
	else:
		if args.pseudoExperiment:
			#rawHistoForPSE = dataFile.Get( hist+('JetHT_Run2016' if args.miniTree else '') )
			#numBkgEvents = rawHistoForPSE.Integral()
			#dataFile = createPseudoExperiment( rawHistoForPSE, '', numBkgEvents, 0, 3000, True )
			newFile = TFile( "pseudoExperiment"+args.decay+".root", 'open' )
			dataFile = newFile.Get( "hbkgPSE" )
		fitParameters = rootFitter( dataFile, 
				hist+( 'JetHT_Run2016' if args.miniTree else ''), 
				1, 
				bkgFunctions,
				minX, 
				maxX, 
				rebinX, 
				False )

	dictDataAndFunc = OrderedDict()
	dictPullResChi2NDF = OrderedDict()
	for func in bkgFunctions:
		dictDataAndFunc[ func[0].GetName() ] = func[0]
		dictPullResChi2NDF[ func[0].GetName() ] = residualAndPulls( fitParameters[2], 
										fitParameters[3], 
										func[0], 
										fitParameters[5],
										minX, maxX, rebinX )

	dictFtest = OrderedDict()
	for key1, key2 in combinations(dictPullResChi2NDF.keys(), r = 2):
		dictFtest[ key1+'_'+key2 ] = doftest( dictPullResChi2NDF[ key1 ][2], 
							dictPullResChi2NDF[ key2 ][2], 
							dictDataAndFunc[ key1 ].GetNpar(),  
							dictDataAndFunc[ key2 ].GetNpar(), 
							len(fitParameters[2]) )

		print '|----> Ftest for ', key1, 'and', key2, ':', dictFtest[ key1+'_'+key2 ], 'Residual1, Residual2, Fvalue, CL, altCL'
		
##########################################################

def drawBiasTest( listMasses, folderRootfiles ):
	"""docstring for plotLimits"""
	
	dictTest = {}
	dictTestHisto = {}
	listMeans = {}
	listMeansErr = {}
	for t in range( 0, args.numBiasTest ):
		listMeans[t] = []
		listMeansErr[t] = []

	for mass in listMasses:
		
		for t in range( 0, args.numBiasTest ):
			dictTest[ mass+t ] = TChain( 'tree_fit_sb' )
			dictTest[ mass+t ].Add( folderRootfiles+'fitDiagnostics_RPVStopStopToJets_'+args.decay+'_M-'+str(mass)+'_Resolved_'+args.cut+'_BiasTest_signal'+str(args.signalInj)+'_'+args.version+'_Index0ToIndex'+str(t)+'.root') 
			dictTestHisto[ mass+t ] = TH1F( 'h'+str(mass)+str(t), 'h'+str(mass)+str(t), 
							(16 if (args.signalInj==1) else 40), 
							(-4 if (args.signalInj==1) else -10), 
							(4 if (args.signalInj==1) else 10))
			#dictTest[ mass+t ].Draw( "(r-"+str(args.signalInj)+")/rErr>>h"+str(mass)+str(t))
			dictTest[ mass+t ].Draw( "(mu-"+str(args.signalInj)+")/muErr>>h"+str(mass)+str(t))
			
			dictTestHisto[ mass+t ].SetMarkerSize(2)
			dictTestHisto[ mass+t ].SetLineWidth(2)
			tmpMax = dictTestHisto[ mass+t ].GetMaximumBin()
			print tmpMax
			#sys.exit(0)
			dictTestHisto[ mass+t ].SetStats(True)
			dictTestHisto[ mass+t ].Fit("gaus")
			dictTestHisto[ mass+t ].SetLineColor(t+1)
			dictTestHisto[ mass+t ].GetFunction("gaus").SetLineColor(t+1)
			listMeans[t].append( dictTestHisto[ mass+t ].GetFunction("gaus").GetParameter( 1 ) )
			listMeansErr[t].append( dictTestHisto[ mass+t ].GetFunction("gaus").GetParError( 1 ) )


		c = TCanvas("c"+str(mass), "c"+str(mass), 800, 600)
		c.cd()
		gStyle.SetOptFit(1)

		dictTestHisto[ mass ].Draw("es")
		for t in range( 1, args.numBiasTest ):
			dictTestHisto[ mass+t ].Draw("e sames")

		c.Update()
		st = dictTestHisto[ mass ].GetListOfFunctions().FindObject("stats")
		st.SetX1NDC(.75)
		st.SetX2NDC(.95)
		st.SetY1NDC(.75)
		st.SetY2NDC(.95)
		st.SetTextColor(1)
		st1 = dictTestHisto[ mass+1 ].GetListOfFunctions().FindObject("stats")
		st1.SetX1NDC(.75)
		st1.SetX2NDC(.95)
		st1.SetY1NDC(.55)
		st1.SetY2NDC(.75)
		st1.SetTextColor(2)
		st2 = dictTestHisto[ mass+2 ].GetListOfFunctions().FindObject("stats")
		st2.SetX1NDC(.75)
		st2.SetX2NDC(.95)
		st2.SetY1NDC(.35)
		st2.SetY2NDC(.55)
		st2.SetTextColor(3)
		st3 = dictTestHisto[ mass+3 ].GetListOfFunctions().FindObject("stats")
		st3.SetX1NDC(.75)
		st3.SetX2NDC(.95)
		st3.SetY1NDC(.15)
		st3.SetY2NDC(.35)
		st3.SetTextColor(4)
		c.Modified()
		c.SaveAs( 'Plots/BiasTest_'+args.decay+'_M-'+str(mass)+'_signal'+str(args.signalInj)+'_ResolvedAnalysis_'+args.version+'.'+args.extension )


	masses = array( 'd', listMasses )	
	massesErr = array( 'd', [0]*len(listMasses) )	

	xline = array('d', [ 100, 2000, 2000, 100 ] )
	yline = array('d', [ -0.5, -.5, .5, .5 ] )
	box = TGraph( 4, xline, yline )
	box.SetFillColorAlpha(18, 0.5)
	box.SetFillStyle(3350)

	graphTest = {}
	for t in range( 0, args.numBiasTest ):
		graphTest[t] = TGraphErrors( len( masses), masses, array( 'd', listMeans[t] ), massesErr, array( 'd', listMeansErr[t] ) )

	c1 = TCanvas("c1"+str(mass), "c1"+str(mass), 800, 600)
	c1.cd()

	legend = TLegend(.50,.70,.90,.88)
	legend.SetTextSize(0.025)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
	#legend.SetHeader('95% CL upper limits')

	graphTest[0].GetXaxis().SetTitle("Resonance mass [GeV]")
	graphTest[0].GetYaxis().SetTitle("Means")
	graphTest[0].GetYaxis().SetTitleOffset(0.9)
	graphTest[0].GetYaxis().SetRangeUser(-2,2)
	graphTest[0].SetMarkerStyle(21)
	graphTest[0].SetMarkerColor(1)
	legend.AddEntry(graphTest[0],"Compare with P3","lp")
	graphTest[0].Draw("AP")
	box.Draw("F")
	gPad.RedrawAxis()
	graphTest[0].Draw("P")
	for t in range( 1, args.numBiasTest ):
		graphTest[t].SetMarkerStyle(21)
		graphTest[t].SetMarkerColor(t+1)
		legend.AddEntry(graphTest[t],"Compare function"+str(t),"lp")
		graphTest[t].Draw("P")
	#graphTest[1].Draw("P")
	#graphTest[2].Draw("P")
	#graphTest[3].Draw("P")
	#legend.AddEntry(graphTest[1],"Compare function CDF3","lp")
	#legend.AddEntry(graphTest[2],"Compare function expoPoli3","lp")
	#legend.AddEntry(graphTest[3],"Compare function altExpo4","lp")

    	legend.Draw()

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(c, 4, 0)

	#c1.SetLogy()
	c1.SaveAs( 'Plots/BiasTest_'+args.decay+'_signal'+str(args.signalInj)+'_MeanValues_ResolvedAnalysis'+args.version+'.'+args.extension)

##########################################################
def tmpPlot( ):
	"""docstring for tmpPlot"""

	masses = array( 'd', [ 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510])
	PSEvalues = array( 'd', [ 2.24678561, 1.966548158, 2.048505244, 2.084535605, 1.874874205, 1.381561921, 1.368623312, 1.242284854, 1.019427999, 0.7358551882, 0.798614449, 0.622256033, 0.4626225542, 0.4708705545, 0.4515574049, 0.2168075023, 0.2667836463, 0.1738658039] ) 
	Datavalues = array( 'd', [ 2.261087917, 1.966712071, 2.046098917, 2.077534816, 1.801417583, 1.34448088, 1.350154925, 1.209888639, 1.024352611, 0.731742915, 0.7915318358, 0.6102611321, 0.4600476406, 0.4581055092, 0.4559566661, 0.2083981714, 0.2500557876, 0.1666223384] )

	xline = array('d', [ 100, 2000 ] )
	yline = array('d', [.5, .5 ] )
	box = TGraph( 2, xline, yline )
	box.SetFillColorAlpha(18, 0.5)
	box.SetLineColor(3)
	box.SetLineWidth(2)
	box.SetFillStyle(3350)

	PSEgraph = TGraph( len( masses), masses, PSEvalues )
	Datagraph = TGraph( len( masses), masses, Datavalues )

	c1 = TCanvas("c1", "c1", 800, 600)
	c1.cd()

	legend = TLegend(.60,.70,.95,.88)
	legend.SetTextSize(0.025)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
	#legend.SetHeader('95% CL upper limits')

	PSEgraph.GetXaxis().SetTitle("Resonance mass [GeV]")
	PSEgraph.GetYaxis().SetTitle("R(m)")
	PSEgraph.GetYaxis().SetTitleOffset(0.9)
	#PSEgraph.GetYaxis().SetRangeUser(-2,2)
	PSEgraph.SetMarkerStyle(22)
	PSEgraph.SetMarkerColor(1)
	legend.AddEntry(PSEgraph,"Pseudoexperiment","p")
	Datagraph.SetMarkerStyle(23)
	Datagraph.SetMarkerColor(2)
	legend.AddEntry(Datagraph,"Data","p")
	PSEgraph.Draw("AP")
	Datagraph.Draw("P")
	box.Draw("L")
	gPad.RedrawAxis()

    	legend.Draw()

	CMS_lumi.extraText = "Preliminary"
	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(c1, 4, 0)

	c1.SaveAs( 'Plots/StartingPointFitTest_'+args.decay+'_ResolvedAnalysis'+args.version+'.'+args.extension)

##########################################################
####  Code below is just for example, DONT USE IT
def rooFitter( dataFile, bkgFile, inFileSignal, hist, scale, P4, minX, maxX, rebinX ):

	myWS = RooWorkspace("myWS")
	''' CONFIG FOR 732 WORKS PERFECT
	x = RooRealVar( "x", "x", 50., 190. )
	myWS.factory("EXPR:bkg_pdf('pow(1-(x/13000.0),p1)/pow(x/13000.0,p2+p3*log(x/13000.))', {x[50,190],p1[2116,2117],p2[-62,63],p3[-5,-4]})")
	myWS.factory("Gaussian:sig_pdf(x, mean[93,94], sigma[5,6])")
	myWS.factory("SUM:model(nsig[0,100000]*sig_pdf, nbkg[0,1000000]*bkg_pdf)")
	myWS.factory("SUM:model_b(nbkg[0,1000000]*bkg_pdf)")

	myWS.factory("EXPR:bkg_pdf('pow(1-(x/13000.0),p1)/pow(x/13000.0,p2+p3*log(x/13000.))', {x[70,600], p1[-4000,4000],p2[-300,300],p3[0,40]})")

	'''

	myWS.factory("x["+str(P4GausParameters[7])+","+str(P4GausParameters[8])+"]")
	bins = (P4GausParameters[8]-P4GausParameters[7])/10
	myWS.var("x").setBins(int(bins))
	myWS.factory("EXPR:bkg_pdf('pow(1-(x/13000.0),p1)/pow(x/13000.0,p2+p3*log(x/13000.))', {x, p1["+str(P4GausParameters[1])+"],p2["+str(P4GausParameters[2])+"],p3["+str(P4GausParameters[3])+"]})")
	#myWS.factory("EXPR:bkg_pdf('pow(1-(x/13000.0),p1)/pow(x/13000.0,p2+p3*log(x/13000.))', {x, p1[-2000,1000],p2[-200,200],p3[-20,20]})")
	#myWS.factory("Gaussian:sig_pdf(x, mean["+str(args.mass)+"], sigma[0,10])")
	myWS.factory("Gaussian:sig_pdf(x, mean["+str(P4GausParameters[5])+"], sigma["+str(P4GausParameters[6])+"])")
	#myWS.factory("Gaussian:sig_pdf(x, mean[93.29], sigma[5.52])")
	#myWS.factory("Gaussian:sig_pdf(x, mean[90,100], sigma[0,10])")
	myWS.factory("SUM:model_bkg( nbkg[0,100000]*bkg_pdf )")
	#myWS.factory("SUM:model_sig( nsig[0,10000]*sig_pdf )")
	myWS.factory("SUM:model( nbkg[0,100000]*bkg_pdf , nsig[0,100000]*sig_pdf)")
	myWS.Print()

	bkg_pdf = myWS.pdf("model_bkg")
	#signal_pdf = myWS.pdf("model_sig")
	#bkg_pdf = myWS.pdf("bkg_"+args.extension)
	#signal_pdf = myWS.pdf("sig_"+args.extension)
	pdf = myWS.pdf("model")

	mass = RooArgList( myWS.var("x") )
	h1 = inFileBkg.Get(hist)
	Bkg = h1.Clone()
	#Bkg.Scale( 1.5 )
	hData = h1.Clone()
	bkg = RooDataHist( 'bkg', 'bkg', mass, Bkg)

	#hSignal.Scale(0.1)
	data_sig = RooDataHist( 'data_sig', 'data_sig', mass, hSignal)
	#getattr( myWS, 'import')(data_sig)
	hData.Add( hSignal )
	data = RooDataHist( 'data', 'data', mass, hData)

	c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	c1.SetLogy()
	xframe = myWS.var("x").frame()
	bkg.plotOn( xframe )
	xframe.Draw()
	xframe.GetXaxis().SetTitle("Average Pruned Mass [GeV]")

	#### MINOS better than MIGRAD http://pprc.qmul.ac.uk/~bevan/yeti/fitting.pdf
	nll = bkg_pdf.createNLL(bkg, RooFit.Offset(1), RooFit.NumCPU(3), RooFit.Optimize(2) )
	m = RooMinuit(nll)
	m.migrad()
	m.hesse()
	m.minos()
	#bkg_pdf.fitTo( bkg, RooFit.Extended(kTRUE), RooFit.SumW2Error(kFALSE) )
	#bkg_pdf.fitTo( bkg, RooFit.Extended(), RooFit.Strategy(2), RooFit.Minos(), RooFit.Save(), RooFit.PrintEvalErrors(-1), RooFit.SumW2Error(kTRUE) ) 
	#bkg_pdf.fitTo( bkg, RooFit.Strategy(2), RooFit.Minos(), RooFit.Save(), RooFit.PrintEvalErrors(-1), RooFit.SumW2Error(kTRUE) ) 
	#bkg_pdf.fitTo( bkg,RooFit.Save(true),RooFit.Minimizer("Minuit2", "Migrad"),RooFit.SumW2Error(kTRUE), RooFit.PrintEvalErrors(-1) )
	bkg_pdf.plotOn( xframe )
	#bkg_pdf.plotOn( xframe, RooFit.Components("bkg_"+args.extension), RooFit.LineStyle(kDashed) )
	bkg_pdf.paramOn( xframe, RooFit.Layout(0.6,0.9,0.94))
	xframe.Draw()
	xframe.SetMaximum(100000)
	xframe.SetMinimum(0.00001)
	c1.SaveAs('Plots/'+hist+"_QCD_FitP4Gaus_rooFit_"+args.version+"."+args.extension)
	del c1

	'''
	c2 = TCanvas('c2', 'c2',  10, 10, 750, 500 )
	c2.SetLogy()
	x2frame = myWS.var("y").frame()
	data_sig1.plotOn( x2frame )
	x2frame.Draw()
	x2frame.GetXaxis().SetTitle("Average Pruned Mass [GeV]")

	#### MINOS better than MIGRAD http://pprc.qmul.ac.uk/~bevan/yeti/fitting.pdf
	nll2 = signal_pdf.createNLL(data_sig1, RooFit.Offset(1), RooFit.NumCPU(3), RooFit.Optimize(2) )
	m2 = RooMinuit(nll2)
	m2.migrad()
	m2.hesse()
	m2.minos()
	signal_pdf.plotOn( x2frame )
	signal_pdf.paramOn( x2frame, RooFit.Layout(0.6,0.9,0.94))
	x2frame.Draw()
	c2.SaveAs('Plots/'+hist+"_RPVSt100tojj_FitP4Gaus_rooFit_"+args.version+"."+args.extension)
	del c2
	'''

	c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	xframe = myWS.var("x").frame()
	data.plotOn( xframe )
	xframe.Draw()
	xframe.GetXaxis().SetTitle("Average Pruned Mass [GeV]")
	#pdf.fitTo( data, RooFit.Save(true), RooFit.PrintEvalErrors(-1) ) # RooFit.Minimizer("Minuit2", "Migrad") )
	#pdf.fitTo( data,RooFit.Save(true),RooFit.Minimizer("Minuit2", "Migrad"),RooFit.SumW2Error(kTRUE), RooFit.PrintEvalErrors(-1) )
	nll_data = pdf.createNLL(data, RooFit.Offset(1), RooFit.NumCPU(3), RooFit.Optimize(2) )
	m_data = RooMinuit(nll_data)
	m_data.migrad()
	m_data.hesse()
	m_data.minos()
	pdf.plotOn( xframe )
	pdf.plotOn( xframe, RooFit.Components("bkg_"+args.extension), RooFit.LineStyle(kDashed) )
	pdf.plotOn( xframe, RooFit.Components("sig_"+args.extension), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed) );
	pdf.paramOn( xframe, RooFit.Layout(0.6,0.9,0.94))
	xframe.Draw()
	c1.SaveAs('Plots/'+hist+'_QCD_RPVSt'+str(args.mass)+'tojj_FitP4Gaus_rooFit.'+args.extension)
	del c1

	#### PSEUDOEXPERIMENT AS DATA
	'''
	numBkg = round( myWS.var("nbkg").getVal() )
	#numBkg = round( myWS.var("nbkg").getVal()+myWS.var("nsig").getVal() )
	numData = random.randint( numBkg-round(sqrt(numBkg)), numBkg+round(sqrt(numBkg)) )
	data_obs = myWS.pdf("model_bkg").generateBinned(RooArgSet(mass),numData, RooFit.Name("data_obs")) 
	#print myWS.var("nbkg").getVal(), myWS.var("nsig").getVal()
	c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	c1.SetLogy()
	xframe = myWS.var("x").frame()
	data_obs.plotOn( xframe )
	xframe.Draw()
	c1.SaveAs('Plots/'+hist+"_PseudoExp_FitP4Gaus_rooFit_"+args.version+"."+args.extension)
	del c1
	'''

	####### QCD AS DATA
	#data_obs = RooDataHist( 'data_obs', 'data_obs', mass, Bkg)
	data_obs = RooDataHist( 'data_obs', 'data_obs', mass, hData)

	getattr( myWS, 'import')(data_obs)
	'''
	modelConfig = RooStats.ModelConfig( 'modelConfig', myWS )
	modelConfig.SetPdf( myWS.pdf("model") )
	#modelConfig.SetPdf( myWS.pdf("model_sig") )
	modelConfig.SetPdf( myWS.pdf("model_bkg") )
	poi = RooArgSet( myWS.var("nsig") )
	modelConfig.SetParametersOfInterest( poi )
	obs = RooArgSet( myWS.var("x") )
	modelConfig.SetObservables( obs )
	myWS.defineSet("nuisParams","p1,p2,p3,nbkg")
	modelConfig.SetNuisanceParameters( myWS.set("nuisParams") )
	getattr( myWS, 'import')(modelConfig)
	'''

	outputRootFile = '/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/Rootfiles/workspace_RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'_Resolved_'+args.cut+'_'+args.version+'.root' 
	myWS.writeToFile(outputRootFile, true )
	myWS.Print()

	####### Creating datacard
	outputfile = open('/afs/cern.ch/work/a/algomez/Substructure/CMSSW_7_1_5/src/MyLimits/datacard_RPVStop'+str(args.mass)+'.txt','w')
	outputfile.write("imax 1 channels\n")
	outputfile.write("jmax 1 backgrounds\n")
	outputfile.write("kmax 0 *\n")
	outputfile.write("-------------------------------\n")
	outputfile.write("shapes * * "+outputRootFile+" myWS:$PROCESS \n")
	outputfile.write("-------------------------------\n")
	outputfile.write("bin           1\n")
	outputfile.write("observation  -1\n")
	outputfile.write("-------------------------------\n")
	outputfile.write("bin           1          1\n")
	outputfile.write("process     sig_pdf bkg_pdf\n")
	outputfile.write("process       0          1\n")
	outputfile.write('rate          '+str( round( myWS.var("nsig").getVal() ) )+' '+str( round( myWS.var("nbkg").getVal() ) )+' \n')
	outputfile.write("-------------------------------\n")
	outputfile.write("# lumi    lnN     1.025         -     \n")
	outputfile.write("# GausSigma  param       45.1220       5.7784  \n")
	outputfile.write("# GausMean  param       1000.0000       10.0000  \n")
	outputfile.write("# SigNormFit   lnN    1.0600       - \n")
	outputfile.write("# SigNormPDF   lnN    1.0300       - \n")
	outputfile.write("# SigNormJES   lnN    1.0400       - \n")
	outputfile.write("# SigNormPU   lnN    1.0300        - \n")
	outputfile.write("# SigNormISR   lnN    1.1000        - \n")
	outputfile.write("# SigNormBtag   lnN    1.0000       - -\n")
	outputfile.write("# b  param       78.8566       33.5562\n")
	outputfile.write("# c  param       -8.4161       10.6159\n")
	outputfile.write("# d  param       -0.8140       1.5043\n")
	outputfile.write("# BkgNorm    lnN     -       2.0000\n")
	outputfile.close()

##########################################################
def rooFitterTree( inFileBkg, inFileSignal, inFileData, hist):
	"""function to run Roofit and save workspace for RooStats"""

	warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='.*class stack<RooAbsArg\*,deque<RooAbsArg\*> >' )
	
	myWS = RooWorkspace("myWS")

	x = RooRealVar( "x", "x", 50., 180. )
	myWS.factory("EXPR:bkg_pdf('pow(1-(x/13000.0),p1)/pow(x/13000.0,p2+p3*log(x/13000.))', {x[50,180],p1[0,1000],p2[0,100],p3[0,10]})")
	myWS.factory("Gaussian:sig_pdf(x, mean[90,110], sigma[0,10])")
	myWS.factory("SUM:model(nsig[0,10000]*sig_pdf, nbkg[0,1000000]*bkg_pdf)")
	myWS.Print()

	#x = myWS.var("x")
	mass = myWS.var("mass")
	pdf = myWS.pdf("model")

	#mass = RooArgList( x )
	massAveForFit = RooRealVar( "massAveForFit", "massAveForFit", 50., 180. )
	bkgTree = inFileBkg.Get('/RUNATree' )
	bkg = RooDataSet( "bkg", "bkg", RooArgSet(massAveForFit), RooFit.Import( bkgTree ) )
	signalTree = inFileSignal.Get('/RUNATree' )
	signal = RooDataSet( "signal", "signal", RooArgSet(massAveForFit), RooFit.Import( signalTree ) )
	dataTree = inFileData.Get('/RUNATree' )
	data = RooDataSet( "data", "data", RooArgSet(massAveForFit), RooFit.Import( dataTree ) )
	
	getattr( myWS, 'import')(data)

	c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	xframe = massAveForFit.frame()
	data.plotOn( xframe )
	xframe.Draw()
	xframe.GetXaxis().SetTitle( histYaxis )
	pdf.fitTo( data, RooFit.Save(true), RooFit.Minimizer("Minuit2", "Migrad") )
	#pdf.fitTo( data, RooFit.Save(true) , RooFit.Minimizer("Minuit2", "Migrad"), RooFit.SumW2Error(kTRUE) )
	pdf.plotOn( xframe )
	pdf.plotOn( xframe, RooFit.Components("bkg_"+args.extension), RooFit.LineStyle(kDashed) )
	pdf.plotOn( xframe, RooFit.Components("sig_"+args.extension), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed) );
	pdf.paramOn( xframe, RooFit.Layout(0.6,0.9,0.94))
	xframe.Draw()
	c1.SaveAs('Plots/'+hist+"_QCD_RPVSt100tojj_FitP4Gaus_rooFitTree_"+args.version+"."+args.extension)
	del c1

	modelConfig = RooStats.ModelConfig( 'modelConfig', myWS )
	modelConfig.SetPdf( myWS.pdf("model") )
	poi = RooArgSet( myWS.var("nsig") )
	modelConfig.SetParametersOfInterest( poi )
	obs = RooArgSet( myWS.var("x") )
	modelConfig.SetObservables( obs )
	myWS.defineSet("nuisParams","p1,p2,p3,nbkg")
	modelConfig.SetNuisanceParameters( myWS.set("nuisParams") )

	getattr( myWS, 'import')(modelConfig)
	myWS.writeToFile("Rootfiles/workspace_QCD_RPVSt100tojj_FitP4Gaus_rooFitTree.root", True )
	myWS.Print()
####################################################################


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--process', action='store', default='Full', help='Type of fit to use.' )
	parser.add_argument('-v', '--version', action='store', default='v00', help='For Boosted of Resolved.' )
	parser.add_argument('-d', '--decay', action='store', default='UDD312', help='Decay, example: jj, bj.' )
	parser.add_argument('-m', '--mass', action='store', type=int, default=300, help='Decay, example: jj, bj.' )
	parser.add_argument('-s', '--signalInj', action='store', type=int, default=1, help='Signal injection for bias test.' )
	parser.add_argument('-q', '--qcd', action='store', default='Pt', dest='qcd', help='Type of QCD binning, example: HT.' )
	parser.add_argument('-l', '--lumi', action='store', type=float, default=1787, help='Luminosity, example: 1.' )
	parser.add_argument('-f', '--final', action='store_true', default=False, help='Final distributions.' )
	parser.add_argument('-b', '--bkgAsData', action='store_true', default=False, help='Background as data.' )
	parser.add_argument('-t', '--miniTree', action='store_true', default=False, help='miniTree: if plots coming from miniTree or RUNAnalysis.' )
	parser.add_argument('-e', '--extension', action='store', default='png', help='Extension of plots.' )
	parser.add_argument('-C', '--cut', action='store', default='delta', help='cut, example: cutDEta' )
	parser.add_argument('-F', '--func', action='store', default='P3', help='Function, example: P3, P4, P5' )
	parser.add_argument('-P', '--pseudoExperiment', action='store_true', default=False, help='Run pseudoexperiment.' )
	parser.add_argument('-a', '--acceptance', action='store_true', default=False, help='Create acceptance plots.' )
	parser.add_argument('-W', '--window', action='store_true', default=False, help='Mass window.' )
	parser.add_argument('-D', '--doubleGaus', action='store_true', default=False, help='Using double gaussian.' )
	parser.add_argument('-n', '--numBiasTest', action='store', type=int, default=3, help='Number of functions for bias test' )
	parser.add_argument('-k', '--comparison', action='store_true', default=False, help='Comparison plot.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)


	filePrefix = 'Rootfiles/RUNMiniResolvedAnalysis'
	hist = 'massAve_'+args.cut+'_'
	scale = args.lumi 
	rebinX = 20 
	if ( args.mass > 0 ): 
		listMass = [ args.mass ]
	else: 
		listMass = range( 400, 1050, 50 ) + range( 1100, 1400, 100 ) 
		#if '312' in args.decay: 
		#	listMass = [ 200, 220, 240 ] + range( 300, 1050, 50 ) + range( 1100, 1400, 100 ) 
		#else: listMass = range( 200, 320, 20 ) + range( 350, 1050, 50 ) #+ ( [1100, 1200 ] if 'CSVv2L' in args.cut else [] )


	QCDSF = 0.255
	outputDir = "Plots/"
	#signalFilename = filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_Moriond17_80X_V2p4_'+args.version+'.root'
	signalFilename = filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_Moriond17_80X_V2p4_v09p1.root'
	signalFile =  TFile.Open(signalFilename)
	#bkgFile = TFile.Open(filePrefix+'_QCD'+args.qcd+'All_Moriond17_80X_V2p4_'+args.version+'.root')
	bkgFile = TFile.Open(filePrefix+'_QCD'+args.qcd+'All_Moriond17_80X_V2p4_v09p1.root')
	#bkgFile2 = TFile.Open(filePrefix+'_QCDHTAll_Moriond17_80X_V2p4_'+args.version+'.root')
	#bkgFile3 = TFile.Open(filePrefix+'_QCDHerwigAll_Moriond17_80X_V2p4_'+args.version+'.root')
	dataFile = TFile.Open(filePrefix+'_JetHT_Run2016_80X_V2p4_'+args.version+'.root')
	#dataFile = TFile.Open(filePrefix+'_JetHT_Run2016_80X_V2p4_v09p10.root')
	

	###### Input parameters
	histYaxis = "Average dijet mass [GeV]"
	minFit = 260
	maxFit = 1600
	CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 1 ) )+" fb^{-1}"
	CMS_lumi.extraText = "Preliminary Simulation"

	######## Fit Functions
	fitFunctions = {}
	fitFunctions['P5'] = [ TF1("P5", "[0]*TMath::Power(1-(x/13000.0),[1])/(TMath::Power(x/13000.0,([2]+([3]*TMath::Log(x/13000.))+([4]*TMath::Power(TMath::Log(x/13000.),2)))))",0,2000), 
			( [1e-10, 47, 14, 1.94, 0.03 ] if 'UDD312' in args.decay else [1e-10, 39, 14, 1.91, 0.03]  ), 
			'pow(1-@0/13000,@1)/pow(@0/13000,@2+@3*log(@0/13000)+@4*pow(log(@0/13000),2))' ]

	fitFunctions['P4'] = [ TF1("P4", "[0]*TMath::Power(1-(x/13000.0),[1])/(TMath::Power(x/13000.0,([2]+([3]*TMath::Log(x/13000.0)))))",0,2000), 
			( [ 1e-22, -21, 26, 3 ] if 'UDD312' in args.decay else [3e-22, -24, 25, 3] ), 
			'pow(1-@0/13000,@1)/pow(@0/13000,@2+@3*log(@0/13000))' ]

	fitFunctions['P3'] = [ TF1("P3", "[0]*TMath::Power(1-(x/13000.0),[1])/(TMath::Power(x/13000.0,[2]))",0,2000), 
			( [ 3619, 114, 0.15] if 'UDD312' in args.decay else [753, 106, -0.02] ), 
			'pow(1-@0/13000,@1)/pow(@0/13000,@2)' ]

	#fitFunctions['P2'] = [ TF1("P2", "[0] / TMath::Power(x/13000.0,[1])",0,2000), 
	#		[ 1, 1 ], 
	#		'pow(@0/13000,-@1)' ]


	fitFunctions['altDijet'] = [ TF1("altDijet", "[0]+ TMath::Power( x/13000.0, [1]+([2]*TMath::Log(x/13000.0)))", 0, 3000 ),
			[], #[ -1, 1, 0.00001 ], #[ 1000, 1, 1 ],
			'exp( (@1*log(@0/13000)) + (@2*pow(log(@0/13000),2)) + ( @3*pow(log(@0/13000),3)) + ( @4*pow(log(@0/13000),4)))']

	fitFunctions['altExpo3'] = [ TF1("altExpo3", "[0]*pow( ( 1 - (x/13000) ), [2] )/pow(x,[1])", 0, 2000 ), 
			[15906., 0.157, 114.], 
			'pow( (1 - (@0/13000) ), @2 )/ pow(@0,@1)']

	fitFunctions['altExpo4'] = [ TF1("altExpo4", "[0]*pow( ( 1 - (x/13000) + [3]*x*x/TMath::Power(13000,2) ), [2] )/pow(x,[1])", 0, 2000 ), 
			[16., -1.27, 179., 1.98], #[ 1000, -0.1, 50, 0 ],
			'pow( (1 - (@0/13000) + (@3*pow(@0,2)/pow(13000,2)) ), @2 )/ pow(@0,@1)']

	fitFunctions['altExpo5'] = [ TF1("altExpo5", "[0]*pow( ( 1 - (x/13000) + [3]*x*x/TMath::Power(13000,2) - [4]*x*x*x/TMath::Power(13000,3) ), [2] )/pow(x,[1])", 0, 2000 ), 
			[ 1000, -0.1, 50, 0, 0 ],
			'pow( (1 - (@0/13000) + (@3*pow(@0,2)/pow(13000,2) - (@4*pow(@0,3)/pow(13000,3)) ), @2 )/ pow(@0,@1)']

	fitFunctions['expoPoli3'] = [ TF1("expoPoli3", "exp([0]+[1]*x+[2]*x*x)", 0, 2000 ), 
			[8., -0.009, -9e-09],
			'exp( (@1*@0) + (@2*pow(@0,2)) )']

	fitFunctions['expoPoli4'] = [ TF1("expoPoli4", "exp([0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3))", 0, 2000 ), 
			[-67.84, -0.061, -4.28e-05, 7.25e-08],
			'exp( (@1*@0) + (@2*pow(@0,2)) + (@3*pow(@0,3)))']

	fitFunctions['expoPoli5'] = [ TF1("expoPoli5", "exp([0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4))", 0, 2000 ), 
			[ 1, -0.001, -0.0000001, 0.000000001, 0.0000000000001 ],
			'exp( (@1*@0) + (@2*pow(@0,2)) + (@3*pow(@0,3)+ (@4*pow(@0,4))))']

	#### http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_477_v10.pdf
	fitFunctions['atlas3'] = [ 
			TF1("atlas3", "[0]* TMath::Power((1-TMath::Power(x/13000.0,0.33)),[1]) / TMath::Power(x/13000.0,[2])",0,2000), 
			[ 10000000000, 40, -1 ], 
			'pow(1-pow(@0/13000,0.33),@1)/pow(@0/13000,@2)' ]

	fitFunctions['atlas4'] = [ 
			TF1("atlas4", "[0]* TMath::Power((1-TMath::Power(x/13000.0,0.33)),[1]) / TMath::Power(x/13000.0,[2]+[3]*(TMath::Power( TMath::Log(x/13000.), 2)))",0,2000), 
			[9.75e-13, 2.84, 12., -0.26],
			'pow(1-pow(@0/13000,0.33),@1)/pow(@0/13000, (@2+@3*pow( log(@0/13000.) ,2)))' ]

	#### https://arxiv.org/pdf/1110.5302.pdf
	fitFunctions['CDF3'] = [ 
			TF1("CDF3", "[0]* TMath::Power( (1-(x/13000.0) ),[2]) / TMath::Power(x,[1])",0,2000), 
			[3791209., 1.21, 87.54],
			'pow(1-@0/13000.0,@2)/pow(@0,@1)' ]

	fitFunctions['CDF4'] = [ 
			TF1("CDF4", "[0]* TMath::Power( (1-(x/13000.0) + [3]*TMath::Power(x/13000.0,2) ),[2]) / TMath::Power(x,[1])",0,2000), 
			#[ 1, -1, 100, 1 ], 
			[13., -1.31, 181., 2.01],
			'pow((1-@0/13000.0+@3*pow(@0/13000,2)),@2)/pow(@0,@1)' ]

	fitFunctions['Landau'] = [ 
			TF1("Landau", "[0]*TMath::Landau(-x,[1],[2])",0,2000), 
			[ ], #10000, 1000, 1000  ], 
			'TMath::Landau(@0,@1,@2)' ]

	fitFunctions['P1'] = [ TF1("P1", "[0] / (TMath::Power(x/13000.0,[1]))",0,2000), [ 0] ]
	if args.doubleGaus: fitFunctions['gaus'] = [ TF1("gaus", "gaus(0)+gaus(3)", 0, 2000), [ ] ]
	else: fitFunctions['gaus'] = [ TF1("gaus", "gaus", 0, 2000), [ ] ]
	fitFunctions['P4Gaus'] = [ TF1("P4Gaus", "[0]*pow(1-(x/13000.0),[1])/pow(x/13000.0,[2]+([3]*log(x/13000.)))+gaus(4)",0,2000), [] ]

	######## Uncertainties
	uncDict = {}
	##      mass   =  btag, 	JES, 	JER, 	PU, 	PDF,	UDD323 JES, 	JER, 	PU, 	PDF	312 fit 	323 fit
	uncDict[ 200 ] = [ 0.001, 	0.031, 	0.073, 	0.005,	0.175,		0.031, 	0.057, 	0.002, 	0.171,	[180, 225, 5],	[ 180, 225, 5] ]
	uncDict[ 220 ] = [ 0.001, 	0.028, 	0.028, 	0.005, 	0.172,		0.037, 	0.013, 	0.005, 	0.174,	[200, 250, 5],	[ 190, 250, 5] ]
	uncDict[ 240 ] = [ 0.001, 	0.019, 	0.007, 	0.011, 	0.17,		0.033,  0.001, 	0.009, 	0.178,	[210, 265, 5],	[ 200, 270, 5] ]
	uncDict[ 260 ] = [ 0.001, 	0.0, 	0.0, 	0.0, 	0.,		0.027, 	0.01, 	0.005, 	0.163,	[220, 280, 5],	[ 220, 290, 5] ]
	uncDict[ 280 ] = [ 0.001, 	0.0, 	0.0, 	0.0, 	0.,		0.022,  0.017, 	0.009, 	0.179,	[250, 310, 5],	[ 240, 305, 5] ]
	uncDict[ 300 ] = [ 0.001, 	0.016, 	0.005, 	0.002, 	0.18,		0.019, 	0.005, 	0.005, 	0.182,	[275, 325, 5],	[ 260, 320, 5] ]
	uncDict[ 350 ] = [ 0.002, 	0.021, 	0.001, 	0.002, 	0.195,		0.024, 	0.014, 	0.005, 	0.194,	[310, 380, 5],	[ 300, 380, 10] ]
	uncDict[ 400 ] = [ 0.002, 	0.018, 	0.014, 	0.001, 	0.208,		0.017, 	0.023, 	0.003, 	0.206,	[340, 460, 20],	[ 320, 460, 20] ]
	uncDict[ 450 ] = [ 0.002, 	0.023, 	0.016, 	0.003, 	0.214,		0.025, 	0.011, 	0.005, 	0.218,	[410, 480, 5],	[ 380, 500, 10] ]
	uncDict[ 500 ] = [ 0.002, 	0.025, 	0.028, 	0.001, 	0.24,		0.037, 	0.031, 	0.016, 	0.238,	[450, 540, 5],	[ 420, 540, 10] ]
	uncDict[ 550 ] = [ 0.002, 	0.023, 	0.015, 	0.007, 	0.253,		0.037, 	0.005, 	0.008, 	0.249,	[500, 590, 10],	[ 460, 590, 10] ]
	uncDict[ 600 ] = [ 0.003, 	0.033, 	0.022, 	0.003, 	0.266,		0.049, 	0.026, 	0.012, 	0.27,	[500, 700, 20],	[ 500, 700, 20] ]
	uncDict[ 650 ] = [ 0.003, 	0.022, 	0.006, 	0.004, 	0.279,		0.032, 	0.04, 	0.007, 	0.29,	[560, 720, 20],	[ 560, 700, 10] ]
	uncDict[ 700 ] = [ 0.003, 	0.033, 	0.041, 	0.006, 	0.302,		0.051, 	0.032, 	0.013, 	0.291,	[540, 840, 20],	[ 540, 820, 20] ]
	uncDict[ 750 ] = [ 0.003, 	0.04, 	0.003, 	0.008, 	0.301,		0.092, 	0.025, 	0.013, 	0.313,	[640, 840, 20],	[ 640, 820, 20] ]
	uncDict[ 800 ] = [ 0.003, 	0.045, 	0.049, 	0.005, 	0.335,		0.085, 	0.001, 	0.003, 	0.336,	[660, 900, 20],	[ 700, 880, 30 ] ]
	uncDict[ 850 ] = [ 0.004, 	0.033, 	0.006,	0.011, 	0.358,		0.073, 	0.018, 	0.005, 	0.334,	[720, 960, 30],	[ 740, 950, 30] ]
	uncDict[ 900 ] = [ 0.003, 	0.024, 	0.024, 	0.019, 	0.365,		0.05, 	0.017, 	0.007, 	0.332,	[750, 990, 30],	[ 750, 1020, 30] ]
	uncDict[ 950 ] = [ 0.003, 	0.032, 	0.001, 	0.006, 	0.387,		0.05, 	0.001, 	0.001, 	0.381,	[800, 1100, 30],	[ 800, 1100, 30] ]
	uncDict[ 1000 ] = [ 0.004, 	0.037, 	0.039, 	0.015, 	0.408,		0.044, 	0.023, 	0.047, 	0.398,	[850, 1150, 30],	[ 850, 1150, 30] ]
	uncDict[ 1100 ] = [ 0.004, 	0.065, 	0.028, 	0.008, 	0.446,		0.095, 	0.001, 	0.032, 	0.435,	[980, 1220, 30],	[ 960, 1280, 40] ]
	uncDict[ 1200 ] = [ 0.008, 	0.062, 	0.019, 	0.017, 	0.499,		0.089, 	0.035, 	0.028, 	0.488,	[1080, 1280, 40],	[ 1040, 1360, 40] ]
	uncDict[ 1300 ] = [ 0.008, 	0.062, 	0.019, 	0.017, 	0.499,		0.089, 	0.035, 	0.028, 	0.488,	[1100, 1400, 50],	[ 1100, 1460, 40] ]
	uncDict[ 1400 ] = [ 0.008, 	0.062, 	0.019, 	0.017, 	0.499,		0.089, 	0.035, 	0.028, 	0.488,	[1200, 1500, 50],	[ 1040, 1560, 40] ]
	uncDict[ 1500 ] = [ 0.008, 	0.062, 	0.019, 	0.017, 	0.499,		0.089, 	0.035, 	0.028, 	0.488,	[1300, 1600, 50],	[ 1040, 1660, 40] ]
	
	############ run process
	if 'Data' in args.process:
		CMS_lumi.extraText = ('Simulation ' if args.bkgAsData or args.pseudoExperiment else '')+"Preliminary"
		p = Process( target=FitterCombination, args=( ( bkgFile if args.bkgAsData else dataFile ), 
				bkgFile, 
				( signalFilename if args.final else ''), 
				hist, 
				scale, 
				( [ fitFunctions[args.func], fitFunctions['P4'], fitFunctions['P5'] ] if args.comparison else [fitFunctions[args.func]] ), 
				#[ fitFunctions[args.func], fitFunctions['CDF3'], fitFunctions['expoPoli3'], fitFunctions['altExpo3'], fitFunctions['atlas4'] ], 
				#160, 1700, 
				#20 
				350, ( 1477 if 'UDD323' in args.decay else 1639),   ##259
				[ len(massBins)-1, "resoBin", massBins ]
				))

	elif 'RPV' in args.process:
		process = 'RPVSt'+str(args.mass)+'tojj'
		p = Process( target=rootFitter, args=( signalFile, hist, scale, [fitFunctions['gaus']], args.mass-50, args.mass+50, rebinX, True, False ) ) 

	elif 'Limit' in args.process:
		p = Process( target=createCards, args=( ( bkgFile if args.bkgAsData else dataFile ),
				bkgFile, 
				signalFilename,
				listMass,
				#(listMass if len(listMass)==1 else range(400, 1600, 100)), 
				hist, 
				scale, 
				[fitFunctions[args.func]], 
				350, maxFit, 1 ) ) 

	elif 'Bias' in args.process:
		p = Process( target=createCards, args=( ( bkgFile if args.bkgAsData else dataFile ), 
			bkgFile, 
			signalFilename,
			(listMass if len(listMass)==1 else range(400, 1600, 100)), 
			hist, 
			scale, 
			[ fitFunctions[args.func], fitFunctions['CDF3'], fitFunctions['expoPoli4'], fitFunctions['altExpo3'], fitFunctions['atlas4'] ], 
			350, maxFit, 1 ) )

	elif 'Fisher' in args.process:
		p = Process( target=FisherTest, args=( ( bkgFile if args.bkgAsData else dataFile ), 
			hist, 
			[ fitFunctions['P3'], fitFunctions['P4'], fitFunctions['P5'] ], 
			350, ( 1477 if 'UDD323' in args.decay else 1639),   ##259
			[ len(massBins)-1, "resoBin", massBins ]
			) )

	elif 'biasDraw' in args.process: 
		p = Process( target=drawBiasTest, args=( (listMass if len(listMass)==1 else  range(400, 1500, 100)),
			'/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/' ) )

	elif 'tmp' in args.process: tmpPlot() 

	p.start()
	p.join()
