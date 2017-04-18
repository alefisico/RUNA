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
TVirtualFitter.SetMaxIterations(50000000)		######### Trick to increase number of iterations

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



def createPseudoExperiment():
	"""docstring for createPseudoExperiment"""
	h1 = inFileBkg.Get(hist)
	binSize = h1.GetBinWidth(1)
	P4PSE = TF1("P4PSE", "[0]*pow(1-x/13000.0,[1])/pow(x/13000.0,[2]+[3]*log(x/13000.))", minX, maxX);
	P4PSE.SetParameter(0,bkgParameters[0])				
	P4PSE.SetParameter(1,bkgParameters[1])
	P4PSE.SetParameter(2,bkgParameters[2])
	P4PSE.SetParameter(3,bkgParameters[3])
	randomNumEventsQCD = random.randint( bkgParameters[4]-round(sqrt(bkgParameters[4])), bkgParameters[4]+round(sqrt(bkgParameters[4])) )
	print "randomNumberOf QCD events", randomNumEventsQCD
	hMainPSE = TH1D("hbkgPSE", "hbkgPSE", int( (maxX-minX)/binSize ) , minX, maxX)
	hMainPSE.FillRandom( "P4PSE", int(randomNumEventsQCD) )
	c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	c1.SetLogy()
	gStyle.SetOptFit()
	gStyle.SetStatY(0.94)
	gStyle.SetStatX(0.9)
	gStyle.SetStatW(0.15)
	gStyle.SetStatH(0.15) 
	hMainPSE.GetXaxis().SetTitle( histYaxis )
	hMainPSE.GetYaxis().SetTitle("Events / 10 GeV" ) # dN/dM_{bbjj} [GeV^{-1}]")
	hMainPSE.GetYaxis().SetTitleOffset(1.2);
	hMainPSE.SetTitle("QCD PseudoExperiments")
	hMainPSE.Sumw2()
	hMainPSE.Draw()
	c1.SaveAs(outputDir+hist+"_PseudoExperiment_Fit_ResolvedAnalysis_"+args.version+"."+args.extension)
	del c1



def rootFitter( inFile, hist, scale, fitFunctions, minX, maxX, rebinX, plot, log=True ):
	"""Simple rootFitter"""

	tmpBinContent = []
	tmpBinError = []
	minBinX = int(minX/(rebinX if args.miniTree else rebinX*10 ))+1
	maxBinX = int(maxX/(rebinX if args.miniTree else rebinX*10 ))+1

	#### Initial histo
	rawHisto = inFile.Get( hist )
	rawHisto.SetBinErrorOption(TH1.kPoisson)
	rawHisto.Scale( scale )
	rawHisto.Rebin( rebinX )

	#### extracting bin contents and dividing by bin size
	print '|----> Creating histograms from bin', minBinX, '(', rawHisto.GetBinLowEdge(minBinX),') to ', maxBinX, '(', rawHisto.GetBinLowEdge(maxBinX),')'
	for ibin in range( minBinX, maxBinX ):
		tmpBinContent.append( rawHisto.GetBinContent(ibin) / rebinX )
		tmpBinError.append( ( rawHisto.GetBinError(ibin) if ( rawHisto.GetBinContent(ibin) > 0 ) else 1.8 ) / rebinX )
		#print rawHisto.GetBinLowEdge( ibin ), rawHisto.GetBinContent(ibin) , rebinX, rawHisto.GetBinError(ibin)

	binContents = np.array(tmpBinContent)
	binError = np.array(tmpBinError)
	print '|----> bins :', binContents
	print '|----> bins error :', binError
	
	#### creating histo for fitting
	numBins = maxBinX - minBinX
	finalHisto = TH1D("finalHisto", "finalHisto", numBins, minX, maxX)
	finalHisto.Sumw2()
	for ibin in range( 0, numBins ):
		finalHisto.SetBinContent( ibin+1, binContents[ibin] )
		finalHisto.SetBinError( ibin+1, binError[ibin] )

	fitParameters = OrderedDict()
	fitParErrors = OrderedDict()
	for fitFunc in fitFunctions:

		#### giving initial values to fit
		if( len(fitFunc[1])>0 ):
			for k in range( len(fitFunc[1]) ): fitFunc[0].SetParameter(k, fitFunc[1][k])

		### Fitting
		tmpFit = 1000000000000000000
		keepFitting = True
		while keepFitting:
			#finalHisto.Fit(fitFunc[0],"MIRS","",minX,maxX)
			finalHisto.Fit( fitFunc[0],"ELLSR","",minX,maxX)
			#### this is just a trick to keep fitting...
			tmpFitFunc = finalHisto.GetFunction( fitFunc[0].GetName() )
			chi2Ndf = tmpFitFunc.GetChisquare() / tmpFitFunc.GetNDF()
			tmpStatus = ( tmpFit - fitFunc[0].GetParameter(1) ) / fitFunc[0].GetParameter(1)
			if ((tmpStatus < 0.001) or (chi2Ndf < 1.5 ) ): 
				keepFitting = False
				print "|----> Final fit for ", fitFunc[0].GetName(), 'with chi2/ndf', tmpFitFunc.GetChisquare() , tmpFitFunc.GetNDF()

			else: 
				tmpFit = fitFunc[0].GetParameter(1)
				print "|----> Keep fitting...."

		numParam = fitFunc[0].GetNpar()
		fitParameters[ fitFunc[0].GetName() ] = [ fitFunc[0].GetParameter(k) for k in range( numParam ) ]
		fitParErrors[ fitFunc[0].GetName() ] =  [ fitFunc[0].GetParError(k) for k in range( numParam ) ]
		print "|----> Fitter parameters for", fitFunc[0].GetName(), fitParameters[ fitFunc[0].GetName() ], fitParErrors[ fitFunc[0].GetName() ]


		######### Plotting Histograms
		if plot:
			c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
			if log: c1.SetLogy()
			finalHisto.GetXaxis().SetTitle( histYaxis )
			finalHisto.GetYaxis().SetTitle('dN/dm_{av} / '+str(round(rebinX))+' GeV' ) 
			finalHisto.GetYaxis().SetTitleOffset(0.9);
			finalHisto.GetXaxis().SetRangeUser( minX-50, maxX+50 )
			finalHisto.SetTitle("")
			finalHisto.Draw()
			c1.SaveAs(outputDir+hist.replace('ResolvedAnalysisPlots/','')+"_"+args.process+"_"+fitFunc[0].GetName()+"Fit_ResolvedAnalysis_"+args.version+"."+args.extension)
			del c1
	
	return [ fitParameters, fitParErrors, binContents, binError ]




def histoFunctionFit( nameHisto, initFitFunction, parameters, parErrors, massBin, massBinErr, minX, maxX ):
	"""docstring for histoFunctionFit"""

	fitFunction = initFitFunction.Clone()
	for i in range(0, initFitFunction.GetNpar() ): 
		fitFunction.SetParameter( i, parameters[fitFunction.GetName()][i] )
		fitFunction.SetParError( i, parErrors[fitFunction.GetName()][i] )

	histoFit = TH1D( nameHisto, nameHisto, len(massBin) , minX, maxX)
	histoFit.Sumw2(False)

	for ibin in range( 0, len(massBin)):
		histoFit.SetBinContent( ibin, massBin[ibin] )
		#histoFit.SetBinError( ibin, ( massBinErr[ibin] if ( massBin[ibin] > 0 ) else 1.8 ) )
		histoFit.SetBinError( ibin, massBinErr[ibin] )

	histoFit.SetBinErrorOption(TH1.kPoisson)
			
	histoFit.Fit( fitFunction, "ELLSR", "", minX, maxX )

	return histoFit, fitFunction



def residualAndPulls(dataPoints, dataErrPoints, function, histo, minX, maxX ):
	"""docstring for residualAndPulls"""

	hPull = TH1D("hpull_"+function.GetName(), "hpull_"+function.GetName(), len(dataPoints) , minX, maxX)
	hPull.Sumw2()
	hResidual = TH1D("hresidual_"+function.GetName(), "hresidual_"+function.GetName(), len(dataPoints) , minX, maxX)
	hResidual.Sumw2()

	######## Calculating Pull and Residual
	chi2 = 0 
	nof = 0
	for ibin in range(0, len(dataPoints) ):
	
		binCont = dataPoints[ibin]
		binErr = dataErrPoints[ibin]
		valIntegral = function.Eval( histo.GetBinCenter(ibin) ) 
		diff = (binCont - valIntegral)/ valIntegral
		#errDiff = diff * TMath.Sqrt( TMath.Power( P4Gaus.GetParError(0) / P4Gaus.GetParameter(0),2 ) + TMath.Power( P4Gaus.GetParError(1)/ P4Gaus.GetParameter(1), 2 )  + TMath.Power( P4Gaus.GetParError(2)/ P4Gaus.GetParameter(2), 2 )  + TMath.Power( P4Gaus.GetParError(3)/ P4Gaus.GetParameter(3), 2 ) )
		#errDiff = diff * TMath.Sqrt( TMath.Power( function.GetParError(0) / function.GetParameter(0),2 ) + TMath.Power( function.GetParError(1)/ function.GetParameter(1), 2 )  + TMath.Power( function.GetParError(2)/ function.GetParameter(2), 2 )  + TMath.Power( function.GetParError(3)/ function.GetParameter(3), 2 ) )
		#print binCont, binErr, valIntegral 

		if (binCont != 0):
			pull = (binCont - valIntegral)/ binErr
			chi2 += TMath.Power(pull,2)
			nof += 1
			
			hPull.SetBinContent(ibin, pull)
			hPull.SetBinError(ibin, 1.0)
	
			hResidual.SetBinContent(ibin, diff)
			hResidual.SetBinError(ibin, binErr/valIntegral )

	NDoF = nof - function.GetNpar() - 1
	print '|----> ############### chi2 and nof: ', chi2, nof

	return hPull, hResidual, chi2, NDoF



def FitterCombination( inFileData, inFileBkg, inFileSignal, hist, scale, bkgFunction, minX, maxX, rebinX, runData ):
	"""docstring for FitterCombination"""

	### Fit QCD MC
	print "|----> Fitting MC QCD"
	BkgParameters = rootFitter( inFileBkg, 
			hist+('QCD'+args.qcd+'All' if args.miniTree else ''), 
			scale*0.85, #( 1 if args.miniTree else 0.75 ), 
			bkgFunction, 
			minX, 
			maxX, 
			rebinX, 
			True )

	bkgParameters = BkgParameters[0]
	bkgParErrors = BkgParameters[1]
	bkgpoints = BkgParameters[2]
	bkgpointsErr = BkgParameters[3]

#	SigParameters = rootFitter( inFileSignal, hist, scale, gaus, [], args.mass-50, args.mass+50, rebinX, False )
	
	legend=TLegend(0.18,0.15,0.50,0.35)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)

	if runData:
		print "|----> Fitting Data"
		DataParameters = rootFitter( inFileData, 
				hist+('JetHT_Run2016' if args.miniTree else ''), 
				1, 
				bkgFunction, 
				minX, 
				maxX, 
				rebinX, 
				False )
		#DataParameters = rootFitter( inFileData, hist+('QCDHerwigAll' if args.miniTree else ''), scale, bkgFunction, minX, maxX, rebinX, False )
		dataParameters = DataParameters[0]
		dataParErrors = DataParameters[1]
		points = DataParameters[2]
		pointsErr = DataParameters[3]
		hMain, mainP4 = histoFunctionFit( 'Data', bkgFunction[0][0], dataParameters, dataParErrors, points, pointsErr, minX, maxX )
		legend.AddEntry( hMain, 'Data', 'ep' )
		legend.AddEntry( mainP4, 'Fit to data', 'l' )
		hMCQCD, qcdMCP4 = histoFunctionFit( 'QCD'+args.qcd+'All', bkgFunction[0][0], bkgParameters, bkgParErrors, bkgpoints, bkgpointsErr, minX, maxX )
		#legend.AddEntry( qcdMCP4, 'Fit to MC QCD pythia', 'l' )

		if isinstance( inFileSignal, TFile):
			Bkg2Parameters = rootFitter( inFileSignal, 
							hist+('QCDHerwigAll' if args.miniTree else ''), 
							scale, 
							bkgFunction, 
							minX, 
							maxX, 
							rebinX, 
							True )
			bkg2Parameters = Bkg2Parameters[0]
			bkg2ParErrors = Bkg2Parameters[1]
			bkg2points = Bkg2Parameters[2]
			bkg2pointsErr = Bkg2Parameters[3]
			hBkg2, qcdHTMCP4 = histoFunctionFit( 'QCDHerwig', bkgFunction[0][0], bkg2Parameters, bkg2ParErrors, bkg2points, bkg2pointsErr, minX, maxX )
			#hBkg2, qcdHTMCP4 = histoFunctionFit( 'QCDHT', bkg2Parameters, bkg2points, bkg2pointsErr, minX, maxX )
			legend.AddEntry( qcdHTMCP4, 'Fit to MC QCD powheg', 'l' )
			#legend.AddEntry( qcdHTMCP4, 'Fit to MC QCD madgraph+pythia', 'l' )

	else:
		hMain, mainP4 = histoFunctionFit( 'QCD'+args.qcd+'All', bkgFunction, bkgParameters, bkgParErrors, bkgpoints, bkgpointsErr, minX, maxX )
		points = bkgpoints
		pointsErr = bkgpointsErr
	
	print '|----> DATA Plotted:', points
	print '|----> DATA Err:', pointsErr

	hPull, hResidual, chi2, NDF = residualAndPulls(points, pointsErr, mainP4, hMain, minX, maxX )

	######### Plotting Histograms
	maxXPlot = maxX+500
	tdrStyle.SetPadRightMargin(0.05)
  	tdrStyle.SetPadLeftMargin(0.15)
	gStyle.SetOptFit()
	gStyle.SetStatY(0.91)
	gStyle.SetStatX(0.95)
	gStyle.SetStatW(0.15)
	gStyle.SetStatH(0.30) 
	#c3 = TCanvas('c1', 'c1',  10, 10, 750, 1000 )
	#pad1 = TPad("pad1", "Fit",0,0.40,1.00,1.00,-1)
	#pad2 = TPad("pad2", "Pull",0,0.25,1.00,0.475,-1);
	#pad3 = TPad("pad3", "Residual",0,0,1.00,0.277,-1);
	c3 = TCanvas('c1', 'c1',  10, 10, 1250, 500 )
	pad1 = TPad("pad1", "Fit",0,0.00,0.50,1.00,-1)
	pad2 = TPad("pad2", "Pull",0.50,0.45,1.00,0.95,-1);
	pad3 = TPad("pad3", "Residual",0.50,0,1.00,0.507,-1);
	pad1.Draw()
	pad2.Draw()
	pad3.Draw()

	pad1.cd()
	pad1.SetLogy()
	hMain.SetMarkerStyle(8)
	hMain.GetYaxis().SetTitle("dN/dm_{av} / "+ str(hMain.GetBinWidth(1)) +" GeV" ) # dN/dM_{bbjj} [GeV^{-1}]")
	hMain.GetXaxis().SetTitle( histYaxis )
	hMain.GetYaxis().SetTitleOffset(1.2);
	hMain.SetTitle("")
	hMain.SetMarkerStyle(20)
	#hMain.SetMaximum( 1.5 * hMain.GetMaximum() )
	hMain.Draw('E0')
	hMain.GetXaxis().SetRangeUser( minX, maxXPlot  )
	mainP4.SetLineColor(kBlack)
	mainP4.Draw("same")
	#gaus2.SetLineColor(kRed-4)
	#gaus2.Draw("same")
	mainP4.SetLineColor(kBlue-4)
	mainP4.Draw("same")
	if runData:
		#qcdMCP4.SetLineColor( kMagenta )
		#qcdMCP4.Draw("same")	
		#if isinstance( inFileSignal, TFile):
		#	qcdHTMCP4.SetLineColor( kViolet )
		#	qcdHTMCP4.Draw("same")	
		legend.Draw("same")
	CMS_lumi.relPosX = 0.13
	CMS_lumi.cmsTextSize = 0.60
	CMS_lumi.lumiTextSize = 0.50
	CMS_lumi.CMS_lumi(pad1, 4, 0)
	#labels( hist, '', '', 0.20, 0.45, 'left' )


	pad2.cd()
	pad2.SetGrid()
	gStyle.SetOptStat(0)
	hPull.GetYaxis().SetTitle("#frac{(Data - Fit)}{#sigma_{Data}}")
	hPull.GetYaxis().SetLabelSize(0.08)
	hPull.GetYaxis().SetTitleSize(0.09)
	hPull.GetYaxis().SetTitleOffset(0.70)
	hPull.GetYaxis().CenterTitle()
	hPull.SetMarkerStyle(7)
	#hPull.SetMaximum(3)
	hPull.GetXaxis().SetRangeUser( minX, maxXPlot )
	hPull.Sumw2()
	hPull.Draw("e")
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
	hResidual.SetMaximum(0.59)
	hResidual.SetMinimum(-0.59)
	hResidual.GetXaxis().SetRangeUser( minX, maxXPlot )
	#hResidual.Sumw2()
	hResidual.Draw("e")
	line.Draw("same")
	c3.SaveAs("Plots/"+hist.replace('ResolvedAnalysisPlots/','')+"_"+args.process+"_"+args.version+"Fit"+bkgFunction[0][0].GetName()+"_ResolvedAnalysis_"+args.version+"."+args.extension)
	del c3



def createCards( dataFile, bkgFile, inFileSignal, listMass, hist, scale, bkgFunctions, minX, maxX, rebinX ):
	"""function to run Roofit and save workspace for RooStats"""
	
	warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='.*class stack<RooAbsArg\*,deque<RooAbsArg\*> >' )

	########## Fitting background and taking parameters
	print '|----> Background'
	'''
	bkgFuncParameters = rootFitter( bkgFile, hist+('QCD'+args.qcd+'All' if args.miniTree else ''), scale, bkgFunctions, minX, maxX, rebinX, True)
	'''
	bkgFuncParameters = rootFitter( dataFile, 
			hist+('JetHT_Run2016' if args.miniTree else ''), 
			1, 
			bkgFunctions, 
			minX, 
			maxX, 
			rebinX, 
			True )

	########## Storing parameters in RooRealVar
	dictPar = OrderedDict()
	for func in bkgFuncParameters[0]:
		for par in range( 1, len( bkgFuncParameters[0][func] ) ):
			tmpNamePar = func+'p'+str(par)
			dictPar[ tmpNamePar ] = RooRealVar( tmpNamePar, 
					tmpNamePar, 
					bkgFuncParameters[0][func][par], 
					-abs(bkgFuncParameters[0][func][par])*100, 
					abs(bkgFuncParameters[0][func][par])*100)

	############ Creating data plot from input in bkgFuncParameters
	hData = TH1D("hData", "hData", len(bkgFuncParameters[2]) , minX, maxX)
	hData.Sumw2()
	for ibin in range(0, len( bkgFuncParameters[2] ) ):
		hData.SetBinContent( ibin+1, bkgFuncParameters[2][ibin] )
		hData.SetBinError( ibin+1, bkgFuncParameters[3][ibin] )


	########### Starting signal fits
	listAcceptance = []
	listAccepError = []
	listMass.sort()
	for imass in listMass:

		######## Fitting signal and extracting parameters
		print '|----> Signal'
		SignalParameters = rootFitter( TFile.Open( inFileSignal.replace( str(args.mass), str(imass) ) ), 
				hist+('RPVStopStopToJets_'+args.decay+'_M-'+str(imass) if args.miniTree else ''), 
				scale, 
				[ fitFunctions['gaus'] ], 
				imass-(100 if (imass < 800 ) else 200 ), 
				imass+(100 if (imass < 800 ) else 200 ), 
				rebinX, 
				True )

		print SignalParameters[0]['gaus'][2]
		massWindow = int(SignalParameters[0]['gaus'][2])*2		### size of search
		mass = RooRealVar( 'mass', 'mass', imass-massWindow, imass+massWindow )
		mypdfs = RooArgList()

		### creating RooStats objects for data
		rooDataHist = RooDataHist('data_obs','data_obs',RooArgList(mass),hData)
		rooDataHist2 = RooDataHist('data_obs','data_obs',RooArgList(mass),hData)
		rooDataHist.Print()

		### creating RooStats objects for signal
		sigMean = RooRealVar('sigMean', 'sigMean', imass )
		sigSigma = RooRealVar('sigSigma', 'sigSigma', SignalParameters[0]['gaus'][2] )
		signal = RooGaussian( 'signal', 'signal', mass, sigMean, sigSigma )
		signal.Print()

		signalGaus = TF1( "RPVStop"+str(imass), "gaus", imass-massWindow, imass+massWindow )
		signalGaus.SetParameter( 0, SignalParameters[0]['gaus'][0] ) 
		signalGaus.SetParameter( 1, imass ) 
		signalGaus.SetParameter( 2, SignalParameters[0]['gaus'][2] ) 
		sigAcc = round( signalGaus.Integral( imass-massWindow, imass+massWindow ), 2 )
		print '|----> Signal events:', sigAcc

		eventsInAccep = sigAcc / ( scaleFactor( 'RPVStopStopToJets_'+args.decay+'_M-'+str( imass ) )* scale ) ### because sigAcc is scaled
		totalEvents = search( dictEvents, 'RPVStopStopToJets_'+args.decay+'_M-'+str( imass ) )[0]
		failedEvents = totalEvents - eventsInAccep
		accXeff = eventsInAccep / totalEvents 
		accXeffErr = TMath.Sqrt( (1/failedEvents) + (1/eventsInAccep) ) * failedEvents * eventsInAccep / TMath.Power( ( totalEvents ), 2 )
		listAcceptance.append( accXeff ) 
		listAccepError.append( accXeffErr )


		'''
		if args.jesUnc:
			print ' |---> Adding JES'
			hSigSyst['JESUp'] = hSignal.Clone()
			hSigSyst['JESDown'] = hSignal.Clone()

		if args.jerUnc:
			print ' |---> Adding JER'
			hSigSyst['JERUp'] = hSignal.Clone()
			hSigSyst['JERDown'] = hSignal.Clone()
		'''
		#hBkg = inFileBkg.Get(hist)
		#hData = hBkg.Clone()
		#hData.Add( hSignal )


		### creating RooStats objects for background
		dictRooPdf = {}
		for bkgFunc in bkgFunctions:

			tmpBkgFunct = bkgFunc[0].Clone()
			for p in range( 0, tmpBkgFunct.GetNpar()): 
				tmpBkgFunct.SetParameter( p, bkgFuncParameters[0][ tmpBkgFunct.GetName() ][p] )
			bkgAcc = round( tmpBkgFunct.Integral( imass-massWindow, imass+massWindow ), 2 )
			print '|----> Background events in mass (', imass, '):', bkgAcc
		
			tmpList = RooArgList() 
			tmpList.add( mass ) 
			for par in dictPar: 
				if tmpBkgFunct.GetName() in par: tmpList.add( dictPar[par] )
			
			if 'Bias' in args.process:
				print '|----> Adding ', bkgFunc[0].GetName()+'_pdf', ' to bias test'
				dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ] = RooGenericPdf( bkgFunc[0].GetName()+'_pdf' , bkgFunc[2], tmpList ) 
				dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ].fitTo( rooDataHist )
				mypdfs.add( dictRooPdf[ bkgFunc[0].GetName()+'_pdf' ] )

			else:
				background = RooGenericPdf('background', bkgFunc[2], tmpList ) 
				background.Print()

		#### S+B model this is to fit in RooFit
		#signal_norm = RooRealVar( 'signal_norm', 'signal_norm', sigAcc, 0, 20.*TMath.Sqrt(sigAcc))
		#signal_norm = RooRealVar( 'signal_norm', 'signal_norm', sigAcc, 0, 10000000000 ) #+(20.*TMath.Sqrt(sigAcc)))
		#signal_norm.setConstant()
		#signal_norm.Print()
		#background_norm = RooRealVar('background_norm','background_norm', bkgAcc, 0., 20.*TMath.Sqrt(bkgAcc))
		#background_norm.Print()
		#model = RooAddPdf("model","s+b",RooArgList(background,signal),RooArgList(background_norm,signal_norm))
		#res = model.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Strategy(0), RooFit.SumW2Error(kFALSE))
		#res.Print()
		###################################################

		##### creating workspace
		outputRootFile = '/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/Rootfiles/workspace_RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'_Resolved_'+args.cut+'_'+('BiasTest_' if 'Bias' in args.process else '')+args.version+'.root' 
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
			myWS.Print()
			myWS.writeToFile(outputRootFile, True)

		print '|----> Workspace created:', outputRootFile

		##### write a datacard
		datacard = open('/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/Datacards/datacard_RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'_Resolved_'+args.cut+'_'+('BiasTest_' if 'Bias' in args.process else '')+args.version+'.txt','w')
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
		datacard.write('------------------------------\n')
		if 'Bias' in args.process: datacard.write('pdf_index\tdiscrete\n')
		for bF in bkgFuncParameters[0]:
			for par in range( 1, len(bkgFuncParameters[0][bF]) ):
				datacard.write( bF+'p'+str(par)+'\tparam\t'+str(bkgFuncParameters[0][bF][par])+'\t'+str(bkgFuncParameters[1][bF][par])+'\n' )
				#datacard.write( bF+'p'+str(par)+'\tflatParam\n')
		datacard.write("lumi\tlnN\t1.062\t-\n")
		datacard.write("trigger\tlnN\t1.05\t-\n")
		datacard.write("sigSigma\tparam\t"+str(SignalParameters[0]['gaus'][1])+'\t'+str(SignalParameters[1]['gaus'][1])+"\n")
		datacard.write("sigMean\tparam\t"+str(SignalParameters[0]['gaus'][2])+'\t'+str(SignalParameters[1]['gaus'][2])+"\n")
		#datacard.write("# SigNormPDF   lnN    1.0300       - \n")
		datacard.write("SigNormJES\tlnN\t1.02\t-\n")
		datacard.write("SigNormJER\tlnN\t1.10\t-\n")
		#datacard.write("# SigNormPU   lnN    1.0300        - \n")
		#datacard.write("# SigNormISR   lnN    1.1000        - \n")
		#datacard.write("# SigNormBtag   lnN    1.0000       - -\n")
		datacard.close()
		print '|----> Datacard created:', datacard
	
	if 'Limit' in args.process:

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
		#accXeffGraph.GetYaxis().SetRangeUser( 0.0001, 0.1  )

		#legend.Draw()
		CMS_lumi.extraText = "Simulation Preliminary"
		CMS_lumi.lumi_13TeV = ''
		CMS_lumi.relPosX = 0.12
		CMS_lumi.CMS_lumi(canAccep, 4, 0)
		canAccep.SaveAs( 'Plots/signalAcceptance_massAve_'+args.cut+'_'+args.decay+'_ResolvedAnalysis_'+args.version+'.'+args.extension )
		del canAccep



def doftest( RSS1, RSS2, NDF1, NDF2, nPar1, nPar2, nBinsFit):
	"""docstring for doftest"""

	Fvalue = ( ( RSS1 - RSS2 ) / ( nPar2 - nPar1 ) ) / ( RSS2 / ( nBinsFit - nPar2 ) )

	Fdist = TF1("Fdistr","TMath::Sqrt( (TMath::Power([0]*x,[0]) * TMath::Power([1],[1])) / (TMath::Power([0]*x + [1],[0]+[1])) / (x*TMath::Beta([0]/2,[1]/2)) )",0,100)
	Fdist.SetParameter( 0, (nPar2-nPar1) )
	Fdist.SetParameter( 1, (nBinsFit-nPar2) )
	CL = 1 - Fdistr.Integral( 0.00000001, Fvalue )
	altCL =  1. - TMath.FDistI( Fvalue, nPar2-nPar1, nBinsFit-nPar2 )

	return [ RSS1, RSS2, Fvalue, CL, altCL, Fdist ]

	

def FisherTest( dataFile, hist, bkgFunctions, minX, maxX, rebinX ):
	"""docstring for FisherTest"""

	fitParameters = rootFitter( dataFile, 
			hist+('JetHT_Run2016' if args.miniTree else ''), 
			1, 
			bkgFunctions,
			minX, 
			maxX, 
			rebinX, 
			False )

	dictDataAndFunc = OrderedDict()
	dictPullResChi2NDF = OrderedDict()
	for func in bkgFunctions:
		dictDataAndFunc[ func[0].GetName() ] = histoFunctionFit( 'Data_'+func[0].GetName(),
									func[0], 
									fitParameters[0], 
									fitParameters[1], 
									fitParameters[2],
									fitParameters[3],
									minX, maxX ) 
		
		dictPullResChi2NDF[ func[0].GetName() ] = residualAndPulls( fitParameters[2], 
									fitParameters[3], 
									dictDataAndFunc[ func[0].GetName() ][1],  
									dictDataAndFunc[ func[0].GetName() ][0], 
									minX, maxX )

	dictFtest = OrderedDict()
	for key1, key2 in combinations(dictPullResChi2NDF.keys(), r = 2):
		dictFtest[ key1+'_'+key2 ] = doftest( dictPullResChi2NDF[ key1 ][2], dictPullResChi2NDF[ key2 ][2], dictPullResChi2NDF[ key1 ][3], dictPullResChi2NDF[ key2 ][3], dictDataAndFunc[ key1 ][1].GetNpar(),  dictDataAndFunc[ key2 ][1].GetNpar(), len(fitParameters[2]) )
		print '|----> Ftest for ', key1, 'and', key2, ':', dictFtest[ key1+'_'+key2 ], 'Residual1, Residual2, Fvalue, CL, altCL, Fdist'
		

def drawBiasTest( listMasses, folderRootfiles ):
	"""docstring for plotLimits"""
	
	dictTest = {}
	dictTestHisto = {}
	listMeans = []
	listMeans1 = []
	listMeans2 = []
	listMeans3 = []
	listMeansErr = []
	listMeansErr1 = []
	listMeansErr2 = []
	listMeansErr3 = []
	for mass in listMasses:
		
		for t in [ 0, 1, 2, 3]:
			dictTest[ mass+t ] = TChain( 'tree_fit_sb' )
			dictTest[ mass+t ].Add( folderRootfiles+'mlfit_RPVStopStopToJets_'+args.decay+'_M-'+str(mass)+'_Resolved_'+args.cut+'_BiasTest_signal'+str(args.signalInj)+'_'+args.version+'_Index0ToIndex'+str(t)+'.root') 
			dictTestHisto[ mass+t ] = TH1F( 'h'+str(mass)+str(t),  'h'+str(mass)+str(t), 32, (-4 if (args.signalInj==1) else -1), (4 if (args.signalInj==1) else 7))
			dictTest[ mass+t ].Draw( "(mu-"+str(args.signalInj)+")/muErr>>h"+str(mass)+str(t))
			
			dictTestHisto[ mass+t ].SetMarkerSize(2)
			dictTestHisto[ mass+t ].SetLineWidth(2)
			dictTestHisto[ mass+t ].SetStats(True)
			dictTestHisto[ mass+t ].Fit("gaus")

		c = TCanvas("c"+str(mass), "c"+str(mass), 800, 600)
		c.cd()
		gStyle.SetOptFit(1)
		dictTestHisto[ mass ].SetLineColor(kBlack)
		dictTestHisto[ mass ].GetFunction("gaus").SetLineColor(1)
		listMeans.append( dictTestHisto[ mass ].GetFunction("gaus").GetParameter( 1 ) )
		listMeansErr.append( dictTestHisto[ mass ].GetFunction("gaus").GetParError( 1 ) )

		dictTestHisto[ mass+1 ].SetLineColor(kBlue)
		dictTestHisto[ mass+1 ].GetFunction("gaus").SetLineColor(kBlue)
		listMeans1.append( dictTestHisto[ mass+1 ].GetFunction("gaus").GetParameter( 1 ) )
		listMeansErr1.append( dictTestHisto[ mass+1 ].GetFunction("gaus").GetParError( 1 ) )

		dictTestHisto[ mass+2 ].SetLineColor(kRed)
		dictTestHisto[ mass+2 ].GetFunction("gaus").SetLineColor(kRed)
		listMeans2.append( dictTestHisto[ mass+2 ].GetFunction("gaus").GetParameter( 1 ) )
		listMeansErr2.append( dictTestHisto[ mass+2 ].GetFunction("gaus").GetParError( 1 ) )

		dictTestHisto[ mass+3 ].SetLineColor(8)
		dictTestHisto[ mass+3 ].GetFunction("gaus").SetLineColor(8)
		listMeans3.append( dictTestHisto[ mass+3 ].GetFunction("gaus").GetParameter( 1 ) )
		listMeansErr3.append( dictTestHisto[ mass+3 ].GetFunction("gaus").GetParError( 1 ) )

		dictTestHisto[ mass ].Draw("es")
		dictTestHisto[ mass+1 ].Draw("e sames")
		dictTestHisto[ mass+2 ].Draw("e sames")
		dictTestHisto[ mass+3 ].Draw("e sames")

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
		st1.SetTextColor(kBlue)
		st2 = dictTestHisto[ mass+2 ].GetListOfFunctions().FindObject("stats")
		st2.SetX1NDC(.75)
		st2.SetX2NDC(.95)
		st2.SetY1NDC(.35)
		st2.SetY2NDC(.55)
		st2.SetTextColor(kRed)
		st3 = dictTestHisto[ mass+3 ].GetListOfFunctions().FindObject("stats")
		st3.SetX1NDC(.75)
		st3.SetX2NDC(.95)
		st3.SetY1NDC(.15)
		st3.SetY2NDC(.35)
		st3.SetTextColor(8)
		c.Modified()
		c.SaveAs( 'Plots/BiasTest_'+args.decay+'_M-'+str(mass)+'_signal'+str(args.signalInj)+'_ResolvedAnalysis_'+args.version+'.'+args.extension )



	masses = array( 'd', listMasses )	
	massesErr = array( 'd', [0]*len(listMasses) )	
	graphTest = TGraphErrors( len( masses), masses, array( 'd', listMeans ), massesErr, array( 'd', listMeansErr ) )
	graphTest1 = TGraphErrors( len( masses), masses, array( 'd', listMeans1 ), massesErr, array( 'd', listMeansErr1 ) )
	graphTest2 = TGraphErrors( len( masses), masses, array( 'd', listMeans2 ), massesErr, array( 'd', listMeansErr2 ) )
	graphTest3 = TGraphErrors( len( masses), masses, array( 'd', listMeans3 ), massesErr, array( 'd', listMeansErr3 ) )

	c1 = TCanvas("c1"+str(mass), "c1"+str(mass), 800, 600)
	c1.cd()

	legend = TLegend(.70,.70,.90,.88)
	legend.SetTextSize(0.025)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
	#legend.SetHeader('95% CL upper limits')

	graphTest.GetXaxis().SetTitle("Resonance mass [GeV]")
	graphTest.GetYaxis().SetTitle("Means")
	graphTest.GetYaxis().SetTitleOffset(1.1)
	graphTest.GetYaxis().SetRangeUser(-0.4,0.4)
	graphTest.SetMarkerStyle(21)
	graphTest.SetMarkerColor(1)
	graphTest1.SetMarkerStyle(21)
	graphTest1.SetMarkerColor(kBlue)
	graphTest2.SetMarkerStyle(21)
	graphTest2.SetMarkerColor(kRed)
	graphTest3.SetMarkerStyle(21)
	graphTest3.SetMarkerColor(8)
	#graph_exp_2sigma.GetXaxis().SetNdivisions(1005)

	graphTest.Draw("AP")
	graphTest1.Draw("P")
	graphTest2.Draw("P")
	graphTest3.Draw("P")

	legend.AddEntry(graphTest,"Compare with itself","lp")
	legend.AddEntry(graphTest1,"Compare function1","lp")
	legend.AddEntry(graphTest2,"Compare function2","lp")
	legend.AddEntry(graphTest3,"Compare function3","lp")
    	legend.Draw()

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(c, 4, 0)

	#c1.SetLogy()
	c1.SaveAs( 'Plots/BiasTest_'+args.decay+'_signal'+str(args.signalInj)+'_MeanValues_ResolvedAnalysis'+args.version+'.'+args.extension)

##########################################################
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
	outputfile.write("# lumi    lnN     1.045         -     \n")
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
	parser.add_argument('-t', '--miniTree', action='store_true', default=False, help='miniTree: if plots coming from miniTree or RUNAnalysis.' )
	parser.add_argument('-e', '--extension', action='store', default='png', help='Extension of plots.' )
	parser.add_argument('-C', '--cut', action='store', default='delta', help='cut, example: cutDEta' )
	parser.add_argument('-F', '--func', action='store', default='P3', help='Function, example: P3, P4, P5' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)


	filePrefix = 'Rootfiles/RUNMiniResolvedAnalysis'
	hist = 'massAve_'+args.cut+'_'
	scale = args.lumi 
	rebinX = 20 
	if ( args.mass > 0 ): listMass = [ args.mass ]
	else: 
		if '312' in args.decay: listMass = [ 240, 300, 350 ] + range( 450, 850, 50 ) + [ 950, 1000, 1100 ] # + range( 1000, 1600, 100 ) 
		else: listMass = [ 240, 280, 300, 350 ] + range( 450, 900, 50 ) + [950, 1000, 1100, 1200, 1300] #+  range( 1100, 1600, 100 )


	outputDir = "Plots/"
	signalFilename = filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_Moriond17_80X_V2p4_'+args.version+'.root'
	signalFile =  TFile.Open(signalFilename)
	bkgFile = TFile.Open(filePrefix+'_QCD'+args.qcd+'All_Moriond17_80X_V2p4_'+args.version+'.root')
	#bkgFile2 = TFile.Open(filePrefix+'_QCDHTAll_Moriond17_80X_V2p4_'+args.version+'.root')
	#bkgFile3 = TFile.Open(filePrefix+'_QCDHerwigAll_Moriond17_80X_V2p4_'+args.version+'.root')
	dataFile = TFile.Open(filePrefix+'_JetHT_Run2016_80X_V2p4_'+args.version+'.root')
	

	###### Input parameters
	histYaxis = "Average dijet mass [GeV]"
	minFit = 200
	maxFit = 1200
	CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 1 ) )+" fb^{-1}"
	CMS_lumi.extraText = "Preliminary Simulation"

	######## Fit Functions
	fitFunctions = {}
	fitFunctions['P5'] = [ TF1("P5", "[0]*TMath::Power(1-(x/13000.0),[1])/(TMath::Power(x/13000.0,[2]+([3]*TMath::Log(x/13000.))+([4]*TMath::Power(TMath::Log(x/13000.),2))))",0,2000), 
			[ 0, 100, 2, 0.1, 0.01], 
			'(pow(1-@0/13000,@1)/pow(@0/13000,@2+@3*log(@0/13000)+@4*pow(log(@0/13000),2))' ]

	fitFunctions['P4'] = [ TF1("P4", "[0]*TMath::Power(1-(x/13000.0),[1])/(TMath::Power(x/13000.0,[2]+([3]*TMath::Log(x/13000.))))",0,2000), 
			[ 0, 100, 2, 1], 
			'pow(1-@0/13000,@1)/pow(@0/13000,@2+@3*log(@0/13000))' ]

	fitFunctions['P3'] = [ TF1("P3", "[0]* TMath::Power(1-(x/13000.0),[1]) / (TMath::Power(x/13000.0,[2]))",0,2000), 
			[ 0, 100, 2], 
			'pow(1-@0/13000,@1)/pow(@0/13000,@2)' ]

	fitFunctions['P2'] = [ TF1("P2", "[0] / TMath::Power(x/13000.0,[1])",0,2000), 
			[ ], 
			'pow(@0/13000,-@1)' ]

	fitFunctions['dijet'] = [ TF1("dijet", "[0]+ TMath::Exp( ([1]*TMath::Log(x/13000.0)) + ([2]*TMath::Power(TMath::Log(x/13000.0),2)) + ([3]*TMath::Power(TMath::Log(x/13000.0),3)))", 0, 2000 ), 
			[ ],
			'exp( (@1*log(@0/13000)) + (@2*pow(log(@0/13000),2)) + ( @3*pow(log(@0/13000),3)))']

	fitFunctions['expoPoli'] = [ TF1("expoPoli", "[0]+exp([1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)", 0, 1000 ), 
			[ ],
			'(exp((@1*@0)+(@2*pow(@0,2))+(@3*pow(@0,3))+(@4*pow(@0,4))))']

	fitFunctions['atlas'] = [ 
			TF1("atlas", "[0]* TMath::Power((1-TMath::Power(x/13000.0,0.33)),[1]) / TMath::Power(x/13000.0,[2])",0,2000), 
			[ 1000, 100, 1 ], 
			'(pow(1-pow(@0/13000,0.33),@1)/pow(@0/13000,@2))' ]

	fitFunctions['P1'] = [ TF1("P1", "[0] / (TMath::Power(x/13000.0,[1]))",0,2000), [ 0] ]
	fitFunctions['landau'] = [ TF1("landau","[0]*TMath::Landau(-x,[1],[2])",0,2000), [ ] ]
	fitFunctions['gaus'] = [ TF1("gaus", "gaus", 0, 2000), [ ] ]
	fitFunctions['P4Gaus'] = [ TF1("P4Gaus", "[0]*pow(1-(x/13000.0),[1])/pow(x/13000.0,[2]+([3]*log(x/13000.)))+gaus(4)",0,2000), [] ]


	if 'full' in args.process:
		CMS_lumi.extraText = "Preliminary"
		p = Process( target=FitterCombination, args=( dataFile, bkgFile, signalFile, hist, scale, [fitFunctions['P4']], minFit, maxFit, rebinX, True ))
	elif 'Data' in args.process:
		CMS_lumi.extraText = "Preliminary"
		p = Process( target=FitterCombination, args=( dataFile, 
			bkgFile, 
			'', #bkgFile3, 
			hist, 
			scale, 
			[fitFunctions[args.func]], 
			minFit, maxFit, rebinX, True ))

	elif 'RPV' in args.process:
		process = 'RPVSt'+str(args.mass)+'tojj'
		p = Process( target=rootFitter, args=( signalFile, hist, scale, [fitFunctions['gaus']], args.mass-50, args.mass+50, rebinX, True, False ) ) 

	elif 'Limit' in args.process:
		p = Process( target=createCards, args=( dataFile, #bkgFile, 
			bkgFile, 
			signalFilename,
			listMass,
			hist, 
			scale, 
			[fitFunctions['P3']], 
			minFit, maxFit, 10 ) )

	elif 'Bias' in args.process:
		p = Process( target=createCards, args=( dataFile, 
			bkgFile, 
			signalFilename,
			listMass,
			hist, 
			scale, 
			[ fitFunctions['P3'], fitFunctions['atlas'], fitFunctions['expoPoli'], fitFunctions['dijet'] ], 
			minFit, maxFit, 1 ) )

	elif 'Fisher' in args.process:
		p = Process( target=FisherTest, args=( dataFile, 
			hist, 
			[ fitFunctions['P2'], fitFunctions['P3'], fitFunctions['P4'], fitFunctions['P5'] ], 
			minFit, maxFit, 10 ) )

	elif 'biasDraw' in args.process: 
		p = Process( target=drawBiasTest, args=( listMass,
			'/afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/' ) )

	p.start()
	p.join()
