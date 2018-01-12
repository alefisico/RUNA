#!/usr/bin/env python

################################
### Creating datacards
################################

from ROOT import *
#from ROOT import RooRealVar, RooDataHist, RooArgList, RooArgSet, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf, RooGaussian
from array import array
from collections import OrderedDict
import argparse
import glob,sys, os
import warnings
import random
import numpy as np
import subprocess
from RUNBkgEstimation import boostedMassAveBins, rebin
try: 
	from RUNA.RUNAnalysis.commonFunctions import *
	import RUNA.RUNAnalysis.CMS_lumi as CMS_lumi 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
except ImportError: 
	sys.path.append('../python') 
	from commonFunctions import *
	import CMS_lumi as CMS_lumi 
	import tdrstyle as tdrstyle


gROOT.Reset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
#gSystem.SetIncludePath('-I$ROOFITSYS/include')
if os.access('RooPowerFunction.cxx', os.R_OK): ROOT.gROOT.ProcessLine('.L RooPowerFunction.cxx+')

#xline = array('d', [0,2000])
#yline = array('d', [0,0])
#line = TGraph(2, xline, yline)
#line.SetLineColor(kRed)

def signalUnc( hSignal, signalMass ):
	"""docstring for signalUnc"""

	######## Signal uncertainties
        hSigSyst = {}
        signalCDF = TGraph(hSignal.GetNbinsX()+1)

        # JES and JER uncertainties
	if args.unc:
		signalCDF.SetPoint(0,0.,0.)
		integral = 0.
		for i in range(1, hSignal.GetNbinsX()+1):
			x = hSignal.GetXaxis().GetBinLowEdge(i+1)
			integral = integral + hSignal.GetBinContent(i)
			signalCDF.SetPoint(i,x,integral)

		hSigSyst['JESUp'] = hSignal.Clone()
		hSigSyst['JESDown'] = hSignal.Clone()

		hSigSyst['JERUp'] = hSignal.Clone()
		hSigSyst['JERDown'] = hSignal.Clone()


        # reset signal histograms
        for key in hSigSyst: hSigSyst[key].Reset()

        # produce JES signal shapes
        if args.unc:
		for q in range(1, hSignal.GetNbinsX()+1):
			xLow = hSignal.GetXaxis().GetBinLowEdge(q)
			xUp = hSignal.GetXaxis().GetBinLowEdge(q+1)
			jes = 1. + jesValue
			xLowPrime = jes*xLow
			xUpPrime = jes*xUp
			hSigSyst['JESDown'].SetBinContent(q, signalCDF.Eval(xUpPrime) - signalCDF.Eval(xLowPrime))
			jes = 1. - jesValue
			xLowPrime = jes*xLow
			xUpPrime = jes*xUp
			hSigSyst['JESUp'].SetBinContent(q, signalCDF.Eval(xUpPrime) - signalCDF.Eval(xLowPrime))

        # produce JER signal shapes
	if args.unc:
		for i in range(1, hSignal.GetNbinsX()+1):
			xLow = hSignal.GetXaxis().GetBinLowEdge(i)
			xUp = hSignal.GetXaxis().GetBinLowEdge(i+1)
			jer = 1. + jerValue
			xLowPrime = jer*(xLow-float(signalMass))+float(signalMass)
			xUpPrime = jer*(xUp-float(signalMass))+float(signalMass)
			hSigSyst['JERDown'].SetBinContent(i, signalCDF.Eval(xUpPrime) - signalCDF.Eval(xLowPrime))
			jer = 1. - jerValue
			xLowPrime = jer*(xLow-float(signalMass))+float(signalMass)
			xUpPrime = jer*(xUp-float(signalMass))+float(signalMass)
			hSigSyst['JERUp'].SetBinContent(i, signalCDF.Eval(xUpPrime) - signalCDF.Eval(xLowPrime))

	return hSigSyst
############################################################################################

def createPseudoExperiment( histo, acceptance ):
	"""docstring for createPseudoExperiment"""

	hPseudo = histo.Clone()
	hPseudo.Reset()
	newNumEvents = random.randint( acceptance-round(TMath.Sqrt(acceptance)), acceptance+round(TMath.Sqrt(acceptance)) )
	print 'Events in sample:', acceptance, ', in PseudoExperiment:', newNumEvents
	hPseudo.FillRandom( histo, newNumEvents ) 
	return hPseudo
############################################################################################
	

def shapeCards( datahistosFile, histosFile, signalFile, signalSample, hist, signalMass, minMass, maxMass, outputName, outputFileTheta ):
	"""function to run Roofit and save workspace for RooStats"""
	warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='.*class stack<RooAbsArg\*,deque<RooAbsArg\*> >' )

	############################################################# DATA
	#hData = histosFile.Get('massAve_deltaEtaDijet_QCDPtAll')
	#hData = datahistosFile.Get('massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_ABCDProj')
	dataFile = TFile( datahistosFile )
	hData = dataFile.Get(hist+'_DATA')
	hData.Rebin( args.reBin )
	#hData = hData.Rebin( len( boostedMassAveBins )-1, hData.GetName(), boostedMassAveBins )
	#hData = histosFile.Get(hist+'_QCDPtAll_A')
	#hData.Add(htmpSignal)
	#hData.Scale(1/hData.Integral())
	#maxMass = boostedMassAveBins[ hData.FindLastBinAbove( 0, 1) ]
	#minMass = signalMass-30
	#maxMass = signalMass+30
	massAve = RooRealVar( 'massAve', 'massAve', minMass, maxMass  )
	#massAveData = RooRealVar( 'massAveData', 'massAveData', minMass, maxMass  )
	rooDataHist = RooDataHist('rooDatahist','rooDatahist',RooArgList(massAve), hData ) 
	rooDataHist.Print()
	############################################################################################

	
	####################### Signal 
	if 'gaus' in args.job: 
		hSignal = TH1F( 'massAve_RPVStop', 'massAve_RPVStop', maxMass/args.reBin, minMass, maxMass)
		for q in range( hSignal.GetNbinsX()+1 ):
			gausEval = signalFile.Eval( hSignal.GetXaxis().GetBinCenter( q ) )
			hSignal.SetBinContent( q, gausEval )

		#meanSig = RooRealVar( 'meanSig', 'mean of signal', sigGaus.GetParameter( 1 ) )
		#sigmaSig = RooRealVar( 'sigmaSig', 'sigma of signal', sigGaus.GetParameter( 2 ) )
		#signalPdf = RooGaussian( 'signal', 'signal', massAve, meanSig, sigmaSig )
		#signalPdf.Print()
	else:
		signalHistosFile = TFile( signalFile )
		hSignal = signalHistosFile.Get(hist+'_'+signalSample)
		hSignal.Rebin( args.reBin )
	#hSignal.Scale ( twoProngSF*antithreeProngSF )

	signalXS = search(dictXS, 'RPVStopStopToJets_UDD312_M-'+str(signalMass) )
	rooSigHist = RooDataHist( 'rooSigHist', 'rooSigHist', RooArgList(massAve), hSignal )
	sigAcc = rooSigHist.sumEntries()  #round(hSignal.Integral( hSignal.GetXaxis().FindBin( minMass ), hSignal.GetXaxis().FindBin( maxMass )), 2)
	rooSigHist.Print()

	#signal = RooHistPdf('signal','signal',RooArgSet(massAve),rooSigHist)
	#signal.Print()
	#signal_norm = RooRealVar('signal_norm','signal_norm',0,-1e+04,1e+04)
	#if args.fitBonly: signal_norm.setConstant()
	#signal_norm.Print()

	#####################################################################
	hSigSyst = signalUnc( hSignal, signalMass ) 
        hSigSystDataHist = {}
        if args.unc:
		hSigSystDataHist['JESUp'] = RooDataHist('hSignalJESUp','hSignalJESUp',RooArgList(massAve),hSigSyst['JESUp'])
		hSigSystDataHist['JESDown'] = RooDataHist('hSignalJESDown','hSignalJESDown',RooArgList(massAve),hSigSyst['JESDown'])

		hSigSystDataHist['JERUp'] = RooDataHist('hSignalJERUp','hSignalJERUp',RooArgList(massAve),hSigSyst['JERUp'])
		hSigSystDataHist['JERDown'] = RooDataHist('hSignalJERDown','hSignalJERDown',RooArgList(massAve),hSigSyst['JERDown'])


	#################################### Background
	if args.withABCDTFactor:
		htmpBkg = dataFile.Get('massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_ABCDProj')
		htmpBkg.Rebin( args.reBin )
	else: 
		newBkgHistoFile = datahistosFile.replace( 'DATA', 'DATA_ABCDBkg' )
		newBkgFile = TFile( newBkgHistoFile )
		htmpBkg = newBkgFile.Get('massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_ABCDProj' )
		if (htmpBkg.GetBinWidth( 1 ) != args.reBin ): 
			print '|----- Bin size in DATA_C histogram is different than rest.'
			sys.exit(0)
	#hBkg = histosFile.Get('massAve_prunedMassAsymVsdeltaEtaDijet_QCDPtAll_ABCDProj')
	#hBkg = histosFile.Get(hist+'_QCDPtAll_BCD')
	#htmpBkg = htmpBkg.Rebin( len( boostedMassAveBins )-1, htmpBkg.GetName(), boostedMassAveBins )
	hBkg = htmpBkg.Clone()
	hBkg.Reset()
	for ibin in range( htmpBkg.GetNbinsX() ):
		binCont = htmpBkg.GetBinContent( ibin )
		binErr = htmpBkg.GetBinError( ibin )
		if binCont == 0:
			hBkg.SetBinContent( ibin, 0 )
			hBkg.SetBinError( ibin, 1.8 )
		else:
			hBkg.SetBinContent( ibin, binCont )
			hBkg.SetBinError( ibin, binErr )
	#hBkg.Scale(1/hBkg.Integral())

	hPseudo = createPseudoExperiment( hBkg, bkgAcc )

	###### Adding statistical uncertanty
	hBkgStatUncUp = hBkg.Clone()
	hBkgStatUncUp.Reset()
	hBkgStatUncDown = hBkg.Clone()
	hBkgStatUncDown.Reset()
	for i in range( hBkg.GetNbinsX()+1 ):
		cont = hBkg.GetBinContent(i) 
		contErr = hBkg.GetBinError(i)
		hBkgStatUncUp.SetBinContent( i,  cont + (1*contErr) )
		hBkgStatUncDown.SetBinContent( i, cont - (1*contErr) )
	hBkgStatUncUpDataHist = RooDataHist( 'hBkgStatUncUp', 'hBkgStatUncUp', RooArgList(massAve), hBkgStatUncUp )
	hBkgStatUncDownDataHist = RooDataHist( 'hBkgStatUncDown', 'hBkgStatUncDown', RooArgList(massAve), hBkgStatUncDown )

	if 'template' in args.job:
		rooBkgHist = RooDataHist( 'rooBkgHist', 'rooBkgHist', RooArgList(massAve), hBkg )
		bkgAcc = rooBkgHist.sumEntries()
		rooBkgHist.Print()
		background = RooHistPdf('background','background',RooArgSet(massAve),rooBkgHist)
		background.Print()

	else:
		massAveBkg = RooRealVar( 'massAveBkg', 'massAveBkg', minMass, maxMass  )
		p1 = RooRealVar('p1','p1', 1 ,0.,100.)
		p2 = RooRealVar('p2','p2', 1 ,0.,60.)
		p3 = RooRealVar('p3','p3', 1 , -10.,10.)

		bkgAcc = round(hBkg.Integral( hBkg.GetXaxis().FindBin( minMass ), hBkg.GetXaxis().FindBin( maxMass )), 2)
		background = RooGenericPdf('background','(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(1300,1300,1300),RooArgList(massAveBkg,p1,p2,p3))
		background.Print()


        hBkgSyst = {}
        hBkgSystDataHist = {}

	print ' |---> Adding bkg unc'
	hBkgSyst['BkgUncUp'] = hBkg.Clone()
	hBkgSyst['BkgUncDown'] = hBkg.Clone()

	for key in hBkgSyst:
		hBkgSyst[key].Reset()
		hBkgSyst[key].SetName(hBkgSyst[key].GetName() + '_' + key)

	for q in range(0, hBkg.GetNbinsX()):
		binCont = hBkg.GetBinContent( q )
		bkgUncUp = 1. + (args.bkgUncValue/100.)
		hBkgSyst['BkgUncUp'].SetBinContent(q, binCont*bkgUncUp )
		bkgUncDown = 1. - (args.bkgUncValue/100.)
		hBkgSyst['BkgUncDown'].SetBinContent(q, binCont*bkgUncDown )
	hBkgSystDataHist['BkgUncUp'] = RooDataHist('hBkgBkgUncUp','hBkgBkgUncUp',RooArgList(massAve),hBkgSyst['BkgUncUp'])
	hBkgSystDataHist['BkgUncDown'] = RooDataHist('hBkgBkgUncDown','hBkgBkgUncDown',RooArgList(massAve),hBkgSyst['BkgUncDown'])

	############################################################################################


	#model = RooAddPdf("model","s+b",RooArgList(background,signal),RooArgList(background_norm,signal_norm))
	#res = model.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Strategy(0))
	#res.Print()


	############################ Create Workspace
	myWS = RooWorkspace("myWS")
        getattr(myWS,'import')(rooBkgHist,RooFit.Rename("background"))
        #getattr(myWS,'import')(background,RooFit.Rename("background"))
        #getattr(myWS,'import')(signal_norm)
        #getattr(myWS,'import')(background_norm)
	'''
	if 'gaus' in args.job: 
		getattr(myWS,'import')(signalPdf,RooFit.Rename("signal")) 
		getattr(myWS,'import')(signalPdfJESUp,RooFit.Rename("signal__JESUp")) 
		getattr(myWS,'import')(signalPdfJESDown,RooFit.Rename("signal__JESDown")) 
	else: 
	'''
	getattr(myWS,'import')(rooSigHist,RooFit.Rename("signal"))
	if args.unc:
		getattr(myWS,'import')(hSigSystDataHist['JESUp'],RooFit.Rename("signal__JESUp"))
		getattr(myWS,'import')(hSigSystDataHist['JESDown'],RooFit.Rename("signal__JESDown"))
		getattr(myWS,'import')(hSigSystDataHist['JERUp'],RooFit.Rename("signal__JERUp"))
		getattr(myWS,'import')(hSigSystDataHist['JERDown'],RooFit.Rename("signal__JERDown"))
	getattr(myWS,'import')(hBkgSystDataHist['BkgUncUp'],RooFit.Rename("background__BkgUncUp"))
	getattr(myWS,'import')(hBkgSystDataHist['BkgUncDown'],RooFit.Rename("background__BkgUncDown"))
	getattr(myWS,'import')(hBkgStatUncUpDataHist, RooFit.Rename("background__BkgStatUncUp") )
	getattr(myWS,'import')(hBkgStatUncDownDataHist, RooFit.Rename("background__BkgStatUncDown") )
        getattr(myWS,'import')(rooDataHist,RooFit.Rename("data_obs"))
        myWS.Print()
	outputRootFile = statDir+'/Rootfiles/workspace_'+outputName+'.root'
        myWS.writeToFile(outputRootFile, True)
	print ' |----> Workspace created in root file:\n', outputRootFile
	'''
	c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
#	c1.SetLogy()
	xframe = myWS.var("massAve").frame()
	signalPdf.plotOn( xframe )
	xframe.Draw()
	c1.SaveAs('test'+args.extension)
	del c1
	'''
	############################################################################################
        
	
	######################### write a datacard

	dataCardName = statDir+'/Datacards/datacard_'+outputName+'.txt'
        datacard = open( dataCardName ,'w')
        datacard.write('imax 1\n')
        datacard.write('jmax 1\n')
        datacard.write('kmax *\n')
        datacard.write('---------------\n')
        if args.unc: datacard.write('shapes * * '+outputRootFile+' myWS:$PROCESS myWS:$PROCESS__$SYSTEMATIC\n')
	else: datacard.write("shapes * * "+outputRootFile+" myWS:$PROCESS \n")
        datacard.write('---------------\n')
        datacard.write('bin '+signalSample+'\n')
        datacard.write('observation -1\n')
        datacard.write('------------------------------\n')
        datacard.write('bin          '+signalSample+'          '+signalSample+'\n')
        datacard.write('process      signal     background\n')
        datacard.write('process      0          1\n')
        #datacard.write('rate         -1         -1\n')
        datacard.write('rate         '+str(sigAcc)+'         '+str(bkgAcc)+'\n')
        datacard.write('------------------------------\n')
	datacard.write('lumi  lnN    %f         -\n'%(lumiValue))
	datacard.write('pu  lnN    %f         -\n'%(puValue))
        datacard.write('JES  shape   1          -\n')
	datacard.write('JER  shape   1          -\n')
        #flat parameters --- flat prior
	datacard.write('BkgUnc  shape   -	   1 \n')
	#NcombineUnc = ( 1 / TMath.Sqrt( args.bkgUncValue / 100. ) ) - 1
	#datacard.write('background_norm  gmN '+str(int(round(NcombineUnc)))+'  -  '+str( round(bkgAcc/NcombineUnc,2) )+'\n')
        #datacard.write('p1  flatParam\n')
	datacard.write('BkgStatUnc  shape   -	   1 \n')
        datacard.close()
	print ' |----> Datacard created:\n', dataCardName
	############################################################################################

	########## Theta
	if args.theta:
		print ' |----> Creating Theta file\n', outputFileTheta
		outFile = TFile( outputFileTheta, 'update')
		tmpName = 'rpvstopjj'+str(signalMass) 
		hSignal.SetName('massAve__'+tmpName)
		hSignal.Write()
		hSigSyst['JESDown'].SetName('massAve__'+tmpName+'__jes__down' )
		hSigSyst['JESDown'].Write()
		hSigSyst['JESUp'].SetName('massAve__'+tmpName+'__jes__up' )
		hSigSyst['JESUp'].Write()
		hSigSyst['JERDown'].SetName('massAve__'+tmpName+'__jer__down' )
		hSigSyst['JERDown'].Write()
		hSigSyst['JERUp'].SetName('massAve__'+tmpName+'__jer__up' )
		hSigSyst['JERUp'].Write()
		if (signalMass == 100): #or (signalMass == 170):
			hBkg.SetName('massAve__background')
			hBkg.Write()
			hBkgSyst['BkgUncDown'].SetName('massAve__background__unc__down')
			hBkgSyst['BkgUncDown'].Write()
			hBkgSyst['BkgUncUp'].SetName('massAve__background__unc__up')
			hBkgSyst['BkgUncUp'].Write()
			hData.SetName('massAve__DATA')
			hData.Write()
		outFile.Close()
############################################################################################

def plotGausQuantities( massList, varList, varErrList, varName, varTitle, func, Log=False, Grid=False ):
	"""docstring for plotGausQuantities"""

	zeroList = [0]*len(massList)
	varGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', varList), array('d', zeroList ), array( 'd', varErrList) )
	can = TCanvas(varName, varName, 10, 10, 750, 500 )
	if Log: can.SetLogy()
	if Grid: can.SetGrid()
	gStyle.SetOptFit(1)
	fitFunc = TF1( "fitFunc"+varName, func, 0, 500)
	for i in range(3): varGraph.Fit( fitFunc, 'MIR', '', 80, 400 )
	varGraph.SetMarkerStyle( 21 )
	varGraph.GetXaxis().SetTitle('Resonance mass [GeV]')
	varGraph.GetYaxis().SetTitle( varTitle )
	varGraph.GetYaxis().SetTitleOffset(0.95)
	varGraph.Draw('APS')
	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.lumi_13TeV = ""
	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(can, 4, 0)
	can.Update()
	st2 = varGraph.GetListOfFunctions().FindObject("stats")
	st2.SetX1NDC( (.70 if Log else .15) )
	st2.SetX2NDC( (.90 if Log else .35) )
	st2.SetY1NDC( .76 )
	st2.SetY2NDC( .91 )
	can.Modified()
	can.SaveAs( 'Plots/'+varName+'_'+args.decay+'_'+args.version+'.'+args.extension )
	del can

	return fitFunc 

############################################################################################

def createGausShapes( massList, newMassList, name, xmin, xmax, rebinX, labX, labY, log, plot=False):
	"""docstring for plot"""

	gausParam = {}
	histos = {}
	newGausFunct = {}

	constList = []
	constErrList = []
	acceptanceList = []
	acceptanceErrList = []
	meanList = []
	meanErrList = []
	sigmaList = []
	sigmaErrList = []

	for xmass in massList:

		inFileSample = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_RPVStopStopToJets_'+args.decay+'_M-'+str(xmass)+'_Moriond17_80X_V2p4_'+args.version+'.root'
		gStyle.SetOptFit(1)
		signalHistosFile = TFile.Open( inFileSample )
		hSignal = signalHistosFile.Get( name+'_RPVStopStopToJets_'+args.decay+'_M-'+str(xmass))
		tmpResolution = int(-1.78 + ( 0.1097 * xmass) + ( -0.0002897 * xmass*xmass ) + ( 3.18e-07 * xmass*xmass*xmass))
		hSignal.Rebin( tmpResolution )
		hSignal.Scale ( args.lumi*twoProngSF*antithreeProngSF )
		htmpSignal = hSignal.Clone()
		htmpSignal.Reset()

		tmpResolution = 4*tmpResolution
		sigGaus = TF1( 'sigGaus', 'gaus', 0, 1000 )
		for i in range(2): 
			sigGaus.SetParameter( 1, xmass )
			hSignal.Fit( sigGaus, 'MIR', '', xmass-tmpResolution, xmass+tmpResolution)

		gausParam[ xmass ] = sigGaus 
		acceptanceList.append( sigGaus.Integral( xmass-tmpResolution, xmass+tmpResolution ) )
		acceptanceErrList.append( sigGaus.IntegralError( xmass-tmpResolution, xmass+tmpResolution ) )
		constList.append( sigGaus.GetParameter( 0 ) )
		constErrList.append( sigGaus.GetParError( 0 ) )
		meanList.append( sigGaus.GetParameter( 1 ) )
		meanErrList.append( sigGaus.GetParError( 1 ) )
		sigmaList.append( sigGaus.GetParameter( 2 ) )
		sigmaErrList.append( sigGaus.GetParError( 2 ) )

		if plot:
			can1 = TCanvas('c'+str(xmass), 'c'+str(xmass),  10, 10, 750, 500 )
			hSignal.GetXaxis().SetRangeUser( xmass-200 , xmass+200 )
			hSignal.Draw()
			can1.SaveAs( 'Plots/test'+str(xmass)+'.'+args.extension )
			del can1
	
	zeroList = [0]*len(massList)
	'''
	numEventsGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', numEventsList), array('d', zeroList ), array( 'd', numEventsErrList) )
	canNumEvents = TCanvas('NumberEvents', 'NumberEvents',  10, 10, 750, 500 )
	gStyle.SetOptFit(1)
	numEventsFit = TF1("numEventsFit", "pol2", 60, 300 )
	for i in range(3): numEventsGraph.Fit( numEventsFit, 'MIR' )
	numEventsGraph.SetMarkerStyle( 21 )
	numEventsGraph.GetXaxis().SetTitle('Average pruned mass [GeV]')
	numEventsGraph.GetYaxis().SetTitle('Number of Events')
	numEventsGraph.GetYaxis().SetTitleOffset(0.95)
	numEventsGraph.Draw('AP')
	canNumEvents.SaveAs( 'Plots/numEvents_'+args.decay+'_'+args.cutTop+'_'+args.version+'.'+args.extension )
	del canNumEvents
	'''
	constFit = plotGausQuantities( massList,
					constList, constErrList,
					'Constant', 'Constant parameter',
					'[0]*TMath::Power(1-(x/13000.0),[1])/(TMath::Power(x/13000.0,([2]+([3]*TMath::Log(x/13000.0))+([4]*TMath::Power(TMath::Log(x/13000.),2)))))',
					Log=True
					)

	meanFit = plotGausQuantities( massList,
					meanList, meanErrList,
					'signalMean', 'Mean value from fit [GeV]',
					'pol1',
					Grid=True
					)

	sigmaFit = plotGausQuantities( massList,
					sigmaList, sigmaErrList,
					'signalSigma', 'Sigma value from fit [GeV]',
					'pol3',
					Log=True
					)

	for xmass in newMassList:
		sigNewGaus = TF1( 'sigNewGaus', 'gaus', 0, 500 )
		sigNewGaus.SetParameters( constFit.Eval( xmass ), 
						meanFit.Eval( xmass ), 
						sigmaFit.Eval( xmass ) )
						#(-1.78 + ( 0.1097 * xmass) + ( -0.0002897 * xmass*xmass ) + ( 3.18e-07 * xmass*xmass*xmass)))
		newGausFunct[ xmass ] = sigNewGaus 
	dummy=1
	can1 = TCanvas('c', 'c',  10, 10, 750, 500 )
	gausParam[ 80 ].Draw()
	gausParam[ 80 ].GetXaxis().SetRangeUser( 50, 400  )
	for x in gausParam:
		if x != 80: 
			gausParam[ x ].SetLineColor(dummy)
			gausParam[ x ].Draw("same")
		dummy=dummy+1
	can1.SaveAs( 'Plots/GaussShapes_'+args.decay+'_'+args.cutTop+'_'+args.version+'.'+args.extension )
	del can1

	dummy2 = 1
	canNewGaus = TCanvas('canNewGaus', 'canNewGaus',  10, 10, 750, 500 )
	#canNewGaus.SetLogy()
	newGausFunct[ 160 ].Draw()
	newGausFunct[ 160 ].GetXaxis().SetRangeUser( 50, 450  )
	#newGausFunct[ 80 ].SetMinimum(0.001)
	for x in newGausFunct:
		#if x != 80: 
		newGausFunct[ x ].SetLineColor(dummy2)
		newGausFunct[ x ].Draw("same")
		dummy2=dummy2+1
	canNewGaus.SaveAs( 'Plots/NewGaussShapes_'+args.decay+'_'+args.cutTop+'_'+args.version+'.'+args.extension)
	del canNewGaus

	return newGausFunct
############################################################################################

def binByBinCards( datahistosFile, bkghistosFile, signalFile, signalSample, hist, signalMass, minMass, maxMass, outputName ):
	"""docstring for binByBinCards"""

	####################### Data
	dataFile = TFile( datahistosFile )
	hData = dataFile.Get( hist+'_JetHT_Run2016')
	hData = rebin( hData, ( args.reBin if (args.reBin>0) else 'reso' ) )

	####################### Signal 
	if not args.ttbarAsSignal:
		if 'gaus' in args.job: 
			sigGaus = signalFile
			hSignal = TH1F( 'massAve_RPVStop', 'massAve_RPVStop', 1000, 0, 1000)
			hSignal = rebin( hSignal, ( args.reBin if (args.reBin>0) else 'reso' ) )
			for q in range( hSignal.GetNbinsX()+1 ):
				gausEval = signalFile.Eval( hSignal.GetXaxis().GetBinCenter( q ) )
				hSignal.SetBinContent( q, gausEval )
				#print hSignal.GetBinLowEdge(q), gausEval
			#c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
			#hSignal.Draw()
			#c1.SaveAs('test.'+args.extension)
			#del c1
		else:
			signalHistosFile = TFile( signalFile )
			hSignal = signalHistosFile.Get(hist+'_'+signalSample)
			hSignal = rebin( hSignal, ( args.reBin if (args.reBin>0) else 'reso' ) )
			hSignal.Scale( args.lumi )
			hSignal.Scale ( twoProngSF*antithreeProngSF )
			sigGaus = TF1( 'sigGaus', 'gaus', 0, 500 )
			tmpResolution = 2*(-1.78 + ( 0.1097 * signalMass) + ( -0.0002897 * signalMass*signalMass ) + ( 3.18e-07 * signalMass*signalMass*signalMass))
			hSignal.Fit( sigGaus, 'MIR', '', (signalMass-tmpResolution), (signalMass+tmpResolution) )
			#c1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
			#hSignal.GetXaxis().SetRangeUser(60, 100)
			#hSignal.Draw()
			#c1.SaveAs('test.'+args.extension)
			#del c1
		signalMassWidth = sigGaus.GetParameter(2)
		signalMaximum = sigGaus.GetParameter(1)

		hSigSyst = signalUnc( hSignal, signalMass ) 
		if args.signalInjec and (signalMass == 100):
			print ' |----> Runnning PseudoData'
			hPseudoData = createPseudoExperiment( hData, hData.GetEntries() )
			hPseudoSignal = createPseudoExperiment( hSignal, hSignal.GetEntries() )
			hPseudoData.Add( hPseudoSignal )
			hData = hPseudoData.Clone()

	##### Bkg estimation
	if args.mcAsData: 
		hDataC = bkghistosFile[ 'QCDPtAll' ][0].Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_QCDPtAll'+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_C')
		hDataC.Scale( bkghistosFile[ 'QCDPtAll' ][1] )
	else: hDataC = dataFile.Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016'+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_C')
	#else: hDataC = dataFile.Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_C')
	hDataC = rebin( hDataC, ( args.reBin if (args.reBin>0) else 'reso' ) )
	#testFile = TFile('test.root', 'recreate')
	#hDataC.Write()
	#hData.Write()
	#testFile.Close()
	#sys.exit(0)

	if args.withABCDTFactor:
		hDataB = dataFile.Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016'+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_B')
		#hDataB = dataFile.Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_B')
		hDataB = rebin( hDataB, ( args.reBin if (args.reBin>0) else 'reso' ) )
		hDataD = dataFile.Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016'+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_D')
		#hDataD = dataFile.Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_D')
		hDataD = rebin( hDataD, ( args.reBin if (args.reBin>0) else 'reso' ) )
	else:
		newBkgHistoFile = datahistosFile.replace( 'V2p4', 'V2p4_combinedBD'+(args.numBtags if 'UDD323' in args.decay else '')+'_'+args.cutTop+'_ABCDEst' )
		#newBkgHistoFile = datahistosFile.replace( 'V2p4', 'V2p4_combinedBD'+'_'+args.cutTop+'_ABCDEst' )
		newBkgFile = TFile( newBkgHistoFile )
		hDataRatioBD = newBkgFile.Get('massAve_prunedMassAsymVsdeltaEtaDijet_DATAMinusResBkg_RatioBD' )
		#hDataC = newBkgFile.Get( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016'+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_C')
		if (hDataRatioBD.GetBinWidth( 15 ) != hDataC.GetBinWidth( 15 ) ): 
			print '|----- Bin size in DATA_C histogram is different than rest.'
			sys.exit(0)


	bkgHistos = OrderedDict()
	for sample in bkghistosFile:
		if 'QCD' in sample and args.mcAsData: 
			hMCAsData = bkghistosFile[ sample ][0].Get( hist+'_'+ sample )
			hMCAsData = rebin( hMCAsData, ( args.reBin if (args.reBin>0) else 'reso' ) )
			hMCAsData.Scale( twoProngSF*antithreeProngSF * bkghistosFile[ sample ][1] )
		else: 
			bkgHistos[ sample ] = bkghistosFile[ sample ][0].Get( hist+'_'+ sample )
			bkgHistos[ sample ] = rebin( bkgHistos[ sample ], ( args.reBin if (args.reBin>0) else 'reso' ) )
			bkgHistos[ sample ].Scale( twoProngSF*antithreeProngSF * bkghistosFile[ sample ][1] )
			if args.ttbarAsSignal and ('TT' in sample): hSignal = bkgHistos[ sample ].Clone()
	
	if args.mcAsData:
		for sam in bkgHistos: hMCAsData.Add( bkgHistos[ sam ] )
		hData = hMCAsData.Clone()


	if not args.ttbarAsSignal:
		if (args.reBin>0) :
			lowEdgeWindow = int(signalMaximum/args.reBin - 2*( int( signalMassWidth )/args.reBin ))+1
			highEdgeWindow = int(signalMaximum/args.reBin + 2*( int( signalMassWidth )/args.reBin ))+2
		else:
			peak = min(range(len(boostedMassAveBins)), key=lambda x:abs(boostedMassAveBins[x]-int(signalMaximum)))
			#print peak, boostedMassAveBins[peak]
			lowEdgeWindow = peak - ( 1 if (signalMass in [200, 300, 350, 400] ) else 2 )
			highEdgeWindow = peak + (2 if (signalMass in [ 100, 120, 200, 300, 400 ] ) else 1 )
	else: 
		signalMassWidht = 20 ## dummy
		lowEdgeWindow = 30
		highEdgeWindow = 40
	print '%'*30, signalMassWidth, lowEdgeWindow, highEdgeWindow, boostedMassAveBins[lowEdgeWindow], boostedMassAveBins[highEdgeWindow], lowEdgeWindow*args.reBin, highEdgeWindow*args.reBin

	#### test
	#for ibin in range(1, hSignal.GetNbinsX()+1):
	#	print ibin, hSignal.GetBinContent( ibin ), hSignal.GetBinCenter( ibin )
	############################

	totalData = 0
	totalBkg = 0
	totalBkgStat = 0
	totalSig = 0
	totalSigStat = 0

	combineCards = 'combineCards.py '
	accDict = OrderedDict()
	for ibin in range( lowEdgeWindow+1, highEdgeWindow+2 ):

		### Signal
		sigAcc = hSignal.GetBinContent( ibin )
		print 'signal', hSignal.GetBinLowEdge( ibin ), sigAcc
		if ( sigAcc == 0 ): continue
		totalSig += sigAcc
		sigAccUnc = hSignal.GetBinError( ibin )
		totalSigStat += sigAccUnc
		sigStatUnc = 1 + sigAccUnc/sigAcc  
		accDict[ 'signal' ] = [ round(sigAcc,3), round(sigStatUnc,3) ]
		if args.unc:
			if hSigSyst['JERDown'].GetBinContent( ibin ) != 0:
				try: sigShapeJERDown = round(1/( hSigSyst['JERDown'].GetBinContent( ibin )/ sigAcc ), 3)
				except ZeroDivisionError: sigShapeJERDown = 1.8 
			else: sigShapeJERDown = 1
			if hSigSyst['JERUp'].GetBinContent( ibin ) != 0:
				try: sigShapeJERUp = round(1/( hSigSyst['JERUp'].GetBinContent( ibin )/ sigAcc ), 3)
				except ZeroDivisionError: sigShapeJERUp = 1.8 
			else: sigShapeJERUp = 1

			### it has to be asymmetical: https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/530/2.html
			if hSigSyst['JESDown'].GetBinContent( ibin ) != 0:
				try: sigShapeJESDown = round(1/( hSigSyst['JESDown'].GetBinContent( ibin )/ sigAcc ), 3)
				except ZeroDivisionError: sigShapeJESDown = 1.8 
			else: sigShapeJESDown = 1
			if hSigSyst['JESUp'].GetBinContent( ibin ) != 0:
				try: sigShapeJESUp = round(( hSigSyst['JESUp'].GetBinContent( ibin )/ sigAcc ), 3)
				except ZeroDivisionError: sigShapeJESUp = 1.8 
			else: sigShapeJESUp = 1
		#print sigShapeJERDown, sigShapeJERUp, sigShapeJESDown, sigShapeJESUp

		### data
		contData = hData.GetBinContent( ibin )
		totalData += contData
		#print 'data', hData.GetBinLowEdge( ibin ), contData

		### bkg
		contDataC = hDataC.GetBinContent( ibin )
		if args.withABCDTFactor:
			contDataB = hDataB.GetBinContent( ibin )
			contDataD = hDataD.GetBinContent( ibin )
			try: tf = contDataB/contDataD
			except ZeroDivisionError: tf = 0
			try: errBD = 1 + (TMath.Sqrt( TMath.Power( TMath.Sqrt( contDataD ) / contDataD, 2 ) + TMath.Power( TMath.Sqrt( contDataB ) / contDataB, 2 ) ) / tf )
			except ZeroDivisionError: errBD = 1
		else:
			tf = round( hDataRatioBD.GetBinContent( ibin ), 3 )
			errBD = 1+ ( hDataRatioBD.GetBinError( ibin ) / tf )
		bkgAcc = tf * contDataC
		totalBkg += bkgAcc
		accDict[ 'qcd' ] = [ round(bkgAcc,3), round(args.bkgUncValue,3) ]

		#### adding MC bkgs
		for sample in bkgHistos:
			mcbkgacc = abs(bkgHistos[ sample ].GetBinContent( ibin ))
			#print bkgHistos[ sample ].GetBinLowEdge( ibin )
			try: mcbkgstatunc = 1 + ( abs(bkgHistos[ sample ].GetBinError( ibin ))/mcbkgacc )
			except ZeroDivisionError: mcbkgstatunc = 1.8
			totalBkg += mcbkgacc
			totalBkgStat += abs(bkgHistos[ sample ].GetBinError( ibin ))
			accDict[ sample.lower() ] = [ round(mcbkgacc,3), round(mcbkgstatunc,3) ]
		if args.ttbarAsSignal: del accDict[ 'TT' ]
		#print ibin, sigAcc, accDict, contData
				

		######################## Creating  datacards
		dataCardName = statDir+'/Datacards/datacard_'+outputName+'_bin'+str(ibin)+'.txt'
		datacard = open( dataCardName ,'w')
		datacard.write('imax 1\n')
		datacard.write('jmax *\n')
		datacard.write('kmax *\n')
		datacard.write(('-'*30)+'\n')
		datacard.write('bin\t\t'+signalSample+'bin'+str(ibin)+'\n')
		datacard.write('observation\t\t'+str(int(contData))+'\n')
		datacard.write(('-'*30)+'\n')
		datacard.write('bin\t\t' + '\t'.join([ signalSample+'bin'+str(ibin) ]*len(accDict)) + '\n' ) 
		datacard.write('process\t\t' + '\t'.join([ sample for sample, info in accDict.items() ]) + '\n' ) 
		datacard.write('process\t\t' + '\t'.join([ str(q) for q in range( len(accDict) ) ]) + '\n' )
		datacard.write('rate\t\t' + '\t'.join([ str(info[0]) for sample, info in accDict.items() ]) + '\n' ) 
		datacard.write(('-'*30)+'\n')
		### ABCD region in combine
		tmpList = ['-']*len(accDict)
		datacard.write('qcdABCDmethod'+str(ibin)+'\tgmN\t'+str(int(contDataC))+'\t'+ '\t'.join( [ str(tf) if i==1 else '-' for i in range(len(accDict)) ] ) +'\n' ) 

		if args.unc: 
			datacard.write('lumi\t\tlnN\t'+ '\t'.join( [ '-' if i==1 else str(round(lumiValue,3)) for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('trigger\t\tlnN\t'+ '\t'.join( [ '-' if i==1 else str(round(triggerValue,3)) for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('twoProngSFSys\t\tlnN\t'+ '\t'.join( [ '-' if i==1 else str(round(twoProngValue,3)) for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('antithreeProngSFSys\t\tlnN\t'+ '\t'.join( [ '-' if i==1 else str(round(antithreeProngValue,3)) for i in range(len(accDict)) ] ) +'\n' ) 

			datacard.write('qcdtfUnc\t\tlnN\t'+ '\t'.join( [ str(round(1+(args.bkgUncValue/100.),3)) if i==1 else '-' for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('qcdtfStatUnc\t\tlnN\t'+ '\t'.join( [ str(errBD) if i==1 else '-' for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('mcUnc\t\tlnN\t'+ '\t'.join( [ str(MCunc) if i>1 else '-' for i in range(len(accDict)) ] ) +'\n' ) 
			### statistical uncertanties
			tmp = 0
			for sample, acc in accDict.items():
				if not 'qcd' in sample:
					datacard.write(sample+'StatUnc'+str(ibin)+'\t\tlnN\t' + '\t'.join([ str(acc[1]) if tmp==i else '-' for i in range(len(accDict)) ]) + '\n' )
				tmp+=1

			datacard.write('JESshape\t\tlnN\t'+ '\t'.join( [ '-' if i!=0 else str(sigShapeJESDown)+'/'+str(sigShapeJESUp) for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('JESaccept\t\tlnN\t'+ '\t'.join( [ '-' if i!=0 else str(uncJES) for i in range(len(accDict)) ]) +'\n' ) 

			datacard.write('JERshape\t\tlnN\t'+ '\t'.join( [ '-' if i!=0 else str(sigShapeJERDown)+'/'+str(sigShapeJERUp) for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('JERaccept\t\tlnN\t'+ '\t'.join( [ '-' if i!=0 else str(uncJER) for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('PUaccept\t\tlnN\t'+ '\t'.join( [ '-' if i!=0 else str(uncPU) for i in range(len(accDict)) ] ) +'\n' ) 
			datacard.write('PDFaccept\t\tlnN\t'+ '\t'.join( [ '-' if i!=0 else str(uncPDF) for i in range(len(accDict)) ] ) +'\n' ) 
			if 'UDD323' in args.decay: datacard.write('Btagaccept\t\tlnN\t'+ '\t'.join( [ '-' if i!=0 else str(uncBtag) for i in range(len(accDict)) ] ) +'\n' ) 


		datacard.close()
		combineCards += 'Bin'+str(ibin)+'='+dataCardName+' '
		print ' |----> Datacard created:\n', dataCardName

	if args.runCombine:
		print ' |----> Running combinedCards.py:\n', statDir+'/Datacards/datacard_'+outputName+'_bins.txt'
		subprocess.call(  combineCards+' > '+statDir+'/Datacards/datacard_'+outputName+'_bins.txt', shell=True  )
		print ' |----> Running combine -M Asymptotic:'
		newOutputName = outputName.replace('pruned', 'Boosted')
		subprocess.call( 'combine -M Asymptotic '+statDir+'/Datacards/datacard_'+outputName+'_bins.txt --rMax 0.5 -n _'+newOutputName, shell=True )
		#subprocess.call( 'combine -M Asymptotic '+statDir+'/Datacards/datacard_'+outputName+'_bins.txt -n _'+outputName+'_2sigma', shell=True )
		subprocess.call( 'mv higgsCombine* '+statDir, shell=True )
		print ' |----> Done. Have a wonderful day. :D'

	print '@'*10, signalMass, ' & [ ', int(boostedMassAveBins[lowEdgeWindow]), ',', int(boostedMassAveBins[highEdgeWindow]), ' ] & ', int(totalData), ' & ', round(totalBkg,1), '$\pm$', round(totalBkgStat,1), ' & ', round(totalSig,1), '$\pm$', round(totalSigStat,1), '\\\\'


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-j', '--job', action='store', dest='job', default='template', help='Process: template or fit.' )
	parser.add_argument('-i', '--injSig', dest='signalInjec', action="store_true", default=False, help='Signal injection test.' )
	parser.add_argument('-t', '--ttbarAsSignal', dest='ttbarAsSignal', action="store_true", default=False, help='ttbar as signal' )
	parser.add_argument('-l', '--lumi', dest='lumi', action="store", default=1, type=int, help='Luminosity, example: 1.' )
	parser.add_argument('-nV', '--bkgUncValue', dest='bkgUncValue', type=int, default=10, help='Value for bkg nomralization uncertainty.' )
	parser.add_argument('-u', '--unc', dest='unc', action="store_true", default=False, help='Luminosity, example: 1.' )
	parser.add_argument('-g', '--grom', action='store', default='pruned', dest='grooming', help='Grooming Algorithm, example: Pruned, Filtered.' )
	parser.add_argument('-d', '--decay', action='store', default='UDD312', dest='decay', help='Decay, example: UDD312, UDD323.' )
	parser.add_argument('-a', '--withABCDTFactor', action="store_true", default=False, help='Regular ABCD (False) or alternative ABCD (true).' )
	parser.add_argument('-R', '--rebin', dest='reBin', type=int, default=5, help='Data: data or pseudoData.' )
    	parser.add_argument('-e', "--theta", dest="theta", action="store_true", default=False, help="Create theta file.")
	parser.add_argument('-v', '--version', action='store', default='v05', dest='version', help='Version of rootfiles: v05.' )
	parser.add_argument('-m', '--mass', dest='massValue', action="store", type=int, default=-1, help='To run in a mass point only' )
	parser.add_argument('-c', '--cutTop', action='store', default='jet1Tau32', dest='cutTop', help='Version of rootfiles: v05.' )
    	parser.add_argument('-M', "--mcAsData", dest="mcAsData", action="store_true", default=False, help="Create theta file.")
	parser.add_argument('-E', '--extension', action='store', default='png', dest='extension', help='Extension of plots: png, pdf' )
	parser.add_argument('-b', '--numBtags', dest='numBtags', default='2btag', help='Number of btags' )
    	parser.add_argument('-r', "--runCombine", dest="runCombine", action="store_true", default=False, help="Create theta file.")
    	parser.add_argument('-T', "--createTable", dest="createTable", action="store_true", default=False, help="Create table.")

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	###### Input parameters
	masses = OrderedDict()
	minMass = 60 
	maxMass = 451 
	lumiValue = 1.025
	twoProngSF = 1.21
	antithreeProngSF = 0.96
	jesValue = 0.0
	jerValue = 0.23
	triggerValue = 1.03
	twoProngValue = 1.23
	antithreeProngValue = 1.01
	#ttbarSF = ( 1.00 if 'jet1Tau32' in args.cutTop else 0.96 )
	MCunc = ( 1.10 if 'jet1Tau32' in args.cutTop else 1.25 )

	uncJES  = 1.012
	uncJER  = 1.018
	uncPU   = 1.01
	uncPDF  = 1.01
	uncBtag = 1.01

	statDir = os.getcwdu()+'/../../RUNStatistics/test/'
	outputFileTheta = ''
	if args.theta:
		outputFileTheta = statDir+'/Rootfiles/theta_histos_Bin'+str(args.reBin)+'_'+args.version+'.root'
		if 'gaus' in args.job: outputFileTheta = outputFileTheta.replace(args.version, args.version+'_GaussShape')
		if args.withABCDTFactor:  outputFileTheta = outputFileTheta.replace(args.version, args.version+'_withABCDTFactor')
		files = glob.glob(outputFileTheta)
		for f in files: os.remove(f)

	if 'gaus' in args.job: 

		massList = range( 80, 200, 20) + range(200, 550, 50 )
		gausFunctList = createGausShapes( 
						#range( 80, 260, 20 )+[300, 350, 400], 
						[ 80, 100, 120, 140, 160, 180, 200, 300, 350, 400 ],
						massList,
						'massAve_'+( 'deltaEtaDijet' if 'UDD312' in args.decay else args.numBtags ), 
						0, 500, args.reBin, 
						0.85, 0.45, 
						False, 
						plot=True )
	else: 
		if 'UDD312' in args.decay: massList = [ 80, 100, 120, 140, 160, 180, 200, 220, 240, 300, 350, 400, 450, 500, 550 ]
		else: massList = [ 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 450, 500, 550 ]

	if args.massValue > 0: massList = [ args.massValue ]
	if args.signalInjec: massList = massList * 1000
	if args.ttbarAsSignal: massList = [ 'TT' ]

	dummy0 = 0
	for mass in massList: #range( len(massList) ):
		if isinstance( mass, int): signalSample = 'RPVStopStopToJets_'+args.decay+'_M-'+str(mass)
		else: signalSample = str(mass)
		dataFileHistos = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_JetHT_Run2016_80X_V2p4_'+args.version+'.root'
		bkgFileHistos = {}

		bkgFileHistos[ 'Dibosons' ] = [ TFile.Open('Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_Dibosons_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*1., 'Dibosons', kMagenta+2 ]
		#bkgFileHistos[ 'ZJetsToQQ' ] = [ TFile.Open('Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_ZJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*1.45, 'Z + Jets', kOrange]
    		bkgFileHistos[ 'DYJetsToQQ' ] = [ TFile.Open('Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_DYJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*1.45, 'DY + Jets', kOrange]
		bkgFileHistos[ 'WJetsToQQ' ] = [ TFile.Open('Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_WJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*1.35, 'W + Jets', 38]
		bkgFileHistos[ 'TT' ] = [ TFile.Open('Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_TT_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*0.98, 't #bar{t} + Jets', kGreen+2 ]
		if args.mcAsData: bkgFileHistos[ 'QCDPtAll' ] = [ TFile.Open('Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_QCDPtAll_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*0.28, 'MC QCD multijets', kBlue ]

		#if args.unc: outputName = signalSample+'_'+args.grooming+'_deltaEtaDijet_'+args.cutTop+'_Bin'+str(args.reBin)+'_'+args.version
		if args.unc: outputName = signalSample+'_'+args.grooming+'_'+args.cutTop+'_Bin'+str(( args.reBin if args.reBin>0 else 'Reso' ))+'_'+args.version
		else: outputName = signalSample+'_'+args.grooming+'_'+args.cutTop+'_NOSys_'+args.version
		if args.signalInjec: 
			outputName = outputName.replace( args.grooming, args.grooming+'_signalInjectionTest'+str(dummy0) )
			dummy0 += 1
		if args.withABCDTFactor: outputName = outputName.replace( args.grooming, args.grooming+'_withABCDTFactor' )
		if args.mcAsData: outputName = outputName.replace( args.grooming, args.grooming+'_MCasDATA' )
		if 'UDD323' in args.decay: outputName = outputName.replace( args.cutTop, args.numBtags+'_'+args.cutTop )
		if 'gaus' in args.job: 
			outputName = outputName+'_GaussShape'
			signalFileHistos = gausFunctList[ mass ] 
		else: 
			signalFileHistos = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming+'_RPVStopStopToJets_'+args.decay+'_M-'+str( mass )+'_Moriond17_80X_V2p4_'+args.version+'.root'

		print '#'*50 
		print ' |----> Creating datacard and workspace for '+( str(mass) if args.ttbarAsSignal else 'RPV St'+str( mass ) )
		print '#'*50  

		if 'bin' in args.job: 
			binByBinCards( dataFileHistos, bkgFileHistos, signalFileHistos, 
					signalSample, 
					'massAve_'+( args.numBtags if 'UDD323' in args.decay else 'deltaEtaDijet' ), 
					#'massAve_deltaEtaDijet', 
					mass, 
					minMass, maxMass, 
					outputName ) 

		else: shapeCards( dataFileHistos, TFile(bkgFileHistos), signalFileHistos, 
					signalSample, 
					'massAve_'+( args.numBtags if 'UDD323' in args.decay else 'deltaEtaDijet' ), 
					mass, 
					minMass, maxMass, 
					outputName, 
					outputFileTheta ) 
