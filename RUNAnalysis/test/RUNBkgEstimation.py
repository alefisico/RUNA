#!/usr/bin/env python
'''
File: DrawHistogram.py
Author: Alejandro Gomez Espinosa
Email: gomez@physics.rutgers.edu
Description: My Draw histograms. Check for options at the end.
'''

#from ROOT import TFile, TH1F, THStack, TCanvas, TMath, gROOT, gPad
from ROOT import *
import time, os, math, sys
from array import array
import argparse
from collections import OrderedDict
from DrawHistogram import Rebin2D
try:
	from RUNA.RUNAnalysis.histoLabels import labels, labelAxis, finalLabels, setSelection
	from RUNA.RUNAnalysis.scaleFactors import * #scaleFactor as SF
	from RUNA.RUNAnalysis.commonFunctions import * 
	import RUNA.RUNAnalysis.CMS_lumi as CMS_lumi 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
except ImportError:
	sys.path.append('../python')
	from histoLabels import labels, labelAxis, finalLabels
	from scaleFactors import * #scaleFactor as SF
	from commonFunctions import * 
	import CMS_lumi as CMS_lumi 
	import tdrstyle as tdrstyle

gROOT.Reset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)

xline = array('d', [0,2000])
yline = array('d', [1, 1])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)

yline11 = array('d', [1.1, 1.1])
line11 = TGraph(2, xline, yline11)
line11.SetLineColor(kBlue)

yline09 = array('d', [0.9, 0.9])
line09 = TGraph(2, xline, yline09)
line09.SetLineColor(kBlue)

yline0 = array('d', [0,0])
line0 = TGraph(2, xline, yline0)
line0.SetLineColor(kRed)


def listOfCont( histo ):
 	"""docstring for listOfCont"""
	tmpListContent = []
	tmpListError = []
	for ibin in range( histo.GetNbinsX() ): 
		tmpListContent.append( histo.GetBinContent( ibin ) )
		tmpListError.append( histo.GetBinError( ibin ) )
	return tmpListContent, tmpListError
######################################################


def BCDHisto( tmpHisto, BHisto, CHisto, DHisto ):
	"""docstring for BCDHisto: simple ABCD, order between B or C does not matter"""

	tmpHisto.Reset()
	for jbin in range( 0, tmpHisto.GetNbinsX()+1 ):
		Nominal_Side = BHisto.GetBinContent( jbin )
		Side_Nominal = CHisto.GetBinContent( jbin )
		NoSignal = DHisto.GetBinContent( jbin )
		if NoSignal != 0: 
			Bkg = Nominal_Side*Side_Nominal/NoSignal
			try: BkgError = Bkg * TMath.Sqrt( 
					TMath.Power(( TMath.Sqrt( Nominal_Side ) / Nominal_Side ), 2) + 
					TMath.Power(( TMath.Sqrt( Side_Nominal ) / Side_Nominal ), 2) + 
					TMath.Power(( TMath.Sqrt( NoSignal ) / NoSignal ), 2) )
			except ZeroDivisionError: BkgError = 0
		else: 
			Bkg = 0
			BkgError = 1.8		### Poisson errors for bins with 0 content
		tmpHisto.SetBinContent( jbin, Bkg )
		tmpHisto.SetBinError( jbin, BkgError )
	return tmpHisto
########################################################

def addSysBand( histo, uncSys, colour, additionalSys='' ):
	"""docstring for addSysBand"""
	
	hclone = histo.Clone()
	hclone.Reset()
	for i in range( 0, histo.GetNbinsX()+1 ):
		cont = histo.GetBinContent( i )
		contError = histo.GetBinError( i )
		hclone.SetBinContent( i, cont )
		if additionalSys:
			addCont = additionalSys.GetBinContent( i )
			addContError = additionalSys.GetBinError( i )
		totalErr = TMath.Sqrt( TMath.Power(((cont*uncSys) -cont), 2 ) +  TMath.Power( contError, 2 ) + ( TMath.Power( addContError, 2 ) if additionalSys else 0 ) )
		hclone.SetBinError( i, totalErr  )
	#hclone.SetFillStyle(3018)
	hclone.SetFillColor( colour )
	return hclone
######################################################################


def makePulls( histo1, histo2 ):
	"""docstring for makePulls"""

	histoPulls = histo1.Clone()
	histoPulls.Reset()
	pullsOnly = TH1F( 'pullsOnly', 'pullsOnly', 14, -3, 3 )
	for ibin in range(0, histo1.GetNbinsX() ):
		try: pull = ( histo1.GetBinContent( ibin ) - histo2.GetBinContent( ibin ) ) / histo1.GetBinError( ibin ) 
		except ZeroDivisionError: pull = 0
		histoPulls.SetBinContent( ibin, pull )
		histoPulls.SetBinError( ibin, 1 )
		pullsOnly.Fill( pull )

	return histoPulls, pullsOnly
######################################################################

def rebin( histo, binning ):
	"""docstring for rebin"""

	oldhisto = histo.Clone()
	if 'reso' in str(binning):  
		boostedMassAveBins = array( 'd', [ 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 300, 310, 320, 330, 340, 350, 500, 510 ] )
		newhisto = oldhisto.Rebin( len( boostedMassAveBins )-1, oldhisto.GetName(), boostedMassAveBins )
	elif 'ratio' in str(binning): 
		boostedMassAveBins = array( 'd', [ 0, 25, 75, 150, 250, 450, 500 ] )
		#boostedMassAveBins = array( 'd', [ 0, 25, 50, 75, 100, 150, 200, 300, 400, 500 ] )
		newhisto = oldhisto.Rebin( len( boostedMassAveBins )-1, oldhisto.GetName(), boostedMassAveBins )
	else:
		newhisto = oldhisto.Rebin( binning )

	return newhisto
######################################################################
	

######### THIS IS AN OLD FUNCTION... I AM USING TGRAPHASYMMERRORS NOW
def ratioPlots( histo1, histo2, graphOnly=False ):
	"""docstring for ratioPlots"""

	chi2 = 0
	ndf = 0
	#h1errFull = histo1.Clone()
	#h1errFull.Reset()
	#h1errh2 = histo1.Clone()
	#h1errh2.Reset()
	ratioList = []
	binCenterList = []
	ratioLogNErrXPlusList = []
	ratioLogNErrXMinusList = []

	for ibin in range( 0, histo1.GetNbinsX() ):
		binCenterList.append( histo1.GetXaxis().GetBinCenter( ibin ) )
		x = histo1.GetBinContent( ibin )
		xErr = histo1.GetBinError( ibin )
		y = histo2.GetBinContent( ibin )
		yErr = histo2.GetBinError( ibin )
		try: 
			ratio = x/y
			ratioErrX = ratio * TMath.Sqrt( TMath.Power( xErr/x, 2) + TMath.Power( yErr/y, 2)  )
			ratioErrY = ratio* yErr / y
			ratioLogNErrXPlus = TMath.Sqrt( TMath.Power( ( (x/(y-yErr)) - ratio ), 2 ) + TMath.Power( ( ((x+xErr)/y) - ratio ) , 2) ) 
			ratioLogNErrXMinus = TMath.Sqrt( TMath.Power( ( (x/(y+yErr)) - ratio ), 2 )  + TMath.Power( ( ((x-xErr)/y) - ratio ) , 2) ) 
		except ZeroDivisionError: 
			ratio = 0
			ratioErrX = 0
			ratioErrY = 0
			ratioLogNErrXPlus = 0
			ratioLogNErrXMinus = 0
		#print x, xErr, y, yErr, ratio, ratioErrX, ratioLogNErrX, ratioErrY, ratioLogNErrY
		#print x, xErr, y, yErr, ratio, ratioLogNErrXPlus, ratioLogNErrXMinus
		#h1errFull.SetBinContent( ibin, ratio )
		#h1errFull.SetBinError( ibin, ratioErrX )
		#h1errh2.SetBinContent( ibin, ratio )
		#h1errh2.SetBinError( ibin, ratioErrY )
		ratioList.append( ratio )
		ratioLogNErrXPlusList.append( ratioLogNErrXPlus )
		ratioLogNErrXMinusList.append( ratioLogNErrXMinus )
		'''
		if ibin < 35 :
			if (y>0):
				try: chi2 += ((y-x)*(y-x))/( (yErr*yErr) + (xErr*xErr) )
				except ZeroDivisionError: chi2 += 0
				ndf += 1
		'''

	'''
	histoAsymErrors = histo1.Clone()
	histoAsymErrors.Reset()
	histoAsymErrors.Sumw2(False)
	for ibin in range( histoAsymErrors.GetNbinsX() ):
		histoAsymErrors.SetBinContent( ibin, ratioList[ibin] )
	histoAsymErrors.SetBinErrorOption(TH1.kPoisson)
	'''

	zeroArray = array( 'd', ( [ 0 ] * (len( ratioList )) ) )
	if graphOnly: asymErrors = TGraph( len(ratioList), array('d', binCenterList), array('d', ratioList) ) 
	else: asymErrors = TGraphAsymmErrors( len(ratioList), array('d', binCenterList), array('d', ratioList), zeroArray, zeroArray, array('d',ratioLogNErrXMinusList), array('d', ratioLogNErrXPlusList) )

	#return h1errFull, h1errh2, chi2, ndf-1, asymErrors
	return asymErrors
######################################################################

def ABCDwithTF( histoB, function, fitResults ):
	"""docstring for ABCDwithTF: creates histogram with ABCD prediction using the transfer function """

	listFitValues = []
	listFitErrors = []
	listBinCenter = []
	histoBCD = histoB.Clone()
	histoBCD.Reset()
	histoRatioCD = histoB.Clone()
	histoRatioCD.Reset()

	for ibin in range( 1, histoBCD.GetNbinsX() ):

		contB = histoB.GetBinContent( ibin )
		errorB = histoB.GetBinError( ibin )
		binCenter = histoB.GetBinCenter( ibin )
		factorCD = function.Eval( binCenter )
		contBCD = contB * factorCD

		err = array( 'd', [0] )   ### error in fit
		fitResults.GetConfidenceIntervals( 1, 1, 1, array('d',[binCenter]), err, 0.683, False ) 
		#print ibin, contB, errorB, binCenter, factorCD, contBCD, err[0]
		listFitValues.append( factorCD )
		listFitErrors.append( err[0] )
		listBinCenter.append( binCenter )

		try: errBCD = contBCD* TMath.Sqrt( TMath.Power( err[0]/factorCD, 2 ) + TMath.Power( errorB/contB, 2 ) )
		except ZeroDivisionError: errBCD = 1.8


		histoBCD.SetBinContent( ibin, contBCD )
		histoBCD.SetBinError( ibin, errBCD )
		histoRatioCD.SetBinContent( ibin, factorCD )
		histoRatioCD.SetBinError( ibin, err[0] )

	returnList = [ histoBCD, histoRatioCD, listFitValues, listFitErrors, listBinCenter ]
	return returnList
######################################################################


def ABCDTFunctionCalculation( nameInRoot, binning, minX, maxX, 
				hDataB, hrawDataC, hrawDataD, 
				hqcdmcB, hqcdmcC, hqcdmcD, 
				hDataMinusResonantBkgB, hDataMinusResonantBkgC, hDataMinusResonantBkgD, 
				hDataBkgMinusResonantBkgB, hqcdmcBtagB, 
				typePlot, plot=True, rootFile=False ):
	"""docstring for ernativeABCDcombined: fits ratio B/D and multiply by C"""
	
	############# Rebinning data for transfer function 
	hDataC = rebin( hrawDataC, binning ) 
	hDataD = rebin( hrawDataD, binning ) 
	hDataCD = hDataC.Clone()
	hDataCD.Reset()
	hDataCD.Divide( hDataC, hDataD ) #, 1., 1., '' )
	########################################

	#### data minus resonant bkgs
	hDataMinusResBkgC = rebin( hDataMinusResonantBkgC, binning ) 
	hDataMinusResBkgD = rebin( hDataMinusResonantBkgD, binning ) 
	hDataMinusResBkgCD = hDataMinusResBkgC.Clone()
	hDataMinusResBkgCD.Reset()
	hDataMinusResBkgCD.Divide( hDataMinusResBkgC, hDataMinusResBkgD) #, 1., 1., '' )
	######################################################

	############## QCD MC only
	hqcdmcC = rebin( hqcdmcC, binning ) 
	hqcdmcD = rebin( hqcdmcD, binning ) 
	hqcdmcCD = hqcdmcC.Clone()
	hqcdmcCD.Reset()
	hqcdmcCD.Divide( hqcdmcC, hqcdmcD) #, 1, 1, '' )
	######################################################

	################ Fit 
	#### for prunedMassAsymVsdeltaEtaDijet
	#fitFunction = '1/([0]+TMath::Exp([1]+([2]*x*x)))'  ## v1
	fitFunction = '1/([0]+TMath::Exp([1]+([2]*x*x*x)))'   ## v2
	#fitFunction = '1/([2]+TMath::Exp(-[0]*(x-[1])))'  ## v3
	#fitFunction = 'pol1'

	print ' |----> Fit to QCD MC'
	fitqcdmcCD = TF1( 'fitqcdmcCD', fitFunction, 0, 500 )
	fitqcdmcCD.SetParameter( 0, 2 )
	#fitqcdmcCD.SetParameter( 1, 1 )  #### 1, 150
	#fitqcdmcCD.SetParameter( 2, 0.55 )
	fitqcdmcCDResult = TFitResultPtr( hqcdmcCD.Fit( fitqcdmcCD, 'ELLSR', '', minX, maxX ) )
	fitqcdmcCDResult = TFitResultPtr( hqcdmcCD.Fit( fitqcdmcCD, 'MISR', '', minX, maxX ) )

	print ' |----> Fit to data'
	fitCD = TF1( 'fitCD', fitFunction, 0, 500 )
	for p in range(fitqcdmcCD.GetNpar() ): fitCD.SetParameter( p, fitqcdmcCD.GetParameter( p ) )
	fitCDResult = TFitResultPtr( hDataCD.Fit( fitCD, 'ELLSR', '', minX, maxX ) )
	fitCDResult = TFitResultPtr( hDataCD.Fit( fitCD, 'MISR', '', minX, maxX ) )

	print ' |----> Fit to data minus resonant bkg'
	fitWOResBkgCD = TF1( 'fitWOResBkgCD', fitFunction, 0, 500 )
	for p in range(fitqcdmcCD.GetNpar() ): fitWOResBkgCD.SetParameter( p, fitqcdmcCD.GetParameter( p ) )
	fitWOResBkgCDResult =  TFitResultPtr(hDataMinusResBkgCD.Fit( fitWOResBkgCD, 'ELLSR', '', minX, maxX ) )
	fitWOResBkgCDResult =  TFitResultPtr(hDataMinusResBkgCD.Fit( fitWOResBkgCD, 'MISR', '', minX, maxX ) )
	######################################################
	
	########## Create histograms with prediction
	### all data
	dataABCDwithTFList = ABCDwithTF( hDataB, fitCD, fitCDResult)  
	hDataBCD = dataABCDwithTFList[0]
	hDataRatioCD = dataABCDwithTFList[1]

	### QCD MC
	hqcdmcBCD = ABCDwithTF( hqcdmcB, fitqcdmcCD, fitqcdmcCDResult)[0]  #### everything from MC 
	hqcdmcBCDwTFWOResBkg = ABCDwithTF( hqcdmcB, fitWOResBkgCD, fitWOResBkgCDResult)[0]    ###### TF from data applied to QCD MC

	### Data minus resonant bkg
	dataMinusResBkgABCDwithTFList = ABCDwithTF( hDataMinusResonantBkgB, fitWOResBkgCD, fitWOResBkgCDResult)  
	hDataMinusResBkgBCD = dataMinusResBkgABCDwithTFList[0]
	hDataMinusResBkgRatioCD = dataMinusResBkgABCDwithTFList[1]
	listFitValues = dataMinusResBkgABCDwithTFList[2]
	listFitErrors = dataMinusResBkgABCDwithTFList[3]
	listBinCenter = dataMinusResBkgABCDwithTFList[4] 

	### data minus resonant bkg with btag requirement
	dataBtagMinusResBkgABCDwithTFList = ABCDwithTF( hDataBkgMinusResonantBkgB, fitWOResBkgCD, fitWOResBkgCDResult)  
	hDataBtagMinusResBkgBCD = dataBtagMinusResBkgABCDwithTFList[0]
	hDataBtagMinusResBkgRatioCD = dataBtagMinusResBkgABCDwithTFList[1]
	hqcdmcBtagBCDwTFWOResBkg = ABCDwithTF( hqcdmcBtagB, fitWOResBkgCD, fitWOResBkgCDResult)[0]    ###### TF from data applied to QCD MC
	
	#### QCD MC bkg with TFactor from data
	hQCDMCHybridTFactorBCD = hqcdmcB.Clone()	
	hQCDMCHybridTFactorBCD.Reset()
	htmpDataC = rebin( hrawDataC, 5 ) 
	htmpDataD = rebin( hrawDataD, 5 ) 
	htmpDataCD = htmpDataC.Clone()
	htmpDataCD.Reset()
	htmpDataCD.Divide( htmpDataC, htmpDataD )

	for ibin in range( 1, hQCDMCHybridTFactorBCD.GetNbinsX() ):
		contB = hqcdmcB.GetBinContent( ibin )
		errorB = hqcdmcB.GetBinError( ibin )
		factorCD = htmpDataCD.GetBinContent( ibin )
		errFactorCD = htmpDataCD.GetBinError( ibin )
		contBCD = contB * factorCD
		#print contB, contBCD, factorCD, hDataCD.GetNbinsX(), hqcdmcB.GetNbinsX()
		try: errBCD = contBCD* TMath.Sqrt( TMath.Power( errFactorCD/factorCD, 2 ) + TMath.Power( errorB/contB, 2 ) )
		except ZeroDivisionError: errBCD = 1.8
		hQCDMCHybridTFactorBCD.SetBinContent( ibin, contBCD )
		hQCDMCHybridTFactorBCD.SetBinError( ibin, errBCD )

	#### Create rootfile for limit setting
	if rootFile:
		tmpFileName = dataFileName.replace('V2p4', 'V2p4_'+typePlot+'_ABCDEst')
		tmpFile = TFile(tmpFileName, 'recreate' )
		#hDataBCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_ABCDProj' )
		#hDataBCD.Write()
		#hDataRatioCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_RatioBD' )
		#hDataRatioCD.Write()
		hDataMinusResBkgBCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_DATAMinusResBkg_ABCDProj' )
		hDataMinusResBkgBCD.Write()
		hDataMinusResBkgRatioCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_DATAMinusResBkg_RatioBD' )
		hDataMinusResBkgRatioCD.Write()
		hDataBtagMinusResBkgBCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_'+args.numBtags+'_DATAMinusResBkg_ABCDProj' )
		hDataBtagMinusResBkgBCD.Write()
		hDataBtagMinusResBkgRatioCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_'+args.numBtags+'_DATAMinusResBkg_RatioBD' )
		hDataBtagMinusResBkgRatioCD.Write()
		hDataB.Write()
		tmpFile.Close()

	##### Plot bkg estimation
	if plot:
		fitWOResBkgCDUp = TGraph( len(listFitValues), array( 'd', listBinCenter ), np.add( listFitValues, listFitErrors ) ) 
		fitWOResBkgCDDown = TGraph( len(listFitValues), array( 'd', listBinCenter ), np.subtract( listFitValues, listFitErrors ) ) 

		canCD = TCanvas('canCD', 'canCD',  10, 10, 750, 500 )
		if not args.final: gStyle.SetOptFit(1)

		### qcd mc
		hqcdmcCD.GetYaxis().SetTitle( 'Ratio B/D' )
		hqcdmcCD.GetYaxis().SetTitleOffset(0.75)
		hqcdmcCD.GetXaxis().SetTitle( 'Average '+args.grooming+' jet mass [GeV]' )
		hqcdmcCD.SetStats( True )
		hqcdmcCD.SetLineColor(kGreen+2)
		hqcdmcCD.GetXaxis().SetRangeUser( minX, maxX )
		hqcdmcCD.GetYaxis().SetRangeUser( 0, 1 )
		hqcdmcCD.Draw()

		fitqcdmcCD.SetLineWidth(1)
		fitqcdmcCD.SetLineColor(kGreen+2)
		fitqcdmcCD.Draw("same")

		### all data 
		hDataCD.SetStats( True )
		hDataCD.SetMarkerStyle(20)
		hDataCD.SetMarkerColor(kBlue)
		hDataCD.SetLineColor(kBlue)
		hDataCD.Draw("sames")

		fitCD.SetLineWidth(2)
		fitCD.SetLineColor(kRed)
		fitCD.Draw("same")

		#### data minus resonant bkg
		hDataMinusResBkgCD.SetStats( True )
		hDataMinusResBkgCD.SetMarkerStyle(21)
		hDataMinusResBkgCD.SetMarkerColor(kRed)
		hDataMinusResBkgCD.Draw("sames")

		fitWOResBkgCD.SetLineWidth(2)
		fitWOResBkgCD.SetLineColor(kBlue)
		fitWOResBkgCD.Draw("same")

		fitWOResBkgCDUp.SetLineColor(kRed)
		fitWOResBkgCDUp.SetLineStyle(2)
		fitWOResBkgCDUp.Draw('same pc')
		fitWOResBkgCDDown.SetLineStyle(2)
		fitWOResBkgCDDown.Draw('same pc')
		fitWOResBkgCDDown.SetLineColor(kRed)

		CMS_lumi.extraText = "Preliminary"
		CMS_lumi.relPosX = 0.13
		CMS_lumi.CMS_lumi(canCD, 4, 0)

		if args.final: 
			legend=TLegend(0.15,0.65,0.55,0.85)
			legend.SetTextSize(0.04)
		else: 
			legend=TLegend(0.50,0.15,0.95,0.35)
			legend.SetTextSize(0.035)
		legend.SetFillStyle(0)
		legend.AddEntry( hDataCD, 'Data (uncorrected)', 'pl' )
		legend.AddEntry( hDataMinusResBkgCD, 'Data minus res. bkg.', 'pl' )
		legend.AddEntry( hqcdmcCD, 'QCD MC', 'lep' )
		legend.AddEntry( fitCD, 'Fit to data minus res. bkg.', 'l' )
		legend.AddEntry( fitWOResBkgCDUp, 'Fit unc. to data minus res. bkg.', 'l' )
		legend.AddEntry( fitqcdmcCD, 'Fit to MC', 'l' )
		legend.Draw("same")

		if not args.final:
			canCD.Update()
			st2 = hqcdmcCD.GetListOfFunctions().FindObject("stats")
			st2.SetX1NDC(.12)
			st2.SetX2NDC(.32)
			st2.SetY1NDC(.76)
			st2.SetY2NDC(.91)
			st2.SetTextColor(kGreen+2)
			st1 = hDataCD.GetListOfFunctions().FindObject("stats")
			st1.SetX1NDC(.32)
			st1.SetX2NDC(.52)
			st1.SetY1NDC(.76)
			st1.SetY2NDC(.91)
			st1.SetTextColor(kBlue)
			st3 = hDataMinusResBkgCD.GetListOfFunctions().FindObject("stats")
			st3.SetX1NDC(.52)
			st3.SetX2NDC(.72)
			st3.SetY1NDC(.76)
			st3.SetY2NDC(.91)
			st3.SetTextColor(kRed)
			canCD.Modified()

		outputFileNameCD = nameInRoot+'_'+typePlot+'_CD_'+args.grooming+'_QCD'+args.qcd+'_bkgEstimationPlots'+args.version+'.'+args.extension
		canCD.SaveAs('Plots/'+outputFileNameCD)

	
	return hDataBCD, hqcdmcBCD, hqcdmcBCDwTFWOResBkg, hDataMinusResBkgBCD, hQCDMCHybridTFactorBCD, hDataBtagMinusResBkgBCD, hqcdmcBtagBCDwTFWOResBkg
######################################################################


def bkgEstimation( dataFile, bkgFiles, signalFiles, cutFinal, xmin, xmax, rebinX, cutTop, rebinTopX, labX, labY, log, Norm=False ):
	"""docstring for bkgEstimation"""

	nameInRoot = 'massAve_'+cutFinal
	##### Opening Bkg histograms, rebin, clone unbinned histo
	MCBkgHistos = OrderedDict()
	unbinnedMCBkgHistos = OrderedDict()
	for bkgSamples in bkgFiles:
		bkgNameHisto = ( nameInRoot+'_'+bkgSamples if args.miniTree else 'BoostedAnalysisPlots/'+nameInRoot )
		scale = bkgFiles[ bkgSamples ][1] 
		for side in [ '_A', '_B', '_C', '_D' ]:
			MCBkgHistos[ bkgSamples+side ] = bkgFiles[ bkgSamples ][0].Get( bkgNameHisto+side )
			unbinnedMCBkgHistos[ bkgSamples+side ] = MCBkgHistos[ bkgSamples+side ].Clone()
			MCBkgHistos[ bkgSamples+side ] = rebin( MCBkgHistos[ bkgSamples+side ], ( rebinX if 'simple' in args.binning else args.binning ) )
			MCBkgHistos[ bkgSamples+side ].Scale( scale * twoProngSF * antiTau32SF )
			unbinnedMCBkgHistos[ bkgSamples+side ].Scale( scale * twoProngSF * antiTau32SF )

			MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ] = bkgFiles[ bkgSamples ][0].Get( bkgNameHisto+'_'+args.numBtags+side )
			unbinnedMCBkgHistos[ bkgSamples+'_'+args.numBtags+side ] = MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ].Clone()
			MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ] = rebin( MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ], ( rebinX if 'simple' in args.binning else args.binning ) )
			MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ].Scale( scale * twoProngSF * antiTau32SF )
			unbinnedMCBkgHistos[ bkgSamples+'_'+args.numBtags+side ].Scale( scale * twoProngSF * antiTau32SF )
		if 'simple' in args.binning: 
			MCBkgHistos[ bkgSamples+'_A' ].SetFillColor( bkgFiles[ bkgSamples ][3] )
			MCBkgHistos[ bkgSamples+'_'+args.numBtags+'_A' ].SetFillColor( bkgFiles[ bkgSamples ][3] )

		MCBkgHistos[ bkgSamples+'_'+cutTop ] = bkgFiles[ bkgSamples ][0].Get( 'massAve_'+cutTop+'_'+bkgSamples )
		MCBkgHistos[ bkgSamples+'_'+cutTop ] = rebin( MCBkgHistos[ bkgSamples+'_'+cutTop ], rebinTopX )
		MCBkgHistos[ bkgSamples+'_'+cutTop ].Scale( scale*10 * antiTwoProngSF * tau32SF )  	#### times 10 because data is 100% for top region
		MCBkgHistos[ bkgSamples+'_'+cutTop ].SetFillColor( bkgFiles[ bkgSamples ][3] )
	#######################################
	
	##### Opening signal histograms, rebin, clone unbinned histo
	sigHistos = OrderedDict()
	dummySig=0
	for signalSamples in signalFiles:
		signalNameHisto = ( nameInRoot+'_RPVStopStopToJets_'+args.decay+'_M-'+str(signalSamples) if args.miniTree else 'BoostedAnalysisPlots/'+nameInRoot )
		sigHistos[ signalSamples ] = signalFiles[ signalSamples ][0].Get( signalNameHisto+'_A' )
		sigHistos[ signalSamples ] = rebin( sigHistos[ signalSamples ], rebinX )
		sigHistos[ signalSamples ].Scale( signalFiles[ signalSamples ][1] * twoProngSF * antiTau32SF )
		sigHistos[ signalSamples ].SetLineColor( signalFiles[ signalSamples ][3] )
		sigHistos[ signalSamples ].SetLineWidth( 3 )
		sigHistos[ signalSamples ].SetLineStyle( 2+dummySig )
		sigHistos[ signalSamples ].GetXaxis().SetRangeUser( int(signalSamples)-30, int(signalSamples)+30 )
		dummySig+=8
	#######################################
	
	##### Opening data histograms, rebin, clone unbinned histo
	dataHistos = OrderedDict() 
	unbinnedDataHistos = OrderedDict()
	dataNameHisto = ( nameInRoot+'_JetHT_Run2016' if args.miniTree else 'BoostedAnalysisPlots/'+nameInRoot )

	for a in [ '_A', '_B', '_C', '_D' ]:
		dataHistos[ 'DATA'+a ] = dataFile.Get( dataNameHisto+a )
		dataHistos[ 'DATA'+a ] = rebin( dataHistos[ 'DATA'+a ], ( rebinX if 'simple' in args.binning else args.binning ) )
		dataHistos[ 'DATA_'+args.numBtags+a ] = dataFile.Get( dataNameHisto+'_'+args.numBtags+a )
		dataHistos[ 'DATA_'+args.numBtags+a ] = rebin( dataHistos[ 'DATA_'+args.numBtags+a ], ( rebinX if 'simple' in args.binning else args.binning ) )
		unbinnedDataHistos[ 'DATA'+a ] = dataFile.Get( dataNameHisto+a )
		unbinnedDataHistos[ 'DATA_'+args.numBtags+a ] = dataFile.Get( dataNameHisto+'_'+args.numBtags+a )

	dataHistos[ 'DATA_'+cutTop ] = dataFile.Get( 'massAve_'+cutTop+'_JetHT_Run2016' )
	dataHistos[ 'DATA_'+cutTop ] = rebin( dataHistos[ 'DATA_'+cutTop ], rebinTopX  )
	#######################################

	#### Top region Bkg estimation
	print '|---> bkg Estimation top region'
	hBkgsMinusTTbarTopRegion = dataHistos[ 'DATA_'+cutTop ].Clone()  #### all bkgs for Top region
	hBkgsMinusTTbarTopRegion.Reset()

	for isamples in MCBkgHistos:
		if cutTop in isamples:
			if not 'TT' in isamples: hBkgsMinusTTbarTopRegion.Add( MCBkgHistos[ isamples ].Clone() )

	stackTopRegion = THStack( 'stackTopRegion', 'stackTopRegion' )
	stackTopRegionNames = []
	for isamples in bkgFiles:
		stackTopRegion.Add( MCBkgHistos[ isamples+'_'+cutTop ].Clone() )
		stackTopRegionNames.append( [ MCBkgHistos[ isamples+'_'+cutTop ].Clone(), bkgFiles[ isamples ][2] ] ) 

	hDataMinusBkgTopRegion = dataHistos[ 'DATA_'+cutTop ].Clone()
	hDataMinusBkgTopRegion.Add( hBkgsMinusTTbarTopRegion, -1 )

	ttbarSF = makePlots( 'massAve_'+cutTop, 
			dataHistos[ 'DATA_'+cutTop ], 'DATA', 
			stackTopRegion, 'Bkg prediction', 
			rebinTopX, xmin, xmax, 
			[ hDataMinusBkgTopRegion, MCBkgHistos[ 'TT_'+cutTop ] ], "(Data-bkgs)/MCttbar", 
			'', 'TopRegion_'+cutTop, 
			True, 
			stackHistos=stackTopRegionNames,
			addHisto=MCBkgHistos[ 'TT_'+cutTop ], 
			addUncBand=False,
			topRegion=True,
			doFit=True,
			)

	### Applying SF to ttbar
	for isamples in MCBkgHistos:
		if 'TT' in isamples: MCBkgHistos[ isamples ].Scale( ttbarSF[0] )
		if cutTop in isamples: del MCBkgHistos[ isamples ]
	#######################################

	##### adding MC bkgs 
	allBkgHistos = OrderedDict()
	for k in [ '_A', '_B', '_C', '_D' ]:
		allBkgHistos[ 'allBkg'+k ] = dataHistos[ 'DATA'+k ].Clone()
		allBkgHistos[ 'allBkg_'+args.numBtags+k ] = dataHistos[ 'DATA_'+args.numBtags+k ].Clone()
	for h in allBkgHistos: allBkgHistos[ h ].Reset()

	hBkgMinusResonantBkgA = dataHistos[ 'DATA_A' ].Clone()	
	hBkgMinusResonantBkgA.Reset()
	hDataMinusResonantBkgB = unbinnedDataHistos[ 'DATA_B' ].Clone()
	hDataMinusResonantBkgC = dataHistos[ 'DATA_C' ].Clone()
	hDataMinusResonantBkgD = unbinnedDataHistos[ 'DATA_D' ].Clone()

	hBkgBtagMinusResonantBkgA = dataHistos[ 'DATA_'+args.numBtags+'_A' ].Clone()	
	hBkgBtagMinusResonantBkgA.Reset()
	hDataBtagMinusResonantBkgB = unbinnedDataHistos[ 'DATA_'+args.numBtags+'_B' ].Clone()
	hDataBtagMinusResonantBkgC = dataHistos[ 'DATA_'+args.numBtags+'_C' ].Clone()
	hDataBtagMinusResonantBkgD = unbinnedDataHistos[ 'DATA_'+args.numBtags+'_D' ].Clone()
	
	for isamples in MCBkgHistos:
		if 'btag' in isamples:
			if '_A' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hBkgBtagMinusResonantBkgA.Add( MCBkgHistos[ isamples ].Clone() )

			if '_B' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_B' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hDataBtagMinusResonantBkgB.Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )

			elif '_C' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_C' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hDataBtagMinusResonantBkgC.Add( MCBkgHistos[ isamples ].Clone(), -1 )

			elif '_D' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_D' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hDataBtagMinusResonantBkgD.Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )

		else:
			if '_A' in isamples: 
				allBkgHistos[ 'allBkg_A' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hBkgMinusResonantBkgA.Add( MCBkgHistos[ isamples ].Clone() )

			if '_B' in isamples: 
				allBkgHistos[ 'allBkg_B' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hDataMinusResonantBkgB.Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )

			elif '_C' in isamples: 
				allBkgHistos[ 'allBkg_C' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hDataMinusResonantBkgC.Add( MCBkgHistos[ isamples ].Clone(), -1 )

			elif '_D' in isamples: 
				allBkgHistos[ 'allBkg_D' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: hDataMinusResonantBkgD.Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )

	#######################################

	if 'simple' in args.binning: binWidth = round(allBkgHistos[ 'allBkg_A' ].GetBinWidth(1))
	else: binWidth = '#sigma_{mass}' 

	##### performing simple ABCD, order doesn't matter
	hDataCR = BCDHisto( dataHistos[ 'DATA_B' ].Clone(), dataHistos[ 'DATA_B' ], dataHistos[ 'DATA_C' ], dataHistos[ 'DATA_D' ] )  
	#######################################

	##### Testing BCD regions
	for q in [ 'B', 'C', 'D']:
		makePlots( nameInRoot, 
				allBkgHistos[ 'allBkg_'+q ], 'All MC Bkgs Region '+q, 
				dataHistos[ 'DATA_'+q ], 'DATA Region '+q, 
				binWidth, xmin, xmax, 
				[ dataHistos[ 'DATA_'+q ], allBkgHistos[ 'allBkg_'+q ] ],
				"DATA/MC", '', q+'_'+cutTop, True)
	#########################################################

	#### Make plot simple ABCD
	makePlots( nameInRoot, 
			allBkgHistos[ 'allBkg_A' ], 'All MC Bkgs SR', 
			hDataCR, 'DATA ABCD Pred', 
			binWidth, xmin, xmax, 
			[ allBkgHistos[ 'allBkg_A' ], hDataCR ], 
			"MC SR/ABCD Pred", '', 'DATA_Bkg'+'_'+cutTop)
	###########################################################

	##### Calculating transfer function
	hDataBCD, hQCDMCBCD, hQCDMCHybridTFunctionBCD, hDataWOResBkgBCD, hQCDMCHybridTFactorBCD, hDataBtagWOResBkgBCD,  hQCDMCBtagHybridTFunctionBCD = ABCDTFunctionCalculation( 
			nameInRoot, 
			25, xmin, xmax, 
			dataHistos[ 'DATA_C' ], 
			unbinnedDataHistos[ 'DATA_B' ], 
			unbinnedDataHistos[ 'DATA_D' ], 
			MCBkgHistos[ 'QCD'+args.qcd+'All_C' ].Clone(), 
			unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_B' ].Clone(), 
			unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_D' ].Clone(), 
			hDataMinusResonantBkgC, 
			hDataMinusResonantBkgB, 
			hDataMinusResonantBkgD, 
			hDataBtagMinusResonantBkgC,  				## for btag
			MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_C' ].Clone(), 	## for btag
			'combined'+'_'+cutTop, 
			rootFile=True )

	#### Plot QCD compared with ABCD Hybrid using Transfer function
	hQCDMCHybridTFunctionBCDtmp = hQCDMCHybridTFunctionBCD.Clone()
	makePlots( nameInRoot, 
			MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], 'MC QCD'+args.qcd+'All', 
			hQCDMCHybridTFunctionBCDtmp, 'Hybrid ABCD Pred. TFunction.', 
			binWidth, xmin, xmax, 
			[ MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], hQCDMCHybridTFunctionBCDtmp ],
			"MC SR/ABCD Pred", '', 
			'QCD'+args.qcd+'All_Log_HybridTFunction_'+cutTop, 
			True, addUncBand=False )

	makePlots( nameInRoot, 
			MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], 'MC QCD'+args.qcd+'All', 
			hQCDMCHybridTFactorBCD, 'Hybrid ABCD Pred. TFactor.', 
			binWidth, xmin, xmax, 
			[ MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], hQCDMCHybridTFactorBCD ],
			"MC SR/ABCD Pred", '', 
			'QCD'+args.qcd+'All_Log_HybridTFactor_'+cutTop, 
			True, addUncBand=False )
	############################################

	######## FINAL ESTIMATION
	### stack plot
	stackABCD = THStack( 'stackABCD', 'stackABCD' )
	stackABCDNames = []
	for isamples in bkgFiles:
		if not 'QCD' in isamples:
			stackABCD.Add( MCBkgHistos[ isamples+'_A' ].Clone() )
			stackABCDNames.append( [ MCBkgHistos[ isamples+'_A' ].Clone(), bkgFiles[ isamples ][2] ] ) 

	### adding ABCD 
	hABCDOnly = hDataWOResBkgBCD.Clone()
	hABCDOnly.SetFillColor( kBlue )
	stackABCD.Add( hABCDOnly )
	addAllBkg = hABCDOnly.Clone()
	addAllBkg.Add( hBkgMinusResonantBkgA )
	stackABCDNames.append( [ hABCDOnly, '#splitline{QCD multijets}{from ABCD}' ] )

	### adding systematics
	hABCDOnlySys = addSysBand( hABCDOnly, 1.10, kBlack, additionalSys=dataHistos['DATA_C'])
	httbarOnlySys = addSysBand( MCBkgHistos[ 'TT_A' ].Clone(), 1.10, kBlack ) 		##### TT systematics
	hABCDOnlySys.Add( httbarOnlySys )
	hwjetsOnlySys = addSysBand( MCBkgHistos[ 'WJetsToQQ_A' ].Clone(), 1.10, kBlack )
	hABCDOnlySys.Add( hwjetsOnlySys )
	hzjetsOnlySys = addSysBand( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone(), 1.10, kBlack )
	hABCDOnlySys.Add( hzjetsOnlySys )
	hdibosonsOnlySys = addSysBand( MCBkgHistos[ 'Dibosons_A' ].Clone(), 1.10, kBlack )
	hABCDOnlySys.Add( hdibosonsOnlySys )
	hABCDOnlySys.SetLineColor( kBlue )
	hABCDOnlySys.SetLineWidth( 2 )
	stackABCDNames.append( [ hABCDOnlySys, 'Bkg. unc.' ] )
	### systematics in ratio plot
	hRatioSys = hABCDOnlySys.Clone()
	hRatioSys.Reset()
	for ibin in range( 0, hRatioSys.GetNbinsX()+1 ): 
		hRatioSys.SetBinContent( ibin, 1 )
		try: ratioErr = hABCDOnlySys.GetBinError(ibin) / hABCDOnlySys.GetBinContent(ibin)  
		except ZeroDivisionError: ratioErr=0
		hRatioSys.SetBinError( ibin, ratioErr )

	makePlots( nameInRoot, 
			dataHistos[ 'DATA_A' ], 'DATA', 
			stackABCD, 'Bkg prediction', 
			binWidth, xmin, xmax, 
			[dataHistos[ 'DATA_A' ], addAllBkg], "Data/Bkg", 
			'', 'Log_BCDPlusMCbkgs_'+cutTop, 
			True, 
			addHisto=addAllBkg, 
			stackHistos=stackABCDNames, 
			addUncBand=[ hABCDOnlySys, hRatioSys ], 
			signalHistos=sigHistos)
	###########################################

	##### Full closure test
	stackMCABCD = THStack( 'stackMCABCD', 'stackMCABCD' )
	stackMCABCD.Add( MCBkgHistos[ 'Dibosons_A' ].Clone() )
	stackMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone() )
	stackMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_A' ].Clone() )
	stackMCABCD.Add( MCBkgHistos[ 'TT_A' ].Clone() )

	hQCDMCBCD.SetFillColor( kBlue )
	stackMCABCD.Add( hQCDMCBCD )
	hQCDMCBCD.Add( MCBkgHistos[ 'TT_A' ].Clone() )
	hQCDMCBCD.Add( MCBkgHistos[ 'WJetsToQQ_A' ].Clone() )
	hQCDMCBCD.Add( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone() )
	hQCDMCBCD.Add( MCBkgHistos[ 'Dibosons_A' ].Clone() )

	makePlots( nameInRoot, 
			allBkgHistos[ 'allBkg_A' ], 'All SM Bkg from MC', 
			stackMCABCD , 'MC Bkg prediction', 
			binWidth, xmin, xmax, 
			[ allBkgHistos[ 'allBkg_A' ], hQCDMCBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
			'', 'MCBkg_Log_BCD_'+cutTop, 
			True, 
			addHisto=hQCDMCBCD, 
			stackHistos=[ [ hQCDMCBCD, 'QCD ABCD from MC'], 
				[ MCBkgHistos[ 'TT_A' ].Clone(), 't #bar{t} + Jets' ], 
				[ MCBkgHistos[ 'WJetsToQQ_A' ].Clone(), 'W + Jets'], 
				[ MCBkgHistos[ 'ZJetsToQQ_A' ].Clone(), 'Z + Jets'], 
				[ MCBkgHistos[ 'Dibosons_A' ].Clone(), 'Dibosons' ] ] 
			)
	###########################################

	##### Closure: QCD MC replace by Hybrid ABCD
	stackMCABCD = THStack( 'stackMCABCD', 'stackMCABCD' )
	stackMCABCD.Add( MCBkgHistos[ 'Dibosons_A' ].Clone() )
	stackMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone() )
	stackMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_A' ].Clone() )
	stackMCABCD.Add( MCBkgHistos[ 'TT_A' ].Clone() )

	hQCDMCHybridTFunctionBCD.SetFillColor( kBlue )
	stackMCABCD.Add( hQCDMCHybridTFunctionBCD )
	hQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'TT_A' ].Clone() )
	hQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'WJetsToQQ_A' ].Clone() )
	hQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone() )
	hQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'Dibosons_A' ].Clone() )

	makePlots( nameInRoot, 
			allBkgHistos[ 'allBkg_A' ], 'All SM Bkg from MC', 
			stackMCABCD , 'MC Bkg prediction', 
			binWidth, xmin, xmax, 
			[ allBkgHistos[ 'allBkg_A' ], hQCDMCHybridTFunctionBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
			'', 'HybridBkg_Log_BCD_'+cutTop, 
			True, 
			addHisto=hQCDMCHybridTFunctionBCD, 
			stackHistos=[ [ hQCDMCHybridTFunctionBCD, 'QCD Hybrid ABCD from MC'], 
				[ MCBkgHistos[ 'TT_A' ].Clone(), 't #bar{t} + Jets' ], 
				[ MCBkgHistos[ 'WJetsToQQ_A' ].Clone(), 'W + Jets'], 
				[ MCBkgHistos[ 'ZJetsToQQ_A' ].Clone(), 'Z + Jets'], 
				[ MCBkgHistos[ 'Dibosons_A' ].Clone(), 'Dibosons' ] ] 
			)
	###########################################


	###########################################
	#### Btag bkg estimation
	##### performing simple ABCD, order doesn't matter
	hDataCRBtag = BCDHisto( dataHistos[ 'DATA_'+args.numBtags+'_B' ].Clone(), dataHistos[ 'DATA_'+args.numBtags+'_B' ], dataHistos[ 'DATA_'+args.numBtags+'_C' ], dataHistos[ 'DATA_'+args.numBtags+'_D' ] )  
	#######################################
	
	##### Testing BCD regions
	for q in [ 'B', 'C', 'D']:
		makePlots( nameInRoot, 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_'+q ], 'All MC Bkgs Region '+q, 
				dataHistos[ 'DATA_'+args.numBtags+'_'+q ], 'DATA Region '+q, 
				binWidth, xmin, xmax, 
				[ dataHistos[ 'DATA_'+args.numBtags+'_'+q ], allBkgHistos[ 'allBkg_'+args.numBtags+'_'+q ] ],
				"DATA/MC", '', q+'_'+args.numBtags+'_'+cutTop, True)
	#########################################################

	#### Make plot simple ABCD
	makePlots( nameInRoot, 
			allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All MC Bkgs SR', 
			hDataCRBtag, 'DATA ABCD Pred', 
			binWidth, xmin, xmax, 
			[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hDataCRBtag ], 
			"MC SR/ABCD Pred", '', 'DATA_Bkg'+'_'+args.numBtags+'_'+cutTop)
	###########################################################

	##### Calculating transfer function for btagging
	hBtagDataBCD, hBtagQCDMCBCD, hBtagQCDMCHybridTFunctionBCD, hBtagDataWOResBkgBCD, hBtagQCDMCHybridTFactorBCD, hBtagDataBtagWOResBkgBCD,  hBtagQCDMCBtagHybridTFunctionBCD = ABCDTFunctionCalculation( 
			nameInRoot, 
			25, xmin, xmax, 
			dataHistos[ 'DATA_'+args.numBtags+'_C' ], 
			unbinnedDataHistos[ 'DATA_'+args.numBtags+'_B' ], 
			unbinnedDataHistos[ 'DATA_'+args.numBtags+'_D' ], 
			MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_C' ].Clone(), 
			unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_B' ].Clone(), 
			unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_D' ].Clone(), 
			hDataBtagMinusResonantBkgC, 
			hDataBtagMinusResonantBkgB, 
			hDataBtagMinusResonantBkgD, 
			hDataBtagMinusResonantBkgC,  				## for btag
			MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_C' ].Clone(), 	## for btag
			'combined'+args.numBtags+'_'+cutTop, 
			rootFile=True )

	#### Plot QCD compared with ABCD Hybrid using Transfer function
	makePlots( nameInRoot, 
			MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], 'MC QCD'+args.qcd+'All', 
			hBtagQCDMCHybridTFunctionBCD, 'Hybrid ABCD Pred. TFunction.', 
			binWidth, xmin, xmax, 
			[ MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], hBtagQCDMCHybridTFunctionBCD ],
			"MC SR/ABCD Pred", '', 
			'QCD'+args.qcd+'All_'+args.numBtags+'_Log_HybridTFunction_'+cutTop, 
			True, addUncBand=False )

	makePlots( nameInRoot, 
			MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], 'MC QCD'+args.qcd+'All', 
			hBtagQCDMCHybridTFactorBCD, 'Hybrid ABCD Pred. TFactor.', 
			binWidth, xmin, xmax, 
			[ MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], hBtagQCDMCHybridTFactorBCD ],
			"MC SR/ABCD Pred", '', 
			'QCD'+args.qcd+'All_'+args.numBtags+'_Log_HybridTFactor_'+cutTop, 
			True, addUncBand=False )
	############################################

	#### FINAL ESTIMATION WITH BTAGGING using TFunction 
	### stack plot
	stackBtagABCD = THStack( 'stackBtagABCD', 'stackBtagABCD' )
	stackBtagABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
	stackBtagABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

	### adding ABCD 
	hBtagABCDOnly = hBtagDataWOResBkgBCD.Clone()
	hBtagABCDOnly.SetFillColor( kBlue )
	stackBtagABCD.Add( hBtagABCDOnly )
	addAllBkgBtag = hBtagABCDOnly.Clone()
	addAllBkgBtag.Add( hBkgBtagMinusResonantBkgA )

	### adding systematics
	hBtagABCDOnlySys = addSysBand( hBtagABCDOnly, 1.10, kBlack, additionalSys=dataHistos['DATA_'+args.numBtags+'_C'])
	hBtagttbarOnlySys = addSysBand( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack ) 		##### TT systematics
	hBtagABCDOnlySys.Add( hBtagttbarOnlySys )
	hBtagwjetsOnlySys = addSysBand( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
	hBtagABCDOnlySys.Add( hBtagwjetsOnlySys )
	hBtagzjetsOnlySys = addSysBand( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
	hBtagABCDOnlySys.Add( hBtagzjetsOnlySys )
	hBtagdibosonsOnlySys = addSysBand( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
	hBtagABCDOnlySys.Add( hBtagdibosonsOnlySys )
	hBtagABCDOnlySys.SetLineColor( kBlue )
	hBtagABCDOnlySys.SetLineWidth( 2 )
	### systematics in ratio plot
	hBtagRatioSys = hBtagABCDOnlySys.Clone()
	hBtagRatioSys.Reset()
	for ibin in range( 0, hBtagRatioSys.GetNbinsX()+1 ): 
		hBtagRatioSys.SetBinContent( ibin, 1 )
		try: ratioErr = hBtagABCDOnlySys.GetBinError(ibin) / hBtagABCDOnlySys.GetBinContent(ibin)  
		except ZeroDivisionError: ratioErr=0
		hBtagRatioSys.SetBinError( ibin, ratioErr )

	makePlots( nameInRoot, 
			dataHistos[ 'DATA_'+args.numBtags+'_A' ], 'DATA', 
			stackBtagABCD, 'Bkg prediction', 
			binWidth, xmin, xmax, 
			[dataHistos[ 'DATA_'+args.numBtags+'_A' ], addAllBkgBtag], "Data/Bkg", 
			'', args.numBtags+'_Log_BCDPlusMCbkgs_'+cutTop, 
			True, 
			addHisto=addAllBkgBtag, 
			stackHistos=[ [ hBtagABCDOnly, 'QCD from ABCD'], 
				[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't #bar{t} + Jets' ], 
				[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'W + Jets'], 
				[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
				[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ], 
				[ hBtagABCDOnlySys, 'Bkg. uncertainty' ] ], 
			addUncBand=[ hBtagABCDOnlySys, hBtagRatioSys ], 
			signalHistos=sigHistos)
	###########################################

	##### Full closure test btagging
	stackBtagMCABCD = THStack( 'stackBtagMCABCD', 'stackBtagMCABCD' )
	stackBtagMCABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
	stackBtagMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagMCABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

	hBtagQCDMCBCD.SetFillColor( kBlue )
	stackBtagMCABCD.Add( hBtagQCDMCBCD )
	hBtagQCDMCBCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )
	hBtagQCDMCBCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	hBtagQCDMCBCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	hBtagQCDMCBCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )

	makePlots( nameInRoot, 
			allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All SM Bkg from MC', 
			stackBtagMCABCD , 'MC Bkg prediction', 
			binWidth, xmin, xmax, 
			[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hBtagQCDMCBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
			'', args.numBtags+'_MCBkg_Log_BCD_'+cutTop, 
			True, 
			addHisto=hBtagQCDMCBCD, 
			stackHistos=[ [ hBtagQCDMCBCD, 'QCD ABCD from MC'], 
				[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't #bar{t} + Jets' ], 
				[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'W + Jets'], 
				[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
				[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ] ] 
			)
	###########################################

	##### Closure: QCD MC replace by Hybrid ABCD for btagging
	stackBtagHybridMCABCD = THStack( 'stackBtagHybridMCABCD', 'stackBtagHybridMCABCD' )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

	hBtagQCDMCHybridTFunctionBCD.SetFillColor( kBlue )
	stackBtagHybridMCABCD.Add( hBtagQCDMCHybridTFunctionBCD )
	hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )
	hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )

	makePlots( nameInRoot, 
			allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All SM Bkg from MC', 
			stackBtagHybridMCABCD , 'MC Bkg prediction', 
			binWidth, xmin, xmax, 
			[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hBtagQCDMCHybridTFunctionBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
			'', args.numBtags+'_HybridBkg_Log_BCD_'+cutTop, 
			True, 
			addHisto=hBtagQCDMCHybridTFunctionBCD, 
			stackHistos=[ [ hBtagQCDMCHybridTFunctionBCD, 'QCD Hybrid ABCD from MC'], 
				[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't #bar{t} + Jets' ], 
				[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'W + Jets'], 
				[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
				[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ] ] 
			)
	###########################################

	
	##### Closure for btag selection: QCD MC replace by Hybrid ABCD FROM INCLUSIVE
	stackBtagHybridMCABCD = THStack( 'stackBtagHybridMCABCD', 'stackBtagHybridMCABCD' )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagHybridMCABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

	hQCDMCBtagHybridTFunctionBCD.SetFillColor( kBlue )
	stackBtagHybridMCABCD.Add( hQCDMCBtagHybridTFunctionBCD )
	hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )
	hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )

	makePlots( nameInRoot, 
			allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All SM Bkg from MC', 
			stackBtagHybridMCABCD , 'MC Bkg prediction', 
			binWidth, xmin, xmax, 
			[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hQCDMCBtagHybridTFunctionBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
			'', args.numBtags+'_HybridInclBkg_Log_BCD_'+cutTop, 
			True, 
			addHisto=hQCDMCBtagHybridTFunctionBCD, 
			stackHistos=[ [ hQCDMCBtagHybridTFunctionBCD, 'QCD Hybrid ABCD incl. from MC'], 
				[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't #bar{t} + Jets' ], 
				[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'W + Jets'], 
				[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
				[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ] ]
			#addUncBand=[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ]], 
			)
	###########################################
	
	#### FINAL ESTIMATION WITH BTAGGING using TFunction from inclusive
	### stack plot
	stackBtagABCD = THStack( 'stackBtagABCD', 'stackBtagABCD' )
	stackBtagABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
	stackBtagABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
	stackBtagABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

	### adding ABCD 
	hBtagABCDOnly = hDataBtagWOResBkgBCD.Clone()
	hBtagABCDOnly.SetFillColor( kBlue )
	stackBtagABCD.Add( hBtagABCDOnly )
	addAllBkgBtag = hBtagABCDOnly.Clone()
	addAllBkgBtag.Add( hBkgBtagMinusResonantBkgA )

	### adding systematics
	hBtagABCDOnlySys = addSysBand( hBtagABCDOnly, 1.10, kBlack, additionalSys=dataHistos['DATA_'+args.numBtags+'_C'])
	hBtagttbarOnlySys = addSysBand( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack ) 		##### TT systematics
	hBtagABCDOnlySys.Add( hBtagttbarOnlySys )
	hBtagwjetsOnlySys = addSysBand( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
	hBtagABCDOnlySys.Add( hBtagwjetsOnlySys )
	hBtagzjetsOnlySys = addSysBand( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
	hBtagABCDOnlySys.Add( hBtagzjetsOnlySys )
	hBtagdibosonsOnlySys = addSysBand( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
	hBtagABCDOnlySys.Add( hBtagdibosonsOnlySys )
	hBtagABCDOnlySys.SetLineColor( kBlue )
	hBtagABCDOnlySys.SetLineWidth( 2 )
	### systematics in ratio plot
	hBtagRatioSys = hBtagABCDOnlySys.Clone()
	hBtagRatioSys.Reset()
	for ibin in range( 0, hBtagRatioSys.GetNbinsX()+1 ): 
		hBtagRatioSys.SetBinContent( ibin, 1 )
		try: ratioErr = hBtagABCDOnlySys.GetBinError(ibin) / hBtagABCDOnlySys.GetBinContent(ibin)  
		except ZeroDivisionError: ratioErr=0
		hBtagRatioSys.SetBinError( ibin, ratioErr )

	makePlots( nameInRoot, 
			dataHistos[ 'DATA_'+args.numBtags+'_A' ], 'DATA', 
			stackBtagABCD, 'Bkg prediction', 
			binWidth, xmin, xmax, 
			[dataHistos[ 'DATA_'+args.numBtags+'_A' ], addAllBkgBtag], "Data/Bkg", 
			'', args.numBtags+'_Log_Incl_BCDPlusMCbkgs_'+cutTop, 
			True, 
			addHisto=addAllBkgBtag, 
			stackHistos=[ [ hBtagABCDOnly, 'QCD from ABCD Incl.'], 
				[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't #bar{t} + Jets' ], 
				[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'W + Jets'], 
				[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
				[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ], 
				[ hBtagABCDOnlySys, 'Bkg. uncertainty' ] ], 
			addUncBand=[ hBtagABCDOnlySys, hBtagRatioSys ], 
			signalHistos=sigHistos)
	###########################################


##############################################################


def makePlots( nameInRoot, 
		tmphisto1, labelh1, 
		tmphisto2, labelh2, 
		binWidth, xmin, xmax, 
		ratioList, labelRatio, ratio2, 
		typePlot, 
		log=False, 
		addUncBand='', 
		addHisto='', 
		stackHistos='', 
		doFit=False, 
		signalHistos='',
		topRegion=False ):
	"""docstring for makePlots"""

	histo1 = tmphisto1.Clone()
	histo2 = tmphisto2.Clone()
	if 'DATA' in labelh1:  
		legend=TLegend(0.45,0.70,0.95,0.89)
		legend.SetNColumns(2)
	else:
		legend=TLegend(0.15,0.80,0.95,0.89)
		legend.SetNColumns(3)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.035)
	if 'DATA' in labelh1: legend.AddEntry( histo1, labelh1, 'ep' )
	else: legend.AddEntry( histo1, labelh1, 'l' )
	if signalHistos: 
		for sample in signalHistos: legend.AddEntry( signalHistos[sample], 'M_{#tilde{t}} = '+str(sample)+' GeV', 'l' )
	if isinstance(tmphisto2, THStack):
		for sh in stackHistos: legend.AddEntry( sh[0], sh[1], 'f' ) 
	else: legend.AddEntry( histo2, labelh2, 'pl' )

	histo1.GetYaxis().SetTitle('Events / '+(str(int(binWidth)) if 'simple' in args.binning else binWidth )+' GeV')
	histo1.GetYaxis().SetTitleOffset( 0.9 )
	histo1.GetXaxis().SetRangeUser( xmin, xmax )
	if 'MC' in labelh1: 
		histo1.SetLineColor(kRed-4)
		histo1.SetLineWidth(2)
		histo1.SetFillColor(0)
	if not isinstance( histo2, THStack ):
		histo2.GetXaxis().SetRangeUser( xmin, xmax )
		histo2.SetLineColor(kBlue)
		histo2.SetLineWidth(2)
		if 'MC' in labelh2: histo2.SetLineStyle(2)
		else: histo2.SetLineStyle(1)

	tdrStyle.SetPadRightMargin(0.05)
	tdrStyle.SetPadLeftMargin(0.15)
	can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
	gStyle.SetOptStat(0)
	pad1 = TPad("pad1", "Main",0,0.30,1.00,1.00,-1)
	pad2 = TPad("pad2", "Ratio",0,0.00,1.00,0.30,-1);
	pad1.Draw()
	pad2.Draw()

	pad1.cd()
	pad1.SetBottomMargin(0)
	if log: 
		pad1.SetLogy() 	
		if 'Region' in labelh1: histo1.SetMaximum( 2.5* max( histo1.GetMaximum(), histo2.GetMaximum() ) )
		else: histo1.SetMaximum( 15000 )
		histo1.SetMinimum( 0.5 )
	else: 
		pad1.SetGrid()
	if 'DATA' in labelh1: 
		histo1.SetMarkerStyle(8)
		histo1.Draw("PE")
	else: histo1.Draw("histe")

	if isinstance( histo2, THStack ): 
		histo2.Draw('hist same')
		histo1.Draw("histe same")
	else: histo2.Draw('hist E0 same')
	
	if addUncBand and ( len(addUncBand)>0 ): 
		addUncBand[0].SetFillStyle(3005)
		addUncBand[0].Draw("same E2")
		addHisto.SetFillColor( 0 )
		addHisto.Draw("hist same")
		histo1.Draw("PE same")
	if signalHistos: 
		for sample in signalHistos: signalHistos[ sample ].Draw("hist same")

	if isinstance( histo2, THStack ): tmpHisto2 = addHisto
	else: tmpHisto2 = histo2.Clone()
	tmpHisto1 = histo1.Clone()
	try: 
		res = array( 'd', ( [ 0 ] * tmpHisto1.GetNbinsX() ) )
		chi2Ndf =  round( tmpHisto1.Chi2Test(tmpHisto2, 'WWCHI2/NDFP', res), 2 )
		chi2 =  round( tmpHisto1.Chi2Test(tmpHisto2, 'WWCHI2'), 2 )
		chi2Test = TLatex( 0.6, 0.70, '#chi^{2}/ndf Test = '+ str( int(chi2) )+'/'+str( int(chi2/chi2Ndf) ) )
		chi2Test.SetNDC()
		chi2Test.SetTextFont(42) ### 62 is bold, 42 is normal
		chi2Test.SetTextSize(0.035)
		if args.final: chi2Test.Draw()
	except ZeroDivisionError: print ' |---> chi2Test failed. ZeroDivisionError'

	if 'DATA' in labelh1: tmpLabel = 'DATA'
	else: tmpLabel = 'SR'
	numEvents = TLatex( 0.5, 0.75, 'events '+tmpLabel+'/Bkg = '+ str( round( tmpHisto1.Integral(),2 ) )+'/'+str( round( tmpHisto2.Integral(),2 ) ) )
	numEvents.SetNDC()
	numEvents.SetTextFont(42) ### 62 is bold, 42 is normal
	numEvents.SetTextSize(0.035)
	if args.final: numEvents.Draw()

	CMS_lumi.extraText = ("Preliminary" if 'DATA' in labelh1 else "Simulation Preliminary")
	CMS_lumi.relPosX = 0.13
	if topRegion: CMS_lumi.lumi_13TeV = str( round( (args.lumi*10/1000.), 1 ) )+" fb^{-1}"
	CMS_lumi.CMS_lumi(pad1, 4, 0)
	if topRegion: CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 1 ) )+" fb^{-1}"
	legend.Draw()
	pad1.RedrawAxis()

	pad2.cd()
	pad2.SetGrid()
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.3)
	tmpPad2= pad2.DrawFrame(xmin,0.1,xmax,1.9)
	tmpPad2.SetXTitle( 'Average '+args.grooming+' jet mass [GeV]' )
	tmpPad2.SetYTitle( labelRatio )
	tmpPad2.SetTitleSize(0.12, "x")
	tmpPad2.SetTitleSize( (0.09 if 'frac' in labelRatio else 0.12  ), 'y')
	tmpPad2.SetLabelSize(0.10, 'x')
	tmpPad2.SetLabelSize(0.12, 'y')
	tmpPad2.SetTitleOffset( (0.6 if 'frac' in labelRatio else 0.5), 'y')
	tmpPad2.SetNdivisions(505, 'x' )
	tmpPad2.SetNdivisions(505, 'y' )
	pad2.Modified()
	
	if isinstance( ratioList, list ):
		ratio = TGraphAsymmErrors()
		ratio.Divide( ratioList[0], ratioList[1], 'pois' )
	else:
		ratio = ratioList.Clone()
	#labelAxis( nameInRoot, ratio, args.grooming )
	#ratio.GetXaxis().SetRangeUser( xmin, xmax )
	ratio.SetMarkerStyle(8)
	ratio.SetLineColor(kBlack)
	#ratio.SetMaximum( 2. )
	#ratio.SetMinimum( 0. )
	#if isinstance( ratio, TGraphAsymmErrors ): ratio.Draw('p')
	ratio.Draw('P')

	if signalHistos: 
		tmpSignalBkg = {}
		for sample in signalHistos: 
			tmpSignalBkg[ sample ] = signalHistos[sample].Clone()
			tmpSignalBkg[ sample ].SetFillColorAlpha( signalFiles[ sample ][3], 0.2 ) 
			tmpSignalBkg[ sample ].Add( addHisto )
			tmpSignalBkg[ sample ].Divide( addHisto )
			tmpSignalBkg[ sample ].Draw("hist same")
		whiteBox = TGraph(4, array('d', [xmin-10, xmax+10, xmax+10, xmin-10]), array('d', [0, 0, 1, 1]))
		whiteBox.SetFillColor(kWhite)
		whiteBox.Draw("F same")
		pad2.RedrawAxis()
		pad2.RedrawAxis('g')


	if addUncBand and ( len(addUncBand) > 1 ):
		addUncBand[1].SetFillStyle(3005)
		addUncBand[1].Draw( 'same E2' )
		ratio.Draw('same P')
	elif not topRegion:
		line11.Draw("same")
		line09.Draw("same")
		line.Draw("same")

	if doFit:
		tmpFit = TF1( 'tmpFit', 'pol0', 100, 250 )
		ratio.Fit( 'tmpFit', '', '', 100, 250 )
		tmpFit.SetLineColor( kGreen )
		tmpFit.SetLineWidth( 2 )
		tmpFit.Draw("same")
		chi2Test = TLatex( 0.7, 0.8, '#splitline{#chi^{2}/ndF = '+ str( round( tmpFit.GetChisquare(), 2 ) )+'/'+str( int( tmpFit.GetNDF() ) )+'}{p0 = '+ str( round( tmpFit.GetParameter( 0 ), 2 ) ) +' #pm '+str(  round( tmpFit.GetParError( 0 ), 2 ) )+'}' )
		chi2Test.SetNDC()
		chi2Test.SetTextFont(42) ### 62 is bold, 42 is normal
		chi2Test.SetTextSize(0.10)
		chi2Test.Draw('same')

	if isinstance( ratio2, TH1 ):
		ratio2.SetFillStyle(3004)
		ratio2.SetFillColor( kRed )
		ratio2.Draw('same E2')

	outputFileName = nameInRoot+'_'+typePlot+'_'+args.grooming+'_QCD'+args.qcd+'_bkgEstimationPlots'+args.version+'.'+args.extension
	if not 'simple' in args.binning: outputFileName = outputFileName.replace( typePlot, typePlot+'_ResoBasedBin' )
	print 'Processing.......', outputFileName
	can.SaveAs( 'Plots/'+ outputFileName )
	del can

	if topRegion: return [ tmpFit.GetParameter(0), tmpFit.GetParError(0) ]

#######################################################################



def tmpPlotBkgEstimation( dataFile, nameInRoot, xmin, xmax, rebinX, labX, labY, log, Norm=False ):
	"""docstring for tmpPlotBkgEstimation"""

	hData_A =  dataFile.Get( 'BoostedAnalysisPlots/'+nameInRoot+'_A' )
	hData_A.Rebin( rebinX )
	hData_B =  dataFile.Get( 'BoostedAnalysisPlots/'+nameInRoot+'_B' )
	hData_B.Rebin( rebinX )
	hData_C =  dataFile.Get( 'BoostedAnalysisPlots/'+nameInRoot+'_C' )
	hData_C.Rebin( rebinX )
	hData_D =  dataFile.Get( 'BoostedAnalysisPlots/'+nameInRoot+'_D' )
	hData_D.Rebin( rebinX )
	
	tmpBC = hData_B.Clone()
	tmpBC.Reset()
	tmpBC.Multiply( hData_B, hData_C, 1, 1, '')
	tmpBCD = hData_B.Clone()
	tmpBCD.Reset()
	tmpBCD.Divide( tmpBC, hData_D, 1, 1, '')

	tmphSignalCR = hData_A.Clone()
	tmphSignalCR.Reset()
	tmphSignalCR.Divide( hData_A, tmpBCD, 1., 1., '' )

	binWidth = hData_A.GetBinWidth(1)

	legend3=TLegend(0.55,0.75,0.90,0.87)
	legend3.SetFillStyle(0)
	legend3.SetTextSize(0.03)
	legend3.AddEntry( hData_A, 'DATA - SR' , 'l' )
	legend3.AddEntry( tmpBCD, 'DATA - ABCD Pred', 'pl' )

	hData_A.SetLineColor(kRed-4)
	hData_A.SetLineWidth(2)
	hData_A.GetYaxis().SetTitle('Events / '+str(binWidth))
	hData_A.GetXaxis().SetRangeUser( xmin, xmax )
	hData_A.SetMaximum( 1.2* max( hData_A.GetMaximum(), tmpBCD.GetMaximum() ) )
	tmpBCD.SetLineColor(kBlue)
	tmpBCD.SetLineWidth(2)
	tmpBCD.SetLineStyle(2)

	tdrStyle.SetPadRightMargin(0.05)
	tdrStyle.SetPadLeftMargin(0.15)
	can = TCanvas('c2', 'c2',  10, 10, 750, 750 )
	pad1 = TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
	pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
	pad1.Draw()
	pad2.Draw()

	pad1.cd()
	pad1.SetGrid()
	#if log: pad1.SetLogy() 	
	hData_A.Draw("histe")
	tmpBCD.Draw('histe same')

	CMS_lumi.extraText = "Preliminary"
	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(pad1, 4, 0)
	legend3.Draw()
	#if not (labX and labY): labels( name, '', '' )
	#labels( name1, '', '' ) #, labX, labY )

	pad2.cd()
	pad2.SetGrid()
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.3)
	
	labelAxis( nameInRoot, hData_A, args.grooming )
	tmphSignalCR.GetXaxis().SetRangeUser( xmin, xmax )
	tmphSignalCR.SetMarkerStyle(8)
	tmphSignalCR.GetXaxis().SetTitleOffset(1.1)
	tmphSignalCR.GetXaxis().SetLabelSize(0.12)
	tmphSignalCR.GetXaxis().SetTitleSize(0.12)
	tmphSignalCR.GetYaxis().SetTitle("SR/ABCD Pred")
	tmphSignalCR.GetYaxis().SetLabelSize(0.12)
	tmphSignalCR.GetYaxis().SetTitleSize(0.12)
	tmphSignalCR.GetYaxis().SetTitleOffset(0.55)
	tmphSignalCR.SetMaximum( 2. )
	tmphSignalCR.SetMinimum( 0. )
	tmphSignalCR.GetYaxis().SetNdivisions(505)
	tmphSignalCR.Draw()
	line.Draw("same")

	outputFileName = nameInRoot+'_BkgPlusRPVSt'+str(args.mass)+'_'+args.grooming+'_QCD'+args.qcd+'_bkgEstimationPlots'+args.version+'.'+args.extension
	print 'Processing.......', outputFileName
	can.SaveAs( 'Plots/'+ outputFileName )
	del can
##############################################################



def plotSimpleBkgEstimation( rootFile, bkg, nameInRoot, xmin, xmax, rebinX, labX, labY, log, Norm=False ):
	"""docstring for plotSimpleBkgEstimation"""

	outputFileName = nameInRoot+'_'+bkg+'_pruned_bkgEstimationPlots'+args.version+'.'+args.extension
	print 'Processing.......', outputFileName

	bkgHistos = OrderedDict()
	bkgHistos[ nameInRoot+'_'+bkg+'_A' ] = rootFile.Get( nameInRoot+'_'+bkg+'_A' )
	bkgHistos[ nameInRoot+'_'+bkg+'_B' ] = rootFile.Get( nameInRoot+'_'+bkg+'_B' )
	bkgHistos[ nameInRoot+'_'+bkg+'_C' ] = rootFile.Get( nameInRoot+'_'+bkg+'_C' )
	bkgHistos[ nameInRoot+'_'+bkg+'_D' ] = rootFile.Get( nameInRoot+'_'+bkg+'_D' )

	histoBC = bkgHistos[ nameInRoot+'_'+bkg+'_A' ].Clone()
	histoBC.Reset()
	histoBC.Multiply( bkgHistos[ nameInRoot+'_'+bkg+'_B' ], bkgHistos[ nameInRoot+'_'+bkg+'_C' ], 1, 1, '')
	histoBCD = bkgHistos[ nameInRoot+'_'+bkg+'_A' ].Clone()
	histoBCD.Reset()
	histoBCD.Divide( histoBC, bkgHistos[ nameInRoot+'_'+bkg+'_D' ], 1, 1, '')

	hRatiohBkg = ratioPlots( bkgHistos[ nameInRoot+'_'+bkg+'_A' ], histoBCD ) 
	makePlots( nameInRoot, bkgHistos[ nameInRoot+'_'+bkg+'_A' ], bkg+' SR', histoBCD, bkg+' MC ABCD Pred.', 5, xmin, xmax, hRatiohBkg, "MC SR/ABCD Pred", '', bkg+'Bkg_Log', True)
##############################################################




def plot2DBkgEstimation( rootFile, dataFile, sample, nameInRoot, scale, titleXAxis, titleXAxis2, Xmin, Xmax, rebinx, Ymin, Ymax, rebiny, legX, legY ):
	"""docstring for plot"""

	outputFileName = nameInRoot+'_'+sample+'_'+(args.decay+'_' if 'UDD323' in args.decay else '' )+args.grooming+'_bkgEstimationPlots'+args.version+'.'+args.extension
	print 'Processing.......', outputFileName

	bkgHistos = OrderedDict()
	if isinstance(rootFile, dict):
		for bkg in rootFile:
			bkgHistos[ bkg+'_A' ] = Rebin2D( rootFile[ bkg ][0].Get( nameInRoot+'_'+bkg+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_A' ), rebinx, rebiny )
			bkgHistos[ bkg+'_B' ] = Rebin2D( rootFile[ bkg ][0].Get( nameInRoot+'_'+bkg+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_B' ), rebinx, rebiny )
			bkgHistos[ bkg+'_C' ] = Rebin2D( rootFile[ bkg ][0].Get( nameInRoot+'_'+bkg+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_C' ), rebinx, rebiny )
			bkgHistos[ bkg+'_D' ] = Rebin2D( rootFile[ bkg ][0].Get( nameInRoot+'_'+bkg+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_D' ), rebinx, rebiny )

		hBkg = bkgHistos[ bkg+'_B' ].Clone()
		hBkg.Reset()
		for samples in bkgHistos:
			bkgHistos[ samples ].Scale( scale )
			hBkg.Add( bkgHistos[ samples ].Clone() )
	else: 
		bkgHistos[ sample+'_A' ] = Rebin2D( rootFile.Get( nameInRoot+'_'+sample+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_A' ), rebinx, rebiny )
		bkgHistos[ sample+'_B' ] = Rebin2D( rootFile.Get( nameInRoot+'_'+sample+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_B' ), rebinx, rebiny )
		bkgHistos[ sample+'_C' ] = Rebin2D( rootFile.Get( nameInRoot+'_'+sample+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_C' ), rebinx, rebiny )
		bkgHistos[ sample+'_D' ] = Rebin2D( rootFile.Get( nameInRoot+'_'+sample+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_D' ), rebinx, rebiny )

		hBkg = bkgHistos[ sample+'_B' ].Clone()
		hBkg.Reset()
		for samples in bkgHistos: 
			bkgHistos[ samples ].Scale( scale )
			hBkg.Add( bkgHistos[ samples ].Clone() )

	'''
	if isinstance( dataFile, TFile ):

		bkgHistos[ 'DATA_A' ] = Rebin2D( dataFile.Get( nameInRoot+'_DATA_A' ), rebinx, rebiny )
		bkgHistos[ 'DATA_B' ] = Rebin2D( dataFile.Get( nameInRoot+'_DATA_B' ), rebinx, rebiny )
		bkgHistos[ 'DATA_C' ] = Rebin2D( dataFile.Get( nameInRoot+'_DATA_C' ), rebinx, rebiny )
		bkgHistos[ 'DATA_D' ] = Rebin2D( dataFile.Get( nameInRoot+'_DATA_D' ), rebinx, rebiny )

		bkgHistos[ 'DATA_A' ].Add( bkgHistos[ 'TT_A' ], -1 )
		bkgHistos[ 'DATA_B' ].Add( bkgHistos[ 'TT_B' ], -1 )
		bkgHistos[ 'DATA_C' ].Add( bkgHistos[ 'TT_C' ], -1 )
		bkgHistos[ 'DATA_D' ].Add( bkgHistos[ 'TT_D' ], -1 )

		bkgHistos[ 'DATA_A' ].Add( bkgHistos[ 'WJetsToQQ_A' ], -1 )
		bkgHistos[ 'DATA_B' ].Add( bkgHistos[ 'WJetsToQQ_B' ], -1 )
		bkgHistos[ 'DATA_C' ].Add( bkgHistos[ 'WJetsToQQ_C' ], -1 )
		bkgHistos[ 'DATA_D' ].Add( bkgHistos[ 'WJetsToQQ_D' ], -1 )

		hBkg = bkgHistos[ 'DATA_A' ].Clone()
		hBkg.Reset()
		for samples in bkgHistos:
			if '_A' not in samples: 
				print samples
				if 'DATA' in samples: hBkg.Add( bkgHistos[ samples ].Clone() )
	'''

	if 'DATA' in sample: CMS_lumi.extraText = "Preliminary"
	else: CMS_lumi.extraText = "Simulation Preliminary"
	hBkg.GetXaxis().SetTitle( titleXAxis )
	hBkg.GetYaxis().SetTitleOffset( 0.9 )
	hBkg.GetYaxis().SetTitle( titleXAxis2 )
	corrFactor = hBkg.GetCorrelationFactor()
	textBox=TLatex()
	textBox.SetNDC()
	textBox.SetTextSize(0.05) 
	textBox.SetTextFont(62) ### 62 is bold, 42 is normal
	textBox.SetTextAlign(31)

	if (Xmax or Ymax):
		hBkg.GetXaxis().SetRangeUser( Xmin, Xmax )
		hBkg.GetYaxis().SetRangeUser( Ymin, Ymax )

	tdrStyle.SetPadRightMargin(0.12)
	hBkg.SetMaximum( ( 100 if 'UDD323' in args.decay else 1500 ) )
	hBkg.SetMinimum( 0.01 )
	can = TCanvas('c3', 'c3',  750, 500 )
	can.SetLogz()
	hBkg.Draw('colz')
	textBox.DrawLatex(0.85, 0.85, ( 'M_{#tilde{t}} = '+args.mass+' GeV' if 'RPV' in sample else sample ) )
	textBox1 = textBox.Clone()
	textBox1.DrawLatex(0.85, 0.8, 'Corr. Factor = '+str(round(corrFactor,2)))
	textBox2 = textBox.Clone()
	textBox2.SetTextSize(0.12)
	textBox2.DrawLatex(0.27, 0.4, 'B  D')
	textBox3 = textBox.Clone()
	textBox3.SetTextSize(0.12)
	textBox3.DrawLatex(0.27, 0.25, 'A  C')

	xline = array('d', [0,1])
	yline = array('d', [1.5, 1.5])
	line = TGraph(2, xline, yline)
	line.SetLineColor(kBlack)
	line.SetLineWidth(5)
	line.Draw("same")
	xline2 = array('d', [0.1,0.1])
	yline2 = array('d', [0, 5])
	line2 = TGraph(2, xline2, yline2)
	line2.SetLineColor(kBlack)
	line2.SetLineWidth(5)
	line2.Draw("same")

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(can, 4, 0)

	can.SaveAs( 'Plots/'+outputFileName )
	#can.SaveAs( 'Plots/'+outputFileName.replace(''+ext, 'gif') )
	del can


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--proc', action='store', default='1D', help='Process to draw, example: 1D, 2D, MC.' )
	parser.add_argument('-d', '--decay', action='store', default='UDD312', dest='decay', help='Decay, example: UDD312, UDD323.' )
	parser.add_argument('-b', '--binning', action='store', default='simple', help='Binning: resoBased or simple' )
	parser.add_argument('-v', '--version', action='store', default='v05', help='Version: v01, v02.' )
	parser.add_argument('-t', '--miniTree', action='store_true', default=False, help='miniTree: if plots coming from miniTree or RUNAnalysis.' )
	parser.add_argument('-g', '--grooming', action='store', default='pruned', help='Grooming Algorithm, example: Pruned, Filtered.' )
	parser.add_argument('-m', '--mass', action='store', default='100', help='Mass of Stop, example: 100' )
	parser.add_argument('-q', '--qcd', action='store', default='Pt', dest='qcd', help='Type of QCD binning, example: HT.' )
	parser.add_argument('-l', '--lumi', action='store', type=float, default=1787, help='Luminosity, example: 1.' )
	parser.add_argument('-e', '--extension', action='store', default='png', help='Extension of plots.' )
	parser.add_argument('-B', '--numBtags', action='store', default='2btag', help='Number of btags: 1btag or 2btag' )
	parser.add_argument('-f', '--final', action='store_true', default=False, help='Final distributions.' )
	parser.add_argument('-c', '--cutTop', action='store', default='jet1Tau32', help='Cut for ttbar SF.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	
	CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 1 ) )+" fb^{-1}"
	
	if 'Pt' in args.qcd: 
		#bkgLabel='(w QCD pythia8)'
		QCDSF = ( 1 if 'Puppi' in args.grooming else 0.366 ) #0.46 )
	else: 
		#bkgLabel='(w QCD madgraphMLM+pythia8)'
		QCDSF = 1

	twoProngSF = 1 #1.10 
	antiTwoProngSF = 1 #0.86 
	tau32SF = 1# 1.07
	antiTau32SF = 1 #0.57 

	if args.miniTree: filePrefix = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming
	else: filePrefix = 'Rootfiles/RUNAnalysis' 

	bkgFiles = OrderedDict() 
	signalFiles = {}
	signalFiles = {}
	dataFileName = filePrefix+'_JetHT_Run2016_80X_V2p4_'+args.version+'.root'
	dataFile = TFile.Open(dataFileName)

	signalFiles[ args.mass ] = [ TFile.Open(filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-'+args.mass+'_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'M_{#tilde{t}} = '+args.mass+' GeV', kRed]
	if args.final:
		signalFiles[ '180' ] = [ TFile.Open(filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-170_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'M_{#tilde{t}} = 170 GeV', kMagenta]
		signalFiles[ '240' ] = [ TFile.Open(filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-240_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'M_{#tilde{t}} = 240 GeV', 28]
	bkgFiles[ 'Dibosons' ] = [ TFile.Open(filePrefix+'_Dibosons_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'Dibosons', kMagenta+2 ]
    	bkgFiles[ 'ZJetsToQQ' ] = [ TFile.Open(filePrefix+'_ZJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'Z + Jets', kOrange]
    	bkgFiles[ 'WJetsToQQ' ] = [ TFile.Open(filePrefix+'_WJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'W + Jets', 38]
	bkgFiles[ 'TT' ] = [ TFile.Open(filePrefix+'_TT_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 't #bar{t} + Jets', kGreen+2 ]
	bkgFiles[ 'QCD'+args.qcd+'All' ] = [ TFile.Open(filePrefix+'_QCD'+args.qcd+'All_Moriond17_80X_V2p4_'+args.version+'.root'), QCDSF*args.lumi, 'MC QCD multijets', kBlue ]


	massMinX = 60
	massMaxX = 350
	jetMassHTlabY = 0.20
	jetMassHTlabX = 0.85

	if '2D' in args.proc: 
		for signal in signalFiles: 
			plot2DBkgEstimation( 
					signalFiles[ signal ][0], '', 'RPVStopStopToJets_'+args.decay+'_M-'+str(signal), 
					('' if args.miniTree else 'BoostedAnalysisPlots/')+'prunedMassAsymVsdeltaEtaDijet', 
					signalFiles[ signal ][1],
					'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 
					0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		for bkg in bkgFiles: 
			plot2DBkgEstimation( 
					bkgFiles[ bkg ][0], '', bkg, 
					('' if args.miniTree else 'BoostedAnalysisPlots/')+'prunedMassAsymVsdeltaEtaDijet', 
					bkgFiles[ bkg ][1],
					'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 
					0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		plot2DBkgEstimation( 
				dataFile, '', 'JetHT_Run2016', 
				('' if args.miniTree else 'BoostedAnalysisPlots/')+'prunedMassAsymVsdeltaEtaDijet', 
				1,
				'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 
				0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		#plot2DBkgEstimation( bkgFiles, dataFile, 'DATAMinusResBkg', 'prunedMassAsymVsdeltaEtaDijet', 'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		#for bkg in bkgFiles: plot2DBkgEstimation( bkgFiles, bkg, 'prunedMassAsymVsdeltaEtaDijet', 'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		#for bkg in signalFiles: plot2DBkgEstimation( signalFiles[ bkg ][0], 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass), 'prunedMassAsymVsdeltaEtaDijet', 'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)

	elif 'simple' in args.proc:
		for bkg in bkgFiles: 
			plotSimpleBkgEstimation( bkgFiles[ bkg ][0], bkg, 'massAve_jet2Tau21VsprunedMassAsym', massMinX, massMaxX, 5, '', '', False )
			plotSimpleBkgEstimation( bkgFiles[ bkg ][0], bkg, 'massAve_jet2Tau21VsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )
			plotSimpleBkgEstimation( bkgFiles[ bkg ][0], bkg, 'massAve_prunedMassAsymVsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )
		plotSimpleBkgEstimation( dataFile, 'DATA', 'massAve_jet2Tau21VsprunedMassAsym', massMinX, massMaxX, 5, '', '', False )
		plotSimpleBkgEstimation( dataFile, 'DATA', 'massAve_jet2Tau21VsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )
		plotSimpleBkgEstimation( dataFile, 'DATA', 'massAve_prunedMassAsymVsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )

	else: 
		bkgEstimation( dataFile, bkgFiles, signalFiles, 
				args.grooming+'MassAsymVsdeltaEtaDijet', 
				massMinX, massMaxX, 5, 
				args.cutTop, 20,
				'', '', False )
