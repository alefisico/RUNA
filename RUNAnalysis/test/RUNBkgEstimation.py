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

boostedMassAveBins = array( 'd', [ 0, 60, 64, 68, 72, 77, 82, 88, 94, 100, 106, 113, 121, 129, 137, 146, 155, 164, 174, 184, 195, 206, 217, 228, 240, 253, 265, 278, 291, 304, 318, 332, 346, 360, 374, 389, 404, 420, 435, 451, 468, 485, 503, 521, 540, 560, 582, 604, 628, 1000 ] )
#boostedMassAveBins = array( 'd', range(0, 500, 10) )

def listOfCont( histo ):
 	"""docstring for listOfCont"""
	tmpListContent = []
	tmpListError = []
	for ibin in range( histo.GetNbinsX() ): 
		tmpListContent.append( histo.GetBinContent( ibin ) )
		tmpListError.append( histo.GetBinError( ibin ) )
	return tmpListContent, tmpListError
######################################################

def truncatedHisto( histo1, maxCont ):
	"""docstring for truncatedHisto"""
	
	newHisto = histo1.Clone()
	newHisto.Reset()
	
	for ibin in range( 0, histo1.GetNbinsX()+1 ):

		if ( histo1.GetBinLowEdge( ibin ) >= maxCont ): 
			newHisto.SetBinContent( ibin, 0 )
			newHisto.SetBinError( ibin, 0 )
		else:
			newHisto.SetBinContent( ibin, histo1.GetBinContent( ibin ) )
			newHisto.SetBinError( ibin, histo1.GetBinError( ibin ) )

	return newHisto
######################################################

def substractTruncatedHistos( histo1, histo2, maxCont, add=1 ):
	"""docstring for substractTruncatedHistos"""
	
	if ( histo1.GetNbinsX() != histo2.GetNbinsX() ): 
		print 'histo1 and histo2 dont have same binning'
		sys.exit(0)
	else:
		subsHisto = histo1.Clone()
		subsHisto.Reset()
		
		histo1SF = histo1.Integral()/histo1.GetEntries()
		histo2SF = histo2.Integral()/histo2.GetEntries()
		for ibin in range( 0, histo1.GetNbinsX()+1 ):
			contHisto1 = histo1.GetBinContent( ibin ) 
			contErrHisto1 = histo1.GetBinError( ibin )/histo1SF
			contHisto2 = histo2.GetBinContent( ibin ) 
			contErrHisto2 = histo2.GetBinError( ibin )/histo2SF 
			newCont = contHisto1 + ( add*contHisto2 )
			newContErr = TMath.Sqrt( (TMath.Abs( histo1SF ) * TMath.Power( contErrHisto1, 2)) + (TMath.Abs( histo2SF ) * TMath.Power( contErrHisto2, 2))) 
			
			if ( histo1.GetBinLowEdge( ibin ) >= maxCont ): 
				subsHisto.SetBinContent( ibin, contHisto1 )
				subsHisto.SetBinError( ibin, contErrHisto1 )
			else:
				subsHisto.SetBinContent( ibin, newCont )
				subsHisto.SetBinError( ibin, newContErr )

	return subsHisto
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
		#print i, cont, contError, ( addContError if additionalSys else 0), totalErr
	#hclone.SetFillStyle(3018)
	hclone.SetFillColor( colour )
	return hclone
######################################################################


def makePulls( histo1, histo2 ):
	"""docstring for makePulls"""

	histoPulls = histo1.Clone()
	histoPulls.Reset()
	pullsOnly = TH1F( 'pullsOnly'+histo1.GetName()+histo2.GetName(), 'pullsOnly'+histo1.GetName()+histo2.GetName(), 14, -3, 3 )
	for ibin in range(0, histo1.GetNbinsX()+1 ):
		try: pull = ( histo1.GetBinContent( ibin ) - histo2.GetBinContent( ibin ) ) / histo2.GetBinError( ibin ) 
		except ZeroDivisionError: pull = 0
		#try: pull = TMath.Sqrt(2)*TMath.ErfInverse(-1+2* Math.poisson_cdf( int(histo1.GetBinContent( ibin )) , histo2.GetBinContent( ibin )))
		#except ValueError: pull = 0
		#print ibin, histo1.GetBinCenter( ibin ), histo1.GetBinContent( ibin ), histo2.GetBinContent( ibin ), histo1.GetBinError( ibin ), pull
		histoPulls.SetBinContent( ibin, pull )
		#histoPulls.SetBinError( ibin, 1 )
		pullsOnly.Fill( pull )

	return histoPulls, pullsOnly
######################################################################

def rebin( histo, binning ):
	"""docstring for rebin"""

	oldhisto = histo.Clone()
	if 'reso' in str(binning):  
		newhisto = oldhisto.Rebin( len( boostedMassAveBins )-1, oldhisto.GetName(), boostedMassAveBins )
	elif 'ratio' in str(binning): 
		#if 'Btag' in str(binning): diffBins = array( 'd', [ 0, 60, 100, 125, 150, 175, 200, 250, 300, 400, 500 ] )
		#else: diffBins = array( 'd', [ 0, 60, 80, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500 ] )
		diffBins = array( 'd', [ 0, 60, 80, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500 ] )
		newhisto = oldhisto.Rebin( len( diffBins )-1, oldhisto.GetName(), diffBins )
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

def ABCDwithTF( histoB, function, fitResults, increaseUnc=False ):
	"""docstring for ABCDwithTF: creates histogram with ABCD prediction using the transfer function """

	listFitValues = []
	listFitErrors = []
	listBinCenter = []
	histoBCD = histoB.Clone()
	#histoBCD.Reset()
	histoRatioCD = histoB.Clone()
	histoRatioCD.Reset()

	btagUnc = [0.06, 0.045, 0.042, 0.039, 0.036, 0.032, 0.028, 0.024, 0.021, 0.018, 0.015, 0.013, 0.013, 0.013, 0.014, 0.015, 0.015, 0.014, 0.009, -0.0, 0.009, 0.024, 0.039, 0.052, 0.063, 0.07, 0.074, 0.076, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077] 	##### it is calculated at lines ~1071
	#btagUnc = [0.067, 0.066, 0.066, 0.064, 0.063, 0.061, 0.059, 0.057, 0.055, 0.053, 0.05, 0.048, 0.045, 0.041, 0.037, 0.034, 0.03, 0.026, 0.023, 0.02, 0.018, 0.016, 0.014, 0.013, 0.013, 0.013, 0.013, 0.013, 0.014, 0.015, 0.015, 0.015, 0.015, 0.014, 0.012, 0.01, 0.006, 0.002, -0.004, 0.005, 0.011, 0.018, 0.025, 0.032, 0.039, 0.045, 0.051, 0.056, 0.06, 0.064, 0.067, 0.069, 0.071, 0.073, 0.074, 0.075, 0.076, 0.076, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077]  ##### 5 bins
	#btagUnc = [0.067, 0.065, 0.062, 0.058, 0.054, 0.049, 0.043, 0.035, 0.028, 0.022, 0.017, 0.014, 0.013, 0.013, 0.014, 0.015, 0.015, 0.011, 0.004, 0.002, 0.015, 0.029, 0.042, 0.053, 0.062, 0.068, 0.072, 0.074, 0.076, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077]    #### 10 bins
	#btagUnc = [0.066, 0.063, 0.057, 0.05, 0.041, 0.03, 0.02, 0.014, 0.013, 0.014, 0.015, 0.012, 0.002, 0.011, 0.032, 0.051, 0.064, 0.071, 0.075, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077]   ##### 15 bins

	for ibin in range( 1, histoBCD.GetNbinsX() ):

		contB = abs(histoB.GetBinContent( ibin ))
		errorB = histoB.GetBinError( ibin )
		binCenter = histoB.GetBinCenter( ibin )
		factorCD = function.Eval( binCenter )
		contBCD = contB * factorCD

		err = array( 'd', [0] )   ### error in fit
		fitResults.GetConfidenceIntervals( 1, 1, 1, array('d',[binCenter]), err, 0.683, False ) 
		errFactorCD = ( (btagUnc[ibin-1]+err[0] ) if increaseUnc else err[0] )
		#errFactorCD = ( err[0] )
		listFitValues.append( factorCD )
		listFitErrors.append( errFactorCD )
		listBinCenter.append( binCenter )

		try: errBCD = contBCD* TMath.Sqrt( TMath.Power( errFactorCD/factorCD, 2 ) + TMath.Power( errorB/contB, 2 ) )
		except ZeroDivisionError: errBCD = 1.8
		#print ibin, contB, errorB, binCenter, factorCD, errFactorCD, contBCD, errBCD

		histoBCD.SetBinContent( ibin, contBCD )
		histoBCD.SetBinError( ibin, errBCD )
		histoRatioCD.SetBinContent( ibin, factorCD )
		histoRatioCD.SetBinError( ibin, errFactorCD )

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
	if (args.ratioABCD == 'BD'): 
		#if 'Btag' in str(binning): fitFunction = '(1/([0]+TMath::Exp([1]+([2]*x*x*x))))'
		#else: fitFunction = '(1/([0]+TMath::Exp([1]+([2]*x*x)-([3]*x*x*x))))' 
		fitFunction = '(1/([0]+TMath::Exp([1]+([2]*x*x)-([3]*x*x*x))))' 
		#if 'Btag' in str(binning): fitFunction = 'pol6' #'(x<250)*(1/([0]+TMath::Exp([1]+([2]*x*x))))+(x>200)*pol3(3)' 
		#else:
#			fitOne = [ '(1/([0]+TMath::Exp([1]+([2]*x*x)-([3]*x*x*x))))', 60, 250 ]
#			fitTwo = [ 'pol2', 250, 500 ]
			#fitFunction = '(x<250)*(1/([0]+TMath::Exp([1]+([2]*x*x)+([3]*x*x*x))))+(x>250)*pol2(4)' 
			#fitFunction = '1/([0]+TMath::Exp([1]+([2]*x*x)+([3]*x*x*x)))' 
#			fitFunction = 'pol4' 
	else:
		if 'Btag' in str(binning): fitFunction = '(x<150)*(pol3)+(x>150)*(1/([4]+TMath::Exp([5]+([6]*x*x*x))))'
		else: fitFunction = '(x<150)*(pol3)+(x>150)*(1/([4]+TMath::Exp([5]+([6]*x*x*x))))' 

	print ' |----> Fit to QCD MC'
	#fitOneFunction = TF1( 'fitOneFunction', fitOne[0], fitOne[1], fitOne[2] )
	#fitOneFunctionNpar = fitOneFunction.GetNpar()
	#fitTwoFunction = TF1( 'fitTwoFunction', fitTwo[0], fitTwo[1], fitTwo[2] )
	#fitTwoFunctionNpar = fitTwoFunction.GetNpar()
	#fitqcdmcCD = TF1( 'fitqcdmcCD', fitOne[0]+'+'+fitTwo[0]+'('+str(fitOneFunctionNpar)+')', 0, 1000 )
	#hqcdmcCD.Fit( fitOneFunction, 'R' )
	#hqcdmcCD.Fit( fitTwoFunction, 'R+' )
	#for p in range(fitOneFunctionNpar): fitqcdmcCD.SetParameter( p, fitOneFunction.GetParameter( p ) )
	#for p in range(fitTwoFunctionNpar): fitqcdmcCD.SetParameter( (p+fitOneFunctionNpar), fitTwoFunction.GetParameter( p ) )
	fitqcdmcCD = TF1( 'fitqcdmcCD', fitFunction, 0, 1000 )
	fitqcdmcCD.SetParameter( 0, 2 )
	#fitqcdmcCD.SetParameter( 1, -0.5 )
	#fitqcdmcCDResult = TFitResultPtr( hqcdmcCD.Fit( fitqcdmcCD, 'ELLSR', '', minX, maxX ) )
	fitqcdmcCDResult = TFitResultPtr( hqcdmcCD.Fit( fitqcdmcCD, 'MISR', '', minX, maxX ) )

	print ' |----> Fit to data'
	fitCD = TF1( 'fitCD', fitFunction, 0, 1000 )
	#fitCD = TF1( 'fitCD', fitOne[0]+'+'+fitTwo[0]+'('+str(fitOneFunctionNpar)+')', 0, 1000 )
	for p in range(fitqcdmcCD.GetNpar() ): fitCD.SetParameter( p, fitqcdmcCD.GetParameter( p ) )
	#fitCDResult = TFitResultPtr( hDataCD.Fit( fitCD, 'ELLSR', '', minX, maxX ) )
	fitCDResult = TFitResultPtr( hDataCD.Fit( fitCD, '0MISR', '', minX, maxX ) )

	print ' |----> Fit to data minus resonant bkg'
	fitWOResBkgCD = TF1( 'fitWOResBkgCD', fitFunction, 0, 1000 )
	#fitWOResBkgCD = TF1( 'fitWOResBkgCD', fitOne[0]+'+'+fitTwo[0]+'('+str(fitOneFunctionNpar)+')', 0, 1000 )
	for p in range(fitqcdmcCD.GetNpar() ): fitWOResBkgCD.SetParameter( p, fitqcdmcCD.GetParameter( p ) )
	hDataMinusResBkgCD.Fit( fitWOResBkgCD, '', '', minX, maxX ) 
	hDataMinusResBkgCD.Fit( fitWOResBkgCD, 'ELLSR', '', minX, maxX ) 
	fitWOResBkgCDResult =  TFitResultPtr(hDataMinusResBkgCD.Fit( fitWOResBkgCD, 'MISR', '', minX, maxX ) )

	######################################################
	
	########## Create histograms with prediction
	### all data
	dataABCDwithTFList = ABCDwithTF( hDataB, fitCD, fitCDResult)  
	hDataBCD = dataABCDwithTFList[0]
	hDataRatioCD = dataABCDwithTFList[1]

	### QCD MC
	hqcdmcBCDTF = ABCDwithTF( hqcdmcB, fitqcdmcCD, fitqcdmcCDResult)  #### everything from MC 
	hqcdmcBCDwTFWOResBkg = ABCDwithTF( hqcdmcB, fitWOResBkgCD, fitWOResBkgCDResult)[0]    ###### TF from data applied to QCD MC
	hqcdmcBCD = hqcdmcBCDTF[0]
	listQCDFitValues = hqcdmcBCDTF[2]
	listQCDFitErrors = hqcdmcBCDTF[3]
	listQCDBinCenter = hqcdmcBCDTF[4] 

	### Data minus resonant bkg
	print '&'*30, 'ABCD'
	dataMinusResBkgABCDwithTFList = ABCDwithTF( hDataMinusResonantBkgB.Clone(), fitWOResBkgCD, fitWOResBkgCDResult, increaseUnc=( True if 'Btag' in str(binning) else False ))  
	hDataMinusResBkgBCD = dataMinusResBkgABCDwithTFList[0]
	hDataMinusResBkgRatioCD = dataMinusResBkgABCDwithTFList[1]
	listFitValues = dataMinusResBkgABCDwithTFList[2]
	listFitErrors = dataMinusResBkgABCDwithTFList[3]
	listBinCenter = dataMinusResBkgABCDwithTFList[4] 

	### data minus resonant bkg with btag requirement ##### FAKE BTAG
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
		hDataMinusResonantBkgB.SetName('massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_woResBkg_C')
		hDataMinusResonantBkgB.SetTitle('massAve_prunedMassAsymVsdeltaEtaDijet_JetHT_Run2016_woResBkg__C')
		hDataMinusResonantBkgB.Write()
		hDataMinusResBkgBCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_DATAMinusResBkg_ABCDProj' )
		hDataMinusResBkgBCD.SetTitle( 'massAve_prunedMassAsymVsdeltaEtaDijet_DATAMinusResBkg_ABCDProj' )
		hDataMinusResBkgBCD.Write()
		hDataMinusResBkgRatioCD.SetName( 'massAve_prunedMassAsymVsdeltaEtaDijet_DATAMinusResBkg_RatioBD' )
		hDataMinusResBkgRatioCD.SetTitle( 'massAve_prunedMassAsymVsdeltaEtaDijet_DATAMinusResBkg_RatioBD' )
		hDataMinusResBkgRatioCD.Write()
		#hDataB.Write()
		tmpFile.Close()

	##### Plot bkg estimation
	if plot:
		#fitQCDUp = TGraph( len(listQCDFitValues), array( 'd', listQCDBinCenter ), np.add( listQCDFitValues, listQCDFitErrors ) ) 
		#fitQCDDown = TGraph( len(listQCDFitValues), array( 'd', listQCDBinCenter ), np.subtract( listQCDFitValues, listQCDFitErrors ) ) 
		fitWOResBkgCDBandValues = array( 'd', np.add( listFitValues, listFitErrors ) ) +  array( 'd', np.subtract( listFitValues, listFitErrors )[::-1])
		listBinCenterBand = array( 'd', listBinCenter ) +  array( 'd', listBinCenter[::-1])
		fitWOResBkgCDBand = TGraph( len(listBinCenterBand), listBinCenterBand, fitWOResBkgCDBandValues )

		btagUnc = [0.06, 0.045, 0.042, 0.039, 0.036, 0.032, 0.028, 0.024, 0.021, 0.018, 0.015, 0.013, 0.013, 0.013, 0.014, 0.015, 0.015, 0.014, 0.009, -0.0, 0.009, 0.024, 0.039, 0.052, 0.063, 0.07, 0.074, 0.076, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077, 0.077] 	##### it is calculated at lines ~1071
		fitWOResBkgCDBandBtagValues = array( 'd', np.add( np.add( listFitValues, listFitErrors ), btagUnc )) +  array( 'd', np.subtract( np.subtract( listFitValues, listFitErrors ), btagUnc )[::-1])
		fitWOResBkgCDBandBtag = TGraph( len(listBinCenterBand), listBinCenterBand, fitWOResBkgCDBandBtagValues )


		tdrStyle.SetPadLeftMargin(0.12)
		canCD = TCanvas('canCD', 'canCD',  10, 10, 750, 500 )
		if not args.final: gStyle.SetOptFit(1)

		### qcd mc
		'''
		hqcdmcCD.GetYaxis().SetTitle( 'Ratio '+ ('B/D' if (args.ratioABCD=='BD') else 'C/D' ) )
		hqcdmcCD.GetYaxis().SetTitleOffset(0.75)
		hqcdmcCD.GetXaxis().SetTitle( '#bar{m} [GeV]' ) #'Average '+args.grooming+' jet mass [GeV]' )
		hqcdmcCD.SetStats( True )
		hqcdmcCD.SetLineColor(kBlue)
		hqcdmcCD.SetMarkerStyle(21)
		hqcdmcCD.SetMarkerColor(kBlue-2)
		hqcdmcCD.GetXaxis().SetRangeUser( minX, maxX )
		hqcdmcCD.GetYaxis().SetRangeUser( ( 0 if (args.ratioABCD=='BD') else 1 ), ( 1.5 if (args.ratioABCD=='BD') else 6.5 ) )
		hqcdmcCD.Draw()
		hqcdmcCD.GetFunction('fitqcdmcCD').SetLineColor(kBlue)
		hqcdmcCD.GetFunction('fitqcdmcCD').SetLineWidth(2)
		hqcdmcCD.GetFunction('fitqcdmcCD').SetLineStyle(2)
		'''

		#### data minus resonant bkg
		hDataMinusResBkgCD.GetYaxis().SetTitle( 'Ratio '+ ('B/D' if (args.ratioABCD=='BD') else 'C/D' ) )
		hDataMinusResBkgCD.GetYaxis().SetTitleOffset(0.75)
		hDataMinusResBkgCD.GetXaxis().SetTitle( '#bar{m} [GeV]' ) #'Average '+args.grooming+' jet mass [GeV]' )
		hDataMinusResBkgCD.GetXaxis().SetRangeUser( minX, maxX )
		hDataMinusResBkgCD.GetYaxis().SetRangeUser( ( 0 if (args.ratioABCD=='BD') else 1 ), ( 1.5 if (args.ratioABCD=='BD') else 6.5 ) )
		hDataMinusResBkgCD.Draw()

		hDataMinusResBkgCD.SetStats( True )
		hDataMinusResBkgCD.SetMarkerStyle(21)
		hDataMinusResBkgCD.SetMarkerColor(kBlack)
		hDataMinusResBkgCD.SetLineColor(kBlack)
		#hDataMinusResBkgCD.Draw("sames")
		hDataMinusResBkgCD.GetFunction('fitWOResBkgCD').SetLineWidth(2)
		hDataMinusResBkgCD.GetFunction('fitWOResBkgCD').SetLineStyle(2)
		hDataMinusResBkgCD.GetFunction('fitWOResBkgCD').SetLineColor(kBlack)


		fitWOResBkgCDBandBtag.SetFillColorAlpha(kRed+2, 0.7)
		#fitWOResBkgCDBandBtag.SetFillColor(kAzure)
		fitWOResBkgCDBandBtag.Draw('F')

		fitWOResBkgCDBand.SetFillColorAlpha(kGray, 0.9)
		fitWOResBkgCDBand.Draw('F')

#		fitQCDUp.SetLineColor(kGreen+2)
#		fitQCDUp.SetLineStyle(2)
#		fitQCDUp.Draw('same pc')
#		fitQCDDown.SetLineColor(kGreen+2)
#		fitQCDDown.SetLineStyle(2)
#		fitQCDDown.Draw('same pc')

		### all data 
		hDataCD.SetStats( True )
		hDataCD.SetMarkerStyle(20)
		hDataCD.SetMarkerColor(kGreen+2)
		hDataCD.SetLineColor(kGreen+2)
		#hDataCD.Draw("sames")
		hDataMinusResBkgCD.Draw("sames")

		CMS_lumi.extraText = ""#"Preliminary"
		CMS_lumi.relPosX = 0.13
		CMS_lumi.CMS_lumi(canCD, 4, 0)

		if args.final: 
			legend=TLegend(0.15,0.70,0.45,0.87)
			legend.SetTextSize(0.04)
		else: 
			legend=TLegend(0.50,0.15,0.95,0.35)
			legend.SetTextSize(0.035)
		legend.SetFillStyle(0)
		legend.AddEntry( hDataMinusResBkgCD, 'Data (corrected)', 'epl' )
		#legend.AddEntry( hDataCD, 'Data (uncorrected)', 'pl' )
		#legend.AddEntry( hqcdmcCD, 'QCD multijets sim.', 'lep' )
		legend.AddEntry( hDataMinusResBkgCD.GetFunction('fitWOResBkgCD'), 'Fit to data', 'l' )
		legend.AddEntry( fitWOResBkgCDBand, 'Fit unc. to data', 'f' )
		legend.AddEntry( fitWOResBkgCDBandBtag, 'Fit unc. to data minus res. bkg. (b-tagged sel.)', 'F' )
		#legend.AddEntry( hqcdmcCD.GetFunction('fitqcdmcCD'), 'Fit to QCD multijets sim.', 'pl' )
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
			st3.SetX1NDC(.12)
			st3.SetX2NDC(.32)
			st3.SetY1NDC(.62)
			st3.SetY2NDC(.76)
			st3.SetTextColor(kRed)
			canCD.Modified()

		outputFileNameCD = nameInRoot+'_'+( typePlot if not args.runEra else typePlot+'_Run2016'+args.runEra )+'_CD_'+args.grooming+'_QCD'+args.qcd+'_bkgEstimationPlots'+args.version+'.'+args.extension
		canCD.SaveAs('Plots/'+outputFileNameCD)
		#canCD.SaveAs('Plots/'+outputFileNameCD.replace(args.extension, 'C'))

	
	return hDataBCD, hqcdmcBCD, hqcdmcBCDwTFWOResBkg, hDataMinusResBkgBCD, hQCDMCHybridTFactorBCD, hDataBtagMinusResBkgBCD, hqcdmcBtagBCDwTFWOResBkg, hDataMinusResBkgRatioCD 
######################################################################


def bkgEstimation( dataFile, bkgFiles, signalFiles, cutFinal, xmin, xmax, rebinX, rebinTopX, labX, labY, log, Norm=False ):
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
			#MCBkgHistos[ bkgSamples+side ] = truncatedHisto( MCBkgHistos[ bkgSamples+side ], 200 )
			unbinnedMCBkgHistos[ bkgSamples+side ] = MCBkgHistos[ bkgSamples+side ].Clone()
			MCBkgHistos[ bkgSamples+side ] = rebin( MCBkgHistos[ bkgSamples+side ], ( rebinX if 'simple' in args.binning else args.binning ) )
			MCBkgHistos[ bkgSamples+side ].Scale( scale * twoProngSF * antiTau32SF , ( '' if 'simple' in args.binning else 'width' ) )

			unbinnedMCBkgHistos[ bkgSamples+side ].Scale( scale * twoProngSF * antiTau32SF )

			MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ] = bkgFiles[ bkgSamples ][0].Get( bkgNameHisto+'_'+args.numBtags+side )
			#MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ] = truncatedHisto( MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ], 200 )
			unbinnedMCBkgHistos[ bkgSamples+'_'+args.numBtags+side ] = MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ].Clone()
			MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ] = rebin( MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ], ( rebinX if 'simple' in args.binning else args.binning ) )
			MCBkgHistos[ bkgSamples+'_'+args.numBtags+side ].Scale( scale * twoProngSF * antiTau32SF, ( '' if 'simple' in args.binning else 'width' ) )
			unbinnedMCBkgHistos[ bkgSamples+'_'+args.numBtags+side ].Scale( scale * twoProngSF * antiTau32SF )


			#MCBkgHistos[ bkgSamples+'_'+args.cutTop+side ] = bkgFiles[ bkgSamples ][0].Get( bkgNameHisto+'_jet2Tau32'+side )
			#MCBkgHistos[ bkgSamples+'_'+args.cutTop+side ] = rebin( MCBkgHistos[ bkgSamples+'_'+args.cutTop+side ], rebinTopX )
			#MCBkgHistos[ bkgSamples+'_'+args.cutTop+side ].Scale( scale * antiTwoProngSF * tau32SF ) #,( '' if 'simple' in args.binning else 'width' ) )

		MCBkgHistos[ bkgSamples+'_A' ].SetFillColor( bkgFiles[ bkgSamples ][3] )
		MCBkgHistos[ bkgSamples+'_'+args.numBtags+'_A' ].SetFillColor( bkgFiles[ bkgSamples ][3] )

		MCBkgHistos[ bkgSamples+'_'+args.cutTop ] = bkgFiles[ bkgSamples ][0].Get( 'massAve_'+args.cutTop+'_'+bkgSamples )
		MCBkgHistos[ bkgSamples+'_'+args.cutTop ] = rebin( MCBkgHistos[ bkgSamples+'_'+args.cutTop ], rebinTopX )
		MCBkgHistos[ bkgSamples+'_'+args.cutTop ].Scale( scale * antiTwoProngSF * tau32SF ) #,( '' if 'simple' in args.binning else 'width' ) )
		MCBkgHistos[ bkgSamples+'_'+args.cutTop ].SetFillColor( bkgFiles[ bkgSamples ][3] )
	#######################################
	
	##### Opening signal histograms, rebin, clone unbinned histo
	sigHistos = OrderedDict()
	dummySig=0
	for signalSamples in signalFiles:
		signalNameHisto = ( nameInRoot+'_RPVStopStopToJets_'+args.decay+'_M-'+str(signalSamples)+('_'+args.numBtags if 'UDD323' in args.decay else '' ) if args.miniTree else 'BoostedAnalysisPlots/'+nameInRoot )
		sigHistos[ signalSamples ] = signalFiles[ signalSamples ][0].Get( signalNameHisto+'_A' )
		sigHistos[ signalSamples ] = rebin( sigHistos[ signalSamples ], ( rebinX if 'simple' in args.binning else args.binning ) )
		sigHistos[ signalSamples ].Scale( signalFiles[ signalSamples ][1] * twoProngSF * antiTau32SF, ( '' if 'simple' in args.binning else 'width' ) )
		sigHistos[ signalSamples ].SetLineColor( signalFiles[ signalSamples ][3] )
		sigHistos[ signalSamples ].SetLineWidth( 3 )
		sigHistos[ signalSamples ].SetLineStyle( 2+dummySig )
		tmpResolution = 2*(-1.78 + ( 0.1097 * int(signalSamples)) + ( -0.0002897 * int(signalSamples)*int(signalSamples) ) + ( 3.18e-07 * int(signalSamples)*int(signalSamples)*int(signalSamples)))
		sigHistos[ signalSamples ].GetXaxis().SetRangeUser( int(signalSamples)-tmpResolution, int(signalSamples)+tmpResolution )
		dummySig+=8
	#######################################
	
	##### Opening data histograms, rebin, clone unbinned histo
	dataHistos = OrderedDict() 
	unbinnedDataHistos = OrderedDict()
	dataMinusBkgsHistos = OrderedDict() 
	unbinnedDataMinusBkgsHistos = OrderedDict()
	dataNameHisto = ( nameInRoot+'_JetHT_Run2016' if args.miniTree else 'BoostedAnalysisPlots/'+nameInRoot )

	for a in [ '_A', '_B', '_C', '_D' ]:
		dataHistos[ 'DATA'+a ] = dataFile.Get( dataNameHisto+a )
		dataHistos[ 'DATA'+a ] = rebin( dataHistos[ 'DATA'+a ], ( rebinX if 'simple' in args.binning else args.binning ) )
		dataHistos[ 'DATA'+a ].Scale( 1, ( '' if 'simple' in args.binning else 'width' ) )
		dataMinusBkgsHistos[ 'DATA'+a ] = dataHistos[ 'DATA'+a ].Clone()
		dataHistos[ 'DATA_'+args.numBtags+a ] = dataFile.Get( dataNameHisto+'_'+args.numBtags+a )
		dataHistos[ 'DATA_'+args.numBtags+a ] = rebin( dataHistos[ 'DATA_'+args.numBtags+a ], ( rebinX if 'simple' in args.binning else args.binning ) )
		dataHistos[ 'DATA_'+args.numBtags+a ].Scale( 1, ( '' if 'simple' in args.binning else 'width' ) )
		dataMinusBkgsHistos[ 'DATA_'+args.numBtags+a ] = dataHistos[ 'DATA_'+args.numBtags+a ].Clone()

		unbinnedDataHistos[ 'DATA'+a ] = dataFile.Get( dataNameHisto+a )
		unbinnedDataMinusBkgsHistos[ 'DATA'+a ] = unbinnedDataHistos[ 'DATA'+a ].Clone()
		unbinnedDataHistos[ 'DATA_'+args.numBtags+a ] = dataFile.Get( dataNameHisto+'_'+args.numBtags+a )
		unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+a ] = unbinnedDataHistos[ 'DATA_'+args.numBtags+a ].Clone()

		#dataHistos[ 'DATA_'+args.cutTop+a ] = dataFile.Get( dataNameHisto+'_jet2Tau32'+a )
		#dataHistos[ 'DATA_'+args.cutTop+a ] = rebin( dataHistos[ 'DATA_'+args.cutTop+a ], rebinTopX  )

	dataHistos[ 'DATA_'+args.cutTop ] = dataFile.Get( 'massAve_'+args.cutTop+'_JetHT_Run2016' )
	dataHistos[ 'DATA_'+args.cutTop ] = rebin( dataHistos[ 'DATA_'+args.cutTop ], rebinTopX  )
	#dataHistos[ 'DATA_'+args.cutTop ].Scale( 1, ( '' if 'simple' in args.binning else 'width' ) )
	#######################################
	'''
	testFile = TFile('test.root', 'recreate')
	for i in MCBkgHistos: MCBkgHistos[i].Write() 
	for i in unbinnedMCBkgHistos: unbinnedMCBkgHistos[i].Write() 
	for i in dataHistos: dataHistos[i].Write() 
	for i in unbinnedDataHistos: unbinnedDataHistos[i].Write() 
	testFile.Write()
	testFile.Close()
	sys.exit(0)
	'''

	#### Top region Bkg estimation
	print '|---> bkg Estimation top region'
	hBkgsMinusTTbarTopRegion = dataHistos[ 'DATA_'+args.cutTop ].Clone()  #### all bkgs for Top region
	hBkgsMinusTTbarTopRegion.Reset()

	for isamples in MCBkgHistos:
		if args.cutTop in isamples:
			if not 'TT' in isamples: hBkgsMinusTTbarTopRegion.Add( MCBkgHistos[ isamples ].Clone() )

	stackTopRegion = THStack( 'stackTopRegion', 'stackTopRegion' )
	for isamples in bkgFiles:
		stackTopRegion.Add( MCBkgHistos[ isamples+'_'+args.cutTop ].Clone() )

	hDataMinusBkgTopRegion = dataHistos[ 'DATA_'+args.cutTop ].Clone()
	hDataMinusBkgTopRegion.Add( hBkgsMinusTTbarTopRegion, -1 )

	hTopOnlySys = addSysBand( MCBkgHistos[ 'QCDPtAll_'+args.cutTop ].Clone(), 1, kBlack )
	httbarTopOnlySys = addSysBand( MCBkgHistos[ 'TT_'+args.cutTop ].Clone(), 1, kBlack ) 
	hTopOnlySys.Add( httbarTopOnlySys )
	hwjetsTopOnlySys = addSysBand( MCBkgHistos[ 'WJetsToQQ_'+args.cutTop ].Clone(), 1, kBlack )
	hTopOnlySys.Add( hwjetsTopOnlySys )
	hzjetsTopOnlySys = addSysBand( MCBkgHistos[ 'DYJetsToQQ_'+args.cutTop ].Clone(), 1, kBlack )
	hTopOnlySys.Add( hzjetsTopOnlySys )
	hdibosonsTopOnlySys = addSysBand( MCBkgHistos[ 'Dibosons_'+args.cutTop ].Clone(), 1, kBlack )
	hTopOnlySys.Add( hdibosonsTopOnlySys )
	hTopOnlySys.SetLineColor( kBlue )
	hTopOnlySys.SetLineWidth( 2 )

	ttbarSF = makePlots( 'massAve_'+args.cutTop, 
			dataHistos[ 'DATA_'+args.cutTop ], 'Data - top sel.', 
			stackTopRegion, 'Bkg prediction', 
			rebinTopX, 100, 260, 
			[ hDataMinusBkgTopRegion, MCBkgHistos[ 'TT_'+args.cutTop ] ], "(Data-bkgs)/MCttbar", 
			'', 'TopRegion_'+args.cutTop, 
			False, 
			stackHistos=[ 
				[ MCBkgHistos[ 'DYJetsToQQ_'+args.cutTop ].Clone(), 'Z(q#bar{q})+jets'], 
				[ MCBkgHistos[ 'QCDPtAll_'+args.cutTop ].Clone(), 'QCD multijets' ], 
				[ MCBkgHistos[ 'Dibosons_'+args.cutTop ].Clone(), 'Dibosons' ], 
				[ MCBkgHistos[ 'TT_'+args.cutTop ].Clone(), 't#bar{t}+jets' ], 
				[ hTopOnlySys, 'Bkg. uncertainty' ], 
				[ MCBkgHistos[ 'WJetsToQQ_'+args.cutTop ].Clone(), "W(q'#bar{q})+jets"]
				],
			addHisto=MCBkgHistos[ 'TT_'+args.cutTop ], 
			addUncBand=[ hTopOnlySys ],
			topRegion=True,
			doFit=True,
			)
	
	### Applying SF to ttbar
	for isamples in MCBkgHistos:
		#if 'TT' in isamples: MCBkgHistos[ isamples ].Scale( ( ( 0.96 if 'p11' in args.version else 1.01 ) if args.runEra else ttbarSF[0] ) )
		if not 'QCD' in isamples: MCBkgHistos[ isamples ].Scale( 1.03 ) #ttbarSF[0] )

	#for ibin in range( 13, 18 ): print ibin, bkgSamples, round(MCBkgHistos[ 'TT_A' ].GetBinContent( ibin )*(MCBkgHistos[ 'TT_A' ].GetBinWidth( ibin )),3)
	#sys.exit(0)


	#### ttbar estimation with TF #### TEST
	'''
	for a in [ '_A', '_B', '_C', '_D' ]: dataHistos[ 'DATA_'+args.cutTop+a ].Add( MCBkgHistos[ 'TT_'+args.cutTop+a ], -1 ) 
	hDataCRTopRegion = BCDHisto( dataHistos[ 'DATA_'+args.cutTop+'_B' ].Clone(), dataHistos[ 'DATA_'+args.cutTop+'_B' ], dataHistos[ 'DATA_'+args.cutTop+'_C' ], dataHistos[ 'DATA_'+args.cutTop+'_D' ])  ### ABCD estimation
	tmphDataCRTopRegion = hDataCRTopRegion.Clone() #### clone ABCD estimation for ratio
	hDataMinusBkgTopRegionABCD = dataHistos[ 'DATA_'+args.cutTop ].Clone() ### data minus estimation
	hDataMinusBkgTopRegionABCD.Add( tmphDataCRTopRegion, -1 )
	hDataCRTopRegion.Add( MCBkgHistos[ 'TT_'+args.cutTop ].Clone()) ### estimation plus MC ttbar

	makePlots( 'massAve_'+args.cutTop, 
			dataHistos[ 'DATA_'+args.cutTop ], 'Data', 
			hDataCRTopRegion, 'Bkg prediction with ABCD', 
			rebinTopX, 100, 240, 
			[ hDataMinusBkgTopRegionABCD, MCBkgHistos[ 'TT_'+args.cutTop ] ], "(Data-bkgs)/MCttbar", 
			'', 'TopRegionABCD_'+args.cutTop, 
			False
				)
	'''
	#######################################

	##### adding MC bkgs 
	allBkgHistos = OrderedDict()
	for k in [ '_A', '_B', '_C', '_D' ]:
		allBkgHistos[ 'allBkg'+k ] = dataHistos[ 'DATA'+k ].Clone()
		allBkgHistos[ 'allBkg_'+args.numBtags+k ] = dataHistos[ 'DATA_'+args.numBtags+k ].Clone()
	for h in allBkgHistos: allBkgHistos[ h ].Reset()

	dataMinusBkgsHistos[ 'DATA_A' ].Reset()
	htmpDataMinusResonantBkgB = dataHistos[ 'DATA_B' ].Clone()
	htmpDataMinusResonantBkgC = dataHistos[ 'DATA_C' ].Clone()
	htmpDataMinusResonantBkgD = dataHistos[ 'DATA_D' ].Clone()

	dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_A' ].Reset()
	htmpDataBtagMinusResonantBkgB = dataHistos[ 'DATA_'+args.numBtags+'_B' ].Clone()
	htmpDataBtagMinusResonantBkgC = dataHistos[ 'DATA_'+args.numBtags+'_C' ].Clone()
	htmpDataBtagMinusResonantBkgD = dataHistos[ 'DATA_'+args.numBtags+'_D' ].Clone()
	
	for isamples in MCBkgHistos:
		if 'btag' in isamples:
			if '_A' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: 
					dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_A' ].Add( MCBkgHistos[ isamples ].Clone() )

			if '_B' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_B' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: 
					dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_B' ].Add( MCBkgHistos[ isamples ].Clone(), -1 )
					unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_B' ].Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )
					#unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_B' ] = substractTruncatedHistos( unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_B' ], unbinnedMCBkgHistos[ isamples ].Clone(), 200, -1 )
					htmpDataBtagMinusResonantBkgB.Add( MCBkgHistos[ isamples ].Clone(), -1 )

			elif '_C' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_C' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: 
					dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_C' ].Add( MCBkgHistos[ isamples ].Clone(), -1 )
					unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_C' ].Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )
					#unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_C' ] = substractTruncatedHistos( unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_C' ], unbinnedMCBkgHistos[ isamples ].Clone(), 250 )
					htmpDataBtagMinusResonantBkgC.Add( MCBkgHistos[ isamples ].Clone(), -1 )

			elif '_D' in isamples: 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_D' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: 
					dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_D' ].Add( MCBkgHistos[ isamples ].Clone(), -1 )
					unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_D' ].Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )
					#unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_D' ] = substractTruncatedHistos( unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_D' ], unbinnedMCBkgHistos[ isamples ].Clone(), 250 )
					htmpDataBtagMinusResonantBkgD.Add( MCBkgHistos[ isamples ].Clone(), -1 )

		else:
			if '_A' in isamples: 
				allBkgHistos[ 'allBkg_A' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: dataMinusBkgsHistos[ 'DATA_A' ].Add( MCBkgHistos[ isamples ].Clone() )

			if '_B' in isamples: 
				allBkgHistos[ 'allBkg_B' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: 
					dataMinusBkgsHistos[ 'DATA_B' ].Add( MCBkgHistos[ isamples ].Clone(), -1 )
					unbinnedDataMinusBkgsHistos[ 'DATA_B' ].Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )
					#unbinnedDataMinusBkgsHistos[ 'DATA_B' ] = substractTruncatedHistos( unbinnedDataMinusBkgsHistos[ 'DATA_B' ], unbinnedMCBkgHistos[ isamples ].Clone(), 250 )
					htmpDataMinusResonantBkgB.Add( MCBkgHistos[ isamples ].Clone(), -1 )

			elif '_C' in isamples: 
				allBkgHistos[ 'allBkg_C' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: 
					dataMinusBkgsHistos[ 'DATA_C' ].Add( MCBkgHistos[ isamples ].Clone(), -1 )
					unbinnedDataMinusBkgsHistos[ 'DATA_C' ].Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )
					#unbinnedDataMinusBkgsHistos[ 'DATA_C' ] = substractTruncatedHistos( unbinnedDataMinusBkgsHistos[ 'DATA_C' ], unbinnedMCBkgHistos[ isamples ].Clone(), 250 )
					htmpDataMinusResonantBkgC.Add( MCBkgHistos[ isamples ].Clone(), -1 )

			elif '_D' in isamples: 
				allBkgHistos[ 'allBkg_D' ].Add( MCBkgHistos[ isamples ].Clone() )
				if not 'QCD' in isamples: 
					dataMinusBkgsHistos[ 'DATA_D' ].Add( MCBkgHistos[ isamples ].Clone(), -1 )
					unbinnedDataMinusBkgsHistos[ 'DATA_D' ].Add( unbinnedMCBkgHistos[ isamples ].Clone(), -1 )
					#unbinnedDataMinusBkgsHistos[ 'DATA_D' ] = substractTruncatedHistos( unbinnedDataMinusBkgsHistos[ 'DATA_D' ], unbinnedMCBkgHistos[ isamples ].Clone(), 250 )
					htmpDataMinusResonantBkgD.Add( MCBkgHistos[ isamples ].Clone(), -1 )
	#######################################

	if 'simple' in args.binning: binWidth = round(allBkgHistos[ 'allBkg_A' ].GetBinWidth(1))
	else: binWidth = '#sigma_{mass}' 

	##### Calculating transfer function
	hDataBCD, hQCDMCBCD, hQCDMCHybridTFunctionBCD, hDataWOResBkgBCD, hQCDMCHybridTFactorBCD, hDataBtagWOResBkgBCD, hQCDMCBtagHybridTFunctionBCD, hDataTFinclusive = ABCDTFunctionCalculation( 
			nameInRoot, 
			'ratio', 
			60, 500, #xmax, 
			dataHistos[ 'DATA_'+scaleFactorABCD ], 
			unbinnedDataHistos[ 'DATA_'+numeratorABCD ], 
			unbinnedDataHistos[ 'DATA_D' ], 
			MCBkgHistos[ 'QCD'+args.qcd+'All_'+scaleFactorABCD ].Clone(), 
			unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_'+numeratorABCD ].Clone(), 
			unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_D' ].Clone(), 
			dataMinusBkgsHistos[ 'DATA_'+scaleFactorABCD ], 
			unbinnedDataMinusBkgsHistos[ 'DATA_'+numeratorABCD ], 
			unbinnedDataMinusBkgsHistos[ 'DATA_D' ], 
			dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_'+scaleFactorABCD ],  				## for btag
			MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_'+scaleFactorABCD ].Clone(), 	## for btag
			'combined'+args.ratioABCD+'_'+args.cutTop, 
			rootFile=True )

	if 'UDD312' in args.decay:

		if not args.final:
			##### performing simple ABCD, order doesn't matter
			#hDataCR = BCDHisto( dataHistos[ 'DATA_B' ].Clone(), dataHistos[ 'DATA_B' ], dataHistos[ 'DATA_C' ], dataHistos[ 'DATA_D' ] )  
			hDataCR = BCDHisto( htmpDataMinusResonantBkgB.Clone(), htmpDataMinusResonantBkgB, htmpDataMinusResonantBkgC, htmpDataMinusResonantBkgD )  
			'''
			testFile = TFile('test.root', 'recreate')
			htmpDataMinusResonantBkgB.Divide(htmpDataMinusResonantBkgD)
			MCBkgHistos[ 'QCDPtAll_B' ].Divide( MCBkgHistos[ 'QCDPtAll_D' ] )
			MCBkgHistos[ 'QCDPtAll_B' ].Write()
			MCBkgHistos[ 'QCDPtAll_C' ].Write()
			htmpDataMinusResonantBkgC.Write()
			htmpDataMinusResonantBkgB.Write()
			testFile.Write()
			testFile.Close()
			#can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
			#htmpDataMinusResonantBkgB.Draw()
			#can.SaveAs('test.png')
			#sys.exit(0)
			'''
			#######################################

			##### Testing BCD regions
			for q in [ 'B', 'C', 'D']:
			#for q in [ 'C']:
				makePlots( nameInRoot, 
						#htmpDataMinusResonantBkgC, '',
						allBkgHistos[ 'allBkg_'+q ], 'All MC Bkgs Region '+q, 
						dataHistos[ 'DATA_'+q ], 'DATA Region '+q, 
						binWidth, xmin, xmax, 
						[ dataHistos[ 'DATA_'+q ], allBkgHistos[ 'allBkg_'+q ] ],
						"DATA/MC", '', q+'_'+args.cutTop, True)
			#########################################################

			#### Make plot simple ABCD
			makePlots( nameInRoot, 
					allBkgHistos[ 'allBkg_A' ], 'All MC Bkgs SR', 
					hDataCR, 'DATA ABCD Pred', 
					binWidth, xmin, xmax, 
					[ allBkgHistos[ 'allBkg_A' ], hDataCR ], 
					"MC SR/ABCD Pred", '', 'DATA_Bkg'+'_'+args.cutTop
					)
			###########################################################

			#### Plot QCD compared with ABCD Hybrid using Transfer function
			'''
			hQCDMCHybridTFunctionBCDtmp = hQCDMCHybridTFunctionBCD.Clone()
			makePlots( nameInRoot, 
					MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], 'MC QCD'+args.qcd+'All', 
					hQCDMCHybridTFunctionBCDtmp, 'Hybrid ABCD Pred. TFunction.', 
					binWidth, xmin, xmax, 
					[ MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], hQCDMCHybridTFunctionBCDtmp ],
					"MC SR/ABCD Pred", '', 
					'QCD'+args.qcd+'All_Log_HybridTFunction_'+args.cutTop, 
					True, addUncBand=False )

			makePlots( nameInRoot, 
					MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], 'MC QCD'+args.qcd+'All', 
					hQCDMCHybridTFactorBCD, 'Hybrid ABCD Pred. TFactor.', 
					binWidth, xmin, xmax, 
					[ MCBkgHistos[ 'QCD'+args.qcd+'All_A' ], hQCDMCHybridTFactorBCD ],
					"MC SR/ABCD Pred", '', 
					'QCD'+args.qcd+'All_Log_HybridTFactor_'+args.cutTop, 
					True, addUncBand=False )
			'''
			############################################

		######## FINAL ESTIMATION
		### stack plot
		stackABCD = THStack( 'stackABCD', 'stackABCD' )
		for isamples in bkgFiles:
			if not 'QCD' in isamples:
				for ibin in range( 18, 22 ): print ibin, isamples, round(MCBkgHistos[ isamples+'_A' ].GetBinContent( ibin )*(MCBkgHistos[ isamples+'_A' ].GetBinWidth( ibin )),3)
				#print MCBkgHistos[ isamples+'_A' ].GetBinContent( 15 ), MCBkgHistos[ isamples+'_A' ].GetBinCenter( 15 )
				stackABCD.Add( MCBkgHistos[ isamples+'_A' ].Clone() )

		### adding ABCD 
		#hABCDOnly = hDataCR.Clone()  
		#hABCDOnly = hDataBCD.Clone()
		hABCDOnly = hDataWOResBkgBCD.Clone()
		hABCDOnly.SetFillColor( kBlue )
		stackABCD.Add( hABCDOnly )
		addAllBkg = hABCDOnly.Clone()
		addAllBkg.Add( dataMinusBkgsHistos[ 'DATA_A' ] )

		### adding systematics
		hABCDOnlySys = addSysBand( hABCDOnly, 1.10, kBlack) #, additionalSys=dataHistos['DATA_C'])
		hABCDOnlySys.SetName('BkgUnc')
		httbarOnlySys = addSysBand( MCBkgHistos[ 'TT_A' ].Clone(), 1.10, kBlack ) 		##### TT systematics
		hABCDOnlySys.Add( httbarOnlySys )
		hwjetsOnlySys = addSysBand( MCBkgHistos[ 'WJetsToQQ_A' ].Clone(), 1.10, kBlack )
		hABCDOnlySys.Add( hwjetsOnlySys )
		#hzjetsOnlySys = addSysBand( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone(), 1.10, kBlack )
		hzjetsOnlySys = addSysBand( MCBkgHistos[ 'DYJetsToQQ_A' ].Clone(), 1.10, kBlack )
		hABCDOnlySys.Add( hzjetsOnlySys )
		hdibosonsOnlySys = addSysBand( MCBkgHistos[ 'Dibosons_A' ].Clone(), 1.10, kBlack )
		hABCDOnlySys.Add( hdibosonsOnlySys )
		hABCDOnlySys.SetLineColor( kBlue )
		hABCDOnlySys.SetLineWidth( 2 )
		### systematics in ratio plot
		hRatioSys = hABCDOnlySys.Clone()
		hRatioSys.Reset()
		for ibin in range( 0, hRatioSys.GetNbinsX()+1 ): 
			hRatioSys.SetBinContent( ibin, 1 )
			try: ratioErr = hABCDOnlySys.GetBinError(ibin) / hABCDOnlySys.GetBinContent(ibin)  
			except ZeroDivisionError: ratioErr=0
			hRatioSys.SetBinError( ibin, ratioErr )

		for ibin in range( 19, 22 ): print ibin, 'data', dataHistos[ 'DATA_A' ].GetBinContent(ibin)*(dataHistos[ 'DATA_A' ].GetBinWidth(ibin)), 'qcd', round(hABCDOnly.GetBinContent(ibin)*(hABCDOnly.GetBinWidth(ibin)), 3), 'signal', round(sigHistos[args.mass].GetBinContent(ibin)*(sigHistos[args.mass].GetBinWidth(ibin)), 3)
		makePlots( nameInRoot, 
				dataHistos[ 'DATA_A' ], 'Data - inclusive sel.', 
				stackABCD, 'Bkg prediction', 
				binWidth, xmin, xmax, 
				[dataHistos[ 'DATA_A' ], hABCDOnlySys], "Data/Bkg", 
				'', 'Log_BCDratio'+args.ratioABCD+'PlusMCbkgs_'+args.cutTop, 
				True, 
				addHisto=addAllBkg, 
				stackHistos=[ 
					[ MCBkgHistos[ 'Dibosons_A' ].Clone(), 'Dibosons' ], 
					[ hABCDOnly, 'QCD multijets'], 
					[ hABCDOnlySys, 'Bkg. uncertainty' ], 
					[ MCBkgHistos[ 'TT_A' ].Clone(), 't#bar{t}+jets' ], 
					int(args.mass),
					[ MCBkgHistos[ 'WJetsToQQ_A' ].Clone(), "W(q'#bar{q})+jets"], 
					200,
					[ MCBkgHistos[ 'DYJetsToQQ_A' ].Clone(), 'Z(q#bar{q})+jets'], 
					],
				addUncBand=[ hABCDOnlySys, hRatioSys ], 
				signalHistos=sigHistos, 
				plotPull=False)
		###########################################

		##### Full closure test
		stackMCABCD = THStack( 'stackMCABCD', 'stackMCABCD' )
		stackMCABCD.Add( MCBkgHistos[ 'Dibosons_A' ].Clone() )
		#stackMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone() )
		stackMCABCD.Add( MCBkgHistos[ 'DYJetsToQQ_A' ].Clone() )
		stackMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_A' ].Clone() )
		stackMCABCD.Add( MCBkgHistos[ 'TT_A' ].Clone() )

		hQCDMCBCD.SetFillColor( kBlue )
		stackMCABCD.Add( hQCDMCBCD )
		hQCDMCBCD.Add( MCBkgHistos[ 'TT_A' ].Clone() )
		hQCDMCBCD.Add( MCBkgHistos[ 'WJetsToQQ_A' ].Clone() )
		#hQCDMCBCD.Add( MCBkgHistos[ 'ZJetsToQQ_A' ].Clone() )
		hQCDMCBCD.Add( MCBkgHistos[ 'DYJetsToQQ_A' ].Clone() )
		hQCDMCBCD.Add( MCBkgHistos[ 'Dibosons_A' ].Clone() )

		makePlots( nameInRoot, 
				allBkgHistos[ 'allBkg_A' ], 'All SM Bkg from MC', 
				stackMCABCD , 'MC Bkg prediction', 
				binWidth, xmin, xmax, 
				[ allBkgHistos[ 'allBkg_A' ], hQCDMCBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
				'', 'MCBkg_Log_BCDratio'+args.ratioABCD+'_'+args.cutTop, 
				True, 
				addHisto=hQCDMCBCD, 
				stackHistos=[ [ hQCDMCBCD, 'QCD ABCD from MC'], 
					[ MCBkgHistos[ 'TT_A' ].Clone(), 't#bar{t}+jets' ], 
					[ MCBkgHistos[ 'WJetsToQQ_A' ].Clone(), "W(q'#bar{q})+jets"], 
					#[ MCBkgHistos[ 'ZJetsToQQ_A' ].Clone(), 'Z + Jets'], 
					[ MCBkgHistos[ 'DYJetsToQQ_A' ].Clone(), 'Z(q#bar{q})+jets'], 
					[ MCBkgHistos[ 'Dibosons_A' ].Clone(), 'Dibosons' ] ] ,
				#addUncBand= False #[ hABCDOnlySys, hRatioSys ], 
				)
		###########################################

		##### Closure: QCD MC replace by Hybrid ABCD
		'''
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
				'', 'HybridBkg_Log_BCD_'+args.cutTop, 
				True, 
				addHisto=hQCDMCHybridTFunctionBCD, 
				stackHistos=[ [ hQCDMCHybridTFunctionBCD, 'QCD Hybrid ABCD from MC'], 
					[ MCBkgHistos[ 'TT_A' ].Clone(), 't#bar{t}+jets' ], 
					[ MCBkgHistos[ 'WJetsToQQ_A' ].Clone(), "W(q'#bar{q})+jets"], 
					[ MCBkgHistos[ 'ZJetsToQQ_A' ].Clone(), 'Z + Jets'], 
					[ MCBkgHistos[ 'Dibosons_A' ].Clone(), 'Dibosons' ] ] 
				)
		'''
		###########################################

	else: 

		##### Calculating transfer function for btagging
		hBtagDataBCD, hBtagQCDMCBCD, hBtagQCDMCHybridTFunctionBCD, hBtagDataWOResBkgBCD, hBtagQCDMCHybridTFactorBCD, hBtagDataBtagWOResBkgBCD, hBtagQCDMCBtagHybridTFunctionBCD, hDataTFbtagged = ABCDTFunctionCalculation( 
				nameInRoot, 
				'ratioBtag', #50, 
				60, 500, #xmax, 
				dataHistos[ 'DATA_'+args.numBtags+'_'+scaleFactorABCD ], 
				#unbinnedDataHistos[ 'DATA_'+args.numBtags+'_'+numeratorABCD ], 
				#unbinnedDataHistos[ 'DATA_'+args.numBtags+'_D' ], 
				unbinnedDataHistos[ 'DATA_'+numeratorABCD ], 
				unbinnedDataHistos[ 'DATA_D' ], 
				MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_'+scaleFactorABCD ].Clone(), 
				#unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_'+numeratorABCD ].Clone(), 
				#unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_D' ].Clone(), 
				unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_'+numeratorABCD ].Clone(), 
				unbinnedMCBkgHistos[ 'QCD'+args.qcd+'All_D' ].Clone(), 
				dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_'+scaleFactorABCD ], 
				#unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_'+numeratorABCD ], 
				#unbinnedDataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_D' ], 
				unbinnedDataMinusBkgsHistos[ 'DATA_'+numeratorABCD ], 
				unbinnedDataMinusBkgsHistos[ 'DATA_D' ], 
				dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_'+scaleFactorABCD ],  				## for btag
				MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_'+scaleFactorABCD ].Clone(), 	## for btag
				'combined'+args.ratioABCD+args.numBtags+'_'+args.cutTop, 
				rootFile=True )

		'''	
		tdrStyle.SetPadLeftMargin(0.12)
		canCD = TCanvas('canCD', 'canCD',  10, 10, 750, 500 )
		if not args.final: gStyle.SetOptFit(1)

		newhDataTFinclusive = hDataTFinclusive.Clone()
		newhDataTFinclusive.Reset()
		tmpList = []
		for ibin in range( 1, hDataTFinclusive.GetNbinsX() ):
			btagValueUp = hDataTFbtagged.GetBinContent( ibin ) + hDataTFbtagged.GetBinError( ibin )
			btagValueDown = hDataTFbtagged.GetBinContent( ibin ) - hDataTFbtagged.GetBinError( ibin )
			inclValueUp = hDataTFinclusive.GetBinContent( ibin ) + hDataTFinclusive.GetBinError( ibin )
			inclValueDown = hDataTFinclusive.GetBinContent( ibin ) - hDataTFinclusive.GetBinError( ibin )
			if (btagValueUp > inclValueUp) and (btagValueDown > inclValueDown): addToUnc = btagValueUp - inclValueUp
			else: addToUnc = inclValueDown - btagValueDown
			#print ibin, btagValueUp, btagValueDown, inclValueUp, inclValueDown, addToUnc, (addToUnc+hDataTFinclusive.GetBinError( ibin )+hDataTFinclusive.GetBinContent( ibin )), hDataTFinclusive.GetBinContent( ibin )-(addToUnc+hDataTFinclusive.GetBinError( ibin ))
			tmpList.append( round(addToUnc, 3 ))
			newhDataTFinclusive.SetBinContent( ibin, hDataTFinclusive.GetBinContent( ibin ) )
			newhDataTFinclusive.SetBinError( ibin, addToUnc+hDataTFinclusive.GetBinError( ibin ) )
		print tmpList

		newhDataTFinclusive.GetYaxis().SetTitle( 'Ratio '+ ('B/D' if (args.ratioABCD=='BD') else 'C/D' ) )
		newhDataTFinclusive.GetYaxis().SetTitleOffset(0.75)
		newhDataTFinclusive.GetXaxis().SetTitle( '#bar{m} [GeV]' ) 
		newhDataTFinclusive.SetLineColor(kBlue-7)
		newhDataTFinclusive.SetFillColorAlpha(kBlue,0.3)
		newhDataTFinclusive.GetXaxis().SetRangeUser( 60, 500 )
		newhDataTFinclusive.GetYaxis().SetRangeUser( ( 0 if (args.ratioABCD=='BD') else 1 ), ( 1.5 if (args.ratioABCD=='BD') else 6.5 ) )
		newhDataTFinclusive.Draw('E4')

		hDataTFinclusive.SetLineWidth(2)
		hDataTFinclusive.SetLineColor(kBlue)
		hDataTFinclusive.Draw("same E1")

		hDataTFbtagged.SetLineWidth(2)
		hDataTFbtagged.SetMarkerColor(kRed)
		hDataTFbtagged.SetLineColor(kRed)
		hDataTFbtagged.Draw("same E1")

		CMS_lumi.extraText = "Preliminary"
		CMS_lumi.relPosX = 0.13
		CMS_lumi.CMS_lumi(canCD, 4, 0)

		legend=TLegend(0.15,0.60,0.55,0.85)
		legend.SetTextSize(0.04)
		legend.SetFillStyle(0)
		legend.AddEntry( hDataTFinclusive, 'Inclusive selection', 'pl' )
		legend.AddEntry( newhDataTFinclusive, 'New uncertainty in Inclusive TF for b-tagged selection', 'f' )
		legend.AddEntry( hDataTFbtagged, 'b-tagged selection', 'pl' )
		legend.Draw("same")

		outputFileNameCD = nameInRoot+'_BDComparison_'+args.grooming+'_bkgEstimationPlots'+args.version+'.'+args.extension
		canCD.SaveAs('Plots/'+outputFileNameCD)
		sys.exit(0)
		'''

		if not args.final:

			###########################################
			#### Btag bkg estimation
			##### performing simple ABCD, order doesn't matter
			#hDataCRBtag = BCDHisto( dataHistos[ 'DATA_'+args.numBtags+'_B' ].Clone(), dataHistos[ 'DATA_'+args.numBtags+'_B' ], dataHistos[ 'DATA_'+args.numBtags+'_C' ], dataHistos[ 'DATA_'+args.numBtags+'_D' ] )  
			hDataCRBtag = BCDHisto( htmpDataBtagMinusResonantBkgB.Clone(), htmpDataBtagMinusResonantBkgB, htmpDataBtagMinusResonantBkgC, htmpDataBtagMinusResonantBkgD )  
			#######################################
			
			##### Testing BCD regions
			for q in [ 'B', 'C', 'D']:
				makePlots( nameInRoot, 
						allBkgHistos[ 'allBkg_'+args.numBtags+'_'+q ], 'All MC Bkgs Region '+q, 
						dataHistos[ 'DATA_'+args.numBtags+'_'+q ], 'DATA Region '+q, 
						binWidth, xmin, xmax, 
						[ dataHistos[ 'DATA_'+args.numBtags+'_'+q ], allBkgHistos[ 'allBkg_'+args.numBtags+'_'+q ] ],
						"DATA/MC", '', q+'_'+args.numBtags+'_'+args.cutTop, True)
			#########################################################

			#### Make plot simple ABCD
			makePlots( nameInRoot, 
					allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All MC Bkgs SR', 
					hDataCRBtag, 'DATA ABCD Pred', 
					binWidth, xmin, xmax, 
					[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hDataCRBtag ], 
					"MC SR/ABCD Pred", '', 'DATA_Bkg'+'_'+args.numBtags+'_'+args.cutTop)
			###########################################################

			#### Plot QCD compared with ABCD Hybrid using Transfer function
			makePlots( nameInRoot, 
					MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], 'MC QCD'+args.qcd+'All', 
					hBtagQCDMCHybridTFunctionBCD, 'Hybrid ABCD Pred. TFunction.', 
					binWidth, xmin, xmax, 
					[ MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], hBtagQCDMCHybridTFunctionBCD ],
					"MC SR/ABCD Pred", '', 
					'QCD'+args.qcd+'All_'+args.numBtags+'_Log_HybridTFunction_ratio'+args.ratioABCD+'_'+args.cutTop, 
					True, addUncBand=False )

			makePlots( nameInRoot, 
					MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], 'MC QCD'+args.qcd+'All', 
					hBtagQCDMCHybridTFactorBCD, 'Hybrid ABCD Pred. TFactor.', 
					binWidth, xmin, xmax, 
					[ MCBkgHistos[ 'QCD'+args.qcd+'All_'+args.numBtags+'_A' ], hBtagQCDMCHybridTFactorBCD ],
					"MC SR/ABCD Pred", '', 
					'QCD'+args.qcd+'All_'+args.numBtags+'_Log_HybridTFactor_ratio'+args.ratioABCD+'_'+args.cutTop, 
					True, addUncBand=False )
			############################################
			'''
			testFile = TFile('test.root', 'recreate')
			hDataBtagMinusResonantBkgB.Divide(hDataBtagMinusResonantBkgD)
			unbinnedMCBkgHistos[ 'QCDPtAll_2btag_B' ].Divide( unbinnedMCBkgHistos[ 'QCDPtAll_2btag_D' ] )
			unbinnedMCBkgHistos[ 'QCDPtAll_2btag_B' ].Write()
			MCBkgHistos[ 'QCDPtAll_2btag_C' ].Write()
			htmpDataBtagMinusResonantBkgC.Write()
			hDataBtagMinusResonantBkgB.Write()
			testFile.Write()
			testFile.Close()
			'''

		for isamples in bkgFiles:
			if not 'QCD' in isamples:
				for ibin in range( 14, 19 ): print ibin, isamples, round(MCBkgHistos[ isamples+'_'+args.numBtags+'_A' ].GetBinContent( ibin )*(MCBkgHistos[ isamples+'_'+args.numBtags+'_A' ].GetBinWidth( ibin )),3)

		#### FINAL ESTIMATION WITH BTAGGING using TFunction 
		### stack plot
		stackBtagABCD = THStack( 'stackBtagABCD', 'stackBtagABCD' )
		stackBtagABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
		stackBtagABCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		stackBtagABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		stackBtagABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

		### adding ABCD 
		hBtagABCDOnly = hBtagDataWOResBkgBCD.Clone()
		#hBtagABCDOnly = hDataCRBtag.Clone() 
		hBtagABCDOnly.SetFillColor( kBlue )
		stackBtagABCD.Add( hBtagABCDOnly )
		addAllBkgBtag = hBtagABCDOnly.Clone()
		addAllBkgBtag.Add( dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_A' ] )

		### adding systematics
		hBtagABCDOnlySys = addSysBand( hBtagABCDOnly, 1.10, kBlack ) #, additionalSys=dataHistos['DATA_'+args.numBtags+'_C'])
		hBtagttbarOnlySys = addSysBand( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack ) 		##### TT systematics
		hBtagABCDOnlySys.Add( hBtagttbarOnlySys )
		hBtagwjetsOnlySys = addSysBand( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
		hBtagABCDOnlySys.Add( hBtagwjetsOnlySys )
		#hBtagzjetsOnlySys = addSysBand( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
		hBtagzjetsOnlySys = addSysBand( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
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

		for ibin in range( 12, 17 ): print ibin, 'data', dataHistos[ 'DATA_'+args.numBtags+'_A' ].GetBinContent(ibin)*(dataHistos[ 'DATA_'+args.numBtags+'_A' ].GetBinWidth(ibin)), 'qcd', round(hBtagABCDOnly.GetBinContent(ibin)*(hBtagABCDOnly.GetBinWidth(ibin)), 3), 'signal', round(sigHistos[args.mass].GetBinContent(ibin)*(sigHistos[args.mass].GetBinWidth(ibin)), 3)

		makePlots( nameInRoot, 
				dataHistos[ 'DATA_'+args.numBtags+'_A' ], 'Data - b-tagged sel.', 
				stackBtagABCD, 'Bkg prediction', 
				binWidth, xmin, xmax, 
				[dataHistos[ 'DATA_'+args.numBtags+'_A' ], hBtagABCDOnlySys], "Data/Bkg", 
				'', args.numBtags+'_Log_BCDratio'+args.ratioABCD+'PlusMCbkgs_'+args.cutTop, 
				True, 
				addHisto=addAllBkgBtag, 
				stackHistos=[ 
					[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ], 
					[ hBtagABCDOnly, 'QCD multijets'], 
					[ hBtagABCDOnlySys, 'Bkg. uncertainty' ], 
					[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't#bar{t}+jets' ], 
					int(args.mass),
					[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), "W(q'#bar{q})+jets"], 
					200,
					[ MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z(q#bar{q})+jets'], 
					],
				addUncBand=[ hBtagABCDOnlySys, hBtagRatioSys ], 
				signalHistos=sigHistos,
				plotPull=False)
		###########################################

		##### Full closure test btagging
		stackBtagMCABCD = THStack( 'stackBtagMCABCD', 'stackBtagMCABCD' )
		stackBtagMCABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
		#stackBtagMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		stackBtagMCABCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		stackBtagMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		stackBtagMCABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

		hBtagQCDMCBCD.SetFillColor( kBlue )
		stackBtagMCABCD.Add( hBtagQCDMCBCD )
		hBtagQCDMCBCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )
		hBtagQCDMCBCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		#hBtagQCDMCBCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		hBtagQCDMCBCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
		hBtagQCDMCBCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )

		makePlots( nameInRoot, 
				allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All SM Bkg from MC', 
				stackBtagMCABCD , 'MC Bkg prediction', 
				binWidth, xmin, xmax, 
				[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hBtagQCDMCBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
				'', args.numBtags+'_MCBkg_Log_BCDratio'+args.ratioABCD+'_'+args.cutTop, 
				True, 
				addHisto=hBtagQCDMCBCD, 
				stackHistos=[ [ hBtagQCDMCBCD, 'QCD ABCD from MC'], 
					[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't#bar{t}+jets' ], 
					[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), "W(q'#bar{q})+jets"], 
					#[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
					[ MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z(q#bar{q})+jets'], 
					[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ] ] 
				)
		###########################################

		if not args.final:

			##### Closure: QCD MC replace by Hybrid ABCD for btagging
			stackBtagHybridMCABCD = THStack( 'stackBtagHybridMCABCD', 'stackBtagHybridMCABCD' )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
			#stackBtagHybridMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

			hBtagQCDMCHybridTFunctionBCD.SetFillColor( kBlue )
			stackBtagHybridMCABCD.Add( hBtagQCDMCHybridTFunctionBCD )
			hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )
			hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			#hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			hBtagQCDMCHybridTFunctionBCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )

			makePlots( nameInRoot, 
					allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All SM Bkg from MC', 
					stackBtagHybridMCABCD , 'MC Bkg prediction', 
					binWidth, xmin, xmax, 
					[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hBtagQCDMCHybridTFunctionBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
					'', args.numBtags+'_HybridBkg_Log_BCDratio'+args.ratioABCD+'_'+args.cutTop, 
					True, 
					addHisto=hBtagQCDMCHybridTFunctionBCD, 
					stackHistos=[ [ hBtagQCDMCHybridTFunctionBCD, 'QCD Hybrid ABCD from MC'], 
						[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't#bar{t}+jets' ], 
						[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), "W(q'#bar{q})+jets"], 
						#[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
						[ MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z(q#bar{q})+jets'], 
						[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ] ] 
					)
			###########################################

			
			##### Closure for btag selection: QCD MC with TF FROM INCLUSIVE
			stackBtagHybridMCABCD = THStack( 'stackBtagHybridMCABCD', 'stackBtagHybridMCABCD' )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
			#stackBtagHybridMCABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagHybridMCABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

			hQCDMCBtagHybridTFunctionBCD.SetFillColor( kBlue )
			stackBtagHybridMCABCD.Add( hQCDMCBtagHybridTFunctionBCD )
			hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )
			hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			#hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			hQCDMCBtagHybridTFunctionBCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )

			makePlots( nameInRoot, 
					allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], 'All SM Bkg from MC', 
					stackBtagHybridMCABCD , 'MC Bkg prediction', 
					binWidth, xmin, xmax, 
					[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], hQCDMCBtagHybridTFunctionBCD ], "#frac{All SM Bkg}{All Bkgs with QCD ABCD}", 
					'', args.numBtags+'_InclBkg_Log_BCDratio'+args.ratioABCD+'_'+args.cutTop, 
					True, 
					addHisto=hQCDMCBtagHybridTFunctionBCD, 
					stackHistos=[ [ hQCDMCBtagHybridTFunctionBCD, 'QCD ABCD from MC (TF from incl.)'], 
						[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't#bar{t}+jets' ], 
						[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), "W(q'#bar{q})+jets"], 
						#[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
						[ MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z(q#bar{q})+jets'], 
						[ MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone(), 'Dibosons' ] ]
					#addUncBand=[ allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ], allBkgHistos[ 'allBkg_'+args.numBtags+'_A' ]], 
					)
			###########################################
			
			#### FINAL ESTIMATION WITH BTAGGING using TFunction from inclusive
			### stack plot
			stackBtagABCD = THStack( 'stackBtagABCD', 'stackBtagABCD' )
			stackBtagABCD.Add( MCBkgHistos[ 'Dibosons_'+args.numBtags+'_A' ].Clone() )
			#stackBtagABCD.Add( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagABCD.Add( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagABCD.Add( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone() )
			stackBtagABCD.Add( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone() )

			### adding ABCD 
			hBtagABCDOnly = hDataBtagWOResBkgBCD.Clone()
			hBtagABCDOnly.SetFillColor( kBlue )
			stackBtagABCD.Add( hBtagABCDOnly )
			addAllBkgBtag = hBtagABCDOnly.Clone()
			addAllBkgBtag.Add( dataMinusBkgsHistos[ 'DATA_'+args.numBtags+'_A' ] )

			### adding systematics
			hBtagABCDOnlySys = addSysBand( hBtagABCDOnly, 1.10, kBlack, additionalSys=dataHistos['DATA_'+args.numBtags+'_C'])
			hBtagttbarOnlySys = addSysBand( MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack ) 		##### TT systematics
			hBtagABCDOnlySys.Add( hBtagttbarOnlySys )
			hBtagwjetsOnlySys = addSysBand( MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
			hBtagABCDOnlySys.Add( hBtagwjetsOnlySys )
			#hBtagzjetsOnlySys = addSysBand( MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
			hBtagzjetsOnlySys = addSysBand( MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone(), 1.10, kBlack )
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
					'', args.numBtags+'_Log_Incl_BCDratio'+args.ratioABCD+'PlusMCbkgs_'+args.cutTop, 
					True, 
					addHisto=addAllBkgBtag, 
					stackHistos=[ [ hBtagABCDOnly, 'QCD from ABCD Incl.'], 
						[ MCBkgHistos[ 'TT_'+args.numBtags+'_A' ].Clone(), 't#bar{t}+jets' ], 
						[ MCBkgHistos[ 'WJetsToQQ_'+args.numBtags+'_A' ].Clone(), "W(q'#bar{q})+jets"], 
						#[ MCBkgHistos[ 'ZJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z + Jets'], 
						[ MCBkgHistos[ 'DYJetsToQQ_'+args.numBtags+'_A' ].Clone(), 'Z(q#bar{q})+jets'], 
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
		topRegion=False,
		plotPull=False ):
	"""docstring for makePlots"""

	histo1 = tmphisto1.Clone()
	histo2 = tmphisto2.Clone()
	if 'Data' in labelh1:  
		legend=TLegend(0.43,0.65,0.93,0.89)
		legend.SetNColumns(2)
	else:
		legend=TLegend(0.15,0.80,0.95,0.89)
		legend.SetNColumns(3)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.035)
	if 'Data' in labelh1: legend.AddEntry( histo1, labelh1, 'elp' )
	else: legend.AddEntry( histo1, labelh1, 'l' )
	if isinstance(tmphisto2, THStack):
		for sh in stackHistos: 
			if isinstance(sh, int): legend.AddEntry( signalHistos[ str(sh) ], 'm_{#tilde{t}} = '+str(sh)+' GeV', 'l' )
			else: legend.AddEntry( sh[0], sh[1], 'f' ) 
	else: legend.AddEntry( histo2, labelh2, 'pl' )
	
	histo1.GetYaxis().SetTitle( ('Events / '+str(int(binWidth))+' GeV' if (( 'simple' in args.binning ) or ( 'jet1Tau32' in nameInRoot)) else '< Events / GeV >' ))
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
	#pad1.SetLogx() 	
	if log: 
		pad1.SetLogy() 	
		if 'Region' in labelh1: histo1.SetMaximum( 2.5* max( histo1.GetMaximum(), histo2.GetMaximum() ) )
		else: histo1.SetMaximum( 15000 )
		#else: histo1.SetMaximum( 2100 )
		#histo1.SetMinimum( 270 )
		histo1.SetMinimum( 0.02 )
	else: 
		pad1.SetGrid()
		histo1.SetMaximum( 1.5* max( histo1.GetMaximum(), histo2.GetMaximum() ) )
	if 'Data' in labelh1: 
		histo1.SetMarkerStyle(8)
		histo1.Draw("PE")
		print histo1.GetBinContent(6), histo1.GetBinError(6), histo1.GetXaxis().GetBinLowEdge(6)
	else: histo1.Draw("histe")

	if isinstance( histo2, THStack ): 
		histo2.Draw('hist same')
		histo1.Draw("PE same")
	else: histo2.Draw('hist E0 same')
	
	if addUncBand and ( len(addUncBand)>0 ): 
		addUncBand[0].SetFillStyle(3005)
		addUncBand[0].Draw("same E2")
		addHisto.SetFillColor( 0 )
		#addHisto.Draw("hist same")
		histo1.Draw("PE same")
	if signalHistos: 
		for sample in signalHistos: signalHistos[ sample ].Draw("hist same")

	if isinstance( histo2, THStack ): tmpHisto2 = addHisto
	else: tmpHisto2 = histo2.Clone()
	tmpHisto1 = histo1.Clone()
	'''
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

	numEvents = TLatex( 0.5, 0.75, 'events '+tmpLabel+'/Bkg = '+ str( round( tmpHisto1.Integral(),2 ) )+'/'+str( round( tmpHisto2.Integral(),2 ) ) )
	numEvents.SetNDC()
	numEvents.SetTextFont(42) ### 62 is bold, 42 is normal
	numEvents.SetTextSize(0.035)
	if args.final: numEvents.Draw()
	'''

	if 'Data' in labelh1: tmpLabel = 'Data'
	else: tmpLabel = 'SR'

	#CMS_lumi.extraText = ("Preliminary" if 'Data' in labelh1 else "Simulation Preliminary")
	CMS_lumi.extraText = ("" if 'Data' in labelh1 else "Simulation")
	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(pad1, 4, 0)
	legend.Draw()
	pad1.RedrawAxis()

	pad2.cd()
	pad2.SetGrid()
	#pad2.SetLogx()
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.3)
	tmpPad2= pad2.DrawFrame(xmin,0.1,xmax,1.9)
	tmpPad2.SetXTitle( '#bar{m} [GeV]' ) #'Average '+args.grooming+' jet mass [GeV]' )
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
		if plotPull:
			gStyle.SetHistMinimumZero()
			ratio, dummy = makePulls( ratioList[0], ratioList[1] )
			ratio.GetXaxis().SetTitle( '#bar{m} [GeV]' ) #'Average '+args.grooming+' jet mass [GeV]' )
			ratio.GetXaxis().SetTitleSize( 0.12 )
			ratio.GetXaxis().SetLabelSize( 0.10 )
			ratio.GetXaxis().SetRangeUser( xmin, xmax )
			ratio.GetXaxis().SetNdivisions(505)
			ratio.GetYaxis().SetTitle( '#frac{Data-Bkg}{Uncertainty}' )
			ratio.GetYaxis().CenterTitle( )
			ratio.GetYaxis().SetTitleOffset( 0.5 )
			ratio.GetYaxis().SetTitleSize( 0.11 )
			ratio.GetYaxis().SetLabelSize( 0.12 )
			ratio.GetYaxis().SetNdivisions(505)
			ratio.SetMaximum( 3. )
			ratio.SetMinimum( -3. )
			ratio.SetFillColor(kRed)
			ratio.SetBarWidth(1)
			ratio.SetBarOffset(0)
			ratio.Draw('bar')
		else: 
			#ratio = TGraphAsymmErrors()
			#ratio.Divide( ratioList[0], ratioList[1], 'pois' )
			ratio = ratioList[0].Clone()
			ratio.Reset()
			for ibin in range( 0, ratioList[0].GetNbinsX() ):
				try: 
					ratio.SetBinContent( ibin, ratioList[0].GetBinContent(ibin)/ratioList[1].GetBinContent(ibin) )
					ratio.SetBinError( ibin, ratioList[0].GetBinError(ibin)/ratioList[1].GetBinContent(ibin) )
				except ZeroDivisionError:
					ratio.SetBinContent( ibin, 0 )
					ratio.SetBinError( ibin, 0 )
			ratio.SetMarkerStyle(8)
			ratio.SetLineColor(kBlack)
			ratio.GetXaxis().SetTitle( '#bar{m} [GeV]' ) #'Average '+args.grooming+' jet mass [GeV]' )
			ratio.GetXaxis().SetTitleSize( 0.12 )
			ratio.GetXaxis().SetLabelSize( 0.10 )
			ratio.GetXaxis().SetRangeUser( xmin, xmax )
			ratio.GetXaxis().SetNdivisions(505)
			ratio.GetYaxis().SetTitle( labelRatio ) #'Data/Bkg' )
			ratio.GetYaxis().CenterTitle( )
			ratio.GetYaxis().SetTitleOffset( 0.5 )
			ratio.GetYaxis().SetTitleSize( 0.11 )
			ratio.GetYaxis().SetLabelSize( 0.12 )
			ratio.GetYaxis().SetNdivisions(505)
			ratio.SetMaximum( 2. )
			ratio.SetMinimum( 0. )
			ratio.Draw('P')
	else:
		ratio = ratioList.Clone()
		ratio.SetMarkerStyle(8)
		ratio.SetLineColor(kBlack)
		ratio.Draw('P')


	if not plotPull:
		##### temp
		if signalHistos: 
			tmpSignalBkg = {}
			for sample in signalHistos: 
				tmpSignalBkg[ sample ] = signalHistos[sample].Clone()
				tmpSignalBkg[ sample ].SetFillColorAlpha( signalFiles[ sample ][3], 0.2 ) 
				tmpSignalBkg[ sample ].Add( addHisto )
				if plotPull:
					tmpSignalBkg[ sample ], dummy = makePulls( tmpSignalBkg[ sample ], addHisto )
				else: 
					tmpSignalBkg[ sample ].Divide( addHisto )
				tmpSignalBkg[ sample ].Draw("hist same")
			if not plotPull:
				whiteBox = TGraph(4, array('d', [xmin-10, xmax+10, xmax+10, xmin-10]), array('d', [0, 0, 1, 1]))
				whiteBox.SetFillColor(kWhite)
				whiteBox.Draw("F same")
				pad2.RedrawAxis()
				pad2.RedrawAxis('g')
		###################

		if addUncBand and ( len(addUncBand) > 1 ):
			addUncBand[1].SetFillStyle(3005)
			addUncBand[1].Draw( 'same E2' )
			ratio.Draw('same P')
		elif not topRegion:
			line11.Draw("same")
			line09.Draw("same")
			line.Draw("same")
	else:
		line0.Draw("same")

	if doFit:
		tmpFit = TF1( 'tmpFit', 'pol0', 140, 220 )
		ratio.Fit( 'tmpFit', '', '', 140, 220 )
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

	outputFileName = nameInRoot+'_'+( typePlot if not args.runEra else typePlot+'_Run2016'+args.runEra )+'_'+args.grooming+'_QCD'+args.qcd+'_bkgEstimationPlots'+args.version+'.'+args.extension
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

	CMS_lumi.extraText = ""#Preliminary"
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
	bkgHistos[ nameInRoot+'_'+bkg+'_A' ] = rebin( bkgHistos[ nameInRoot+'_'+bkg+'_A' ], 'reso' )
	bkgHistos[ nameInRoot+'_'+bkg+'_B' ] = rebin( bkgHistos[ nameInRoot+'_'+bkg+'_B' ], 'reso' )
	bkgHistos[ nameInRoot+'_'+bkg+'_C' ] = rebin( bkgHistos[ nameInRoot+'_'+bkg+'_C' ], 'reso' )
	bkgHistos[ nameInRoot+'_'+bkg+'_D' ] = rebin( bkgHistos[ nameInRoot+'_'+bkg+'_D' ], 'reso' )

	histoBC = bkgHistos[ nameInRoot+'_'+bkg+'_A' ].Clone()
	histoBC.Reset()
	histoBC.Multiply( bkgHistos[ nameInRoot+'_'+bkg+'_B' ], bkgHistos[ nameInRoot+'_'+bkg+'_C' ], 1, 1, '')
	histoBCD = bkgHistos[ nameInRoot+'_'+bkg+'_A' ].Clone()
	histoBCD.Reset()
	histoBCD.Divide( histoBC, bkgHistos[ nameInRoot+'_'+bkg+'_D' ], 1, 1, '')

	hRatiohBkg = ratioPlots( bkgHistos[ nameInRoot+'_'+bkg+'_A' ], histoBCD ) 
	makePlots( nameInRoot, 
			bkgHistos[ nameInRoot+'_'+bkg+'_A' ], 
			bkg+' SR', histoBCD, bkg+' MC ABCD Pred.', 5, xmin, xmax, hRatiohBkg, "MC SR/ABCD Pred", '', bkg+'Bkg_Log', True)
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
		print nameInRoot+'_'+sample+('_'+args.numBtags if 'UDD323' in args.decay else '')+'_A' 
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
	quantileHistos = {}
	for q in range( 1, 10 ):
		quantileHistos[ 'X'+str(q) ] = hBkg.QuantilesX( q/10., hBkg.GetName()+'_qX0p'+str(q) )
		quantileHistos[ 'Y'+str(q) ] = hBkg.QuantilesY( q/10., hBkg.GetName()+'_qY0p'+str(q) )
		quantileHistos[ 'Y'+str(q) ] = convertTH1toTGraph( quantileHistos[ 'Y'+str(q) ] )

	if 'JetHT' in sample: CMS_lumi.extraText = ""#"Preliminary"
	else: CMS_lumi.extraText = "Simulation" # Preliminary"
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
	#for k in quantileHistos:  
	#	quantileHistos[k].SetLineWidth(2)
	#	quantileHistos[k].Draw( ( "C" if ('Y' in k) else "hist same" ) )

	textBox.DrawLatex(0.85, 0.85, ( 'm_{#tilde{t}} = '+args.mass+' GeV' if 'RPV' in sample else 'Data - incl. selection' ) )
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
###########################################

def convertTH1toTGraph( histo ):
	"""docstring for convertTH1toTGraph"""
	
	xValues = []
	yValues = []
	for ibin in range( 0, histo.GetNbinsX()+1 ): 
		#print histo.GetBinCenter(ibin), histo.GetBinContent(ibin)
		xValues.append( histo.GetBinCenter(ibin) )
		yValues.append( histo.GetBinContent(ibin) )
	
	newGraph = TGraph( len(yValues),  array('d', yValues), array('d', xValues) )   ### flipped on purpouse 

	return newGraph

###########################################
##### OLD FUNCTION, USE QUANTILEX INSTEAD
def GetQuantileProfiles(Th2f, cut, name): 
	'''2D histogram, value of correlation point, bookkeeping name'''

	q1 = []
	nxbins = Th2f.GetYaxis().GetNbins();
	xlo = Th2f.GetXaxis().GetBinLowEdge(1);
	xhi = Th2f.GetXaxis().GetBinUpEdge(Th2f.GetXaxis().GetNbins() )
	for i in range(nxbins):
		H = Th2f.ProjectionY("ProjY"+str(i),i+1,i+1)
		probSum = array('d', [cut])
		q = array('d', [0.0]*len(probSum))
		H.GetQuantiles(len(probSum), q, probSum)
		print H.GetQuantiles(len(probSum), q, probSum), probSum, q
		q1.append(q[0])

	H1 = TH1F(name, "", nxbins,xlo,xhi)
	for i in range(nxbins): 
		H1.SetBinContent(i+1,q1[i])

	return H1


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
	parser.add_argument('-l', '--lumi', action='store', type=float, default=35864, help='Luminosity, example: 1.' )
	parser.add_argument('-e', '--extension', action='store', default='png', help='Extension of plots.' )
	parser.add_argument('-B', '--numBtags', action='store', default='2btag', help='Number of btags: 1btag or 2btag' )
	parser.add_argument('-f', '--final', action='store_true', default=False, help='Final distributions.' )
	parser.add_argument('-c', '--cutTop', action='store', default='jet1Tau32', help='Cut for ttbar SF.' )
	parser.add_argument('-r', '--runEra', action='store', default='', help='Run data era.' )
	parser.add_argument('-a', '--ratioABCD', action='store', default='BD', help='Run data era.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	
	CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 1 ) )+" fb^{-1}"
	
	if 'Pt' in args.qcd: QCDSF = ( 1 if 'Puppi' in args.grooming else 0.28 )  ### v09 0.28 
	else:  QCDSF = 1

	twoProngSF = 1.21
	antiTwoProngSF = 0.76
	tau32SF = 1.15
	antiTau32SF = 0.96

	if args.miniTree: filePrefix = 'Rootfiles/RUNMiniBoostedAnalysis_'+args.grooming
	else: filePrefix = 'Rootfiles/RUNAnalysis' 

	bkgFiles = OrderedDict() 
	signalFiles = {}
	dataFileName = filePrefix+'_JetHT_Run2016'+args.runEra+'_80X_V2p4_'+args.version+'.root'
	#dataFileName = filePrefix+'_JetHT_Run2016_80X_V2p4_v09p10.root'
	dataFile = TFile.Open(dataFileName)

	signalFiles[ args.mass ] = [ TFile.Open(filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-'+args.mass+'_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'm_{#tilde{t}} = '+args.mass+' GeV', kRed]
	if args.final:
		signalFiles[ '200' ] = [ TFile.Open(filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-200_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'm_{#tilde{t}} = 200 GeV', kMagenta]
		#signalFiles[ '240' ] = [ TFile.Open(filePrefix+'_RPVStopStopToJets_'+args.decay+'_M-240_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'm_{#tilde{t}} = 240 GeV', kMagenta-4]
	bkgFiles[ 'Dibosons' ] = [ TFile.Open(filePrefix+'_Dibosons_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'Dibosons', kMagenta+2 ]
    	bkgFiles[ 'DYJetsToQQ' ] = [ TFile.Open(filePrefix+'_DYJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*1.45, 'Z(q#bar{q})+jets', kOrange]
    	bkgFiles[ 'WJetsToQQ' ] = [ TFile.Open(filePrefix+'_WJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*1.35, "W(q'#bar{q})+jets", 38]
	bkgFiles[ 'QCD'+args.qcd+'All' ] = [ TFile.Open(filePrefix+'_QCD'+args.qcd+'All_Moriond17_80X_V2p4_'+args.version+'.root'), QCDSF*args.lumi, 'QCD multijets sim.', kBlue ]
	bkgFiles[ 'TT' ] = [ TFile.Open(filePrefix+'_TT_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 't#bar{t}+jets', kGreen+2 ]


	massMinX = 65
	massMaxX = 500
	jetMassHTlabY = 0.20
	jetMassHTlabX = 0.85

	if (args.ratioABCD == 'BD'): 
		numeratorABCD = 'B'
		scaleFactorABCD = 'C'
	else:
		numeratorABCD = 'C'
		scaleFactorABCD = 'B'

	if '2D' in args.proc: 
		for signal in signalFiles: 
			plot2DBkgEstimation( 
					signalFiles[ signal ][0], '', 'RPVStopStopToJets_'+args.decay+'_M-'+str(signal), 
					('' if args.miniTree else 'BoostedAnalysisPlots/')+'prunedMassAsymVsdeltaEtaDijet', 
					signalFiles[ signal ][1],
					'm_{asym}', '#Delta #eta', 
					0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		'''
		for bkg in bkgFiles: 
			plot2DBkgEstimation( 
					bkgFiles[ bkg ][0], '', bkg, 
					('' if args.miniTree else 'BoostedAnalysisPlots/')+'prunedMassAsymVsdeltaEtaDijet', 
					bkgFiles[ bkg ][1],
					'm_{asym}', '#Delta #eta', 
					0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		'''
		plot2DBkgEstimation( 
				dataFile, '', 'JetHT_Run2016', 
				('' if args.miniTree else 'BoostedAnalysisPlots/')+'prunedMassAsymVsdeltaEtaDijet', 
				1,
				'm_{asym}', '#Delta #eta', 
				0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		#plot2DBkgEstimation( bkgFiles, dataFile, 'DATAMinusResBkg', 'prunedMassAsymVsdeltaEtaDijet', 'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		#for bkg in bkgFiles: plot2DBkgEstimation( bkgFiles, bkg, 'prunedMassAsymVsdeltaEtaDijet', 'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)
		#for bkg in signalFiles: plot2DBkgEstimation( signalFiles[ bkg ][0], 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass), 'prunedMassAsymVsdeltaEtaDijet', 'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY)

	elif 'simple' in args.proc:
		#for bkg in bkgFiles: 
		#	plotSimpleBkgEstimation( bkgFiles[ bkg ][0], bkg, 'massAve_jet2Tau21VsprunedMassAsym', massMinX, massMaxX, 5, '', '', False )
		#	plotSimpleBkgEstimation( bkgFiles[ bkg ][0], bkg, 'massAve_jet2Tau21VsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )
		#	plotSimpleBkgEstimation( bkgFiles[ bkg ][0], bkg, 'massAve_prunedMassAsymVsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )
		#plotSimpleBkgEstimation( dataFile, 'DATA', 'massAve_jet2Tau21VsprunedMassAsym', massMinX, massMaxX, 5, '', '', False )
		#plotSimpleBkgEstimation( dataFile, 'DATA', 'massAve_jet2Tau21VsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )
		plotSimpleBkgEstimation( dataFile, 'DATA', 'massAve_prunedMassAsymVsdeltaEtaDijet', massMinX, massMaxX, 5, '', '', False )

	else: 
		bkgEstimation( dataFile, bkgFiles, signalFiles, 
				args.grooming+'MassAsymVsdeltaEtaDijet', 
				massMinX, massMaxX, args.binning, 
				#113, 240, 5,
				10,
				'', '', False )
