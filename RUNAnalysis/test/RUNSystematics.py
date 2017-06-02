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
try:
	from RUNA.RUNAnalysis.histoLabels import labels, labelAxis, finalLabels, setSelection
	from RUNA.RUNAnalysis.scaleFactors import * #scaleFactor as SF
	from RUNA.RUNAnalysis.cuts import selection 
	import RUNA.RUNAnalysis.CMS_lumi as CMS_lumi 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
	from RUNA.RUNAnalysis.commonFunctions import *
except ImportError:
	sys.path.append('../python')
	from histoLabels import labels, labelAxis, finalLabels
	from scaleFactors import * #scaleFactor as SF
	from cuts import selection 
	import CMS_lumi as CMS_lumi 
	import tdrstyle as tdrstyle
	from commonFunctions import *




gROOT.Reset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)

xline = array('d', [0,2000])
yline = array('d', [1,1])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)

jetMassHTlabY = 0.20
jetMassHTlabX = 0.85


def plotSystematics( name, xmin, xmax, labX, labY, log):
	"""docstring for plot"""

	### enlisting masses available
	if 'Boosted' in args.boosted:
		massList = [ 80, 90, 100, 110, 120, 130, 140, 150, 170, 180, 190, 210, 220, 230, 240, 300 ] 
		massWindow = 30 
	else:
		if '312' in args.decay: 
			massList = [ 200, 220, 240 ] + range( 300, 1050, 50 ) + range( 1100, 1300, 100 ) 
			massList.remove( 850 )
		else: massList = [ 200, 220, 240, 280, 300, 350 ] + range( 450, 900, 50 ) + [950, 1000, 1100, 1200, 1300] #+  range( 1100, 1600, 100 )

	nomArray = []
	nomArrayErr = []
	upArray = []
	upArrayErr = []
	downArray = []
	downArrayErr = []
	upOverNomArray = []
	errUpOverNomArray = []
	downOverNomArray = []
	errDownOverNomArray = []

	if args.batchSys: folder = '/cms/gomez/archiveEOS/Archive/v8020/Analysis/'+args.version+'/'
	else: folder = 'Rootfiles/'

	## looping over the masses
	for xmass in massList:

		gStyle.SetOptFit(1)
		rootFile = TFile.Open( folder+'RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+str(xmass)+'_80X_V2p4_'+args.version+'.root' )  ## opening root file
		outputFileName = name+'_'+args.decay+'RPVSt'+str(xmass)+'_'+args.unc+args.boosted+'_'+args.version+'.'+args.ext 
		print 'Processing.......', outputFileName

		### opening different histograms
		histos = {}
		histos[ 'Up' ] = rootFile.Get( 'ResolvedAnalysisPlots/'+name ) 
		histos[ 'Down' ] = rootFile.Get( 'ResolvedAnalysisPlots'+args.unc+'Up/'+name )
		histos[ 'Nominal' ] = rootFile.Get( 'ResolvedAnalysisPlots'+args.unc+'Down/'+name )

		scale = 1 / (scaleFactor( 'RPVStopStopToJets_UDD312_M-'+str(xmass) ) )  ## removing scaling of histogram
		for k in histos: 
			histos[ k ].Scale( scale )
			histos[ k ] = histos[ k ].Rebin( args.reBin )

		#histos[ 'Up' ] = histos[ 'Up' ].Rebin( len( boostedMassAveBins )-1, histos[ 'Up' ].GetName(), boostedMassAveBins )
		#histos[ 'Down' ] = histos[ 'Down' ].Rebin( len( boostedMassAveBins )-1, histos[ 'Down' ].GetName(), boostedMassAveBins )
		#histos[ 'Nominal' ] = histos[ 'Nominal' ].Rebin( len( boostedMassAveBins )-1, histos[ 'Nominal' ].GetName(), boostedMassAveBins )

		#if args.reBin > 1: 
		#for k in histos: histos[ k ] = histos[ k ].Rebin( len( boostedMassAveBins )-1, histos[ k ].GetName(), boostedMassAveBins )

		##### window for acceptance
		if 'Resolved' in args.boosted: massWindow = 9.73 + ( 0.029 * xmass ) 	## from generator level 
		else: massWindow = 30 
		lowEdgeWindow = int(xmass/args.reBin -  2*( massWindow )/args.reBin )
		highEdgeWindow = int(xmass/args.reBin + 2*( massWindow )/args.reBin )

		'''
		####### Fits for Nominal/Up/Down
		gausNom = TF1("gausNom", "gaus", 0, 400 )
		gausNom.SetParameter(1, xmass)
		for i in range(0,5): histos[ 'Nominal' ].Fit(gausNom, 'MIR', '', int(xmass)-massWindow, int(xmass)+massWindow )
		gausNomIntegralError = gausNom.IntegralError(0, 400)
		gausUp = TF1("gausUp", "gaus", 0, 400 )
		gausUp.SetParameter(1, xmass)
		for i in range(0,5): histos[ 'Up' ].Fit(gausUp, 'MIR', '', int(xmass)-massWindow, int(xmass)+massWindow )
		gausUpIntegralError = gausUp.IntegralError(0, 400)
		gausDown = TF1("gausDown", "gaus", 0, 400 )
		gausDown.SetParameter(1, xmass)
		for i in range(0,5): histos[ 'Down' ].Fit(gausDown, 'MIR', '', int(xmass)-massWindow, int(xmass)+massWindow )
		gausDownIntegralError = gausDown.IntegralError(0, 400)
		'''

		binWidth = histos['Nominal'].GetBinWidth(1)
		#totalNumber = scaleFactor( 'RPVStopStopToJets_UDD312_M-'+str(xmass) ) * args.lumi
		totalNumber = search( dictEvents, 'RPVStopStopToJets_UDD312_M-'+str(xmass) )[0]
		errTotalNumber = TMath.Sqrt( totalNumber )
		'''
		nomArray.append( gausNom.Integral( 0, 400 ) / binWidth / totalNumber )
		nomArrayErr.append( gausNomIntegralError / binWidth  / totalNumber )
		upArray.append( gausUp.Integral( 0, 400 ) / binWidth / totalNumber )
		upArrayErr.append( gausUpIntegralError / binWidth / totalNumber )
		downArray.append( gausDown.Integral( 0, 400 ) / binWidth / totalNumber )
		downArrayErr.append( gausDownIntegralError / binWidth / totalNumber )
		upOverNomArray.append( ( gausUp.Integral( 0, 400 ) - gausNom.Integral( 0, 400 ) ) / gausNom.Integral( 0, 400 ) )
		downOverNomArray.append( ( gausDown.Integral( 0, 400 ) - gausNom.Integral( 0, 400 ) ) / gausNom.Integral( 0, 400 ) )
		'''
		errorIntNom = Double(0)
		eventsInWindow = histos['Nominal'].IntegralAndError( lowEdgeWindow, highEdgeWindow, errorIntNom )
		failedEvents = totalNumber - eventsInWindow
		accXeffErr = sqrt( (1/failedEvents) + (1/eventsInWindow) ) * failedEvents * eventsInWindow / pow( ( totalNumber ), 2 )
		nomArray.append( eventsInWindow / totalNumber ) 
		nomArrayErr.append( accXeffErr )

		errorIntUp = Double(0)
		eventsUpInWindow = histos['Up'].IntegralAndError( lowEdgeWindow, highEdgeWindow, errorIntUp )
		failedEvents = totalNumber - eventsUpInWindow
		accXeffErr = sqrt( (1/failedEvents) + (1/eventsUpInWindow) ) * failedEvents * eventsUpInWindow / pow( ( totalNumber ), 2 )
		upArray.append( eventsUpInWindow / totalNumber ) 
		upArrayErr.append( accXeffErr )

		errorIntDown = Double(0)
		eventsDownInWindow = histos['Down'].IntegralAndError( lowEdgeWindow, highEdgeWindow, errorIntDown )
		failedEvents = totalNumber - eventsDownInWindow
		accXeffErr = sqrt( (1/failedEvents) + (1/eventsDownInWindow) ) * failedEvents * eventsDownInWindow / pow( ( totalNumber ), 2 )
		downArray.append( eventsDownInWindow / totalNumber ) 
		downArrayErr.append( accXeffErr )

		upOverNom = ( eventsUpInWindow - eventsInWindow ) / eventsInWindow
		upOverNomArray.append( upOverNom )
		### errors
		errorSubs = TMath.Sqrt( TMath.Power(errorIntUp,2) + TMath.Power( TMath.Sqrt( eventsInWindow ), 2 ) )
		errorUpOverNom = TMath.Abs( upOverNom ) * TMath.Sqrt( TMath.Power( errorIntUp / eventsUpInWindow, 2 ) + TMath.Power( TMath.Sqrt( eventsInWindow )/ eventsInWindow , 2 ) )
		errUpOverNomArray.append( errorUpOverNom )

		downOverNom = ( eventsDownInWindow - eventsInWindow ) / eventsInWindow
		downOverNomArray.append( downOverNom )
		### errors
		errorSubs = TMath.Sqrt( TMath.Power(errorIntDown,2) + TMath.Power( TMath.Sqrt( eventsInWindow ), 2 ) )
		errorDownOverNom = TMath.Abs( downOverNom ) * TMath.Sqrt( TMath.Power( errorIntDown / eventsDownInWindow, 2 ) + TMath.Power( TMath.Sqrt( eventsInWindow )/ eventsInWindow , 2 ) )
		errDownOverNomArray.append( errorDownOverNom )


		legend=TLegend(0.70,0.75,0.90,0.90)
		legend.SetFillStyle(0)
		legend.SetTextSize(0.03)
		legend.AddEntry( histos[ 'Nominal' ], 'Nominal', 'lp' )
		legend.AddEntry( histos[ 'Up' ], args.unc+' Up', 'lp' )
		legend.AddEntry( histos[ 'Down' ], args.unc+' Down', 'lp' )

		histos[ 'Nominal' ].SetLineWidth(2)
		histos[ 'Up' ].SetLineWidth(2)
		histos[ 'Down' ].SetLineWidth(2)
		histos[ 'Nominal' ].SetLineColor(kBlack)
		histos[ 'Up' ].SetLineColor(kBlue)
		histos[ 'Down' ].SetLineColor(kRed)
		histos[ 'Nominal' ].SetMaximum( 1.2* max( histos[ 'Nominal' ].GetMaximum(), histos[ 'Up' ].GetMaximum(), histos[ 'Down' ].GetMaximum() ) ) 
		if xmax: histos[ 'Nominal' ].GetXaxis().SetRangeUser( xmin, xmax )

		if 'PDF' in args.unc:
			tdrStyle.SetPadRightMargin(0.05)
			tdrStyle.SetPadLeftMargin(0.15)
			can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
			pad1 = TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
			pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
			pad1.Draw()
			pad2.Draw()
			pad1.cd()
		else:
			can = TCanvas('c'+str(xmass), 'c'+str(xmass),  10, 10, 750, 500 )
			if log: can.SetLogy()
		#histos['Sample1'].SetMinimum(10)
		histos['Nominal'].GetXaxis().SetRangeUser( xmass-70, xmass+70   )
		histos['Nominal'].Draw('histes') 
		histos['Up'].Draw('histe sames')
		histos['Down'].Draw('histe sames')
		histos['Nominal'].GetYaxis().SetTitleOffset(0.9)
		histos['Nominal'].GetYaxis().SetTitle( 'Events / '+str(binWidth) )
		'''
		gausNom.SetLineColor(kBlack)
		gausNom.Draw('sames')
		gausUp.SetLineColor(kBlue)
		gausUp.Draw('sames')
		gausDown.SetLineColor(kRed)
		gausDown.Draw('sames')

		can1.Update()
		st1 = histos['Nominal'].GetListOfFunctions().FindObject("stats")
		st1.SetX1NDC(.12);
		st1.SetX2NDC(.32);
		st1.SetY1NDC(.76);
		st1.SetY2NDC(.91);
		#st1.SetTextColor(4);
		st2 = histos['Up'].GetListOfFunctions().FindObject("stats")
		st2.SetX1NDC(.12);
		st2.SetX2NDC(.32);
		st2.SetY1NDC(.61);
		st2.SetY2NDC(.76);
		st2.SetTextColor(kBlue);
		st3 = histos['Down'].GetListOfFunctions().FindObject("stats")
		st3.SetX1NDC(.12);
		st3.SetX2NDC(.32);
		st3.SetY1NDC(.46);
		st3.SetY2NDC(.61);
		st3.SetTextColor(kRed);
		can1.Modified()
		'''

		labelAxis( name, histos['Nominal'], ('pruned' if 'Boosted' in args.grooming else '') )
		legend.Draw()
		CMS_lumi.extraText = "Simulation Preliminary"
		CMS_lumi.relPosX = 0.12
		#CMS_lumi.CMS_lumi( (pad1 if 'PDF' in args.unc else 'can' ), 4, 0)
		if not (labX and labY): labels( name, '' )
		else: labels( name, '', labX, labY )

		if 'PDF' in args.unc:
			pad2.cd()
			pad2.SetGrid()
			pad2.SetTopMargin(0)
			pad2.SetBottomMargin(0.3)

			hRatioUp = histos['Nominal'].Clone()
			hRatioUp.Divide( histos['Up'].Clone() )
			hRatioUp.SetMaximum( 1.3 )
			hRatioUp.SetMinimum( 0.7 )
			hRatioUp.SetLineColor(kRed)
			hRatioUp.GetYaxis().SetNdivisions(505)
			hRatioUp.GetXaxis().SetLabelSize(0.12)
			hRatioUp.GetXaxis().SetTitleSize(0.12)
			hRatioUp.GetYaxis().SetTitle( 'Up(Down)/Nominal' )
			hRatioUp.GetYaxis().SetTitleOffset(1.5)
			hRatioUp.GetYaxis().SetLabelSize(0.12)
			hRatioUp.GetYaxis().SetTitleSize(0.12)
			hRatioUp.GetYaxis().SetTitleOffset(0.45)
			hRatioUp.GetYaxis().CenterTitle()
			hRatioUp.Draw()
			hRatioDown = histos['Nominal'].Clone()
			hRatioDown.Divide( histos['Down'].Clone() )
			hRatioDown.SetLineColor(kBlue)
			hRatioDown.Draw("same")

		can.SaveAs( 'Plots/'+outputFileName )
		del can

	tdrStyle.SetPadRightMargin(0.05)
	tdrStyle.SetPadLeftMargin(0.15)
	can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
	pad1 = TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
	pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
	pad1.Draw()
	pad2.Draw()

	pad1.cd()
	if log: pad1.SetLogy()
	#PT = TText(0.1, 0.1, sample )
	multiGraph = TMultiGraph()
	legend=TLegend(0.70,0.70,0.90,0.90)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)

	zeroList = [0]*len(nomArray)
	#nomGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', nomArray), array('d', zeroList ), array( 'd', nomArrayErr) )
	nomGraph = TGraph( len( massList ), array( 'd', massList), array( 'd', nomArray) )		
	nomGraph.SetMarkerStyle( 21 )
	nomGraph.SetLineColor( kBlack )
	nomGraph.SetLineWidth( 2 )
	multiGraph.Add( nomGraph )
	legend.AddEntry( nomGraph, 'Nominal', 'l' )

	#upGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', upArray), array('d', zeroList ), array( 'd', upArrayErr) )
	upGraph = TGraph( len( massList ), array( 'd', massList), array( 'd', upArray) )		
	upGraph.SetMarkerStyle( 22 )
	upGraph.SetLineColor( kBlue )
	upGraph.SetLineWidth( 2 )
	multiGraph.Add( upGraph )
	legend.AddEntry( upGraph, args.unc+'Up', 'l' )

	#downGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', downArray), array('d', zeroList ), array( 'd', downArrayErr) )
	downGraph = TGraph( len( massList ), array( 'd', massList), array( 'd', downArray) )		
	downGraph.SetMarkerStyle( 23 )
	downGraph.SetLineColor( kRed )
	downGraph.SetLineWidth( 2 )
	multiGraph.Add( downGraph )
	legend.AddEntry( downGraph, args.unc+'Down', 'l' )

	multiGraph.Draw("AP")
	multiGraph.GetXaxis().SetTitle('Average pruned mass [GeV]')
	multiGraph.GetYaxis().SetTitle('Acceptance')
	multiGraph.GetYaxis().SetTitleOffset(0.95)
	legend.Draw()
	CMS_lumi.extraText = "Simulation Preliminary"
	#CMS_lumi.lumi_13TeV = ''
	CMS_lumi.relPosX = 0.14
	CMS_lumi.CMS_lumi(pad1, 4, 0)

	pad2.cd()
	pad2.SetGrid()
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.3)
	multiGraphRatio = TMultiGraph()

	uncAcceptanceError = []
	for m in range(len(massList)): 	
		print massList[m], round(upOverNomArray[m],3), round(downOverNomArray[m],3), round( max( abs( upOverNomArray[m] ), abs( downOverNomArray[m] ) ), 3 )
		uncAcceptanceError.append( round( max( abs( upOverNomArray[m] ), abs( downOverNomArray[m] ) ), 3 ) )
	print '-'*10, args.unc+' list:', uncAcceptanceError
	upOverNomGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', upOverNomArray), array( 'd', [0]*len(massList) ), array( 'd', errUpOverNomArray) )		
	upOverNomGraph.SetMarkerStyle( 20 )
	upOverNomGraph.SetMarkerColor( kBlue )
	multiGraphRatio.Add( upOverNomGraph )
	downOverNomGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', downOverNomArray), array( 'd', [0]*len(massList) ), array( 'd', errDownOverNomArray) )		
	downOverNomGraph.SetMarkerStyle( 24 )
	downOverNomGraph.SetMarkerColor( kRed)
	multiGraphRatio.Add( downOverNomGraph )

	multiGraphRatio.Draw("AP")
	if 'PDF' in args.unc: multiGraphRatio.GetYaxis().SetRangeUser(-0.2, .2)
	else: multiGraphRatio.GetYaxis().SetRangeUser(-0.1,0.1)
	multiGraphRatio.GetYaxis().SetNdivisions(505)
	multiGraphRatio.GetXaxis().SetTitle('Average pruned mass [GeV]')
	multiGraphRatio.GetXaxis().SetLabelSize(0.12)
	multiGraphRatio.GetXaxis().SetTitleSize(0.12)
	multiGraphRatio.GetYaxis().SetTitle('Rel. Error')
	#multiGraphRatio.GetYaxis().SetTitle('(Up(Down)-Nom)/Nom')
	multiGraphRatio.GetYaxis().SetTitleOffset(1.5)
	multiGraphRatio.GetYaxis().SetLabelSize(0.12)
	multiGraphRatio.GetYaxis().SetTitleSize(0.12)
	multiGraphRatio.GetYaxis().SetTitleOffset(0.45)
	multiGraphRatio.GetYaxis().CenterTitle()
	can.SaveAs('Plots/'+name+'_'+args.decay+'RPVSt_'+args.unc+args.boosted+'_'+args.version+'.'+args.ext)
	del can
	

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--proc', action='store', default='1D', dest='process', help='Process to draw, example: 1D, 2D, MC.' )
	parser.add_argument('-d', '--decay', action='store', default='UDD312', dest='decay', help='Decay, example: UDD312, UDD323.' )
	parser.add_argument('-b', '--boosted', action='store', default='Boosted', dest='boosted', help='Boosted or Resolved boosted, example: Boosted' )
	parser.add_argument('-v', '--version', action='store', default='v05', dest='version', help='Version of files: v05.' )
	parser.add_argument('-g', '--grom', action='store', default='pruned', dest='grooming', help='Grooming Algorithm, example: Pruned, Filtered.' )
	parser.add_argument('-C', '--cut', action='store', default='_deltaEtaDijet', dest='cut', help='cut, example: cutDEta' )
	parser.add_argument('-e', '--extension', action='store', default='png', dest='ext', help='Extension of plots.' )
	parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
	parser.add_argument('-R', '--reBin', action='store', default=5, type=float, dest='reBin', help='Rebin number.' )
	parser.add_argument('-B', '--batchSys', action='store_true',  dest='batchSys', default=False, help='Process: all or single.' )

	try: args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)


	CMS_lumi.lumi_13TeV = ''#str( round( ( args.lumi / 1000 ), 2 ) )+" fb^{-1}"
	

	plotSystematics( 'massAve_'+args.cut, 0, 400, 0.85, 0.45, False )
			
