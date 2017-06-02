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
except ImportError:
	sys.path.append('../python')
	from histoLabels import labels, labelAxis, finalLabels
	from scaleFactors import * #scaleFactor as SF
	from cuts import selection 
	import CMS_lumi as CMS_lumi 
	import tdrstyle as tdrstyle




gROOT.Reset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()

xline = array('d', [0,2000])
yline = array('d', [1,1])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)

jetMassHTlabY = 0.20
jetMassHTlabX = 0.85


def plotResolutionCalculator( inFileSample, Groom, name, xmin, xmax, rebinX, labX, labY, log):
	"""docstring for plot"""

	if 'Boosted' in args.boosted: massList = [ 80, 90, 100, 110, 120, 130, 140, 150, 170, 180, 190, 210, 220, 230, 240, 300 ] #, 350 ] 
	else: massList = [ 240, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1100, 1200, 1300 ]
	resolArray = []
	resolErrArray = []
	tmpArray = []

	for xmass in massList:

		gStyle.SetOptFit()
		gStyle.SetStatY(0.91)
		gStyle.SetStatX(0.90)
		gStyle.SetStatW(0.15)
		gStyle.SetStatH(0.15) 
		gStyle.SetTextSize(0.5)
		massAveFile = TFile.Open( inFileSample.replace( str(args.mass), str(xmass) ) )
		outputFileName = name+'_'+args.decay+'RPVSt'+str(xmass)+'_ResolCalc.'+args.extension
		print 'Processing.......', outputFileName

		histos = {}
		histos[ 'massAve' ] = massAveFile.Get( args.boosted+'AnalysisPlots/'+name )
		scale = scaleFactor( 'RPVStopStopToJets_UDD312_M-'+str(xmass) ) 
		#histos[ 'massAve' ].Scale( 1/scale )
		histos[ 'massAve' ].Scale( args.lumi )

		#if rebinX > 1:  histos[ 'massAve' ].Rebin( rebinX )
		histos[ 'massAve' ].Rebin( 5 )

		massWindow = 20
		gausNom = TF1("gaus", "gaus", 0, 500)
		#gausNom.SetParameter(0, histos['massAve'].GetMaximum() )
		gausNom.SetParameter(1, xmass)
		#gausNom.SetParameter(1, 10 )
		for i in range(0,3): histos[ 'massAve' ].Fit(gausNom, 'MIR', '', int(xmass)-massWindow, int(xmass)+massWindow )

		meanGaus = gausNom.GetParameter( 1 )
		resol = gausNom.GetParameter( 2 ) #/ xmass
		resolError = gausNom.GetParError( 2 ) #/xmass 
		#resol = 2.355 * gausNom.GetParameter( 2 ) / meanGaus 
		#resolError = 2.355 * gausNom.GetParError( 2 ) 
		#resol = gausNom.GetParameter( 2 )# / meanGaus 
		resolError = resol * TMath.Sqrt( TMath.Power( (gausNom.GetParError(1)/ meanGaus) , 2) + TMath.Power( (gausNom.GetParError(2)/ gausNom.GetParameter(2)) , 2) )
		resolArray.append( resol )
		resolErrArray.append( resolError )
		tmpArray.append( 0 )
		binWidth = histos['massAve'].GetBinWidth(1)

		#histos[ 'massAve' ].SetLineWidth(2)
		#histos[ 'massAve' ].SetLineColor(kBlack)
		#histos[ 'massAve' ].SetMaximum( 1.2* max( histos[ 'massAve' ].GetMaximum(), histos[ 'Up' ].GetMaximum(), histos[ 'Down' ].GetMaximum() ) ) 
		if xmax: histos[ 'massAve' ].GetXaxis().SetRangeUser( xmin, xmax )

		can = TCanvas('c'+str(xmass), 'c'+str(xmass),  10, 10, 750, 500 )
		if log: can.SetLogy()
		#histos['Sample1'].SetMinimum(10)
		histos['massAve'].GetXaxis().SetRangeUser( xmass-70, xmass+70   )
		histos['massAve'].Draw('histes') 
		gausNom.SetLineColor(kRed)
		gausNom.SetLineWidth(2)
		gausNom.Draw('sames')
		histos['massAve'].GetYaxis().SetTitleOffset(0.9)
		histos['massAve'].GetYaxis().SetTitle( 'Events / '+str(binWidth) )

		labelAxis( name, histos['massAve'], Groom )
		CMS_lumi.relPosX = 0.12
		CMS_lumi.CMS_lumi(can, 4, 0)
		#if not (labX and labY): labels( name, '' )
		#else: labels( name, '', labX, labY )

		can.SaveAs( 'Plots/'+outputFileName )
		del can

	tdrStyle.SetPadRightMargin(0.05)
	tdrStyle.SetPadLeftMargin(0.15)
	can1 = TCanvas('c2', 'c2', 10, 10, 750, 500 )
	#PT = TText(0.1, 0.1, sample )

	print resolArray
	resolGraph = TGraphErrors( len( massList ), array( 'd', massList), array( 'd', resolArray), array( 'd', tmpArray), array( 'd', resolErrArray ) )
	#resolGraph = TGraph( len( massList ), array( 'd', massList), array( 'd', resolArray) )
	resolGraph.SetLineColor( kBlack )
	resolGraph.SetLineWidth( 2 )
	exp = TF1("pol1", "pol1", 80, 250)
	for i in range(0,3): resolGraph.Fit("pol1")

	resolGraph.GetXaxis().SetTitle('Average pruned mass [GeV]')
	resolGraph.GetYaxis().SetTitle('#sigma_{mass}')
	resolGraph.GetYaxis().SetTitleOffset(0.95)
	resolGraph.SetMarkerStyle(20)
	resolGraph.Draw('AP')
	exp.Draw("same")
	#legend.Draw()
	can1.SaveAs('Plots/ResolutionCalculation_'+args.decay+'RPVStop.'+args.extension)
	del can1
	


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--proc', action='store', default='1D', help='Process to draw, example: 1D, 2D, MC.' )
	parser.add_argument('-d', '--decay', action='store', default='UDD312', dest='decay', help='Decay, example: UDD312, UDD323.' )
	parser.add_argument('-b', '--boosted', action='store', default='Boosted', help='Boosted or non boosted, example: Boosted' )
	parser.add_argument('-v', '--version', action='store', default='v00', help='Version: v01, v02.' )
	parser.add_argument('-g', '--grom', action='store', default='pruned', dest='grooming', help='Grooming Algorithm, example: Pruned, Filtered.' )
	parser.add_argument('-m', '--mass', action='store', default='100', help='Mass of Stop, example: 100' )
	parser.add_argument('-C', '--cut', action='store', default='_cutMassAsym', help='cut, example: cutDEta' )
	parser.add_argument('-pu', '--PU', action='store', default='Asympt25ns', help='PU, example: PU40bx25.' )
	parser.add_argument('-s', '--single', action='store', default='all', help='single histogram, example: massAve_cutDijet.' )
	parser.add_argument('-q', '--QCD', action='store', default='Pt', help='Type of QCD binning, example: HT.' )
	parser.add_argument('-c', '--campaign', action='store', default='RunIISpring15MiniAODv2-74X', help='Campaign, example: PHYS14.' )
	parser.add_argument('-l', '--lumi', action='store', type=float, default=2606, help='Luminosity, example: 1.' )
	parser.add_argument('-r', '--range', action='store', default='low', dest='RANGE', help='Trigger used, example PFHT800.' )
	parser.add_argument('-e', '--extension', action='store', default='png', help='Extension of plots.' )
	parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	
	CMS_lumi.lumi_13TeV = '' #str( round((args.lumi/1000.),2) )+" fb^{-1}"
	CMS_lumi.extraText = "Simulation Preliminary"
	
	#plotResolutionCalculator( 'Rootfiles/RUNResolutionCalc_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_RunIIFall15MiniAODv2_v76x_v2p1.root', 'pruned', 'massAve_cutHT', '', '', 1, '', '', False )
	plotResolutionCalculator( 'Rootfiles/RUNResolvedCalc_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_80X_V2p3_v01.root', 'pruned', 'massAve_cutBestPair', '', '', 1, '', '', False )
