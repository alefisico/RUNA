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
import argparse
try:
	from RUNA.RUNAnalysis.histoLabels import labels, labelAxis, finalLabels, setSelection
	import RUNA.RUNAnalysis.CMS_lumi as CMS_lumi 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
	from RUNA.RUNAnalysis.commonFunctions import *
except ImportError:
	sys.path.append('../../RUNAnalysis/python')
	from histoLabels import labels, labelAxis, finalLabels
	import CMS_lumi as CMS_lumi 
	import tdrstyle as tdrstyle
	from commonFunctions import *


gROOT.Reset()
gROOT.SetBatch()
#setTDRStyle()
#gROOT.SetStyle('tdrStyle')
#set the tdr style
gROOT.ForceStyle()
tdrstyle.setTDRStyle()

gStyle.SetOptStat(0)


def plotTriggerEfficiency( inFileSample, sample, triggerDenom, name, cut, xmin, xmax, rebin, labX, labY, log, PU ):
	"""docstring for plot"""

	outputFileName = name+'_'+cut+'_'+triggerDenom+"_"+args.trigger+'_'+sample+'_'+args.version+'_TriggerEfficiency.'+args.extension
	print 'Processing.......', outputFileName

	DenomOnly = inFileSample.Get( args.version+'TriggerEfficiency'+args.trigger+'/'+name+'Denom_'+cut ) #cutDijet' ) #+cut )
	DenomOnly.Rebin(rebin)
	Denom = DenomOnly.Clone()
	PassingOnly = inFileSample.Get( args.version+'TriggerEfficiency'+args.trigger+'/'+name+'Passing_'+cut ) #cutHT' ) #+cut )
	PassingOnly.Rebin(rebin)
	Passing = PassingOnly.Clone()
	print Denom, Passing
	Efficiency = TGraphAsymmErrors( Passing, Denom, 'cp'  )
	#Efficiency = TEfficiency( Passing, Denom )

	binWidth = DenomOnly.GetBinWidth(1)

	legend=TLegend(0.50,0.75,0.90,0.90)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.04)

	DenomOnly.SetLineWidth(2)
	DenomOnly.SetLineColor(kRed-4)
	PassingOnly.SetLineWidth(2)
	PassingOnly.SetLineColor(kBlue-4)

	can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
	pad1 = TPad("pad1", "Histo",0,0.46,1.00,1.00,-1)
	pad2 = TPad("pad2", "Efficiency",0,0.00,1.00,0.531,-1);
	pad1.Draw()
	pad2.Draw()

	pad1.cd()
	if log: pad1.SetLogy()

	legend.AddEntry( DenomOnly, triggerDenom+' (baseline trigger)', 'l' )
	legend.AddEntry( PassingOnly, args.trigger, 'l' )
	#DenomOnly.SetMinimum(10)
	DenomOnly.GetXaxis().SetRangeUser( xmin, xmax )
	DenomOnly.Draw('histe')
	DenomOnly.GetYaxis().SetTitleSize(0.06)
	DenomOnly.GetYaxis().SetTitleOffset(0.8)
	DenomOnly.GetYaxis().SetLabelSize(0.06)
	DenomOnly.GetXaxis().SetTitleOffset(0.8)
	DenomOnly.GetXaxis().SetTitleSize(0.06)
	DenomOnly.GetXaxis().SetLabelSize(0.05)
	PassingOnly.Draw('histe same')
	DenomOnly.GetYaxis().SetTitle( 'Events / '+str(binWidth) )

	CMS_lumi.CMS_lumi(pad1, 4, 0)
	legend.Draw()

	pad2.cd()
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.3)
	Efficiency.SetMarkerStyle(8)
	Efficiency.SetLineWidth(2)
	Efficiency.SetLineColor(kBlue-4)
	#Efficiency.SetFillStyle(1001)
	Efficiency.GetYaxis().SetTitle("Efficiency")
	Efficiency.GetYaxis().SetLabelSize(0.06)
	Efficiency.GetXaxis().SetLabelSize(0.06)
	Efficiency.GetYaxis().SetTitleSize(0.06)
	Efficiency.GetYaxis().SetTitleOffset(0.8)
	Efficiency.SetMinimum(-0.1)
	Efficiency.SetMaximum(1.1)
	Efficiency.GetXaxis().SetLimits( xmin, xmax )
	labelAxis( name, Efficiency, 'Pruned')
	Efficiency.Draw()
	#if not (labX and labY): labels( '', '', '')
	#else: labels( '', '', '', labX, labY, 'left' ) #, sel1='AK8PFHT700TrimMass50' )

	can.SaveAs( 'Plots/'+outputFileName.replace('.','Extended.') )
	del can

	#### Fitting
	#errF = TF1('errF', '0.5*(1+TMath::Erf((x-[0])/[1]))', 500, 1500 )
	#errF = TF1('errF', '0.5*(1+TMath::Erf(([0]*x-[1])/[2]))', 400, 1000 )  ## HT
	#errF = TF1('errF', '0.5*(1+TMath::Erf(([0]*x-[1])/[2]))', 0, 100 )  ## Mass
	#Efficiency.SetStatisticOption(TEfficiency.kFWilson)
	#for i in range(5): eff.Fit(errF, '+')
	#for i in range(5): Efficiency.Fit('errF', 'MIR')
	#print '&'*10, '900', errF.Eval(900)
	#print '&'*10, '1000', errF.Eval(1000)
	gStyle.SetOptFit(1)
	can1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	Efficiency.SetMarkerStyle(8)
	Efficiency.SetMarkerColor(kGray)
	Efficiency.SetMinimum(-0.15)
	#Efficiency.SetMinimum(0.7)
	Efficiency.SetMaximum(1.15)
	Efficiency.GetYaxis().SetLabelSize(0.05)
	Efficiency.GetXaxis().SetLabelSize(0.05)
	Efficiency.GetYaxis().SetTitleSize(0.06)
	Efficiency.GetYaxis().SetTitleOffset(0.8)
	Efficiency.GetXaxis().SetTitleOffset(0.8)
	#Efficiency.GetXaxis().SetLimits( 400, 1200 )
	#Efficiency.GetXaxis().SetLimits( 700, 1050 )
	Efficiency.GetXaxis().SetLimits( xmin, xmax )
	Efficiency.Draw('AP')
	'''
	errF.SetLineColor(kRed)
	errF.SetLineWidth(2)
	errF.Draw('sames')
	can1.Update()
	st1 = Efficiency.GetListOfFunctions().FindObject("stats")
	st1.SetX1NDC(.60);
	st1.SetX2NDC(.90);
	st1.SetY1NDC(.20);
	st1.SetY2NDC(.50);
#	#eff.Draw("same")
	can1.Modified()
	'''
	
	'''
	rightmax = 1.2*PassingOnly.GetMaximum()
	rightmin = PassingOnly.GetMinimum()
	scale = gPad.GetUymax()/rightmax
	PassingOnly.SetLineColor(kBlue-5)
	PassingOnly.Scale( scale )
	PassingOnly.Draw( 'hist same' )
	#axis = TGaxis( gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),-3,rightmax,710,"+L")
	axis = TGaxis( gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),rightmin,rightmax,10,"+L")
	axis.SetTitle('Events / '+str(binWidth) )
	axis.SetTitleColor(kBlue-5)
	axis.SetTitleSize(0.06)
	axis.SetLabelSize(0.05)
	axis.SetTitleFont(42)
	axis.SetLabelFont(42)
	axis.SetLineColor(kBlue-5)
	axis.SetLabelColor(kBlue-5)
	axis.SetTitleOffset(0.7)
	axis.Draw()
	'''
	labelAxis( name, Efficiency, 'Pruned')
	CMS_lumi.relPosX = 0.11
	CMS_lumi.cmsTextSize = 0.7
	CMS_lumi.extraOverCmsTextSize = 0.6
	CMS_lumi.CMS_lumi(can1, 4, 0)
	#if not (labX and labY): labels( name, '', PU, camp,  )
	#else: labels( name, '', PU, camp, labX, labY-0.05 ) #, sel1= [ 'AK8PFHT700TrimMass50' ] )

	can1.SaveAs( 'Plots/'+outputFileName )
	del can1

	return Efficiency

def plotDiffEff( listOfEff, name ):
	"""docstring for plotDiffEff"""

	can = TCanvas('c1', 'c1',  10, 10, 750, 500 )

	legend=TLegend(0.60,0.25,0.90,0.40)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.04)

	dummy = 1
	for sample in listOfEff: 
		legend.AddEntry( listOfEff[ sample ], sample, 'l' )

		listOfEff[ sample ].SetMarkerStyle(8)
		listOfEff[ sample ].SetLineWidth(2)
		listOfEff[ sample ].SetLineColor(dummy)
		listOfEff[ sample ].GetYaxis().SetTitle("Efficiency")
		listOfEff[ sample ].GetYaxis().SetLabelSize(0.06)
		listOfEff[ sample ].GetXaxis().SetLabelSize(0.06)
		listOfEff[ sample ].GetYaxis().SetTitleSize(0.06)
		listOfEff[ sample ].GetYaxis().SetTitleOffset(0.8)
		listOfEff[ sample ].SetMinimum(0.8)
		listOfEff[ sample ].SetMaximum(1.05)
		#listOfEff[ sample ].GetXaxis().SetLimits( 850, 950 )
		listOfEff[ sample ].GetXaxis().SetLimits( 0, 200 )
		if dummy == 1:
			labelAxis( name, listOfEff[ sample ], 'Pruned')
			listOfEff[ sample ].Draw()
		else: 
			listOfEff[ sample ].Draw('same')
		dummy+=1

	legend.Draw('same')
	CMS_lumi.lumi_13TeV = ""
	CMS_lumi.relPosX = 0.11
	CMS_lumi.cmsTextSize = 0.7
	CMS_lumi.extraOverCmsTextSize = 0.6
	CMS_lumi.CMS_lumi(can, 4, 0)
	can.SaveAs( 'Plots/'+name+'_DiffEfficiencies.'+args.extension )
	del can


def plot2DTriggerEfficiency( inFileSample, dataset, triggerDenom, name, cut, xlabel, ylabel, Xmin, Xmax, rebinx, Ymin, Ymax, rebiny, labX, labY, PU ):
	"""docstring for plot"""

	outputFileName = name+'_'+cut+'_'+triggerDenom+"_"+args.trigger+'_'+dataset+'_'+args.version+'TriggerEfficiency.'+args.extension
	print 'Processing.......', outputFileName

	print args.version+'TriggerEfficiency'+args.trigger+'/'+name+'Denom_'+cut
	rawDenom = inFileSample.Get( args.version+'TriggerEfficiency'+args.trigger+'/'+name+'Denom_'+cut )
	Denom = Rebin2D( rawDenom, rebinx, rebiny )
	rawPassing = inFileSample.Get( args.version+'TriggerEfficiency'+args.trigger+'/'+name+'Passing_'+cut )
	Passing = Rebin2D( rawPassing, rebinx, rebiny )
	
	'''
	if ( TEfficiency.CheckConsistency( Passing, Denom ) ): Efficiency = TEfficiency( Passing, Denom )
	else: 
		print '--- Passing and Denom are inconsistent.'
		#sys.exit(0)
	'''

	Efficiency = Denom.Clone() 
	Efficiency.Reset()
	Efficiency.Divide( Passing, Denom, 1, 1, 'B' )

	'''
	for i in range( Efficiency.GetNbinsX() ):
		for j in range( Efficiency.GetNbinsY() ):
			if Efficiency.GetXaxis().GetBinLowEdge(i) == 400: print '400: ', round( Efficiency.GetBinContent( i, j ),2), '\pm', round( Efficiency.GetBinError( i , j), 4 )
			if Efficiency.GetXaxis().GetBinLowEdge(i) == 500: print '500: ', round( Efficiency.GetBinContent( i, j ),2), '\pm', round( Efficiency.GetBinError( i , j), 4 )
			if Efficiency.GetXaxis().GetBinLowEdge(i) == 550: print '550: ', round( Efficiency.GetBinContent( i, j ),2), '\pm', round( Efficiency.GetBinError( i , j), 4 )
	'''
	#Efficiency.SetTitle( ';'+xlabel+';'+ ylabel )
	#eff = Efficiency.CreateHistogram()
	
	tdrStyle.SetPadRightMargin(0.12)
	can = TCanvas('c1', 'c1',  10, 10, 1000, 750 )
	gStyle.SetPaintTextFormat("4.2f")
	Efficiency.SetMaximum(1)
	Efficiency.SetMinimum(0)
	Efficiency.Draw('colz')
	Efficiency.Draw('same text')
	#gPad.Update()
	Efficiency.GetYaxis().SetTitleOffset(1.0)
	Efficiency.SetMarkerSize(2)
	Efficiency.GetXaxis().SetRange( int(Xmin/(10.*rebinx)), int(Xmax/(10.*rebinx)) )
	Efficiency.GetXaxis().SetTitle( xlabel )
	Efficiency.GetYaxis().SetTitle( ylabel )
	Efficiency.GetYaxis().SetRange( int(Ymin/(10.*rebiny)), int(Ymax/(10.*rebiny)) )
	#gPad.Update()

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(can, 4, 0)
	#if not (labX and labY): labels( name, '', '', ''  )
	#else: labels( name, '', '', '', labX, labY ) #, sel1= [ 'AK8PFHT700TrimMass50' ] )

	can.SaveAs( 'Plots/'+outputFileName )
	del can




if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--proc', action='store', default='1D', help='Process to draw, example: 1D, 2D, MC.' )
	parser.add_argument('-d', '--dataset', action='store', default='JetHT', help='Dataset: JetHT, SingleMuon, etc.' )
	parser.add_argument('-v', '--version', action='store', default='Boosted', help='Boosted or non version, example: Boosted' )
	parser.add_argument('-g', '--grom', action='store', default='Pruned', help='Grooming Algorithm, example: Pruned, Filtered.' )
	parser.add_argument('-m', '--mass', action='store', default='100', help='Mass of Stop, example: 100' )
	parser.add_argument('-C', '--cut', action='store', default='_cutDEta', help='cut, example: cutDEta' )
	parser.add_argument('-pu', '--PU', action='store', default='Asympt25ns', help='PU, example: PU40bx25.' )
	parser.add_argument('-s', '--single', action='store', default='all', help='single histogram, example: massAve_cutDijet.' )
	parser.add_argument('-c', '--campaign', action='store', default='RunIISpring15DR74', help='Campaign, example: PHYS14.' )
	parser.add_argument('-l', '--lumi', action='store', default='15.5', help='Luminosity, example: 1.' )
	parser.add_argument('-t', '--trigger', action='store', default='AK8PFHT700TrimMass50', help='Trigger used, example PFHT800.' )
	parser.add_argument('-e', '--extension', action='store', default='png', help='Extension of plots.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	#process = args.proc
	PU = args.PU
	camp = args.campaign
	histo = args.single
	mass = args.mass
	cut = args.cut
	grom = args.grom
	single = args.single
	
	triggerlabX = 0.15
	triggerlabY = 1.0
	jetMassHTlabX = 0.87
	jetMassHTlabY = 0.20

	massMinX = 0
	massMaxX = 3*(int(mass))
	HTMinX = 300
	HTMaxX = 1500
	ptMinX = 100
	ptMaxX = 800


	plotList = [ 

		[ '1D', 'HT', HTMinX, HTMaxX, 1, triggerlabX, triggerlabY, True],
		[ '1D', 'prunedMassAve', 0, massMaxX, 2, triggerlabX, triggerlabY, True],
#		[ '1D', 'trimmedMass', massMinX, massMaxX, 2, triggerlabX, triggerlabY, True],
#		#[ '1D', 'ak4HT', HTMinX, HTMaxX, 5, triggerlabX, triggerlabY, True],
		[ '1D', 'jet1Pt', ptMinX, ptMaxX, 2, triggerlabX, triggerlabY, True],
#		[ '1D', 'jet2Pt', ptMinX, ptMaxX, 2, triggerlabX, triggerlabY, True],
		[ '1D', 'jet1PrunedMass', 0, massMaxX, 1, triggerlabX, triggerlabY, True],
		[ '1D', 'jet1SoftDropMass', 0, massMaxX, 1, triggerlabX, triggerlabY, True],
		#[ '1D', 'jet3Pt', 0, 200, 1, triggerlabX, triggerlabY, True],
		#[ '1D', 'jet4Pt', 0, 200, 1, triggerlabX, triggerlabY, True],

		#[ '2D', 'jetMassHTDenom_noTrigger', 'Leading Trimmed Jet Mass [GeV]', 'H_{T} [GeV]', 0, 200, 2, 100, HTMaxX, 5, 0.85, 0.2],
		#[ '2D', 'jetTrimmedMassHTDenom_noTrigger', 'Leading Trimmed Jet Mass [GeV]', 'H_{T} [GeV]', 0, 200, 2, 100, HTMaxX, 5, 0.85, 0.2],
		#[ '2D', 'jetMassHTDenom_triggerOne', 'Leading Trimmed Jet Mass [GeV]', 'H_{T} [GeV]', 0, 200, 2, 100, HTMaxX, 5, 0.85, 0.2],
		#[ '2D', 'jetMassHTDenom_triggerTwo', 'Leading Trimmed Jet Mass [GeV]', 'H_{T} [GeV]', 0, 200, 2, 100, HTMaxX, 5, 0.85, 0.2],
		#[ '2D', 'jetMassHTDenom_triggerOneAndTwo', 'Leading Trimmed Jet Mass [GeV]', 'H_{T} [GeV]', 0, 200, 2, 100, HTMaxX, 5, 0.85, 0.25],

		[ '2D', 'jet1PrunedMassHT', 'Leading Jet Pruned Mass [GeV]', 'HT [GeV]', 20, 300, 2, HTMinX, 1500, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet1SoftDropMassHT', 'Leading Jet SoftDrop Mass [GeV]', 'HT [GeV]', 20, 200, 2, HTMinX, 1200, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet1PtHT', 'Leading Jet Pt [GeV]', 'HT [GeV]', 150, 700, 5, HTMinX, 1200, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet1PrunedMassjet1Pt', 'Leading Jet Pruned Mass [GeV]', 'Leading Jet Pt [GeV]', 20, 300, 2, 150, 700, 5, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet1SoftDropMassjet1Pt', 'Leading Jet SoftDrop Mass [GeV]', 'Leading Jet Pt [GeV]', 20, 300, 2, 150, 700, 5, jetMassHTlabX, jetMassHTlabY],

		[ '2D', 'jet2PrunedMassHT', '2nd leading Jet Pruned Mass [GeV]', 'HT [GeV]', 20, 300, 2, HTMinX, 1500, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet2SoftDropMassHT', '2nd leading Jet SoftDrop Mass [GeV]', 'HT [GeV]', 20, 200, 2, HTMinX, 1200, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet2PtHT', '2nd leading Jet Pt [GeV]', 'HT [GeV]', 150, 700, 5, HTMinX, 1200, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet2PrunedMassjet2Pt', '2nd leading Jet Pruned Mass [GeV]', '2nd leading Jet Pt [GeV]', 20, 300, 2, 150, 700, 5, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet2SoftDropMassjet2Pt', '2nd leading Jet SoftDrop Mass [GeV]', '2nd leading Jet Pt [GeV]', 20, 300, 2, 150, 700, 5, jetMassHTlabX, jetMassHTlabY],

		[ '2D', 'prunedMassAveHT', 'Leading Jet Pruned Mass [GeV]', 'HT [GeV]', 20, 200, 2, HTMinX, 1200, 10, jetMassHTlabX, jetMassHTlabY],
		#[ '2D', 'jetMassHT', 20, 200, 2, HTMinX, 1200, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'jet4PtHT', '4th jet Pt [GeV]', 'HT [GeV]',  20, 200, 2, HTMinX, 1200, 10, jetMassHTlabX, jetMassHTlabY],
		]

	if 'all' in single: Plots = [ x[1:] for x in plotList if x[0] in args.proc ]
	else: Plots = [ y[1:] for y in plotList if ( ( y[0] in args.proc ) and ( y[1] in single ) )  ]

	if 'all' in grom: Groomers = [ '', 'Trimmed', 'Pruned', 'Filtered' ]
	else: Grommers = [ grom ]


	effList = {}
	#for process in [ 'SingleMu', 'MET', 'JetHT' ]:
	bkgFiles = {}
	signalFiles = {}
	CMS_lumi.extraText = "Preliminary"

	if 'SingleMu' in args.dataset:
		inputTrigger = TFile.Open('Rootfiles/RUNTriggerEfficiencies_SingleMuon_Run2016F_V2p1_v00.root')
		BASEDTrigger = 'Mu50'
		CMS_lumi.lumi_13TeV = "2.6 fb^{-1}"
	elif 'MET' in args.dataset:
		inputTrigger = TFile.Open('Rootfiles/RUNTriggerStudies_MET_Run2015D-16Dec2015-v1_v76x_v1p0_v03p1.root')
		BASEDTrigger = 'PFMET170_HBHENoiseCleaned'
	elif 'Vector' in args.dataset:
		inputTrigger = TFile.Open('Rootfiles/RUNTriggerStudies_VectorDiJet1Jet_M50RunIIFall15MiniAODv2_v04.root')
		CMS_lumi.extraText = "Simulation Preliminary"
		CMS_lumi.lumi_13TeV = ""
		if 'MET' in args.dataset: 
			BASEDTrigger = 'PFMET170_HBHENoiseCleaned'
		elif 'SingleMu' in args.dataset:
			BASEDTrigger = 'Mu50'
		else:
			BASEDTrigger = 'PFHT475'
	else:
		#inputTrigger = TFile.Open('RUNFullTriggerEfficiency_JetHT_Run2015C.root' ) 
		#inputTrigger = TFile.Open('RUNTriggerValAndEff_JetHT_Run2016C.root' ) 
		#inputTrigger = TFile.Open('Rootfiles/RUNTriggerEfficiencies_JetHT_Run2016C_V2p1_v00.root')
		inputTrigger = TFile.Open('Rootfiles/RUNTriggerEfficiencies_JetHT_Run2016C_V2p1_v00p3.root')
		#inputTrigger = TFile.Open('Rootfiles/RUNFullTriggerEfficiency_JetHTRun2016C_v00_1.root')
		BASEDTrigger = 'PFHT475'
		CMS_lumi.lumi_13TeV = "1.9 pb^{-1}"

	for i in Plots:
		if '1D' in args.proc:
			effList[ args.dataset ] = plotTriggerEfficiency( inputTrigger, args.dataset, BASEDTrigger, i[0], cut, i[1], i[2], i[3], i[4], i[5], i[6], PU )
		elif '2D' in args.proc:
			plot2DTriggerEfficiency( inputTrigger, args.dataset, BASEDTrigger, i[0], cut, i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], PU )

	#if '1D' in args.proc: plotDiffEff( effList, Plots[0][0] )

