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


def plotSignalBkg( signalFiles, bkgFiles, dataFile, nameInRoot, name, xmin, xmax, rebinX, labX, labY, log, posLegend, addRatioFit=False, Norm=False ):
	"""docstring for plot"""

	if 'prunedMassAsymVsdeltaEtaDijet' in args.cut: tmpRegion = '_D'
	else: tmpRegion = ''

	if 'DATA' in args.process: outputFileName = name+tmpRegion+'_'+args.grooming+'_DATA_PlusBkgQCD'+args.qcd+'_'+args.boosted+'AnalysisPlots'+args.version+'.'+args.ext 
	else: outputFileName = name+tmpRegion+'_'+args.grooming+'_'+args.decay+'RPVSt'+str(args.mass)+'_PlusBkgQCD'+args.qcd+'_'+args.boosted+'AnalysisPlots'+args.version+'.'+args.ext 
	if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
	print 'Processing.......', outputFileName
	
	#if args.final:  legend=TLegend(0.50,0.70,0.95,0.89)
	if args.final:  legend=TLegend(0.15,0.70,0.55,0.89)
	elif ( 'DATA' in args.process ) and ('massAve' in nameInRoot): legend=TLegend(0.55,0.70,0.95,0.89)
	else:
		if posLegend: legend=TLegend(0.15, 0.60, 0.4, 0.87 )
		else: legend=TLegend(0.65,0.55,0.90,0.87)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)

	if 'DATA' in args.process:
		dataHistos = {}
		dataHistos[ 'DATA' ] = dataFile.Get( nameInRoot+'_JetHT_Run2016'+tmpRegion if args.miniTree else args.boosted+'AnalysisPlots'+('' if 'pruned' in args.grooming else args.grooming)+'/'+nameInRoot  )
		#if 'mass' in nameInRoot: dataHistos[ 'DATA' ] = dataHistos[ 'DATA' ].Rebin( len( boostedMassAveBins )-1, dataHistos[ 'DATA' ].GetName(), boostedMassAveBins )
		#elif rebinX > 1: dataHistos[ 'DATA' ] = dataHistos[ 'DATA' ].Rebin( rebinX )
		if rebinX > 1: dataHistos[ 'DATA' ] = dataHistos[ 'DATA' ].Rebin( rebinX )
		legend.AddEntry( dataHistos[ 'DATA' ], 'Data', 'ep' )
		if Norm: dataHistos[ 'DATA' ].Scale( 1 /dataHistos['DATA'].Integral() )

	signalHistos = OrderedDict()
	binWidth = 0
	maxList = []
	if len(signalFiles) > 0:
		dummySig=0
		for sigSamples in signalFiles:
			signalHistos[ sigSamples ] = signalFiles[ sigSamples ][0].Get( nameInRoot+'_RPVStopStopToJets_'+args.decay+'_M-'+str(sigSamples)+tmpRegion if args.miniTree else args.boosted+'AnalysisPlots'+('' if 'pruned' in args.grooming else args.grooming )+'/'+nameInRoot )
			#signalHistos[ sigSamples ] = signalFiles[ sigSamples ][0].Get(  nameInRoot+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_A' )
			#if 'mass' in nameInRoot: signalHistos[ sigSamples ] = signalHistos[ sigSamples ].Rebin( len( boostedMassAveBins )-1, signalHistos[ sigSamples ].GetName(), boostedMassAveBins )
			if rebinX > 1: signalHistos[ sigSamples ] = signalHistos[ sigSamples ].Rebin( rebinX )
			if signalFiles[ sigSamples ][1] != 1: signalHistos[ sigSamples ].Scale( signalFiles[ sigSamples ][1] ) 
			if not 'Tau32' in args.cut: legend.AddEntry( signalHistos[ sigSamples ], signalFiles[ sigSamples ][2], 'l' if Norm else 'f' )
			#if not 'DATA' in args.process: legend.AddEntry( signalHistos[ sigSamples ], signalFiles[ sigSamples ][2], 'l' if Norm else 'f' )
			print sigSamples, round( signalHistos[ sigSamples ].Integral(), 2 )
			print sigSamples, 'in mass window', round( signalHistos[ sigSamples ].Integral((args.mass-40)/rebinX, (args.mass+40)/rebinX), 2 )
			if Norm:
				signalHistos[ sigSamples ].SetLineColor( signalFiles[ sigSamples ][3] )
				signalHistos[ sigSamples ].SetLineWidth( 3 )
				signalHistos[ sigSamples ].Scale( 1 / signalHistos[ sigSamples ].Integral() )
				maxList.append( signalHistos[ sigSamples ].GetMaximum() )
			else:
				if 'DATA' in args.process:
					signalHistos[ sigSamples ].SetLineColor( signalFiles[ sigSamples ][3] )
					signalHistos[ sigSamples ].SetFillColor(0)
					signalHistos[ sigSamples ].SetLineWidth(3)
					signalHistos[ sigSamples ].SetLineStyle(2+dummySig)
				else:
					signalHistos[ sigSamples ].SetFillStyle( 1001 )
					signalHistos[ sigSamples ].SetFillColor( signalFiles[ sigSamples ][3] )
					signalHistos[ sigSamples ].SetLineColor( signalFiles[ sigSamples ][3] )
			#if 'mass' in nameInRoot: binWidth = '#sigma_{mass}'
			binWidth = int(signalHistos[ sigSamples ].GetBinWidth( 1 ))
			dummySig+=8

	bkgHistos = OrderedDict()
	if len(bkgFiles) > 0:
		for bkgSamples in bkgFiles:
			bkgHistos[ bkgSamples ] = bkgFiles[ bkgSamples ][0].Get( nameInRoot+'_'+bkgSamples+tmpRegion if args.miniTree else args.boosted+'AnalysisPlots'+('' if 'pruned' in args.grooming else args.grooming )+'/'+nameInRoot )
			#bkgHistos[ bkgSamples ] = bkgFiles[ bkgSamples ][0].Get( nameInRoot+'_'+bkgSamples+'_A' )
			#if 'mass' in nameInRoot: bkgHistos[ bkgSamples ] = bkgHistos[ bkgSamples ].Rebin( len( boostedMassAveBins )-1, bkgHistos[ bkgSamples ].GetName(), boostedMassAveBins )
			if rebinX > 1: bkgHistos[ bkgSamples ] = bkgHistos[ bkgSamples ].Rebin( rebinX )
			if bkgFiles[ bkgSamples ][1] != 1: bkgHistos[ bkgSamples ].Scale( bkgFiles[ bkgSamples ][1] ) 
			legend.AddEntry( bkgHistos[ bkgSamples ], bkgFiles[ bkgSamples ][2], 'l' if Norm else 'f' )
			print bkgSamples, round( bkgHistos[ bkgSamples ].Integral(), 2 )
			print bkgSamples, 'in mass window', round( bkgHistos[ bkgSamples ].Integral((args.mass-40)/rebinX, (args.mass+40)/rebinX), 2 )
			if Norm:
				bkgHistos[ bkgSamples ].SetLineColor( bkgFiles[ bkgSamples ][3] )
				bkgHistos[ bkgSamples ].SetLineWidth( 3 )
				bkgHistos[ bkgSamples ].SetLineStyle( 2 )
				bkgHistos[ bkgSamples ].Scale( 1 / bkgHistos[ bkgSamples ].Integral() )
				maxList.append( bkgHistos[ bkgSamples ].GetMaximum() )
			else:
				bkgHistos[ bkgSamples ].SetFillStyle( 1001 )
				bkgHistos[ bkgSamples ].SetFillColor( bkgFiles[ bkgSamples ][3] )
			

	CMS_lumi.extraText = "Simulation Preliminary"
	hBkg = signalHistos[ args.mass ].Clone()
	hBkg.Reset()
	print 'Total Bkg :', round(hBkg.Integral(), 2 )
	#print 'Total Bkg + Signal :', signalHistos[ args.mass ].Integral() + hBkg.Integral(), signalHistos[ args.mass ].Integral() 
	#print 'Contamination :', signalHistos[ args.mass ].Integral() / ( signalHistos[ args.mass ].Integral() + hBkg.Integral() )

	if not Norm:

		for samples in bkgHistos:
			hBkg.Add( bkgHistos[ samples ].Clone() )
		#for samples in bkgHistos:
		#	bkgHistos[ samples ].Scale(1/hBkg.Integral())

		stackHisto = THStack('stackHisto', 'stack')
		for BkgSamples in bkgHistos: 
			stackHisto.Add( bkgHistos[ BkgSamples ] )

		if 'massAve' in nameInRoot:
			if not ('DATA' in args.process) and not args.final: 
				for SigSamples in signalHistos: stackHisto.Add( signalHistos[ SigSamples ] )

  		tdrStyle.SetPadRightMargin(0.05)
  		tdrStyle.SetPadLeftMargin(0.15)
		can = TCanvas('c1', 'c1',  10, 10, 750, ( 600 if args.final else 750 ) )
		if not args.final:
			pad1 = TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
			pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
			pad1.Draw()
			pad2.Draw()

			pad1.cd()
		if log and not args.final: pad1.SetLogy()
		elif log: can.SetLogy()
		#stackHisto.SetMaximum(0.5)
		stackHisto.Draw('hist')
		if args.final: labelAxis( name, stackHisto, args.grooming )
		if 'massAve' in nameInRoot: 
			stackHisto.SetMinimum( 1 )
			if 'massAve_deltaEtaDijet' in nameInRoot: stackHisto.SetMaximum( 1000 ) 
			elif 'massAve_jet2Tau32WOTau21' in nameInRoot: stackHisto.SetMaximum( 1000 ) 
			elif ('massAve_btag' in nameInRoot) or ('massAve_jet2Tau32' in nameInRoot): stackHisto.SetMaximum( 1000 ) 
			elif 'massAve_jet1Tau32' in nameInRoot: stackHisto.SetMaximum( 1000 ) 
		if xmax: stackHisto.GetXaxis().SetRangeUser( xmin, xmax )
		#stackHisto.SetMaximum( 10 ) 
		#stackHisto.SetMinimum( 0.1 ) 
		#hBkg.SetFillStyle(0)
		hBkg.SetLineColor(kBlack)
		hBkg.SetLineStyle(1)
		hBkg.SetLineWidth(1)
		#hBkg.SetFillStyle(3004)
		#hBkg.SetFillColor( kRed )
		#hBkg.Draw("same")
		stackHisto.GetYaxis().SetTitle( 'Events / '+str(binWidth)+' GeV' )
		#stackHisto.GetYaxis().SetTitle( 'Normalized' ) 
		if 'DATA' in args.process: 
			dataHistos[ 'DATA' ].SetMarkerStyle(8)
			dataHistos[ 'DATA' ].Draw('same')
			CMS_lumi.extraText = "Preliminary"
			legend.SetNColumns(2)
			if not 'Tau32' in args.cut:
				for sample in signalHistos: 
					if 'massAve' in nameInRoot: 
						#lowEdgeWindow = int(int(sample) - ( int( massWidthList[int(sample)])*3 ))
						#highEdgeWindow = int(int(sample) + ( int( massWidthList[int(sample)])*3 ))
						lowEdgeWindow = int(int(sample) - 15 )
						highEdgeWindow = int(int(sample) + 15 )
						signalHistos[ sample ].GetXaxis().SetRangeUser( lowEdgeWindow, highEdgeWindow )
					signalHistos[ sample ].Draw("hist same")
		else:  
			tmpHisto = {}
			for sample in signalHistos: 
				tmpHisto[ sample ] = signalHistos[ sample ].Clone()
				tmpHisto[ sample ].SetFillColor(0)
				tmpHisto[ sample ].SetLineStyle(2)
				tmpHisto[ sample ].SetLineWidth(3)
				tmpHisto[ sample ].Draw("hist same")

		#CMS_lumi.lumi_13TeV = ''
		CMS_lumi.relPosX = 0.14
		CMS_lumi.CMS_lumi( (can if args.final else pad1), 4, 0)
		legend.Draw()
		if not args.final: 
			if not (labX and labY): 
				if 'mini' in args.process: finalLabels( 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) ) 
				elif '1D' in args.process:  labels( args.cut, '' )
				else: labels( name, args.camp )
			else: 
				if 'mini' in args.process: finalLabels( 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass), labX, labY ) 
				elif '1D' in args.process:  labels( args.cut, '', labX, labY )
				else: labels( name, args.camp, labX, labY )

		if not args.final:
			pad2.cd()
			pad2.SetGrid()
			pad2.SetTopMargin(0)
			pad2.SetBottomMargin(0.3)
			
			if 'DATA' in args.process:
				tmpPad2= pad2.DrawFrame(xmin,0.5,xmax,1.5)
				labelAxis( name.replace( args.cut, ''), tmpPad2, ( 'softDrop' if 'Puppi' in args.grooming else args.grooming ) )
				tmpPad2.GetYaxis().SetTitle( "Data/Bkg" )
				tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
				tmpPad2.GetYaxis().CenterTitle()
				tmpPad2.SetLabelSize(0.12, 'x')
				tmpPad2.SetTitleSize(0.12, 'x')
				tmpPad2.SetLabelSize(0.12, 'y')
				tmpPad2.SetTitleSize(0.12, 'y')
				tmpPad2.SetNdivisions(505, 'x')
				tmpPad2.SetNdivisions(505, 'y')
				pad2.Modified()
				hRatio = TGraphAsymmErrors()
				hRatio.Divide( dataHistos[ 'DATA' ], hBkg, 'pois' )
				hRatio.SetMarkerStyle(8)
				hRatio.Draw('P')
				
			else:
				hRatio = signalHistos[ args.mass ].Clone()
				hRatio.Reset()
				allBkgWindow = 0
				allSigWindow = 0
				for ibin in range((args.mass-40)/rebinX, (args.mass+40)/rebinX+1 ): 
					binContSignal = signalHistos[ args.mass ].GetBinContent(ibin)
					allSigWindow += binContSignal
					binContBkg = hBkg.GetBinContent(ibin)
					allBkgWindow += binContBkg
					try: value = binContSignal / TMath.Sqrt( binContBkg )
					#try: value = binContSignal / TMath.Sqrt( binContSignal + binContBkg )
					#try: value = binContSignal / ( binContSignal + binContBkg )
					except ZeroDivisionError: continue
					hRatio.SetBinContent( ibin, value )
				ratioLabel = "S / #sqrt{B}"
				print 's/B ', allSigWindow/TMath.Sqrt(allBkgWindow), allSigWindow, allBkgWindow, allSigWindow/allBkgWindow
			
				labelAxis( name, hRatio, ( 'softDrop' if 'Puppi' in args.grooming else args.grooming) )
				hRatio.GetYaxis().SetTitleOffset(1.2)
				hRatio.GetXaxis().SetLabelSize(0.12)
				hRatio.GetXaxis().SetTitleSize(0.12)
				hRatio.GetYaxis().SetTitle( ratioLabel )
				hRatio.GetYaxis().SetLabelSize(0.12)
				hRatio.GetYaxis().SetTitleSize(0.12)
				hRatio.GetYaxis().SetTitleOffset(0.45)
				hRatio.GetYaxis().CenterTitle()
				#hRatio.SetMaximum(0.7)
				if xmax: hRatio.GetXaxis().SetRangeUser( xmin, xmax )
				hRatio.Draw( ("PES" if 'DATA' in args.process else "hist" ) )

			if addRatioFit:
				tmpFit = TF1( 'tmpFit', 'pol0', 120, 240 )
				hRatio.Fit( 'tmpFit', '', '', 120, 240 )
				tmpFit.SetLineColor( kGreen )
				tmpFit.SetLineWidth( 2 )
				tmpFit.Draw("sames")
				chi2Test = TLatex( 0.7, 0.8, '#splitline{#chi^{2}/ndF = '+ str( round( tmpFit.GetChisquare(), 2 ) )+'/'+str( int( tmpFit.GetNDF() ) )+'}{p0 = '+ str( round( tmpFit.GetParameter( 0 ), 2 ) ) +' #pm '+str(  round( tmpFit.GetParError( 0 ), 2 ) )+'}' )
				chi2Test.SetNDC()
				chi2Test.SetTextFont(42) ### 62 is bold, 42 is normal
				chi2Test.SetTextSize(0.10)
				chi2Test.Draw('same')

			if 'DATA' in args.process: 
				hRatio.GetYaxis().SetNdivisions(505)
				line.Draw('same')

		can.SaveAs( 'Plots/'+outputFileName )
		del can

	else:

  		tdrStyle.SetPadRightMargin(0.05)
		can = TCanvas('c1', 'c1', 750, 500 )
		if log: can.SetLogy()
		signalHistos[args.mass].GetYaxis().SetTitleOffset(1.0)
		signalHistos[args.mass].GetYaxis().SetTitle( ( 'Normalized / '+str(int(binWidth))+' GeV' if name in [ 'massAve', 'HT', 'jet1Pt', 'jet2Pt', 'MET' ] else 'Normalized' ) )
		if xmax: signalHistos[args.mass].GetXaxis().SetRangeUser( xmin, xmax )
		labelAxis( name, signalHistos[args.mass], args.grooming )
		signalHistos[args.mass].Draw('hist')
		if len(signalHistos) > 1: signalHistos[ signalHistos.keys()[1] ].Draw('hist same')
		for bkgSamples in bkgHistos: bkgHistos[ bkgSamples ].Draw('hist same')
		if 'DATA' in args.process: 
			dataHistos[ 'DATA' ].SetMarkerStyle(8)
			dataHistos[ 'DATA' ].Draw('same')
			CMS_lumi.extraText = "Preliminary"
		signalHistos[args.mass].SetMaximum( 1.1 * max( maxList ) )

		if not 'DATA' in args.process: CMS_lumi.lumi_13TeV = ''
		CMS_lumi.relPosX = 0.11
		CMS_lumi.CMS_lumi(can, 4, 0)
		legend.Draw()
		if not (labX and labY): labels( ( '' if 'n-1' in nameInRoot else 'presel'), args.camp )
		else: labels( ( '' if 'n-1' in nameInRoot else 'presel'), args.camp, labX, labY )

		can.SaveAs( 'Plots/'+outputFileName )
		del can

def plot2DSignalBkg( bkgFiles, Groom, nameInRoot, name, titleXAxis, titleXAxis2, Xmin, Xmax, rebinx, Ymin, Ymax, rebiny, legX, legY ):
	"""docstring for plot"""

	outputFileName = name+'_'+Groom+'_'+args.decay+'RPVSt'+args.mass+'_QCD'+args.qcd+'_PlusBkg_'+args.boosted+'AnalysisPlots'+args.version+'.'+args.ext 
	print 'Processing.......', outputFileName

	bkgHistos = OrderedDict()
	bkgHistosProfile = OrderedDict()
	tmpText = ''
	if len(bkgFiles) > 0:
		for bkgSamples in bkgFiles:
			tmpText = bkgSamples
			bkgHistos[ bkgSamples ] =  bkgFiles[ bkgSamples ][0].Get( nameInRoot+'_'+bkgSamples+'_Bkg' )
			bkgHistos[ bkgSamples ] = Rebin2D( bkgHistos[ bkgSamples ], rebinx, rebiny )
			bkgHistosProfile[ bkgSamples ] = bkgHistos[ bkgSamples ].ProfileY( bkgSamples, Ymin, Ymax ) 
			if bkgFiles[ bkgSamples ][1] != 1: bkgHistos[ bkgSamples ].Scale( bkgFiles[ bkgSamples ][1] ) 

	CMS_lumi.extraText = "Simulation Preliminary"
	#if 'QCD' in tmpText: 
	#	hBkg = bkgHistos[ 'QCDHT500to700' ].Clone()
	#	for samples in bkgHistos:
	#		if 'QCDHT500to700' not in samples: hBkg.Add( bkgHistos[ samples ].Clone() )
	#else: hBkg = bkgHistos[ tmpText ].Clone()

	hBkg.GetXaxis().SetTitle( titleXAxis )
	hBkg.GetYaxis().SetTitleOffset( 0.9 )
	hBkg.GetYaxis().SetTitle( titleXAxis2 )

	if (Xmax or Ymax):
		hBkg.GetXaxis().SetRangeUser( Xmin, Xmax )
		hBkg.GetYaxis().SetRangeUser( Ymin, Ymax )

	tdrStyle.SetPadRightMargin(0.12)
	can = TCanvas('c1', 'c1',  750, 500 )
	can.SetLogz()
	hBkg.Draw('colz')

	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(can, 4, 0)
	if not (legX and legY): labels( name, args.camp )
	else: labels( name, args.camp, legX, legY )

	can.SaveAs( 'Plots/'+outputFileName )
	#can.SaveAs( 'Plots/'+outputFileName.replace(''+args.ext, 'gif') )
	del can

def plot2D( inFiles, sample, scale, Groom, nameInRoot, name, titleXAxis, titleXAxis2, Xmin, Xmax, rebinx, Ymin, Ymax, rebiny, legX, legY ):
	"""docstring for plot"""

	outputFileName = nameInRoot+'_'+Groom+'_'+sample+'_'+args.boosted+'AnalysisPlots'+args.version+'.'+args.ext 
	print 'Processing.......', outputFileName
	#for samples in inFiles:
	#h1 = inFiles[ samples ][0].Get( nameInRoot+'_'+sample if 'RPV' in sample else nameInRoot+'_'+samples )
	print nameInRoot+'_'+sample
	h1 = inFiles.Get( nameInRoot+'_'+sample if args.miniTree else args.boosted+'AnalysisPlots'+('' if 'pruned' in args.grooming else args.grooming)+'/'+nameInRoot )
	#h1 = inFiles.Get( nameInRoot+'_'+sample  )
	#h1 = inFile.Get( 'AnalysisPlots'+Groom+'/'+name )
	#h1 = inFile.Get( 'TriggerEfficiency'+Groom+'/'+name )
	tmph1 = h1.Clone()
	
	### Rebinning
	nbinsx = h1.GetXaxis().GetNbins()
	nbinsy = h1.GetYaxis().GetNbins()
	xmin  = h1.GetXaxis().GetXmin()
	xmax  = h1.GetXaxis().GetXmax()
	ymin  = h1.GetYaxis().GetXmin()
	ymax  = h1.GetYaxis().GetXmax()
	nx = nbinsx/rebinx
	ny = nbinsy/rebiny
	h1.SetBins( nx, xmin, xmax, ny, ymin, ymax )

	for biny in range( 1, nbinsy):
		for binx in range(1, nbinsx):
			ibin1 = h1.GetBin(binx,biny)
			h1.SetBinContent( ibin1, 0 )
		
	for biny in range( 1, nbinsy):
		by = tmph1.GetYaxis().GetBinCenter( biny )
		iy = h1.GetYaxis().FindBin(by)
		for binx in range(1, nbinsx):
			bx = tmph1.GetXaxis().GetBinCenter(binx)
			ix  = h1.GetXaxis().FindBin(bx)
			bin = tmph1.GetBin(binx,biny)
			ibin= h1.GetBin(ix,iy)
			cu  = tmph1.GetBinContent(bin)
			h1.AddBinContent(ibin,cu)

	h1.Scale( scale )
	h1.GetXaxis().SetTitle( titleXAxis )
	h1.GetYaxis().SetTitleOffset( 1.0 )
	h1.GetYaxis().SetTitle( titleXAxis2 )

	if (Xmax or Ymax):
		h1.GetXaxis().SetRangeUser( Xmin, Xmax )
		h1.GetYaxis().SetRangeUser( Ymin, Ymax )

	tdrStyle.SetPadRightMargin(0.12)
	can = TCanvas('c1', 'c1',  750, 500 )
	can.SetLogz()
	if 'Boosted' in args.boosted: h1.SetMaximum(5000)
	elif 'Resolved' in args.boosted: h1.SetMaximum(30000)
	#	h1.SetMinimum(0.0001)
	#h1.SetMaximum(100)
	#h1.SetMinimum(0.1)
	h1.Draw('colz')

	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(can, 4, 0)
	if not (legX and legY): labels( name, args.camp )
	else: labels( name, args.camp, legX, legY )

	can.SaveAs( 'Plots/'+outputFileName )
	#can.SaveAs( 'Plots/'+outputFileName.replace(''+args.ext, 'gif') )
	del can


def plotCutFlow( signalFiles, bkgFiles, listOfCuts, name, xmax, log, Norm=False ):
	"""docstring for plot"""

	#outputFileName = name+'_'+args.grooming+'_'+args.decay+'RPVSt'+args.mass+'_Bkg_AnalysisPlots'+args.version+'.'+args.ext 
	#print 'Processing.......', outputFileName

	histos = {}
	histosErr = {}
	signalCF = OrderedDict()
	bkgCF = OrderedDict()
	if len(signalFiles) > 0:
		for iSignal in signalFiles:
			line = 'RPV St '+str(iSignal)+" GeV "
			preFile = TFile( folder+'/RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_80X_V2p4_'+(args.version).replace('p1','')+'.root')
			print preFile
			preCF = preFile.Get('hcutflow')
			preCF.Scale( args.lumi )
			for i in range(1,4): line = line + ' & '+str( round(preCF.GetBinContent(i),2) )+' \pm '+ str( round(preCF.GetBinError(i),2) )
			for icut in listOfCuts:
				histos[ iSignal ] = signalFiles[ iSignal ][0].Get(name+'_'+icut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(iSignal))
				if signalFiles[ iSignal ][1] != 1: histos[ iSignal ].Scale( signalFiles[ iSignal ][1] ) 
				signalIntErr = Double(0)
				signalInt =  histos[ iSignal ].IntegralAndError( 0, xmax, signalIntErr )
				line = line + '& $'+str( round(signalInt,2) )+' \pm '+str( round(signalIntErr,2) )+'$ '
			print line 

	if len(bkgFiles) > 0:
		for iBkg in bkgFiles:
			line = iBkg+' ' 
			for icut in listOfCuts:
				histos[ iBkg ] = bkgFiles[ iBkg ][0].Get(name+'_'+icut+'_'+iBkg)
				if bkgFiles[ iBkg ][1] != 1: histos[ iBkg ].Scale( bkgFiles[ iBkg ][1] ) 
				bkgIntErr = Double(0)
				bkgInt =  histos[ iBkg ].IntegralAndError( 0, xmax, bkgIntErr )
				bkgCF[ str(iBkg)+'_'+icut ] =  '$'+str( round(bkgInt,2) )+' \pm '+str( round(bkgIntErr,2) )+'$'
				line = line + '& $'+str( round(bkgInt,2) )+' \pm '+str( round(bkgIntErr,2) )+'$ '
			print line 
		


def plotSignalCutFlow( runaFile, miniRunaFile, xmax, log, Norm=False ):
	"""docstring for plot"""

	outputFileName = 'signalCutFlow_'+args.grooming+'_'+args.decay+'RPVSt'+args.mass+'_Bkg_AnalysisPlots'+args.version+'.'+args.ext 
	print 'Processing.......', outputFileName

	if 'low' in args.RANGE: massList = [ 90, 100, 110, 120, 130, 140, 150 ] 
	else: massList = [ 170, 180, 190, 210, 220, 230, 240 ] 

	cutFlowValues = OrderedDict()
	histos = OrderedDict()
	for m in massList:
		RUNAFile = TFile( runaFile.replace( '100', str(m) ) )
		cutFlowRUNA = RUNAFile.Get('BoostedAnalysisPlots/cutflow')
		tmpCutflow = []
		cutLabels = []
		for i in range( 1, cutFlowRUNA.GetNbinsX()+1 ): 
			tmpCutflow.append( cutFlowRUNA.GetBinContent( i ) )
			cutLabels.append( cutFlowRUNA.GetXaxis().GetBinLabel( i ) )

		miniRUNAFile = TFile( miniRunaFile.replace( '100', str(m) ) )
		cutFlowMiniRUNA = miniRUNAFile.Get('cutflow_RPVStopStopToJets_'+args.decay+'_M-'+str(m))
		for j in range( 1, cutFlowMiniRUNA.GetNbinsX()+1 ): 
			tmpCutflow.append( cutFlowMiniRUNA.GetBinContent( j ) )
			cutLabels.append( cutFlowMiniRUNA.GetXaxis().GetBinLabel( j ) )
		cutFlowValues[ m ] = tmpCutflow

	legend=TLegend(0.60,0.67,0.90,0.87)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)
	dummy=1
	for p in cutFlowValues:
		histos[ p ] = TH1F( 'cutflow'+str(p), 'cutflow'+str(p), len(cutFlowValues[p]), 0, len(cutFlowValues[p]) )
		for q in range( 1, len(cutFlowValues[p])+1 ):
			histos[p].SetBinContent( q, cutFlowValues[p][q-1]/cutFlowValues[p][0] )
			histos[p].GetXaxis().SetBinLabel( q, cutLabels[q-1] )
		legend.AddEntry( histos[p], args.decay+' RPV #tilde{t} '+str(p)+' GeV' , 'l' )
		histos[ p ].SetLineWidth(2)
		histos[ p ].SetLineColor( dummy )
		dummy+= 1


	can = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	#pad1 = TPad("pad1", "Fit",0,0.20,1.00,1.00,-1)
	#pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
	#pad1.Draw()
	#pad2.Draw()

	can.SetLogy()

	can.SetGridx()
	histos[ massList[0] ].GetYaxis().SetTitle( 'Percentage' )
	histos[ massList[0] ].GetYaxis().SetTitleOffset( 0.8 )
	histos[ massList[0] ].GetXaxis().SetRangeUser( 0, xmax )

	histos[ massList[0] ].Draw()
	for k in histos: 
		if (k != massList[0]): histos[ k ].Draw("same")

	legend.Draw()
	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.lumi_13TeV = ''
	CMS_lumi.relPosX = 0.12
	CMS_lumi.CMS_lumi(can, 4, 0)
	can.SaveAs( 'Plots/'+outputFileName )
	del can

def plotSignalShape( nameInRoot, rebinX, massList, massWidthList, log ):
	"""docstring for plot"""

	outputFileName = 'signalShape_'+nameInRoot+'_'+args.grooming+'_'+args.decay+'RPVSt_'+args.boosted+'AnalysisPlots'+args.version+'.'+args.ext 
	print 'Processing.......', outputFileName

	legend=TLegend(0.65,(0.4 if 'Tau21' in nameInRoot else 0.50),0.90,0.87)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.04)
	histos = {}
	functs = {}
	files = {}
	dummy=1
	
	for imass in massList:
		fileName = folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_RPVStopStopToJets_'+args.decay+'_M-'+str(imass)+'_Moriond17_80X_V2p4_'+args.version+'.root'
		files[ imass ] = TFile( fileName )
	
	maxList = []
	multiGraph = TMultiGraph()
	for m in range( len(massList) ): 
		histos[ massList[m] ] = files[ massList[m] ].Get( nameInRoot+'_RPVStopStopToJets_'+args.decay+'_M-'+str( massList[m] ) )
		#histos[ massList[m] ].Scale( 1/(histos[ massList[m] ].Integral()))
		histos[ massList[m] ].Scale( args.lumi )
		if rebinX > 1: histos[ massList[m] ] = histos[ massList[m] ].Rebin( rebinX )
		if 'Boosted' in args.boosted:
			if 'massAve' in args.single: histos[ massList[m] ].GetXaxis().SetRangeUser( massList[m]-(int(massWidthList[m])*3), massList[m]+(int(massWidthList[m])*3) )
			binWidth = histos[ massList[m] ].GetBinWidth(1)
			maxList.append( histos[ massList[m] ].GetMaximum() )
			histos[ massList[m] ].SetLineWidth(2)
			histos[ massList[m] ].SetLineColor( dummy )
			legend.AddEntry( histos[ massList[m] ], 'M_{#tilde{t}} = '+str( massList[m] )+' GeV' , 'l' )
		else:
			tmpResolution = 9.73 + ( 0.029 * massList[m])
			histos[ massList[m] ].Fit( "gaus", 'ELLSR', '', massList[m]-(5*tmpResolution), massList[m]+(5*tmpResolution) )
			massWindow = histos[ massList[m] ].GetFunction("gaus").GetParameter( 2 )* 3
			functs[ massList[m] ] = TF1( "RPVStop"+str(massList[m]), "gaus", massList[m]-massWindow, massList[m]+massWindow )
			functs[ massList[m] ].SetParameter( 0, histos[ massList[m] ].GetFunction("gaus").GetParameter( 0 ) )
			functs[ massList[m] ].SetParameter( 1, massList[m] )
			functs[ massList[m] ].SetParameter( 2, histos[ massList[m] ].GetFunction("gaus").GetParameter( 2 ) )
			tmpTGraph = TGraph( functs[ massList[m] ] )
			tmpTGraph.SetLineColor( dummy )
			tmpTGraph.SetLineWidth( 2 )
			legend.AddEntry( tmpTGraph, 'M_{#tilde{t}} = '+str( massList[m] )+' GeV' , 'l' )
			multiGraph.Add( tmpTGraph )
			del tmpTGraph
		dummy+= 1
		if dummy == 5: dummy = 6
		if dummy == 9: dummy = 40

	can = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	if log: can.SetLogy()

	if 'Boosted' in args.boosted: 
		if 'mass' in nameInRoot: histos[ massList[-1] ].GetXaxis().SetRangeUser( 50, 350 ) 
		elif 'Pt' in nameInRoot: histos[ massList[-1] ].GetXaxis().SetRangeUser( 200, 2000 ) 
		histos[ massList[-1] ].GetYaxis().SetTitle( 'Normalized' ) #Events / '+str(binWidth)+' GeV' )
		#histos[ massList[-1] ].GetXaxis().SetTitle( 'Average puned jet mass [GeV]' ) 
		histos[ massList[-1] ].GetYaxis().SetTitleOffset( 0.95 )
		#histos[ massList[-1] ].SetMinimum( (0.1 if 'high' in args.RANGE else 0.00001 ) )
		histos[ massList[-1] ].SetMaximum(  1.2*max(maxList) )
		histos[ massList[-1] ].Draw( ('hist' if 'Boosted' in args.boosted else 'l') )
		for k in histos: histos[ k ].Draw( ('hist' if 'Boosted' in args.boosted else 'l')+" same")
		labelAxis( nameInRoot, histos[ massList[-1] ], 'pruned' )
	else: 
		can.SetLogy()
		multiGraph.Draw("al")
		multiGraph.GetYaxis().SetTitle( 'Events' )
		multiGraph.GetXaxis().SetTitle( "Stop mass [GeV]" )
		multiGraph.GetYaxis().SetTitleOffset( 0.8 )


	legend.Draw()
	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.lumi_13TeV = ''
	CMS_lumi.relPosX = 0.12
	CMS_lumi.CMS_lumi(can, 4, 0)
	can.SaveAs( 'Plots/'+outputFileName )
	del can


def plotSignalAcceptance( miniRunaFile, nameInRoot, massList, massWidthList, log ):
	"""docstring for plot acceptance"""

	outputFileName = 'signalAcceptance_'+nameInRoot+'_'+args.grooming+'_'+args.decay+'RPVSt_AnalysisPlots_diffVersions.'+args.ext 
	print 'Processing.......', outputFileName

	legend=TLegend(0.60,(0.5 if 'Tau21' in nameInRoot else 0.67),0.90,0.87)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)
	histos = {}
	files = {}
	dummy=1
	
#	if 'Tau21' in nameInRoot: 
	#massList = [80, 100, 120, 140, 150, 170, 190, 210, 230, 300, 350 ]
	#outputFileName = outputFileName.replace( '_low', '' )
#	else:
	accXeffGraph= {}
	massesList = []
	accXeffList = []
	accXeffErrList = []

	for M in massList:
		fileName = miniRunaFile.replace( str(args.mass), str(M) )
		files[ M ] = TFile( fileName )
	
	for m in range( len(massList) ): 
		NAME = 'RPVStopStopToJets_'+args.decay+'_M-'+str( massList[m] )
		events = search( dictEvents, NAME )[0]
		histos[ massList[m] ] = files[ massList[m] ].Get( nameInRoot+'_'+NAME )
		#histos[ massList[m] ].Scale( 1.10 ) ### due to two prong tagger
		eventsInWindow = histos[ massList[m] ].Integral( ) #massList[m]-(int(2*massWidthList[m])), massList[m]+(int(2*massWidthList[m])) )
		failedEvents = events - eventsInWindow
		#acceptance = eventsInHisto/events
		#efficiency = eventsInWindow/eventsInHisto
		#accXeff = acceptance * efficiency
		accXeff = eventsInWindow / events 
		accXeffErr = sqrt( (1/failedEvents) + (1/eventsInWindow) ) * failedEvents * eventsInWindow / pow( ( events ), 2 )
		print m, eventsInWindow, events, accXeff, accXeffErr
		massesList.append( massList[m] )
		accXeffList.append( accXeff )
		accXeffErrList.append( accXeffErr )
		#legend.AddEntry( accXeffGraph[ args.boosted ], '#tau_{21} < 0.45', 'l' )

	accXeffGraph[ args.boosted ] = TGraphErrors(len(massesList), array( 'd', massesList), array( 'd', accXeffList), array('d', [0]*len(massesList)), array('d', accXeffErrList) ) 

	multiGraph = TMultiGraph()
	can = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	#if log: can.SetLogy()
	can.SetLogy()

	accXeffGraph[ args.boosted ].SetLineColor(kRed)
	accXeffGraph[ args.boosted ].SetLineWidth(2)
	accXeffGraph[ args.boosted ].SetMarkerStyle(8)
	multiGraph.Add( accXeffGraph[ args.boosted ] )

	multiGraph.Draw("ap")
	multiGraph.GetYaxis().SetTitle( 'Acceptance #times efficiency' )
	multiGraph.GetXaxis().SetTitle( "Stop mass [GeV]" )
	multiGraph.GetYaxis().SetTitleOffset( 0.8 )
	#histos[ massList[0] ].GetXaxis().SetRangeUser( (50 if 'low' in args.RANGE else 100 ) , ( 250 if 'low' in args.RANGE else 400 ) )

	#histos[ massList[0] ].SetMinimum( (0.1 if 'high' in args.RANGE else 0.00001 ) )
	#histos[ massList[0] ].SetMaximum( 1.2*max(maxList) )
	#histos[ massList[0] ].Draw('hist')

	#legend.Draw()
	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.lumi_13TeV = ''
	CMS_lumi.relPosX = 0.12
	CMS_lumi.CMS_lumi(can, 4, 0)
	can.SaveAs( 'Plots/'+outputFileName )
	del can

def plotSimple( inFile, sample, Groom, name, xmax, labX, labY, log, Norm=False ):
	"""docstring for plot"""

	outputFileName = name+'_'+sample+'_AnalysisPlots'+args.version+'.'+args.ext 
	print 'Processing.......', outputFileName

	histo = inFile.Get( args.boosted+'AnalysisPlots'+Groom+'/'+name )

#	histos.values()[0].SetMaximum( 2* max( listMax ) ) 
#	histos.values()[0].GetXaxis().SetRangeUser( 0, xmax )
	binWidth = histo.GetBinWidth(1)

	legend=TLegend(0.60,0.75,0.90,0.90)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)

	#histo.SetFillColor(48)
	histo.SetFillStyle(1001)

	tdrStyle.SetPadRightMargin(0.05)
	can = TCanvas('c1', 'c1',  10, 10, 750, 500 )

	if log: 
		can.SetLogy()
		outName = outputFileName.replace('_AnalysisPlots','_Log_AnalysisPlots')
	else:
		outName = outputFileName 

	legend.AddEntry( histo, sample, 'f' )
	histo.GetYaxis().SetTitleOffset(0.90)
	histo.Draw('hist')
	histo.GetYaxis().SetTitle( 'Events / '+str(binWidth) )

	labelAxis( name, histo, '' )
	legend.Draw()
	if not (labX and labY): labels( '', sample )
	else: labels( '', 'MC Truth', labX, labY )
	can.SaveAs( 'Plots/'+outName )
	del can

def plotDiffSample( inFileSample1, inFileSample2, sample1, sample2, name, xmax, labX, labY, log, Diff , Norm=False):
	"""docstring for plot"""

	outputFileName = name+'_'+args.grooming+'_'+args.decay+'_Diff'+Diff+'.'+args.ext 
	#outputFileName = name+'_'+args.grooming+'_'+args.decay+'RPVSt'+args.mass+'_Diff'+Diff+'.'+args.ext 
	print 'Processing.......', outputFileName

	histos = {}
	histos[ 'Sample1' ] = inFileSample1.Get( args.boosted+'AnalysisPlots'+(args.grooming)+'/'+name )
	histos[ 'Sample2' ] = inFileSample2.Get( args.boosted+'AnalysisPlots'+(args.grooming)+'/'+name )

	hSample1 = histos[ 'Sample2' ].Clone()
	hSample1.Scale( 1/hSample1.Integral() ) 
	hSample2 = histos[ 'Sample1' ].Clone()
	hSample2.Scale( 1/hSample2.Integral() ) 
	hRatio = TGraphAsymmErrors()
	hRatio.Divide( hSample1, hSample2, 'pois' )

	binWidth = histos['Sample1'].GetBinWidth(1)

	legend=TLegend(0.60,0.75,0.90,0.90)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)

	if not Norm:
		histos[ 'Sample1' ].SetLineWidth(2)
		histos[ 'Sample1' ].SetLineColor(48)
		histos[ 'Sample2' ].SetLineColor(38)
		histos[ 'Sample2' ].SetLineWidth(2)
		#histos[ 'Sample1' ].SetMaximum( 1.2* max( histos[ 'Sample1' ].GetMaximum(), histos[ 'Sample2' ].GetMaximum() ) ) 
		#histos.values()[0].GetXaxis().SetRangeUser( 0, xmax )

		can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
		pad1 = TPad("pad1", "Fit",0,0.25,1.00,1.00,-1)
		pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.25,-1);
		pad1.Draw()
		pad2.Draw()

		pad1.cd()
		if log: 
			pad1.SetLogy()
			outName = outputFileName.replace('_Diff','_Log_Diff')
		else:
			outName = outputFileName 

		legend.AddEntry( histos[ 'Sample1' ], sample1, 'l' )
		legend.AddEntry( histos[ 'Sample2' ], sample2, 'l' )
		histos['Sample1'].SetMinimum(10)
		histos['Sample1'].Draw('hist')
		histos['Sample1'].GetYaxis().SetTitleOffset(1.2)
		histos['Sample2'].Draw('hist same')
		histos['Sample1'].GetYaxis().SetTitle( 'Events / '+str(binWidth) )

		labelAxis( name, histos['Sample1'], args.grooming )
		legend.Draw()
		#if not (labX and labY): labels( name, '13 TeV - Scaled to '+str(lumi)+' fb^{-1}', '' )
		#else: labels( name, '13 TeV - Scaled to '+str(lumi)+' fb^{-1}', '', labX, labY )

		pad2.cd()
		hRatio.SetLineColor(48)
		#hRatio.SetFillStyle(1001)
		hRatio.GetYaxis().SetTitle("Ratio")
		hRatio.GetYaxis().SetLabelSize(0.12)
		hRatio.GetXaxis().SetLabelSize(0.12)
		hRatio.GetYaxis().SetTitleSize(0.12)
		hRatio.GetYaxis().SetTitleOffset(0.45)
		#hRatio.SetMaximum(1.0)
		hRatio.Sumw2()
		hRatio.Draw("histe")

		can.SaveAs( 'Plots/'+outName )
		del can
	else:
		histos[ 'Sample1' ].SetLineWidth(2)
		histos[ 'Sample1' ].SetLineColor(48)
		histos[ 'Sample2' ].SetLineColor(38)
		histos[ 'Sample2' ].SetLineWidth(2)

		can = TCanvas('c1', 'c1',  10, 10, 750, 500 )
		if log: 
			can.SetLogy()
			outName = outputFileName.replace('_Diff','_Log_Norm_Diff')
			histos[ 'Sample1' ].GetYaxis().SetTitleOffset(1.2)
		else:
			outName = outputFileName.replace('_Diff','_Norm_Diff')

		legend.AddEntry( histos[ 'Sample1' ], sample1 , 'l' )
		legend.AddEntry( histos[ 'Sample2' ], sample2 , 'l' )
		histos['Sample1'].GetYaxis().SetTitle( 'Normalized / '+str(binWidth) )

		histos['Sample1'].DrawNormalized()
		histos['Sample2'].DrawNormalized('same')

		legend.Draw()
		labelAxis( name, histos['Sample1'], args.grooming )
		#if not (labX and labY): labels( name, '13 TeV - Scaled to '+lumi+' fb^{-1}', '' )
		#else: labels( name, '13 TeV - Scaled to '+lumi+' fb^{-1}', '', labX, labY )

		can.SaveAs( 'Plots/'+outName )
		del can


def plotQuality( dataFile, bkgFiles, Groom, nameInRoot, name, xmin, xmax, rebinX, labX, labY, log, moveCMSlogo=False, fitRatio=False ):
	"""docstring for plot"""

	outputFileName = name+'_'+Groom+'_QCD'+args.qcd+'_dataQuality'+args.boosted+'Plots'+args.version+'.'+args.ext
	print 'Processing.......', outputFileName

	histos = {}
	histos[ 'Data' ] = dataFile.Get( nameInRoot+'_JetHT_Run2016' if args.miniTree else nameInRoot )
	if rebinX > 1: histos[ 'Data' ].Rebin( rebinX )
	hBkg = histos[ 'Data'].Clone()
	hBkg.Reset()

	for samples in bkgFiles:
		#print nameInRoot+'_'+samples
		histos[ samples ] = bkgFiles[ samples ][0].Get( nameInRoot+'_'+samples if args.miniTree else nameInRoot )
		if bkgFiles[ samples ][1] != 1: histos[ samples ].Scale( bkgFiles[ samples ][1] ) 
		if rebinX > 1: histos[ samples ].Rebin( rebinX )
		hBkg.Add( histos[ samples ].Clone() )

	hData = histos[ 'Data' ].Clone()
	#hData.Scale( 1/hData.Integral() )
	#hBkg.Scale( 1/hBkg.Integral() )
	hRatio = TGraphAsymmErrors()
	hRatio.Divide( hData, hBkg, 'pois' )
	

	binWidth = histos['Data'].GetBinWidth(1)

	if (labY < 0.5) and ( labX < 0.5 ): legend=TLegend(0.20,0.50,0.50,0.62)
	elif (labX < 0.5): legend=TLegend(0.20,0.75,0.50,0.87)
	else: legend=TLegend(0.70,0.75,0.90,0.87)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.04)
	legend.AddEntry( hData, 'DATA' , 'ep' )
	legend.AddEntry( hBkg, 'All MC Bkgs', 'lp' )

	hBkg.SetLineColor(kRed-4)
	hBkg.SetLineWidth(2)
	hData.SetMarkerStyle(8)

	tdrStyle.SetPadRightMargin(0.05)
	tdrStyle.SetPadLeftMargin(0.15)
	can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
	pad1 = TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
	pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
	pad1.Draw()
	pad2.Draw()

	pad1.cd()
	if log: pad1.SetLogy() 	
	hData.Draw("E")
	hBkg.Draw('hist same')
	hData.SetMaximum( 1.2* max( hData.GetMaximum(), hBkg.GetMaximum() )  )
	#hData.GetYaxis().SetTitleOffset(1.2)
	if xmax: hData.GetXaxis().SetRangeUser( xmin, xmax )
	#hData.GetYaxis().SetTitle( 'Normalized / '+str(int(binWidth))+' GeV' )
	hData.GetYaxis().SetTitle( ( 'Events / '+str(int(binWidth))+' GeV' if name in [ 'massAve', 'HT', 'jet1Pt', 'jet2Pt', 'MET' ] else 'Events' ) )

	#CMS_lumi.relPosX = 0.13
	if moveCMSlogo: 
		CMS_lumi.cmsTextOffset = 0.1
		CMS_lumi.relPosX = 0.15
	else: 
		CMS_lumi.cmsTextOffset = 0.0
		CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(pad1, 4, 0)
	labelAxis( name, hData, Groom )
	legend.Draw()
	if 'deltaEtaDijet' in args.cut: finalLabels( 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass), labX, labY )
	else: 
		if 'n-1' in args.cut: label = 'n-1 selection'
		else: label = 'Preselection'
		if not (labX and labY): setSelection( [ label ], '', ''  )
		else: setSelection( [ label ], labX, labY )

	pad2.cd()
	gStyle.SetOptFit(1)
	pad2.SetGrid()
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.3)
	tmpPad2= pad2.DrawFrame(xmin, ( 0 if fitRatio else 0.5), xmax,1.5)
	labelAxis( name.replace( args.cut, ''), tmpPad2, ( 'softDrop' if 'Puppi' in args.grooming else Groom ) )
	tmpPad2.GetYaxis().SetTitle( "Data/Bkg" )
	tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
	tmpPad2.GetYaxis().CenterTitle()
	tmpPad2.SetLabelSize(0.12, 'x')
	tmpPad2.SetTitleSize(0.12, 'x')
	tmpPad2.SetLabelSize(0.12, 'y')
	tmpPad2.SetTitleSize(0.12, 'y')
	tmpPad2.SetNdivisions(505, 'x')
	tmpPad2.SetNdivisions(505, 'y')
	pad2.Modified()
	hRatio.SetMarkerStyle(8)
	hRatio.Draw('P')
	if fitRatio:
		fitLine = TF1( 'fitLine', 'pol0', 800, 5000)
		hRatio.Fit( 'fitLine', 'MIR')
		fitLine.Draw("same")
		pad2.Update()
		st1 = hRatio.GetListOfFunctions().FindObject("stats")
		st1.SetX1NDC(.65)
		st1.SetX2NDC(.95)
		st1.SetY1NDC(.75)
		st1.SetY2NDC(.95)
		#st1.SetTextColor(kRed)
		pad2.Modified()

	can.SaveAs( 'Plots/'+ outputFileName.replace('Plots', ( 'Fit' if fitRatio else '') ) )
	del can


def tmpplotDiffSample( qcdFile, ttbarFile, signalFile, name, xmin, xmax,reBin,  labX, labY, log, Norm=False):
	"""docstring for plot"""


	histos = OrderedDict()
	#histos[ 'DeltaR' ] = qcdFile.Get( 'ResolvedAnalysisPlots/'+name+'_'+args.cut )
	#histos[ 'Mass' ] = qcdFile.Get( 'ResolvedAnalysisPlotsMassPairing/'+name+'_'+args.cut )
	#histos[ 'KinFit' ] = qcdFile.Get( 'ResolvedAnalysisPlotsChi2Pairing/'+name+'_'+args.cut )
	#histos[ 'DeltaR' ] = signalFile.Get( 'ResolvedAnalysisPlots/'+name+'_'+args.cut )
	#histos[ 'Mass' ] = signalFile.Get( 'ResolvedAnalysisPlotsMassPairing/'+name+'_'+args.cut )
	#histos[ 'KinFit' ] = signalFile.Get( 'ResolvedAnalysisPlotsChi2Pairing/'+name+'_'+args.cut )

	#histos[ '=4jets' ] = qcdFile.Get( name+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	#histos[ '>4jets' ] = ttbarFile.Get( name+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	#histos[ '=4jets' ] = qcdFile.Get( name+args.cut+'_QCDPtAll' )
	#histos[ '>4jets' ] = ttbarFile.Get( name+args.cut+'_QCDPtAll' )

	'''
	if ( args.mass > 0 ):
		outputFileName = name+'_'+args.boosted+'_'+args.decay+'RPV'+str(args.mass)+'_'+args.cut+'_4jets_DiffPairingMethods.'+args.ext 
		#outputFileName = name+'_'+args.boosted+'_'+args.decay+'RPV'+str(args.mass)+args.cut+'_DiffNumJets.'+args.ext 
		histos[ 'DeltaR' ] = qcdFile.Get( name+'_minDeltaR_'+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
		histos[ 'DeltaR + MassAsym cut' ] = qcdFile.Get( name+'_minDeltaR_cutMassAsym_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
		histos[ 'Mass' ] = qcdFile.Get( name+'_minMass_'+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
		histos[ 'KinFit' ] = qcdFile.Get( name+'_minChi_'+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
		'''
	'''
		histos[ 'DeltaR + >4 jets' ] = ttbarFile.Get( name+'_minDeltaR_'+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
		histos[ 'Mass + >4 jets' ] = ttbarFile.Get( name+'_minMass_'+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
		histos[ 'KinFit + >4 jets' ] = ttbarFile.Get( name+'_minChi_'+args.cut+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
		'''
	'''
	else: 
		outputFileName = name+'_'+args.boosted+'_QCD_'+args.cut+'_DiffPairingMethods.'+args.ext 
		#outputFileName = name+'_'+args.boosted+'_QCD_'+args.cut+'_DiffDeltaChi.'+args.ext 
		#outputFileName = name+'_'+args.boosted+'_QCD_'+args.cut+'_DiffNumJets.'+args.ext 
		histos[ 'DeltaR' ] = qcdFile.Get( name+'_minDeltaR_'+args.cut+'_QCDPtAll' )
		histos[ 'DeltaR + MassAsym cut' ] = qcdFile.Get( name+'_minDeltaR_cutMassAsym_QCDPtAll' )
		histos[ 'Mass' ] = qcdFile.Get( name+'_minMass_'+args.cut+'_QCDPtAll' )
		#histos[ 'Mass + delta>75' ] = qcdFile.Get( name+'_minMass_'+args.cut+'75_QCDPtAll' )
		#histos[ 'Mass + delta>100' ] = qcdFile.Get( name+'_minMass_'+args.cut+'100_QCDPtAll' )
		#histos[ 'Mass + delta>150' ] = qcdFile.Get( name+'_minMass_'+args.cut+'150_QCDPtAll' )
		#histos[ 'Mass + delta>200' ] = qcdFile.Get( name+'_minMass_'+args.cut+'_QCDPtAll' )
		histos[ 'KinFit' ] = qcdFile.Get( name+'_minChi_'+args.cut+'_QCDPtAll' )
		#histos[ 'KinFit + delta>100' ] = qcdFile.Get( name+'_minChi_'+args.cut+'100_QCDPtAll' )
		#histos[ 'KinFit + delta>150' ] = qcdFile.Get( name+'_minChi_'+args.cut+'150_QCDPtAll' )
		#histos[ 'KinFit + delta>200' ] = qcdFile.Get( name+'_minChi_'+args.cut+'_QCDPtAll' )
		#histos[ 'DeltaR + >4 jets' ] = ttbarFile.Get( name+'_minDeltaR_'+args.cut+'_QCDPtAll' )
		#histos[ 'Mass + >4 jets' ] = ttbarFile.Get( name+'_minMass_'+args.cut+'_QCDPtAll' )
		#histos[ 'KinFit + >4 jets' ] = ttbarFile.Get( name+'_minChi_'+args.cut+'_QCDPtAll' )
	'''
	outputFileName = name+'_'+args.boosted+'_QCD_DiffMethods2Btag.'+args.ext 
	histos[ 'delta200' ] = ttbarFile.Get( name+'_delta_QCDPtAll' )
	histos[ 'delta200+2CSVL' ] = ttbarFile.Get( name+'_delta_2CSVv2L_QCDPtAll' )
	histos[ 'delta200+2CSVM' ] = ttbarFile.Get( name+'_delta_2CSVv2M_QCDPtAll' )
	histos[ 'delta200+2CSVT' ] = ttbarFile.Get( name+'_delta_2CSVv2T_QCDPtAll' )
	'''
	histos[ 'delta200' ] = qcdFile.Get( name+'_delta_QCDPtAll' )
	histos[ 'delta200+1gql' ] = qcdFile.Get( name+'_delta_1qgl_QCDPtAll' )
	histos[ 'delta200+2qgl' ] = qcdFile.Get( name+'_delta_2qgl_QCDPtAll' )
	histos[ 'delta200+4qgl' ] = qcdFile.Get( name+'_delta_4qgl_QCDPtAll' )
	'''
	'''
	outputFileName = name+'_'+args.boosted+'_QCD_'+args.cut+'_DiffdeltaMethods.'+args.ext 
	histos[ 'delta50' ] = ttbarFile.Get( name+'_delta50_QCDPtAll' )
	histos[ 'delta100' ] = ttbarFile.Get( name+'_delta100_QCDPtAll' )
	histos[ 'delta150' ] = ttbarFile.Get( name+'_delta150_QCDPtAll' )
	histos[ 'delta200' ] = ttbarFile.Get( name+'_delta200_QCDPtAll' )
	histos[ 'delta250' ] = ttbarFile.Get( name+'_delta250_QCDPtAll' )
	histos[ 'delta300' ] = ttbarFile.Get( name+'_delta300_QCDPtAll' )
	histos[ 'delta350' ] = ttbarFile.Get( name+'_delta350_QCDPtAll' )
	histos[ 'delta400' ] = ttbarFile.Get( name+'_delta400_QCDPtAll' )
	histos[ 'delta450' ] = ttbarFile.Get( name+'_delta450_QCDPtAll' )
	histos[ 'delta500' ] = ttbarFile.Get( name+'_delta500_QCDPtAll' )
	#print qcdFile, name+'_maxDelta50_QCDPtAll', histos
	'''
	'''
	histos[ 'delta200' ] = signalFile.Get( name+'_delta_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200+2CSVL' ] = signalFile.Get( name+'_delta_2CSVv2L_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200+2CSVM' ] = signalFile.Get( name+'_delta_2CSVv2M_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200+2CSVT' ] = signalFile.Get( name+'_delta_2CSVv2T_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200' ] = signalFile.Get( name+'_delta_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200+1gql' ] = signalFile.Get( name+'_delta_1qgl_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200+2qgl' ] = signalFile.Get( name+'_delta_2qgl_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200+4qgl' ] = signalFile.Get( name+'_delta_4qgl_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	'''
	'''
	outputFileName = name+'_'+args.boosted+'_'+args.decay+'RPV'+str(args.mass)+args.cut+'_DiffdeltaMethods.'+args.ext 
	histos[ 'delta50' ] = signalFile.Get( name+'_delta50_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta100' ] = signalFile.Get( name+'_delta100_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta150' ] = signalFile.Get( name+'_delta150_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta200' ] = signalFile.Get( name+'_delta200_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta250' ] = signalFile.Get( name+'_delta250_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta300' ] = signalFile.Get( name+'_delta300_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta350' ] = signalFile.Get( name+'_delta350_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta400' ] = signalFile.Get( name+'_delta400_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta450' ] = signalFile.Get( name+'_delta450_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	histos[ 'delta500' ] = signalFile.Get( name+'_delta500_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass) )
	'''
	#inFileSample1 = TFile.Open( inFile ) 
	#histos[ 'Sample1' ] = inFileSample1.Get( 'massAve_deltaEtaDijet_TTJets' )
	#histos[ 'QCD CHS+Pruned' ] = qcdFile.Get( 'BoostedAnalysisPlots/'+name )
	#histos[ 'QCD Puppi+SoftDrop' ] = qcdFile.Get( 'BoostedAnalysisPlotsPuppi/'+name )
	#histos[ 'ttbar CHS+Pruned' ] = ttbarFile.Get( 'BoostedAnalysisPlots/'+name )
	#histos[ 'ttbar Puppi+SoftDrop' ] = ttbarFile.Get( 'BoostedAnalysisPlotsPuppi/'+name )
	#histos[ 'signal CHS+Pruned' ] = signalFile.Get( 'BoostedAnalysisPlots/'+name )
	#histos[ 'signal Puppi+SoftDrop' ] = signalFile.Get( 'BoostedAnalysisPlotsPuppi/'+name )

	#binWidth = histos['Sample1'].GetBinWidth(1)
	print 'Processing.......', outputFileName

	legend=TLegend(0.60,0.60,0.90,0.90)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.03)
	dummy = 1
	maxList = []
	for h in histos:
		histos[h].Rebin( reBin )
		histos[h].Scale( args.lumi ) #1/ histos[h].Integral() )
		histos[h].SetLineWidth(2)
		histos[h].SetLineColor( dummy )
		#histos[h].SetLineColor( (kBlue if 'CHS' in h else kRed) )
		#if 'QCD' in h: histos[h].SetLineStyle( 1 )
		#elif 'ttbar' in h: histos[h].SetLineStyle( 2 )
		#else: histos[h].SetLineStyle( 3 )
		maxList.append( histos[ h ].GetMaximum() )
		legend.AddEntry( histos[h], h, 'l' )
		dummy+=1
		if dummy == 5: dummy = 6
		if dummy == 9: dummy = 40
	
	histos.values()[0].SetMaximum( 1.2* max(maxList) )
	#histos.values()[0].SetMinimum( 0.01 )
	if ( args.mass > 0 ): histos.values()[0].GetXaxis().SetRangeUser( int(args.mass)-200, int(args.mass)+200 ) 
	else: histos.values()[0].GetXaxis().SetRangeUser( 0, xmax )

	can = TCanvas('c1', 'c1',  10, 10, 750, 500 )
	if log: can.SetLogy()
	histos.values()[0].Draw('histe')
	histos.values()[0].GetYaxis().SetTitle( 'Events' )
	histos.values()[0].GetYaxis().SetTitleOffset(0.80)
	for q in range(1,len(histos)): histos.values()[q].Draw('hist same')

	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.lumi_13TeV = ''
	CMS_lumi.relPosX = 0.14
	CMS_lumi.CMS_lumi(can, 4, 0)
	labelAxis( name, histos.values()[0], '' )
	legend.Draw()
	#if not (labX and labY): labels( name, '13 TeV - Scaled to '+lumi+' fb^{-1}', '' )
	#else: 
	#labels( 'presel', '', 0.82, 0.7 )

	can.SaveAs( 'Plots/'+outputFileName )
	del can


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--proc', action='store', default='1D', dest='process', help='Process to draw, example: 1D, 2D, MC.' )
	parser.add_argument('-d', '--decay', action='store', default='UDD312', dest='decay', help='Decay, example: UDD312, UDD323.' )
	parser.add_argument('-b', '--boosted', action='store', default='Boosted', help='Boosted or non boosted, example: Boosted' )
	parser.add_argument('-v', '--version', action='store', default='v00', help='Version: v01, v02.' )
	parser.add_argument('-g', '--grom', action='store', default='pruned', dest='grooming', help='Grooming Algorithm, example: Pruned, Filtered.' )
	parser.add_argument('-m', '--mass', type=int, action='store', default='100', help='Mass of Stop, example: 100' )
	parser.add_argument('-C', '--cut', action='store', default='_deltaEtaDijet', help='cut, example: cutDEta' )
	parser.add_argument('-s', '--single', action='store', default='all', help='single histogram, example: massAve_cutDijet.' )
	parser.add_argument('-q', '--qcd', action='store', default='Pt', dest='qcd', help='Type of QCD binning, example: HT.' )
	parser.add_argument('-c', '--camp', action='store', default='RunIISpring15MiniAODv2-74X', help='Campaign, example: PHYS14.' )
	parser.add_argument('-l', '--lumi', action='store', type=float, default=35870, help='Luminosity, example: 1.' )
	parser.add_argument('-r', '--range', action='store', default='low', dest='RANGE', help='Trigger used, example PFHT800.' )
	parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
	parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
	parser.add_argument('-f', '--final', action='store_true', default=False, dest='final',  help='If plot is final' )
	parser.add_argument('-F', '--addFit', action='store_true', default=False, dest='addFit',  help='Plot fit in ratio plot.' )
	parser.add_argument('-t', '--miniTree', action='store_true', default=False, help='miniTree: if plots coming from miniTree or RUNAnalysis.' )
	parser.add_argument('-B', '--batchSys', action='store_true',  dest='batchSys', default=False, help='Process: all or single.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	
	bkgFiles = OrderedDict() 
	signalFiles = OrderedDict()
	CMS_lumi.extraText = "Preliminary"
	CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 1 ) )+" fb^{-1}"
	
	if 'Pt' in args.qcd: 
		bkgLabel='(w QCD pythia8)'
		##QCDSF = ( 0.68 if 'Resolved' in args.boosted else ( 1 if 'Puppi' in args.grooming else 0.46) )  ### 0.54 with all data for resolved ### v08
		QCDSF = ( 0.255 if 'Resolved' in args.boosted else ( 1 if 'Puppi' in args.grooming else 0.366 ) )  ### v09
	else: 
		bkgLabel='(w QCD madgraphMLM+pythia8)'
		QCDSF = 0.66

	if args.batchSys: folder = '/cms/gomez/archiveEOS/Archive/v8020/Analysis/'+args.version+'/'
	else: folder = 'Rootfiles/'

	if 'Resolved' in args.boosted: args.grooming = ''

	if args.miniTree:
		dataFile = TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_JetHT_Run2016_80X_V2p4_'+args.version+'.root')
		signalFiles[ args.mass ] = [ TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_Moriond17_80X_V2p4_'+args.version+'.root'), 
				args.lumi, 
				'M_{#tilde{t}} = '+str(args.mass)+' GeV', 
				kRed]
		if ( 'Norm' in args.process ) or ( 'DATA' in args.process ) or ( 'CF' in args.process ): 
			otherMass = ( '180' if args.boosted == 'Boosted' else '700' )
			signalFiles[ otherMass ] = [ TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_RPVStopStopToJets_'+args.decay+'_M-'+otherMass+'_Moriond17_80X_V2p4_'+args.version+'.root'), 
					args.lumi, 
					'M_{#tilde{t}} = '+otherMass+' GeV', 
					kRed+2]

		if 'Boosted' in args.boosted: 
			bkgFiles[ 'TT' ] = [ TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_TT_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 't #bar{t} + Jets', kGreen+2 ]
			bkgFiles[ 'WJetsToQQ' ] = [ TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_WJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'W + Jets', 38 ]
			bkgFiles[ 'Dibosons' ] = [ TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_Dibosons_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'Dibosons', kMagenta+2 ]
			bkgFiles[ 'ZJetsToQQ' ] = [ TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_ZJetsToQQ_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi, 'Z + Jets', kOrange ]
		bkgFiles[ 'QCD'+args.qcd+'All' ] = [ TFile.Open(folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_QCD'+args.qcd+'All_Moriond17_80X_V2p4_'+args.version+'.root'), args.lumi*QCDSF, 'QCD'+args.qcd+'', kBlue-4 ]
	else:
		dataFile = TFile.Open(folder+'/RUNAnalysis_JetHT_Run2016_80X_V2p4_'+args.version+'.root')
		signalFiles[ args.mass ] = [ TFile.Open(folder+'/RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_80X_V2p4_'+args.version+'.root'), args.lumi, 'M_{#tilde{t}} = '+str(args.mass)+' GeV', kRed]
		#signalFiles[ args.mass ] = [ TFile.Open('~/mySpace/archiveEOS/Archive/v7414/RUNAnalysis_RPVSt350tojj_13TeV_pythia8RunIISpring15MiniAODv2-74X_Asympt25ns_v09_v01.root'), args.lumi, 'M_{#tilde{t}} = '+str(args.mass)+' GeV', kRed]
		if ( 'Norm' in args.process ) or ( 'DATA' in args.process ): 
			otherMass = ( '180' if args.boosted == 'Boosted' else '700' )
			signalFiles[ otherMass ] = [ TFile.Open(folder+'/RUNAnalysis_RPVStopStopToJets_'+args.decay+'_M-'+str(otherMass)+'_80X_V2p4_'+args.version+'.root'), args.lumi, 'M_{#tilde{t}} = '+str(otherMass)+' GeV', kRed]
			#signalFiles[ '800' ] = [ TFile.Open('~/mySpace/archiveEOS/Archive/v7414/RUNAnalysis_RPVStopStopToJets_UDD312_M-800-madgraph_RunIISpring15MiniAODv2-74X_Asympt25ns_v09_v03.root'), args.lumi, 'M_{#tilde{t}} = 800 GeV', kRed]
		if 'Boosted' in args.boosted: 
			bkgFiles[ 'TTJets' ] = [ TFile.Open(folder+'/RUNAnalysis_TT_80X_V2p4_'+args.version+'.root'), args.lumi, 't #bar{t} + Jets', kGreen+2 ]
			bkgFiles[ 'WJetsToQQ' ] = [ TFile.Open(folder+'/RUNAnalysis_WJetsToQQ_80X_V2p4_'+args.version+'.root'), args.lumi , 'W + Jets', 38 ]
			bkgFiles[ 'Dibosons' ] = [ TFile.Open(folder+'/RUNAnalysis_Dibosons_80X_V2p4_'+args.version+'.root'), args.lumi , 'Dibosons', kMagenta+2 ]
			bkgFiles[ 'ZJetsToQQ' ] = [ TFile.Open(folder+'/RUNAnalysis_ZJetsToQQ_80X_V2p4_'+args.version+'.root'), args.lumi, 'Z + Jets', kOrange ]
			#bkgFiles[ 'WWTo4Q' ] = [ TFile.Open(folder+'/RUNAnalysis_WWTo4Q_80X_V2p4_'+args.version+'.root'), args.lumi , 'WW (had)', kMagenta+2 ]
			#bkgFiles[ 'ZZTo4Q' ] = [ TFile.Open(folder+'/RUNAnalysis_ZZTo4Q_80X_V2p4_'+args.version+'.root'), args.lumi, 'ZZ (had)', kOrange+2 ]
			#bkgFiles[ 'WZ' ] = [ TFile.Open(folder+'/RUNAnalysis_WZ_80X_V2p4_'+args.version+'.root'), args.lumi, 'WZ', kCyan ]
		bkgFiles[ 'QCD'+args.qcd+'All' ] = [ TFile.Open(folder+'/RUNAnalysis_QCD'+args.qcd+'All_80X_V2p4_'+args.version+'.root'), args.lumi*QCDSF, 'QCD'+args.qcd+'', kBlue-4 ]
		bkgFiles[ 'QCDPtAll' ] = [ TFile.Open(folder+'/RUNAnalysis_QCDHTAll_80X_V2p4_'+args.version+'.root'), args.lumi*QCDSF, 'QCD'+args.qcd+'', kBlue-4 ]
		#bkgFiles[ 'QCD'+args.qcd+'All' ] = [ TFile.Open('~/mySpace/archiveEOS/Archive/v7414/RUNAnalysis_QCDPtAll_RunIISpring15MiniAODv2-74X_Asympt25ns_v09_v01.root'), args.lumi*QCDSF, 'QCD'+args.qcd+'', kBlue-4 ]


	dijetlabX = 0.85
	dijetlabY = 0.55
	subjet112vs212labX = 0.7
	subjet112vs212labY = 0.88
	polAnglabX = 0.2
	polAnglabY = 0.88
	taulabX = 0.90
	taulabY = 0.85
	cosPhilabX = 0.15
	cosPhilabY = 0.45

	massMinX = 0
	massMaxX = 400
	polAngXmin = 0.7
	polAngXmax = 1.0
	HTMinX = 300
	HTMaxX = 1300
	ptMinX = 100
	ptMaxX = 800


	plotList = [ 
		[ '2D', 'Boosted', 'leadMassHT', 'Leading Jet Mass [GeV]', 'HT [GeV]', 0, massMaxX, 1, 100, HTMaxX, 1, jetMassHTlabX, jetMassHTlabY],

		[ '2D', 'Boosted', 'jet1Tau21VsRhoDDT', 'Leading jet #tau_{21}', 'Leading jet #rho\'', 0, 1, 1, -6, 10, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'jet2Tau21VsRhoDDT', '2nd Leading jet #tau_{21}', '2nd Leading jet #rho\'', 0, 1, 1, -6, 10, 1, jetMassHTlabX, jetMassHTlabY],

		[ '2D', 'Boosted', 'massAsymVsdeltaEtaDijet', 'Mass Asymmetry', '| #eta_{j1} - #eta_{j2} |', 0, 1, 1, 0, 5, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'massAsymVsjet1Tau21', 'Mass Asymmetry', 'Leading jet #tau_{21} |', 0, 1, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'massAsymVsjet2Tau21', 'Mass Asymmetry', '2nd Leading jet #tau_{21} |', 0, 1, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'massAsymVsjet1Tau31', 'Mass Asymmetry', 'Leading jet #tau_{31} |', 0, 1, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'massAsymVsjet2Tau31', 'Mass Asymmetry', '2nd Leading jet #tau_{31} |', 0, 1, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'deltaEtaDijetVsjet1Tau21', '| #eta_{j1} - #eta_{j2} |', 'Leading jet #tau_{21} |', 0, 5, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'deltaEtaDijetVsjet2Tau21', '| #eta_{j1} - #eta_{j2} |', '2nd Leading jet #tau_{21} |', 0, 5, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'deltaEtaDijetVsjet1Tau31', '| #eta_{j1} - #eta_{j2} |', 'Leading jet #tau_{31} |', 0, 5, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'deltaEtaDijetVsjet2Tau31', '| #eta_{j1} - #eta_{j2} |', '2nd Leading jet #tau_{31} |', 0, 5, 1, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'jet1Tau21VsRhoDDT', 'Leading jet #tau_{21}', 'Leading jet #rho\'', 0, 1, 1, -6, 10, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Boosted', 'jet2Tau21VsRhoDDT', '2nd Leading jet #tau_{21}', '2nd Leading jet #rho\'', 0, 1, 1, -6, 10, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'deltavsMassAve', 'Average dijet mass [GeV]', 'Delta', 100, 700, 20, -200, 1000, 20, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAveVsMinChi2', 'Average dijet mass [GeV]', 'Min( #chi^{2} )', 0, 1000, 10, 0, 5, 2, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAveVsMinDeltaR', 'Average dijet mass [GeV]', 'Min( #Delta R )', 0, 1000, 10, 0, 5, 2, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAveVsMinMassAsym', 'Average dijet mass [GeV]', 'Min( Mass Asym )', 0, 1000, 10, 0, 1, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAveVscosThetaStar', 'Average dijet mass [GeV]', '#cos #theta^{*}', 0, 1000, 10, 0, 2, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAveVsXi', 'Average dijet mass [GeV]', '#chi', 0, 1000, 10, 0, 1.5, 1, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAveminChi2vsminDeltaR', 'Average dijet mass (#chi^{2}) [GeV]', 'Average dijet mass (#Delta R)', 0, 1000, 20, 0, 1000, 20, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'minChi2vsminDeltaR', 'Min( #chi^{2} )', 'Min #Delta R', 0, 5, 2, 0, 5, 2, jetMassHTlabX, jetMassHTlabY],
		#[ '2D', 'Resolved', 'delta', 'Average dijet mass [GeV]', 'Delta', 100, 700, 10, 0, 1000, 10, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'jet12vsjet34Mass', 'dijet 1 mass [GeV]', 'dijet 2 mass [GeV]', 0, 1000, 20, 0, 1000, 20, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAvevsjet34Mass', 'average dijet mass [GeV]', 'dijet 2 mass [GeV]', 0, 1000, 20, 0, 1000, 20, jetMassHTlabX, jetMassHTlabY],
		[ '2D', 'Resolved', 'massAvevsjet12Mass', 'average dijet mass [GeV]', 'dijet 1 mass [GeV]', 0, 1000, 20, 0, 1000, 20, jetMassHTlabX, jetMassHTlabY],

		[ '1D', 'Boosted', 'jet1Pt', 100, 1500, 1, '', '', False],
		[ '1D', 'Boosted', 'jet1Eta', -3, 3, 1, '', '', False],
		[ '1D', 'Boosted', 'jet1Mass', 0, massMaxX, 10, '', '', True, False],
		[ '1D', 'Boosted', 'jet1PrunedMass', '', '', 10, '', '', True, False],
		[ '1D', 'Boosted', 'HT', 700, 2000, 5, '', '', False],
		[ '1D', 'Boosted', 'jet2Pt', 100, 1500, 1, '', '', False],
		[ '1D', 'Boosted', 'jet2Eta', -3, 3, 1, '', '', False],
		[ '1D', 'Boosted', 'jet2Mass', 0, massMaxX, 1, '', '', False],
		[ '1D', 'Boosted', 'massAve', 60, 350, 5, 0.92, 0.85, True, False],
		#[ '1DDATA', 'Boosted', 'massAve', 60, 350, (5 if 'deltaEtaDijet' in args.cut else 10 ), 0.92, 0.85, True, False],
		#[ '1DData', 'Boosted', 'massAve', 60, 350, (1 if args.miniTree else 5), 0.92, 0.85, True, False],
		#[ '1DDATA', 'Resolved', 'massAve', 0, 1000, 20, 0.92, 0.85, False, False],

		[ '1D', 'Resolved', 'HT', 700, 5000, 2, '', '', True, False],
		[ '1D', 'Resolved', 'jet1Pt', 100, 1500, 2, '', '', True, False],
		[ '1D', 'Resolved', 'jet2Pt', 0, 1500, 2, '', '', True, False],
		[ '1D', 'Resolved', 'jet3Pt', 0, 500, 2, '', '', True, False],
		[ '1D', 'Resolved', 'jet4Pt', 50, 500, 1, '', '', True, False],
		#[ '1D', 'Resolved', 'massAve', 0, 1000, 2, '', '', True, False],
		[ '1D', 'Resolved', 'massAve', 0, 1500, 10, '', '', True, False],
		[ '1D', 'Boosted', 'deltaEtaDijet', '', '', 5,  0.90, 0.70, False, False],
		[ '1D', 'Boosted', 'prunedMassAsym', '', '', 1, (0.90 if 'n-1' in args.cut else 0.40), (0.70 if 'n-1' in args.cut else 0.4), False, ( False if 'n-1' in args.cut else True )],	
		[ '1D', 'Boosted', 'jet1Tau21', '', '', 1, 0.9, 0.70, False, True],
		[ '1D', 'Boosted', 'jet2Tau21', '', '', 1, 0.9, 0.70, False, True],

		[ 'qual', 'Boosted', 'deltaEtaDijet', 0, 5, 1,  0.90, 0.70, True, False],	### n- 1
		[ 'qual', 'Boosted', 'prunedMassAsym', 0, 1, 1, (0.90 if 'n-1' in args.cut else 0.40), (0.70 if 'n-1' in args.cut else 0.4), False, ( False if 'n-1' in args.cut else True )],	
		[ 'qual', 'Boosted', 'jet1Tau21', 0, 1, 1, (0.90 if 'n-1' in args.cut else 0.40), (0.70 if 'n-1' in args.cut else 0.4), False, True],
		[ 'qual', 'Boosted', 'jet2Tau21', 0, 1, 1, (0.90 if 'n-1' in args.cut else 0.40), (0.70 if 'n-1' in args.cut else 0.4), False, True],
		[ 'qual', 'Boosted', 'jet1Tau32', 0, 1, 1, (0.90 if 'n-1' in args.cut else 0.40), (0.70 if 'n-1' in args.cut else 0.4), False, True],
		[ 'qual', 'Boosted', 'jet2Tau32', 0, 1, 1, (0.90 if 'n-1' in args.cut else 0.40), (0.70 if 'n-1' in args.cut else 0.4), False, True],
		[ 'qual', 'Boosted', 'massAve', 0, 400, 10, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'HT', 700, 2000, 5, 0.90, 0.70, True, False],
		[ 'qual', args.boosted, 'NPV', 0, 50, 1, 0.90, 0.70, False, True],
		[ 'qual', 'Resolved', 'jet1Pt', 100, 1500, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet2Pt', 0, 1500, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet3Pt', 0, 500, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet4Pt', 0, 300, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'massAve', 0, 1000, 20, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'massAsym', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', 'Resolved', 'jetsBtag', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', 'Resolved', 'deltaEta', 0, 5, 1,  0.90, 0.70, True, False],	### n- 1
		[ 'tmp', 'Resolved', 'massAve', 0, 1500, 20, 0.90, 0.70, (False if (args.mass > 0 ) else True), False],
		[ 'tmp', 'Resolved', 'massAsym', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'tmp', 'Resolved', 'deltaEta', 0, 5, 1,  0.90, 0.70, True, False],	### n- 1

		[ 'qual', args.boosted, 'jetNum', 0, 5, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1NeutralHadronEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1NeutralEmEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1ChargedHadronEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1ChargedEmEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1ChargedMultiplicity', 0, 0.5, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1NumConst', 0, 50, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1Pt', 400, 1500, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1Eta', -3, 3, 5, 0.90, 0.70, False, True],
		[ 'qual', 'Boosted', 'jet1PrunedMass', 0, 1000, 10, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2NeutralHadronEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2NeutralEmEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2ChargedHadronEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2ChargedEmEnergyFrac', 0, 1, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2ChargedMultiplicity', 0, 0.5, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2NumConst', 0, 50, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2Pt', 400, 1500, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet2Eta', -3, 3, 5, 0.90, 0.70, False, True],
		[ 'qual', 'Boosted', 'jet2PrunedMass', 0, 1000, 10, 0.90, 0.70, True, False],
		[ 'qual', 'Boosted', 'jet1btagCSVv2', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', 'Boosted', 'jet2btagCSVv2', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', 'Resolved', 'jet1Btag', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', 'Resolved', 'jet2Btag', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', 'Resolved', 'jet3Btag', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', 'Resolved', 'jet4Btag', 0, 1, 1, 0.90, 0.70, False, False],
		[ 'qual', args.boosted, 'MET', 0, 200, 1, 0.90, 0.70, False, True],
		[ 'qual', args.boosted, 'METHT', 0, 200, 1, 0.90, 0.70, True, False],
		#[ 'qual', args.boosted, 'NPV_NOPUWeight', 0, 50, 2, 0.90, 0.70, False, False],
		[ 'qual', 'Resolved', 'neutralHadronEnergyFrac', '', '', 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'neutralEmEnergyFrac', '', '', 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'chargedHadronEnergyFrac', '', '', 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'chargedEmEnergyFrac', '', '', 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'chargedMultiplicity', 0, 0.5, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'numConst', '', '', 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet1Pt', 0, 1500, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet2Pt', 0, 1500, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet3Pt', 0, 1500, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet4Pt', 0, 1500, 5, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet1QGL', 0, 1, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet2QGL', 0, 1, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet3QGL', 0, 1, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'jet4QGL', 0, 1, 1, 0.90, 0.70, True, False],
		[ 'qual', 'Resolved', 'HT', 500, 3000, 10, 0.90, 0.70, True, False],

		#[ 'Norm', 'Boosted', 'jet1Tau1', '', '', 1, taulabX, taulabY, False, False],
		#[ 'Norm', 'Boosted', 'jet1Tau2', '', '', 1, taulabX, taulabY, False, False],
		#[ 'Norm', 'Boosted', 'jet1Tau3', '', '', 1, taulabX, taulabY, False, False],
		#[ 'Norm', 'Boosted', 'jet1RhoDDT', -6, 10, 5, taulabX, taulabY, False, False],
		#[ 'Norm', 'Boosted', 'jet2RhoDDT', -6, 10, 5, taulabX, taulabY, False, False],
		#[ 'Norm', 'Boosted', 'jet1Tau21DDT', 0, 1.3, 1, taulabX, taulabY, False, False],
		#[ 'Norm', 'Boosted', 'jet2Tau21DDT', 0, 1.3, 1, taulabX, taulabY, False, False],
		[ 'Norm', 'Boosted', 'jet1Tau21', '', '', 1, 0.3, 0.55, False, True ], # (True if 'n-1' in args.cut else False)],
		[ 'Norm', 'Boosted', 'jet2Tau21', '', '', 1, 0.3, 0.55, False, True ],
		[ 'Norm', 'Boosted', 'jet2Tau31', '', '', 1, taulabX, taulabY, False, True ], #(True if 'n-1' in args.cut else False)],
		[ 'Norm', 'Boosted', 'jet1Tau31', '', '', 1, taulabX, taulabY, False, True ], #(True if 'n-1' in args.cut else False)],
		[ 'Norm', 'Boosted', 'prunedMassAsym', '', '', 1, 0.40, 0.80, False, False],
		[ 'Norm', 'Boosted', 'deltaEtaDijet', '', '', 5, '', '', False, False],
		[ 'Norm', 'Boosted', 'jet1Tau32', '', '', 1, taulabX, taulabY, False, True],
		[ 'Norm', 'Boosted', 'jet2Tau32', '', '', 1, taulabX, taulabY, False, True],
		[ 'Norm', 'Boosted', 'jet1SubjetPtRatio', '', '', 1, '', '', True, False],
		[ 'Norm', 'Boosted', 'jet2SubjetPtRatio', '', '', 1, '', '', True, False],
		[ 'Norm', 'Boosted', 'subjetPtRatio', '', '', 1, '', '', True, False],
		[ 'Norm', 'Boosted', 'jet1CosThetaStar', '', '', 1, '', '', False, False],
		[ 'NormDATA', 'Resolved', 'massAsym', '', '', 1, '', '', False, False],
		[ 'NormDATA', 'Resolved', 'deltaEta', '', '', 1, '', '', False, False],
		[ 'Norm', 'Resolved', 'minDeltaR', '', '', 1, '', '', False, False],
		[ 'Norm', 'Resolved', 'minChi2', 0, 1, 1, '', '', False, False],
		[ 'Norm', 'Resolved', 'minMassAsym', '', '', 1, '', '', False, False],
		[ 'Norm', 'Resolved', 'deltaR', '', '', 1, '', '', False, False],
		[ 'Norm', 'Resolved', 'jet4Pt', 0, 300, 1, '', '', False, False],
		[ 'Norm', 'Resolved', 'cosThetaStar', '', '', 1, '', '', False, False],
		[ 'Norm', 'Resolved', 'xi', 0, 1.5, 1, '', '', False, False],


		[ 'simple', 'HT',  1000, '', '', False],
		[ 'simple', 'HT',  1000, '', '', True],
		[ 'simple', 'massAve_cutDijet',  massMaxX, '', '', False],
		[ 'simple', 'massAve_cutAsym',  massMaxX, '', '', False],
		[ 'simple', 'massAve_cutCosTheta',  massMaxX, '', '', False],
		[ 'simple', 'massAve_cutSubjetPtRatio',  massMaxX, '', '', False],
		#[ 'simple', 'massAve_cutSubjetPtRatio',  massMaxX, '', '', True ],
		[ 'simple', 'massAve_cutTau31',  massMaxX, '', '', False],
		[ 'simple', 'massAve_cutTau21',  massMaxX, '', '', False],
		

		]

	if 'all' in args.single: Plots = [ x[2:] for x in plotList if ( ( args.process in x[0] ) and ( x[1] in args.boosted ) )  ]
	else: Plots = [ y[2:] for y in plotList if ( ( args.process in y[0] ) and ( y[1] in args.boosted ) and ( y[2] in args.single ) )  ]

	if 'UDD312' in args.decay:
		massList = [ 80, 100, 120, 140, 160, 180, 200, 220, 240, 300, 350] 
		massWidthList = [10]*len(massList)
	else:
		massList = [ 80, 120, 180, 200, 220, 240, 260, 280, 300, 350] 
		massWidthList = [10]*len(massList)

	if 'test' in args.process:
		#plotDiffSample( bkgFiles['QCDPtAll'][0], bkgFiles['QCDPtAll'][0], 'CHS', 'Puppi', 'massAve_cutEffTrigger', '', '', '', '', 'PUMethod')
		plotDiffSample( bkgFiles['QCDPtAll'][0], bkgFiles['QCDPtAll'][0], 'QCDpythia', 'QCDmadgraph', args.single+'_'+args.cut , '', '', '', True, 'QCDSample', True)

	if 'CF' in args.process:
		plotCutFlow( signalFiles, bkgFiles, [ 'cutBestPair', 'massAsym', 'deltaEta', 'delta' ], 'massAve', 1000, True, True )

	if 'Scf' in args.process:
		plotSignalCutFlow(folder+'/RUNAnalysis_RPVStopStopToJets_UDD312_M-100_RunIIFall15MiniAODv2_v76x_v2p0_'+args.version+'.root', folder+'/RUNMiniBoostedAnalysis_'+args.grooming+'_RPVStopStopToJets_UDD312_M-100_'+args.version+'.root', (10 if 'high' in args.RANGE else 12), True, True )

	if 'signal' in args.process:
		plotSignalShape( args.single+'_'+args.cut, 
				( ( 10 if 'Resolved' in args.boosted else 5 ) if 'massAve' in args.single else 20), 
				( massList if 'Boosted' in args.boosted else range( 200, ( 1100 if '2CSVv2M' in args.cut else 1300), 100 )),
				massWidthList,
				False)

	if 'acc' in args.process:
		plotSignalAcceptance( 
				folder+'/RUNMini'+args.boosted+'Analysis'+( '' if 'Resolved' in args.boosted else '_'+args.grooming )+'_RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass)+'_Moriond17_80X_V2p4_'+args.version+'.root', 
				'massAve_'+args.cut, 
				( massList if 'Boosted' in args.boosted else range( 200, 1300, 100 )),
				massWidthList,
				False)


	for i in Plots:
		if args.process in '2D': 
			plot2D( dataFile, 'JetHT_Run2016', 1, args.grooming, i[0]+'_'+args.cut, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )
			for sig in signalFiles: plot2D( signalFiles[ sig ][0], 'RPVStopStopToJets_'+args.decay+'_M-'+str(args.mass), signalFiles[ sig ][1], args.grooming, i[0]+'_'+args.cut, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )
			for bkg in bkgFiles: plot2D( bkgFiles[ bkg ][0], bkg, bkgFiles[ bkg ][1], args.grooming, i[0]+'_'+args.cut, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )
			#plot2D( inputFileTTJets, 'TTJets', args.grooming, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )
			#plot2D( inputFileWJetsToQQ, 'WJets', args.grooming, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )
			#plot2D( inputFileZJetsToQQ, 'ZJets', args.grooming, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )

		elif '1D' in args.process:
			plotSignalBkg( signalFiles, bkgFiles, dataFile, 
					i[0]+'_'+args.cut, 
					i[0]+'_'+args.cut, 
					i[1], i[2], i[3], i[4], i[5], i[6], i[7], 
					addRatioFit=args.addFit )
		
		elif ( 'qual' in args.process ):
			plotQuality( dataFile, bkgFiles, 
					args.grooming, 
					( '' if args.miniTree else (args.boosted+'AnalysisPlots'+('' if 'pruned' in args.grooming else args.grooming)+'/'))+i[0]+'_'+args.cut, 
					i[0]+'_'+args.cut, 
					i[1], i[2], i[3], i[4], i[5], i[6], i[7], 
					fitRatio=args.addFit )
		
		elif 'Norm' in args.process:
			plotSignalBkg( signalFiles, bkgFiles, dataFile, 
					i[0]+('_'+args.cut if args.cut else ''), 
					i[0]+('_'+args.cut if args.cut else ''), 
					i[1], i[2], i[3], i[4], i[5], i[6], i[7], Norm=True )


		elif 'simple' in args.process:
			plotSimple( inputFileTTJets, 'TTJets', args.grooming, i[0], i[1], i[2], i[3], i[4] )
			plotSimple( inputFileWJetsToQQ, 'WJets', args.grooming, i[0], i[1], i[2], i[3], i[4] )
			plotSimple( inputFileZJetsToQQ, 'ZJets', args.grooming, i[0], i[1], i[2], i[3], i[4] )

		elif 'tmp' in args.process:
			tmpplotDiffSample( TFile(folder+'/RUNMiniResolvedAnalysis_'+( 'RPVStopStopToJets_UDD312_M-'+str(args.mass) if (args.mass > 0 ) else 'QCDPtAll' )+'_Moriond17_80X_V2p4_'+args.version+'_4jets.root'), 
				TFile(folder+'/RUNMiniResolvedAnalysis_'+( 'RPVStopStopToJets_UDD312_M-'+str(args.mass) if (args.mass > 0 ) else 'QCDPtAll' )+'_Moriond17_80X_V2p4_'+args.version+'.root'), 
				signalFiles[ args.mass ][0], i[0], i[1], i[2], i[3], i[4], i[5], i[6] ) 

	#if 'tmp' in args.process:
#		tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'QCD'+args.qcd+'All' ][0], signalFiles[ args.mass ][0], 'massAve', 10, 1500, '', '', True)
#		tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'TTJets' ][0], signalFiles[ args.mass ][0], 'massAve'+args.cut, 5, 400, '', '', True)
#		tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'TTJets' ][0], signalFiles[ args.mass ][0], 'jet1Tau21'+args.cut, 1, 1, 0.3, .55, False )
#		tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'TTJets' ][0], signalFiles[ args.mass ][0], 'jet2Tau21'+args.cut, 1, 1, .3, .55, False )
#		tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'TTJets' ][0], signalFiles[ args.mass ][0], 'prunedMassAsym'+args.cut, 1, 1, .3, .55, False )
#		tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'TTJets' ][0], signalFiles[ args.mass ][0], 'deltaEtaDijet'+args.cut, 5, 5, '', '', False )
#		tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'TTJets' ][0], signalFiles[ args.mass ][0], 'massAsym', 1, 1, .3, .55, False )
		#tmpplotDiffSample( bkgFiles[ 'QCD'+args.qcd+'All' ][0], bkgFiles[ 'QCD'+args.qcd+'All' ][0], signalFiles[ args.mass ][0], 'deltaEta', 1, 5, '', '', False )
		#tmpplotDiffSample( bkgFiles[ bkg ][0], bkg, 'jet1PrunedMass'+args.cut, 5, 400, '', '', True )
		#tmpplotDiffSample( bkgFiles[ bkg ][0], bkg, 'jet2PrunedMass'+args.cut, 5, 400, '', '', True )
