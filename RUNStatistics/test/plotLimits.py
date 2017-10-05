#!/usr/bin/env python
'''
File: DrawHistogram.py
Author: Alejandro Gomez Espinosa
Email: gomez@physics.rutgers.edu
Description: My Draw histograms. Check for options at the end.
'''

from ROOT import *
import time, os, math, sys
import argparse
from collections import OrderedDict
try:
	from RUNA.RUNAnalysis.commonFunctions import *
	from RUNA.RUNAnalysis.scaleFactors import *
	import RUNA.RUNAnalysis.CMS_lumi as CMS_lumi 
	import RUNA.RUNAnalysis.tdrstyle as tdrstyle
except ImportError:
	sys.path.append('../python')
	from commonFunctions import *
	from scaleFactors import *
	import CMS_lumi as CMS_lumi 
	import tdrstyle as tdrstyle

gROOT.SetBatch(kTRUE);
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
#gStyle.SetTitleFont(42, "XYZ")
#gStyle.SetTitleSize(0.06, "XYZ")
#gStyle.SetLabelFont(42, "XYZ")
#gStyle.SetLabelSize(0.05, "XYZ")
#gStyle.SetCanvasBorderMode(0)
#gStyle.SetFrameBorderMode(0)
#gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.13)
gStyle.SetPadRightMargin(0.05)
#gStyle.SetPadTopMargin(0.05)
#gStyle.SetPadBottomMargin(0.15)
#gROOT.ForceStyle()


def plotLimits( listMasses  ):
	"""docstring for plotLimits"""
	
	masses = array('d')
	masses_exp = array('d')
	xs_theory = array('d')
	xs_obs_limits = array('d')
	xs_exp_limits = array('d')
	xs_exp_limits_1sigma = array('d')
	xs_exp_limits_1sigma_up = array('d')
	xs_exp_limits_2sigma = array('d')
	xs_exp_limits_2sigma_up = array('d')

	for mass in listMasses:
		masses.append( mass )
		masses_exp.append( mass )

	if args.theta:
		thetaExpectedFile = open('thetaFiles/theta_expected_'+args.version+'.txt','r')
		thetaObservedFile = open('thetaFiles/theta_observed_'+args.version+'.txt','r')
		for line in thetaExpectedFile:
			li=line.strip()
			if not li.startswith("#"):
				mass = float(line.split()[0])
				XS =  search( dictXS, 'RPVStopStopToJets_'+args.decay+'_M-'+str(mass) )
				xs_theory.append( XS )
				xs_exp_limits.append(float(line.split()[1]) *XS )
				xs_exp_limits_2sigma.append(float(line.split()[2]) *XS )
				xs_exp_limits_2sigma_up.append(float(line.split()[3]) *XS )
				xs_exp_limits_1sigma.append(float(line.split()[4]) *XS )
				xs_exp_limits_1sigma_up.append(float(line.split()[5]) *XS )
		for line in thetaObservedFile:
			li=line.strip()
			if not li.startswith("#"):
				mass = float(line.split()[0])
				XS =  search( dictXS, 'RPVStopStopToJets_'+args.decay+'_M-'+str(mass) )
				xs_obs_limits.append(float(line.split()[1]) *XS )

	else:
		for mass in listMasses:

			XS = search( dictXS, 'RPVStopStopToJets_'+args.decay+'_M-'+str(mass) )
			xs_theory.append( XS )

			if 'final' in args.boosted: combineFile = "higgsCombine_RPVStopStopToJets_"+args.decay+"_M-"+str(mass)+'_final_'+args.version+'.'+args.method+".mH120.root"
			else: combineFile = "higgsCombine_RPVStopStopToJets_"+args.decay+"_M-"+str(mass)+'_'+args.boosted+'_'+args.sys+'_'+args.version+'.'+args.method+".mH120.root"
			tmpFile, tmpTree, tmpEntries = getTree( combineFile, "limit" )
			for i in xrange(tmpEntries):
				tmpTree.GetEntry(i)
				tmp = round( tmpTree.quantileExpected, 2)
				if tmp == 0.03: xs_exp_limits_2sigma.append( tmpTree.limit * XS )
				if tmp == 0.16: xs_exp_limits_1sigma.append( tmpTree.limit * XS )
				if tmp == 0.5: 
					xs_exp_limits.append( tmpTree.limit * XS )
					print mass, round( tmpTree.limit * XS, 2)
				if tmp == 0.84: xs_exp_limits_1sigma_up.append( tmpTree.limit * XS )
				if tmp == 0.98: xs_exp_limits_2sigma_up.append( tmpTree.limit * XS ) 
				if tmp == -1: xs_obs_limits.append( tmpTree.limit * XS )

	for i in range(0,len(masses)):
		masses_exp.append( masses[len(masses)-i-1] )
		xs_exp_limits_1sigma.append( xs_exp_limits_1sigma_up[len(masses)-i-1] )
		xs_exp_limits_2sigma.append( xs_exp_limits_2sigma_up[len(masses)-i-1] )

	if args.addComparison: 
		#xs_cdf_masses = array( 'd', [ 50, 100, 100, 50 ] )
		#xs_cdf = array( 'd',  [ 0.1, 0.1, 100000, 100000 ] )
		xs_cdf_masses = array( 'd', [ 70, 90, 100, 125, 150 ] )
		updateLumiCDF = TMath.Sqrt(6.6/lumi)
		xs_cdf = array( 'd',  [ 75*updateLumiCDF*(6430/100),  8.2*updateLumiCDF*(2369/26), 11*updateLumiCDF*(1521.11/15), 14*updateLumiCDF*(574.98/4.4), 37*updateLumiCDF*(249.40/1.5) ] )
		graph_xs_cdf = TGraph(len(xs_cdf_masses),xs_cdf_masses,xs_cdf)
		graph_xs_cdf.SetLineWidth(3)
		#graph_xs_cdf.SetLineStyle(8)
		graph_xs_cdf.SetLineColor(kCyan)
		graph_xs_cdf.SetFillColorAlpha(kCyan, 0.1)

		xs_atlasbj8TeV_masses = array( 'd', [ 100, 125, 150, 175, 200, 225, 250, 275, 300  ] )
		updateLumiAtlas8TeV = TMath.Sqrt(17.4/lumi)
		xs_atlasbj8TeV = array( 'd',  [ 140*updateLumiAtlas8TeV*(1521.11/559.757),  50*updateLumiAtlas8TeV*(574.98/197.12),  28*updateLumiAtlas8TeV*(249.40/80.26),  10.2*updateLumiAtlas8TeV*(121.41/36.79),  7*updateLumiAtlas8TeV*(64.50/18.52),  5*updateLumiAtlas8TeV*(36.38/9.90),  3*updateLumiAtlas8TeV*(21.59/5.57),  2*updateLumiAtlas8TeV*(13.32/3.27),  1.5*updateLumiAtlas8TeV*(8.51/1.99) ] )
		graph_xs_atlasbj8TeV = TGraph(len(xs_atlasbj8TeV_masses),xs_atlasbj8TeV_masses,xs_atlasbj8TeV)
		graph_xs_atlasbj8TeV.SetLineWidth(3)
		#graph_xs_atlasbj8TeV.SetLineStyle(8)
		graph_xs_atlasbj8TeV.SetLineColor(kViolet-5)

		#xs_atlas_masses = array( 'd', [ 250, 350, 350, 250  ] )
		#xs_atlas = array('d', [ 0.1, 0.1, 100000, 100000 ])
		'''
		updateLumi = TMath.Sqrt(15.4/lumi)
		if 'Boosted' in args.boosted:
			tmpAtlasmass = [ 250, 275, 300  ]
			tmpAtlas = [ 6*updateLumi, 4*updateLumi, 3*updateLumi ]
		else:
			tmpAtlasmass = range( 250, 650, 100 ) 
			tmpAtlas = [ 6*updateLumi, 1.08*updateLumi, 0.9*updateLumi, 0.6*updateLumi, 0.3*updateLumi ] 
		'''

		xs_atlas_masses = array( 'd', range( 100, 900, 100 ))
		xs_atlas = array('d', [ 600, 18, 4, 1.5, 0.6, 0.4, 0.25, 0.08])
		graph_xs_atlas = TGraph(len(xs_atlas_masses),xs_atlas_masses,xs_atlas)
		graph_xs_atlas.SetLineWidth(3)
		#graph_xs_atlas.SetLineStyle(8)
		graph_xs_atlas.SetLineColor(kMagenta-2)
		graph_xs_atlas.SetFillColorAlpha(kGreen-2, 0.1)

		xs_atlasbj13TeV_masses = array( 'd', range( 100, 900, 100) )
		xs_atlasbj13TeV = array( 'd',  [  400, 6, 1.5, 0.4, 0.4, 0.2, 0.1, 0.1 ] )
		graph_xs_atlasbj13TeV = TGraph(len(xs_atlasbj13TeV_masses),xs_atlasbj13TeV_masses,xs_atlasbj13TeV)
		graph_xs_atlasbj13TeV.SetLineWidth(3)
		#graph_xs_atlasbj13TeV.SetLineStyle(8)
		graph_xs_atlasbj13TeV.SetLineColor(kRed-5)

		updateLumi8TeV = TMath.Sqrt(12.4/lumi)
		if 'Boosted' in args.boosted:
			tmpCMS8TeVmass = [ 200, 250, 300 ] 
			tmpCMS8TeV = [ 4.*updateLumi8TeV*(64.5/18.52), 2*updateLumi8TeV*(21.59/5.57), 1.1*updateLumi8TeV*(8.51/1.99) ]
		elif 'Resolved' in args.boosted:
			tmpCMS8TeVmass = range( 200, 1100, 100 )
			if '312' in args.decay: 
				tmpCMS8TeV = [ 4.*updateLumi8TeV*(64.5/18.52), 
						1.1*updateLumi8TeV*(8.51/1.99),
						0.4*updateLumi8TeV*(1.83/0.356),
						0.27*updateLumi8TeV*(0.51/0.0855),
						0.11*updateLumi8TeV*(0.17/0.024),
						0.1*updateLumi8TeV*(0.06/0.00811),
						0.08*updateLumi8TeV*(0.028/0.00289),
						0.06*updateLumi8TeV*(0.012/0.0010950),
						0.05*updateLumi8TeV*(0.00615/0.0004354),
						]
			else: 
				tmpCMS8TeV = [ 6.*updateLumi8TeV*(64.5/18.52), 
						1.*updateLumi8TeV*(8.51/1.99),
						0.3*updateLumi8TeV*(1.83/0.356),
						0.2*updateLumi8TeV*(0.51/0.0855),
						0.15*updateLumi8TeV*(0.17/0.024),
						0.12*updateLumi8TeV*(0.06/0.00811),
						0.10*updateLumi8TeV*(0.028/0.00289),
						0.08*updateLumi8TeV*(0.012/0.0010950),
						0.07*updateLumi8TeV*(0.00615/0.0004354),
						]

		xs_cms8TeV_masses = array( 'd', tmpCMS8TeVmass )
		xs_cms8TeV = array('d', tmpCMS8TeV )
		graph_xs_cms8TeV = TGraph(len(xs_cms8TeV_masses),xs_cms8TeV_masses,xs_cms8TeV)
		graph_xs_cms8TeV.SetLineWidth(3)
		graph_xs_cms8TeV.SetLineColor(kBlue)
		graph_xs_cms8TeV.SetFillColorAlpha(kBlue, 0.1)
		#graph_xs_cms8TeV.SetLineStyle(8)
		#graph_xs_cms8TeV.SetLineWidth(-2002)
		#graph_xs_cms8TeV.SetFillStyle(3004)

		xs_cms13TeV2015_masses = array( 'd', [ 80, 90, 100, 110, 120, 130, 140, 150, 170, 180, 190, 210, 220, 230, 240, 300 ] )
		xs_cms13TeV2015 = array( 'd',  [ 1470.6, 994.79, 568.93, 481.14, 327.39, 202.15, 153.35, 155.39, 76.97, 58.28, 60.89, 33.72, 34.41, 25.2, 23.68, 103.52 ] )
		graph_xs_cms13TeV2015 = TGraph(len(xs_cms13TeV2015_masses),xs_cms13TeV2015_masses,xs_cms13TeV2015)
		graph_xs_cms13TeV2015.SetLineWidth(3)
		#graph_xs_cms13TeV2015.SetLineStyle(8)
		graph_xs_cms13TeV2015.SetLineColor(kBlack-2)
		graph_xs_cms13TeV2015.SetFillColorAlpha(kBlack-2, 0.1)

	graph_xs_th = TGraph(len(masses),masses,xs_theory)
	shadow_graph_xs_th = graph_xs_th.Clone()
	graph_xs_th.SetLineWidth(3)
	graph_xs_th.SetLineStyle(2)
	graph_xs_th.SetLineColor(kMagenta)
	shadow_graph_xs_th.SetLineWidth(10)
	shadow_graph_xs_th.SetLineColorAlpha(kMagenta, 0.15);


	graph_exp_2sigma = TGraph(len(masses_exp),masses_exp,xs_exp_limits_2sigma)
	graph_exp_2sigma.SetFillColor(kYellow)

	graph_exp_1sigma = TGraph(len(masses_exp),masses_exp,xs_exp_limits_1sigma)
	graph_exp_1sigma.SetFillColor(kGreen+1)

	graph_exp = TGraph(len(masses),masses,xs_exp_limits) 
	#graph_exp.SetMarkerStyle(24)
	graph_exp.SetLineWidth(3)
	graph_exp.SetLineStyle(2)
	graph_exp.SetLineColor(4)

	graph_obs = TGraph(len(masses),masses,xs_obs_limits)
	graph_obs.SetMarkerStyle(20)
	graph_obs.SetLineWidth(3)
	graph_obs.SetLineStyle(1)
	graph_obs.SetLineColor(1)

	c = TCanvas("c", "",800,600)
	c.cd()

	if args.addComparison: 
		legend = TLegend(.40,.60,.90,.88)
		legend.SetTextSize(0.025)
	else: 
		legend = TLegend(.45,.60,.90,.88)
		legend.SetTextSize(0.03)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
	legend.SetHeader('95% CL upper limits')

	graph_exp_2sigma.GetXaxis().SetTitle("Resonance mass [GeV]")
	graph_exp_2sigma.GetYaxis().SetTitle("(pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow "+("qq" if 'UDD312' in args.decay else "bq" )+") #sigma #times #it{B} [pb] ")
	graph_exp_2sigma.GetYaxis().SetTitleOffset(1.1)
	if 'Boosted' in args.boosted: graph_exp_2sigma.GetYaxis().SetRangeUser(3,1e+04)
	elif 'Resolved' in args.boosted: graph_exp_2sigma.GetYaxis().SetRangeUser(0.01,100)
	else: graph_exp_2sigma.GetYaxis().SetRangeUser(0.01,10000)
	#graph_exp_2sigma.GetXaxis().SetNdivisions(1005)

	graph_exp_2sigma.Draw("AF")
	graph_exp_1sigma.Draw("F")
	graph_exp.Draw("L")
        #shadow_graph_xs_th.Draw("L")
        graph_xs_th.Draw("L")
	#graph_xs_cms13TeV2015.Draw("L")
	if not args.addComparison: graph_obs.Draw("LP")

        #legend.AddEntry(graph_xs_th,"RPV #lambda_{312}^{''} (#tilde{t} #rightarrow qq) cross section","l")
        legend.AddEntry(graph_xs_th,"Top squark pair production #lambda_{"+("312" if '312' in args.decay else '323')+"}^{''} (#tilde{t} #rightarrow "+("qq)" if '312' in args.decay else 'bq)' ),"l")
	if not args.addComparison: legend.AddEntry(graph_obs,"Observed limit","lp")
	legend.AddEntry(graph_exp,"Expected limit","lp")
	legend.AddEntry(graph_exp_1sigma,"Expected #pm 1#sigma","F")
	legend.AddEntry(graph_exp_2sigma,"Expected #pm 2#sigma","F")
	#legend.AddEntry(graph_xs_cms13TeV2015,"Expected limit 2015 data","L")
	if args.addComparison:
		if 'Boosted' in args.boosted:
			if '312' in args.decay:
				graph_xs_cms8TeV.Draw("L")
				graph_xs_atlas.Draw("L")
				graph_xs_cdf.Draw("L")
				graph_xs_cms13TeV2015.Draw("L")
				legend.AddEntry(graph_xs_atlas,"Expected limit, Atlas 13 TeV (resolved) - 36.7 fb^{-1}","L")
				legend.AddEntry(graph_xs_cms8TeV,"Expected limit, CMS 8 TeV (resolved) - 12.4 fb^{-1}","L")
				legend.AddEntry(graph_xs_cdf,"Expected limit, CDF 1.96 TeV (resolved) - 6.6 fb^{-1}","L")
				legend.AddEntry(graph_xs_cms13TeV2015,"Expected limit CMS 2015 data","L")
			else:
				graph_xs_atlasbj8TeV.Draw("L")
				legend.AddEntry(graph_xs_atlasbj8TeV,"Expected limit, Atlas 8 TeV (boosted - #lambda_{323}) - 17.4 fb^{-1}","l")
				graph_xs_atlasbj13TeV.Draw("L")
				legend.AddEntry(graph_xs_atlasbj13TeV,"Expected limit, Atlas 13 TeV (resolved - #lambda_{323}) - 36.7 fb^{-1}","l")

		else: 
			graph_xs_cms8TeV.Draw("L")
			legend.AddEntry(graph_xs_cms8TeV,"Excluded mass regions by CMS 8 TeV (resolved) - 19.4 fb^{-1}","L")
			if '312' in args.decay:
				graph_xs_atlas.Draw("L")
				legend.AddEntry(graph_xs_atlas,"Excluded mass regions by Atlas 13 TeV (resolved) - 36.7 fb^{-1}","L")
			else:
				graph_xs_atlasbj8TeV.Draw("L")
				legend.AddEntry(graph_xs_atlasbj8TeV,"Expected limit, Atlas 8 TeV (boosted - #lambda_{323}) - 17.4 fb^{-1}","l")
				graph_xs_atlasbj13TeV.Draw("L")
				legend.AddEntry(graph_xs_atlasbj13TeV,"Expected limit, Atlas 13 TeV (resolved - #lambda_{323}) - 36.7 fb^{-1}","l")
    	legend.Draw()

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(c, 4, 0)

	c.SetLogy()
	if 'final' in args.boosted: c.SetLogx()
	#fileName = 'xs_limit_%s_%s.%s'%(args.method,args.final_state + ( ('_' + args.postfix) if args.postfix != '' else '' ), args.fileFormat.lower())
	fileName = 'xs_limit_RPVStop_'+args.decay+''+args.boosted+'_'+args.sys+'_'+args.method+'_'+args.version+'.'+args.ext
	if args.theta: fileName = fileName.replace('limit', 'limit_theta')
	if 'gaus' in args.process: fileName = fileName.replace('limit', 'limit_gaus')
	if args.addComparison: fileName = fileName.replace('limit', 'limit_comparison')
	print 'Processing.......', fileName
	c.SaveAs( 'Plots/'+fileName )

def compareLimits( listMasses, diffVersions ):
	"""docstring for compareLimits"""
	
	masses = array('d', listMasses )
	xs_theory = array('d')

	diffLimitsDict = OrderedDict()
	for ver in diffVersions:
		diffLimitsDict[ ver+'_' ] = array('d')
	
	for mass in listMasses:

		XS = search( dictXS, 'RPVStopStopToJets_'+args.decay+'_M-'+str(mass) )
		xs_theory.append( XS )

		for ver in diffVersions:
			if 'Resolved' in args.boosted: 
				if not ( "doubleGaus" in ver ):
					combineFile = "higgsCombine_RPVStopStopToJets_"+args.decay+"_M-"+str(mass)+args.boosted+'_'+ver+'_'+args.version+'.'+args.method+".mH120.root"
				else:
					if mass in [600, 700, 800]:
						combineFile = "higgsCombine_RPVStopStopToJets_"+args.decay+"_M-"+str(mass)+args.boosted+'_'+ver+'_'+args.version+'.'+args.method+".mH120.root"
					else: break
			else: combineFile = "higgsCombine_RPVStopStopToJets_"+args.decay+"_M-"+str(mass)+args.boosted+'_'+args.sys+'_'+args.version+ver+'.'+args.method+".mH120.root"
			tmpFile, tmpTree, tmpEntries = getTree( combineFile, "limit" )
			for i in xrange(tmpEntries):
				tmpTree.GetEntry(i)
				tmp = round( tmpTree.quantileExpected, 2)
				#if tmp == 0.03: diffLimits[ ver+'_' ]_2sigma.append( tmpTree.limit * XS )
				#if tmp == 0.16: diffLimits[ ver+'_' ]_1sigma.append( tmpTree.limit * XS )
				if tmp == 0.5: diffLimitsDict[ ver+'_' ].append( tmpTree.limit * XS )
				#if tmp == 0.84: diffLimits[ ver+'_' ]_1sigma_up.append( tmpTree.limit * XS )
				#if tmp == 0.98: diffLimits[ ver+'_' ]_2sigma_up.append( tmpTree.limit * XS ) 
				#if tmp == -1: xs_obs_limits.append( tmpTree.limit * XS )

	legend = TLegend(.45,.60,.90,.88)
	legend.SetTextSize(0.03)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)

	graph_xs_th = TGraph(len(masses),masses,xs_theory)
	shadow_graph_xs_th = graph_xs_th.Clone()
	graph_xs_th.SetLineWidth(3)
	graph_xs_th.SetLineStyle(2)
	graph_xs_th.SetLineColor(kMagenta)
	shadow_graph_xs_th.SetLineWidth(10)
	shadow_graph_xs_th.SetLineColorAlpha(kMagenta, 0.15);
        legend.AddEntry(graph_xs_th,"Top squark pair production #lambda_{"+("312" if '312' in args.decay else '323')+"}^{''} (#tilde{t} #rightarrow "+("qq)" if '312' in args.decay else 'bq)' ),"l")

	diffLimitsGraphDict = OrderedDict()
	dummy=1
	for l in diffLimitsDict:
		if 'doubleGaus' in l: tmpMasses = array( 'd', [ 600, 700, 800 ])
		else: tmpMasses = masses
		diffLimitsGraphDict[ l ] = TGraph( len(tmpMasses), tmpMasses, diffLimitsDict[l] ) 
		diffLimitsGraphDict[ l ].SetLineWidth(3)
		diffLimitsGraphDict[ l ].SetLineStyle(2)
		diffLimitsGraphDict[ l ].SetLineColor(dummy)
		legend.AddEntry(diffLimitsGraphDict[ l ], ( 'Nominal' if '_' == l else l.replace('_', '') ),"lp")
		dummy+=1

	c = TCanvas("c", "",800,600)
	c.cd()

	diffLimitsGraphDict[ diffVersions[-2]+'_' ].GetXaxis().SetTitle("Resonance mass [GeV]")
	diffLimitsGraphDict[ diffVersions[-2]+'_' ].GetYaxis().SetTitle("#sigma #times #it{B} [pb]")
	diffLimitsGraphDict[ diffVersions[-2]+'_' ].GetYaxis().SetTitleOffset(1.1)
	if 'Boosted' in args.boosted: diffLimitsGraphDict[ diffVersions[-2]+'_' ].GetYaxis().SetRangeUser(3,1e+04)
	elif 'Resolved' in args.boosted: diffLimitsGraphDict[ diffVersions[-2]+'_' ].GetYaxis().SetRangeUser(0.01,100)
	else: diffLimitsGraphDict[ diffVersions[-2]+'_' ].GetYaxis().SetRangeUser(0.01,10000)

	diffLimitsGraphDict[ diffVersions[-2]+'_' ].Draw("AL")
	for l in diffLimitsDict: diffLimitsGraphDict[l].Draw("L")
	graph_xs_th.Draw("L")
    	legend.Draw()

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(c, 4, 0)

	c.SetLogy()
	fileName = 'xs_limit_RPVStop_'+args.decay+'_'+args.boosted+'_'+args.sys+'_'+args.method+'_diff_'+(''.join(diffVersions))+'_'+args.version+'.'+args.ext
	print 'Processing.......', fileName
	c.SaveAs( 'Plots/'+fileName )

def compareFinalLimits( listMasses, categories ):
	"""docstring for compareFinalLimits"""
	
	masses = array('d', listMasses )
	xs_theory = array('d')

	for mass in listMasses:
		XS = search( dictXS, 'RPVStopStopToJets_'+args.decay+'_M-'+str(mass) )
		xs_theory.append( XS )

	legend = TLegend(.45,.60,.90,.88)
	legend.SetTextSize(0.03)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)

	diffLimitsGraphDict = OrderedDict()
	dummy=1
	for cat in [ 'Boosted', 'Resolved', 'Combined' ]: 
		signalStrengh = []
		if 'Boosted' in cat:
			name=categories[0]
			if '312' in args.decay:  listMass = [ 80, 100, 120, 140, 160, 180, 200, 220, 240, 300, 350 ]
			else: listMass = [ 80, 100, 120, 140, 160, 180, 200, 220, 240, 280, 300 ]
		elif 'Resolved' in cat:
			name=categories[1]
			if '312' in args.decay: listMass = [ 200, 220, 240 ] + range( 300, 1050, 50 ) + range( 1100, 1400, 100 ) 
			else: listMass = [ 200, 220, 240, 260, 280 ] + range( 300, 1050, 50 ) #+ range( 1100, 1400, 100 ) 
		else:
			name=categories[1]
			if '312' in args.decay: listMass = range( 80, 260, 20 ) + range( 300, 1050, 50 ) + range( 1100, 1400, 100 ) 
			else: listMass = range( 80, 300, 20 ) + range( 300, 1050, 50 ) #+ range( 1100, 1400, 100 ) 
		
		for imass in listMass:
			combineFile = "higgsCombine_RPVStopStopToJets_"+args.decay+"_M-"+str(imass)+'_'+name+'.'+args.method+".mH120.root"
			XS = search( dictXS, 'RPVStopStopToJets_'+args.decay+'_M-'+str(imass) )

			tmpFile, tmpTree, tmpEntries = getTree( combineFile, "limit" )
			for i in xrange(tmpEntries):
				tmpTree.GetEntry(i)
				tmp = round( tmpTree.quantileExpected, 2)
				#if tmp == 0.03: signalStrengh_2sigma.append( tmpTree.limit * XS )
				#if tmp == 0.16: signalStrengh_1sigma.append( tmpTree.limit * XS )
				if tmp == 0.5: signalStrengh.append( tmpTree.limit * XS )
				#if tmp == 0.84: signalStrengh_1sigma_up.append( tmpTree.limit * XS )
				#if tmp == 0.98: signalStrengh_2sigma_up.append( tmpTree.limit * XS ) 
				#if tmp == -1: xs_obs_limits.append( tmpTree.limit * XS )

		diffLimitsGraphDict[ cat ] = TGraph( len(listMass), array( 'd', listMass), array( 'd', signalStrengh) ) 
		diffLimitsGraphDict[ cat ].SetLineWidth(3)
		diffLimitsGraphDict[ cat ].SetLineStyle(2)
		diffLimitsGraphDict[ cat ].SetLineColor(dummy)
		legend.AddEntry(diffLimitsGraphDict[ cat ], cat,"lp")
		dummy+=1

	graph_xs_th = TGraph(len(masses),masses,xs_theory)
	shadow_graph_xs_th = graph_xs_th.Clone()
	graph_xs_th.SetLineWidth(3)
	graph_xs_th.SetLineStyle(2)
	graph_xs_th.SetLineColor(kMagenta)
	shadow_graph_xs_th.SetLineWidth(10)
	shadow_graph_xs_th.SetLineColorAlpha(kMagenta, 0.15);
        legend.AddEntry(graph_xs_th,"Top squark pair production #lambda_{"+("312" if '312' in args.decay else '323')+"}^{''} (#tilde{t} #rightarrow "+("qq)" if '312' in args.decay else 'bq)' ),"l")

	c = TCanvas("c", "",800,600)
	c.cd()

	graph_xs_th.GetXaxis().SetTitle("Resonance mass [GeV]")
	graph_xs_th.GetYaxis().SetTitle("(pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow "+("qq" if 'UDD312' in args.decay else "bq" )+") #sigma #times #it{B} [pb] ")
	graph_xs_th.GetYaxis().SetTitleOffset(1.1)
	graph_xs_th.GetYaxis().SetRangeUser(0.01,10000)
	graph_xs_th.Draw("AL")
	for l in diffLimitsGraphDict: diffLimitsGraphDict[l].Draw("L")
    	legend.Draw()

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(c, 4, 0)

	c.SetLogy()
	fileName = 'xs_limit_RPVStop_'+args.decay+'_'+args.boosted+'_'+args.method+'_diff_'+(''.join(categories))+'.'+args.ext
	print 'Processing.......', fileName
	c.SaveAs( 'Plots/'+fileName )

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--proc', dest='process', action='store', default='1D', help='Process to draw, example: 1D, 2D, MC.' )
	parser.add_argument('-d', '--decay', dest='decay', action='store', default='UDD312', help='Decay, example: UDD312, UDD323.' )
	parser.add_argument('-b', '--boosted', dest='boosted', action='store', default='Boosted', help='Boosted or non version, example: Boosted' )
	parser.add_argument('-v', '--version', dest='version', action='store', default='v05', help='Version of the root files' )
	parser.add_argument('-t', '--theta', dest='theta', action='store', type=float, default=False, help='Input from theta or not.' )
	parser.add_argument('-l', '--lumi', dest='lumi', action='store', type=float, default=149.9, help='Luminosity, example: 1.' )
	parser.add_argument('-e', '--extension', dest='ext', action='store', default='png', help='Extension of plots.' )
	parser.add_argument('-s', '--sys', dest='sys', action='store', default='NOSys', help='With systematics or not.' )
	parser.add_argument('-m', '--method', dest='method', action='store', default='AsymptoticLimits', help='Limit method: Asymptotic, HybridNew' )
	parser.add_argument('-a', '--addComparison', dest='addComparison', action='store', type=bool, default=False, help='Adding comparisons to plots.' )

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	
	CMS_lumi.extraText = "Preliminary"
	CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 1 ) )+" fb^{-1}"
	#if args.addComparison: 
	lumi = args.lumi/1000
	#if 'gaus' in args.process: listMass = range( 80, 360, 10 )
	if 'Resolved' in args.boosted: 
		if '312' in args.decay:  
			listMass = [ 200, 220, 240 ] + range( 300, 1050, 50 ) + range( 1100, 1400, 100 ) 
		else: listMass = [ 200, 220, 240, 260, 280 ] + range( 300, 1050, 50 ) #+ [1100, 1200 ] 


	elif "Boosted" in args.boosted:
		if '312' in args.decay:  listMass = [ 80, 100, 120, 140, 160, 180, 200, 220, 240, 300, 350 ]
		else: listMass = [ 80, 100, 120, 140, 160, 180, 200, 220, 240, 280, 300 ]

	else:
		if '312' in args.decay:  
			listMass = range( 80, 260, 20) + range( 300, 1050, 50 ) + [ 1100, 1200 ]
		else: 
			listMass = range( 80, 300, 20) + range( 300, 1050, 50 ) + [ 1100, 1200 ]

	if 'compare' in args.process: compareLimits( listMass, (['_1sigma', '_2sigma', '_4sigma', '_10sigma' ] if 'Boosted' in args.boosted else [ 'delta', 'delta_massWindow', 'delta_doubleGaus' ] ) ) 
	elif 'full' in args.process: compareFinalLimits( listMass, [ 'Boosted_jet1Tau32_Bin5_v09p2', 'Resolved_delta_'+( '2CSVv2L_' if 'UDD323' in args.decay else '' )+'massWindow_v09p1', 'Combined' ] ) 
	else: plotLimits( listMass  )

