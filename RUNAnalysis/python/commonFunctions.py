#!/usr/bin/env python

import sys,os,time, re
from math import *
from string import *
from array import array
from ROOT import * 
import numpy as np

##### Support functions
def checkLumi( Run, Lumi, NumEvent):
	"""docstring for checkLumi"""
	result = False
	allEvents = 'Run: '+str(Run)+' LumiSection: '+str(Lumi)+' Event: '+str(NumEvent)
	with open('boostedEventsRPV100tojj.txt') as f:
		lines = f.readlines()
		for i in lines: 
			if allEvents == i: result = True

	return result

def find_nearest(array,value):
	idx = (np.abs(array-value)).argmin()
	return idx

def getTree(filename, treename):
	hfile = TFile(filename)
	if not hfile.IsOpen():
		print "** can't open file %s" % filename
		sys.exit()
	tree = hfile.Get(treename)
	if tree == None:
		print "** can't find tree %s" % treename
		sys.exit()
	entries = tree.GetEntriesFast()
	return (hfile, tree, entries)

#boostedMassAveBins = array( 'd', [ 0, 3, 6, 9, 12, 16, 19, 23, 26, 30, 34, 39, 43, 47, 52, 57, 62, 67, 72, 78, 83, 89, 95, 102, 108, 115, 122, 129, 137, 144, 153, 161, 170, 179, 188, 197, 207, 218, 228, 240, 251, 263, 275, 288, 301, 315, 329, 344, 359, 375, 391, 408, 425, 443, 462, 482, 502 ] )
#boostedMassAveBinSize = array( 'd', [ 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 19, 19, 20, 21] )
boostedMassAveBins = array( 'd', [ 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360, 375, 395, 415, 435, 455, 475, 495, 515 ] )
boostedMassAveBinSize = array( 'd', [ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 20, 20, 20, 20, 20, 20, 20, 20] )

##### from MC true, sigma from fit.
massWidthList = [8.56280909196305, 8.445039648677378, 8.950556420141245, 9.860254339542022, 8.814786972730516, 10.021433248818914, 10.392360104091987, 9.435770844457956, 10.268425520508536, 10.45176971177987, 12.86644189449206, 10.084924444431165, 12.431737065699405, 10.809084324420656, 12.94592267653858, 15.762703291273564]

###### Rebin 2D plots
def Rebin2D( h1, rebinx, rebiny ):
	"""docstring for Rebin2D"""

	tmph1 = h1.Clone()
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
	return h1

def getHistoFromTree( fileName, treeName, plotVar, cuts, histo, numentries=False, skipEvents=0 ):
	"""docstring for getHistoFromTree"""

	chain = TChain( treeName )
	chain.Add( fileName ) 
	if not numentries: numEntries = chain.GetEntries()
	else: numEntries = int( chain.GetEntries())
	print '|---> Plotting: '+plotVar+'>>'+str(histo.GetName()), numEntries, chain.GetEntries(), cuts 
	chain.Draw( plotVar+'>>'+str(histo.GetName()), cuts, 'goff', numEntries, skipEvents ) ### goff no graphics generated

	return histo

def get2DHistoFromTree( fileName, treeName, plotVar1, plotVar2, cuts, histo, numEntries=False, skipEvents=0 ):
	"""docstring for getHistoFromTree"""

	chain = TChain( treeName )
	chain.Add( fileName ) 
	if not numEntries: numEntries = chain.GetEntriesFast()
	else: numEntries = int( chain.GetEntriesFast())
	print '|---> Plotting: '+plotVar1+':'+plotVar2+'>>'+str(histo.GetName()), cuts
	chain.Draw( plotVar2+':'+plotVar1+'>>'+str(histo.GetName()), cuts, 'goff', numEntries, skipEvents  )

	return histo
