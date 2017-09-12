#!/usr/bin/env python
################################################################################################
### Modify script from 
### https://gitlab.cern.ch/CMSDIJET/StatisticalTools/blob/master/scripts/checkWorkspace.py
#################################################################################################

import sys, os, copy, re
from ROOT import * 
from array import array

gROOT.SetBatch()

def main():

     #inputFileWorkspace = TFile("Rootfiles/workspace_RPVStopStopToJets_UDD312_M-300_Resolved_delta_BiasTest_v08p0.root") 
     inputFileWorkspace = TFile("Rootfiles/workspace_RPVStopStopToJets_UDD312_M-400_Resolved_delta_v09p1.root") 
     #inputFileWorkspace = TFile("test.root") 

     workspace = inputFileWorkspace.Get("myWS")
     #workspace = inputFileWorkspace.Get("bkgWS")
     workspace.Print()

     mjj = workspace.var("mass")
     #mjj = workspace.var("massAve")
     mjj.Print()
     '''
     mjj1 = workspace.var("massAveBkg")
     mjj1.Print()
     mjj2 = workspace.var("massAveData")
     mjj2.Print()
     '''

     background = workspace.data("background")
#     backgroundUp = workspace.data("background__BkgUncUp")
#     backgroundDown = workspace.data("background__BkgUncDown")
#     backgroundStatUp = workspace.data("background__BkgStatUncUp")
#     backgroundStatDown = workspace.data("background__BkgStatUncDown")
     data = workspace.data("data_obs")
     signal = workspace.data("signal")
#     signalUp = workspace.data("signal__JESUp")
#     signalDown = workspace.data("signal__JESDown")
#     data.Print()
     #signal.Print()

     data_TH1_fineBinning = data.createHistogram("data_TH1_fineBinning",mjj)    
     data_TH1_fineBinning.SetLineColor(1)
#     signal_TH1_fineBinning = signal.createHistogram("signal_TH1_fineBinning",mjj)    
#     signalUp_TH1_fineBinning = signalUp.createHistogram("signalUp_TH1_fineBinning",mjj)    
#     signalUp_TH1_fineBinning.SetLineColor(kRed)
#     signalDown_TH1_fineBinning = signalDown.createHistogram("signalDown_TH1_fineBinning",mjj)    
#     signalDown_TH1_fineBinning.SetLineColor(kGreen)
#     backgroundUp_TH1_fineBinning = backgroundUp.createHistogram("backgroundUp_TH1_fineBinning",mjj)    
#     background_TH1_fineBinning = background.createHistogram("background_TH1_fineBinning",mjj)    
#     backgroundUp_TH1_fineBinning.SetLineColor(2)
#     backgroundDown_TH1_fineBinning = backgroundDown.createHistogram("backgroundDown_TH1_fineBinning",mjj)    
#     backgroundStatUp_TH1_fineBinning = backgroundStatUp.createHistogram("backgroundStatUp_TH1_fineBinning",mjj)    
#     backgroundStatDown_TH1_fineBinning = backgroundStatDown.createHistogram("backgroundStatDown_TH1_fineBinning",mjj)    
#     print "data_TH1_fineBinning integral = ", data_TH1_fineBinning.Integral()
#     #print "signal_TH1_fineBinning integral = ", signal_TH1_fineBinning.Integral()

     canvas = TCanvas()
     canvas.SetLogy()
     xframe = workspace.var("mass").frame()
     data.plotOn(xframe)
     workspace.pdf("background").plotOn(xframe)
     #workspace.pdf("landau_pdf").plotOn(xframe)
     workspace.pdf("signal").plotOn(xframe)
     xframe.Draw()

     #data_TH1_fineBinning.Rebin(10)    
     #data_TH1_fineBinning.Draw()    
     #signal_TH1_fineBinning.Rebin(10)
     #signal_TH1_fineBinning.Draw('')    
     #signalUp_TH1_fineBinning.Draw("hist same")    
     #signalDown_TH1_fineBinning.Draw("hist same")    
     #backgroundUp_TH1_fineBinning.Draw("")    
     #background_TH1_fineBinning.Rebin(10)
     #background_TH1_fineBinning.Draw()    
     #data_TH1_fineBinning.Draw("same")    
     #backgroundDown_TH1_fineBinning.Draw("")    
     #backgroundStatUp_TH1_fineBinning.Rebin(10)
     #backgroundStatUp_TH1_fineBinning.Draw("")    
     #backgroundStatDown_TH1_fineBinning.Rebin(10)
     #backgroundStatDown_TH1_fineBinning.Draw("same")    
     canvas.SaveAs("test.png")

if __name__ == '__main__':
    main()
