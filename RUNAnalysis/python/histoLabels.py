#!/usr/bin/env python

#from ROOT import TFile, TH1F, THStack, TCanvas, TMath, gROOT, gPad
from ROOT import *
import ROOT as rt
from cuts import selection

varLabels = {}
varLabels[ 'prunedMassAsym' ] = 'Asym'
varLabels[ 'deltaEtaDijet' ] = '| #eta_{j1} - #eta_{j2} |'
varLabels[ 'jet1Tau21' ] = '#tau_{21}^{j1}'
varLabels[ 'jet2Tau21' ] = '#tau_{21}^{j2}'
varLabels[ 'jet1Tau31' ] = '#tau_{31}^{j1}'
varLabels[ 'jet2Tau31' ] = '#tau_{31}^{j2}'

def setSelection( listSel, xMin=0.65, yMax=0.65, align='right' ):

	for i in range( len( listSel ) ):
		textBox=rt.TLatex()
		textBox.SetNDC()
		textBox.SetTextSize(0.04) 
		if 'right' in align: textBox.SetTextAlign(31)
		textBox.SetTextFont(62) ### 62 is bold, 42 is normal
		textBox.DrawLatex(xMin, yMax, listSel[i])
		yMax = yMax -0.05
	
def finalLabels( signal, X=0.92, Y=0.50, align='right' ):
	"""docstring for finalLabels"""

	tmpListSel = selection[signal]

	listSel = [ 'Preselection' ]
	for sel in tmpListSel:
		for lab in varLabels: 
			if sel[0] in lab: listSel.append( varLabels[lab]+' < '+str(sel[1]) )

	setSelection( listSel, X, Y, align) 


def labels( name, camp, X=0.92, Y=0.50, align='right', listSel=[] ):

	if 'cutDijet' in name: listSel = [ 'jet p_{T} > 150 GeV', 'jet |#eta| < 2.4', 'numJets > 1' ] 
	elif 'presel' in name: listSel = [ 'Preselection' ] 
	elif 'cutEffTrigger' in name: listSel = [ 'Preselection' ] 
	elif 'cutHT' in name: listSel = [ 'numJets > 1', 'p_{T}^{j1} = 500 GeV',  'p_{T}^{j2} = 450 GeV' ] 
	elif 'cutTau21' in name: listSel = [ 'numJets > 1', 'p_{T}^{j1} = 500 GeV',  'p_{T}^{j2} = 450 GeV', '#tau_{21} < 0.6']
	elif 'cutMassAsym' in name: listSel = [ 'numJets > 1', 'p_{T}^{j1} = 500 GeV',  'p_{T}^{j2} = 450 GeV', '#tau_{21} < 0.6', '|cos(#theta*)| < 0.2', 'A < 0.1' ]
	elif 'cutDeltaEtaDijet' in name: listSel = [ 'Preselection', '#tau_{21}^{j1,j2} < 0.45', 'Mass_{asym} < 0.1', '|#eta_{j1}-#eta_{j2}|<1.5' ]
	elif 'PFHT800' in name: listSel = [ camp, 'PFHT800', 'HT > 800 GeV', 'numJets > 1', 'A < 0.1', '|cos(#theta*)| < 0.3', 'subjet pt ratio > 0.3']
	elif 'Brock' in name: listSel = [ camp, 'HT > 1600 TeV', 'HT > 800 GeV', 'numJets > 1', 'A < 0.1', '|cos(#theta*)| < 0.3', 'subjet pt ratio > 0.3']
	elif 'cutTau31' in name: listSel = [ camp, 'HT > 800 GeV', 'numJets > 1', 'A < 0.1', '|cos(#theta*)| < 0.3', '#tau_{31} < 0.4']
	#elif 'triggerDATA' in name: listSel =  [  'Eff. (AK8 H_{T} > 700 GeV, trimmed AK8 jet mass > 50 GeV)']
	#elif 'triggerSignal' in name: listSel =  [ 'RPV Stop 100 GeV', triggerUsed.replace('650', '650TrimMass50')]
	#elif 'triggerSignal' in name: listSel =  [ 'pp #rightarrow #tilde{t}(jj) #tilde{t}(jj), M(#tilde{t}) = 100 GeV', 'Eff. (AK8 H_{T} > 650 GeV, trimmed AK8 jet mass > 50 GeV)']
	#elif 'trigger' in name: listSel = [ 'pp #rightarrow #tilde{t}(jj) #tilde{t}(jj), M(#tilde{t}) = 100 GeV', 'AK4 H_{T} > 800 GeV AND', 'AK8 H_{T} > 650 GeV, trimmed AK8 jet mass > 50 GeV' ] 
	#elif 'trigger' in name: listSel = [ 'pp #rightarrow #tilde{t}(jj) #tilde{t}(jj), M(#tilde{t}) = 100 GeV', 'AK4 H_{T} > 800 GeV' ] 
	#elif 'trigger' in name: listSel = [ 'pp #rightarrow #tilde{t}(jj) #tilde{t}(jj), M(#tilde{t}) = 100 GeV', 'AK8 H_{T} > 700 GeV, trimmed AK8 jet mass > 50 GeV' ] 
	#elif 'rigger' in name: listSel = [ 'pp #rightarrow #tilde{t}(jj) #tilde{t}(jj), M(#tilde{t}) = 100 GeV', 'NO Trigger' ] 
	elif 'cut4Jet' in name: listSel = [ camp, 'jet p_{T} > 50 GeV', 'jet |#eta| < 2.4', 'numJets > 3' ] 
	elif 'cut4JetPt' in name: listSel = [ camp, 'jet p_{T} > 50 GeV', 'jet |#eta| < 2.4', 'numJets > 3', 'HT > 800 GeV' ] 
	elif 'cutBestPair' in name: listSel = [ 'Preselection', 'dR_{dijet}' ] 
	elif 'cutDelta' in name: listSel = [ 'Preselection', 'A < 0.1', '| #Delta #eta_{dijet}| < 1', '#Delta > 200 GeV' ]
	elif 'cutDEta' in name: listSel = [ camp,  'HT > 800 GeV', 'numJets > 3','#delta R < 1.5', '|#cos(#theta*)| <  0.7' '#Delta #eta (j^{1},j^{2}) > 1']
	elif '' in name: listSel = [ '' ] 
	else: listSel = [ 'NO LABEL' ]

	setSelection( listSel, X, Y, align) 

def labelAxis(name, histo, Grom ):

	if ('massAve' in name) or ('MassAve' in name): 
		if 'trimmed' in Grom: histo.GetXaxis().SetTitle( 'Average Trimmed Mass [GeV]' )
		elif 'pruned' in Grom: histo.GetXaxis().SetTitle( 'Average Pruned Mass [GeV]' )
		elif 'softDrop' in Grom: histo.GetXaxis().SetTitle( 'Average Soft Drop Mass [GeV]' )
		elif 'filtered' in Grom: histo.GetXaxis().SetTitle( 'Average Filtered Mass [GeV]' )
		else: histo.GetXaxis().SetTitle( 'Average Mass [GeV]' )
	elif 'jet1TrimmedMass' in name: histo.GetXaxis().SetTitle( 'Leading Trimmed Jet Mass [GeV]' )
	elif 'jet1PrunedMass' in name: histo.GetXaxis().SetTitle( 'Leading Pruned Jet Mass [GeV]' )
	elif 'jet1FilteredMass' in name: histo.GetXaxis().SetTitle( 'Leading Filtered Jet Mass [GeV]' )
	elif 'jet1SoftDropMass' in name: histo.GetXaxis().SetTitle( 'Leading Soft Drop Jet Mass [GeV]' )
	elif 'jet2TrimmedMass' in name: histo.GetXaxis().SetTitle( '2nd Leading Trimmed Jet Mass [GeV]' )
	elif 'jet2PrunedMass' in name: histo.GetXaxis().SetTitle( '2nd Leading Pruned Jet Mass [GeV]' )
	elif 'jet2FilteredMass' in name: histo.GetXaxis().SetTitle( '2nd Leading Filtered Jet Mass [GeV]' )
	elif 'jet2SoftDropMass' in name: histo.GetXaxis().SetTitle( '2nd Leading Soft Drop Jet Mass [GeV]' )
	elif 'jet1Rho' in name: histo.GetXaxis().SetTitle( 'Leading Jet #rho [GeV]' )
	elif 'jet2Rho' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #rho [GeV]' )
	elif 'prunedMassAsym' in name: histo.GetXaxis().SetTitle( 'M_{Asym}' )
	elif 'massAsym' in name: histo.GetXaxis().SetTitle( 'M_{Asym}' )
	elif ('jet1CosThetaStar' in name) or ('cosThetaStar1' in name ): histo.GetXaxis().SetTitle( 'cos(#theta *)_{1}' )
	elif ('jet2CosThetaStar' in name) or ('cosThetaStar2' in name ): histo.GetXaxis().SetTitle( 'cos(#theta *)_{2}' )
	elif 'jetEta' in name: histo.GetXaxis().SetTitle( 'Jet #eta' )
	elif 'Tau1_' in name: histo.GetXaxis().SetTitle( '#tau_{1}' )
	elif 'Tau2_' in name: histo.GetXaxis().SetTitle( '#tau_{2}' )
	elif 'Tau3_' in name: histo.GetXaxis().SetTitle( '#tau_{3}' )
	elif 'jet1RhoDDT' in name: histo.GetXaxis().SetTitle( 'Leading Jet #rho\'' )
	elif 'jet2RhoDDT' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #rho\'' )
	elif 'jet1Tau21DDT' in name: histo.GetXaxis().SetTitle( 'Leading Jet #tau_{21}^{DDT} ' )
	elif 'jet2Tau21DDT' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #tau_{21}^{DDT} ' )
	elif 'jet1Tau21_' in name: histo.GetXaxis().SetTitle( 'Leading Jet #tau_{21} ' )
	elif 'jet2Tau21_' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #tau_{21} ' )
	elif 'jet1Tau31_' in name: histo.GetXaxis().SetTitle( 'Leading Jet #tau_{31} ' )
	elif 'jet2Tau31_' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #tau_{31} ' )
	elif 'jet1Tau21' in name: histo.GetXaxis().SetTitle( 'Leading Jet #tau_{21} ' )
	elif 'jet2Tau21' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #tau_{21} ' )
	elif 'jet1Tau31' in name: histo.GetXaxis().SetTitle( 'Leading Jet #tau_{31} ' )
	elif 'jet2Tau31' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #tau_{31} ' )
	elif 'jet1Tau32' in name: histo.GetXaxis().SetTitle( 'Leading Jet #tau_{32} ' )
	elif 'jet2Tau32' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #tau_{32} ' )
	elif 'Tau31_' in name: histo.GetXaxis().SetTitle( '#tau_{31} ' )
	elif 'Tau32_' in name: histo.GetXaxis().SetTitle( '#tau_{32} ' )
	elif 'jet1SubjetPtRatio' in name: histo.GetXaxis().SetTitle( 'Leading jet Subjet Pt_{2}/Pt_{1}' )
	elif 'jet2SubjetPtRatio' in name: histo.GetXaxis().SetTitle( '2nd Leading jet Subjet Pt_{2}/Pt_{1}' )
	elif 'PtRatio_' in name: histo.GetXaxis().SetTitle( 'Subjet Pt_{2}/Pt_{1}' )
	elif 'Mass21' in name: histo.GetXaxis().SetTitle( 'Subjet m_{2}/m_{1}' )
	elif '112Mass' in name: histo.GetXaxis().SetTitle( 'Subjet m_{1}/m_{12}' )
	elif '212Mass' in name: histo.GetXaxis().SetTitle( 'Subjet m_{2}/m_{12}' )
	elif 'PolAngle13412_' in name: histo.GetXaxis().SetTitle( 'cos #psi_{1(34)}^{[12]}' )
	elif 'PolAngle31234_' in name: histo.GetXaxis().SetTitle( 'cos #psi_{3(12)}^{[34]}' )
	elif 'jetMass' in name: histo.GetXaxis().SetTitle( 'Jet Mass [GeV]' )
	elif 'jetPt_' in name: histo.GetXaxis().SetTitle( 'Jet p_{T} [GeV]' )
	elif 'jetsBtag' in name: histo.GetXaxis().SetTitle( 'Jets CSVv2 btagging' )
	elif 'jet1Pt' in name: histo.GetXaxis().SetTitle( 'Leading Jet p_{T} [GeV]' )
	elif 'jet1Mass' in name: histo.GetXaxis().SetTitle( 'Leading Jet Mass [GeV]' )
	elif 'jet1Eta' in name: histo.GetXaxis().SetTitle( 'Leading Jet #eta' )
	elif 'jet2Pt' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet p_{T} [GeV]' )
	elif 'jet3Pt' in name: histo.GetXaxis().SetTitle( '3rd Leading Jet p_{T} [GeV]' )
	elif 'jet4Pt' in name: histo.GetXaxis().SetTitle( '4th Leading Jet p_{T} [GeV]' )
	elif 'jet2Mass' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet Mass [GeV]' )
	elif 'jet2Eta' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet #eta' )
	elif ('jetNum' in name) or ('numJets' in name): histo.GetXaxis().SetTitle( 'Number of Jets' )
	elif 'deltaEta' in name: histo.GetXaxis().SetTitle( '#Delta #eta (j^{1}, j^{2})' )
	elif 'deltaR' in name: histo.GetXaxis().SetTitle( '#Delta R dijet' )
	elif 'NPV' in name: histo.GetXaxis().SetTitle( 'Number of Primary Vertex' )
	elif 'massRes' in name: histo.GetXaxis().SetTitle( 'Fractional mass difference' )
	elif 'neutralHadronEnergy' in name: histo.GetXaxis().SetTitle( 'Neutral hadron energy' )
	elif 'neutralEmEnergy' in name: histo.GetXaxis().SetTitle( 'Neutral EM energy' )
	elif 'chargedHadronEnergy' in name: histo.GetXaxis().SetTitle( 'Charged hadron energy' )
	elif 'chargedEmEnergy' in name: histo.GetXaxis().SetTitle( 'Charged EM energy' )
	elif 'chargedMultiplicity' in name: histo.GetXaxis().SetTitle( 'Charged multiplicity' )
	elif 'numConst' in name: histo.GetXaxis().SetTitle( 'Number of Constituents' )
	elif 'jet1NeutralHadronEnergy' in name: histo.GetXaxis().SetTitle( 'Leading Jet Neutral hadron energy' )
	elif 'jet1NeutralEmEnergy' in name: histo.GetXaxis().SetTitle( 'Leading Jet Neutral EM energy' )
	elif 'jet1ChargedHadronEnergy' in name: histo.GetXaxis().SetTitle( 'Leading Jet Charged hadron energy' )
	elif 'jet1ChargedEmEnergy' in name: histo.GetXaxis().SetTitle( 'Leading Jet Charged EM energy' )
	elif 'jet1ChargedMultiplicity' in name: histo.GetXaxis().SetTitle( 'Leading Jet Charged multiplicity' )
	elif 'jet1NumConst' in name: histo.GetXaxis().SetTitle( 'Leading Jet Number of Constituents' )
	elif 'jet2NeutralHadronEnergy' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet Neutral hadron energy' )
	elif 'jet2NeutralEmEnergy' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet Neutral EM energy' )
	elif 'jet2ChargedHadronEnergy' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet Charged hadron energy' )
	elif 'jet2ChargedEmEnergy' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet Charged EM energy' )
	elif 'jet2ChargedMultiplicity' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet Charged multiplicity' )
	elif 'jet2NumConst' in name: histo.GetXaxis().SetTitle( '2nd Leading Jet Number of Constituents' )
	elif 'METHT' in name: histo.GetXaxis().SetTitle( 'MET/HT' )
	elif 'HT' in name: histo.GetXaxis().SetTitle( 'HT [GeV]' )
	elif 'MET' in name: histo.GetXaxis().SetTitle( 'MET [GeV]' )
	elif 'minDeltaR' in name: histo.GetXaxis().SetTitle( 'min #Delta R' )
	elif 'minChi2' in name: histo.GetXaxis().SetTitle( 'min #Delta #xi^{2}' )
	elif 'minMassAsym' in name: histo.GetXaxis().SetTitle( 'min Mass' )
	#elif '' in name: histo.GetXaxis().SetTitle( '' )
	else: histo.GetXaxis().SetTitle( 'NO LABEL' )
