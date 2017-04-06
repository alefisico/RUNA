// system include files
#include <memory>
#include <vector>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/Framework/interface/GetterOfProducts.h"
#include "FWCore/Framework/interface/ProcessMatch.h"

// Jet Corrections
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

using namespace edm;
using namespace std;
using namespace pat;



typedef struct Jet_struc {
	TLorentzVector p4;
	TLorentzVector subjet0;
	TLorentzVector subjet1;
	double mass;
	double tau1;
	double tau2;
	double tau3;
	double btagCSV;
	double nhf;
	double nEMf;
	double chf;
	double cEMf;
	int numConst;
	double chm;
} JETtype;

class myJet {
	public:
		TLorentzVector p4;
		TLorentzVector subjet0;
		double subjet0BtagCSVv2;
		double subjet0BtagCMVAv2;
		TLorentzVector subjet1;
		double subjet1BtagCSVv2;
		double subjet1BtagCMVAv2;
		double mass;
		double trimmedMass;
		double prunedMass;
		double filteredMass;
		double softDropMass;
		double qgl;
		double tau1;
		double tau2;
		double tau3;
		double btagCSVv2;
		double btagCMVAv2;
		double btagDoubleB;
		double hadronFlavour;
		double nhf;
		double nEMf;
		double chf;
		double cEMf;
		int numConst;
		double chm;
};

inline float massAverage( float m1, float m2 ){ return ( (m1 + m2)/2 ); }

inline float massAsymmetry( float m1, float m2 ){ return abs( ( m1 - m2 )/( m1 + m2 ) ); }

inline float deltaValue( float p1, float p2 ){ return abs( p1 - p2 ); }

inline float calculateCosThetaStar( TLorentzVector jet1, TLorentzVector jet2 ){

	TLorentzVector tmpCM = jet1 + jet2;
	jet1.Boost( -tmpCM.BoostVector() );
	jet2.Boost( -tmpCM.BoostVector() );
	float valueCosThetaStar = TMath::Abs( ( jet1.Px() * tmpCM.Px() +  jet1.Py() * tmpCM.Py() + jet1.Pz() * tmpCM.Pz() ) / (jet1.E() * tmpCM.E() ) ) ;

	return valueCosThetaStar;
}

inline float cosThetaStar( TLorentzVector jet1, TLorentzVector jet2 ){

	TLorentzVector tmpCM = jet1 + jet2;
	TLorentzVector tmpJ1, tmpJ2;
	tmpJ1.SetPtEtaPhiE( jet1.Pt(), jet1.Eta(), jet1.Phi(), jet1.E() );
	tmpJ2.SetPtEtaPhiE( jet2.Pt(), jet2.Eta(), jet2.Phi(), jet2.E() );
	tmpJ1.Boost( -tmpCM.BoostVector() );
	tmpJ2.Boost( -tmpCM.BoostVector() );
	TVector3 tmpV1( tmpJ1.X(), tmpJ1.Y(), tmpJ1.Z() );
	TVector3 tmpV2( tmpJ2.X(), tmpJ2.Y(), tmpJ2.Z() );
	float valueCosThetaStar = TMath::Abs( tmpV1.CosTheta() ) ;

	return valueCosThetaStar;
}

inline bool jetID( double jetEta, double jetE, double jecFactor, double neutralHadronEnergy, double neutralEmEnergy, double chargedHadronEnergy, double muonEnergy, double chargedEmEnergy, double chargedMultiplicity, double neutralMultiplicity, string jetIDtype ){ 

	double jec = 1. / ( jecFactor );
	double NHF = neutralHadronEnergy * jec;
	double NEMF = neutralEmEnergy * jec;
	double CHF = chargedHadronEnergy * jec;
	double MUF = muonEnergy * jec;
	double CEMF = chargedEmEnergy * jec;
	int NumConst = chargedMultiplicity + neutralMultiplicity ; 
	double CHM = chargedMultiplicity * jec;

	bool id = 0;
	if ( jetIDtype == "looseJetID" ) id = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(jetEta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jetEta)>2.4) && abs(jetEta)<=2.7;
	else if ( jetIDtype == "tightJetID" ) id = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(jetEta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jetEta)>2.4) && abs(jetEta)<=2.7;
	else if ( jetIDtype == "tightLepVetoJetID" ) id = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(jetEta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(jetEta)>2.4) && abs(jetEta)<=2.7;
	else LogError("jetID") << "jetID function only takes the following jetID names: looseJetID, tightJetID, tightLepVetoJetID.";

	return id;
}

inline bool checkTriggerBits( vector<string> triggerNames, Handle<vector<float>> triggerBits, Handle<vector<int>> triggerPrescale, TString HLTtrigger, bool baselineTrigger  ){

	float triggerFired = 0;
	for (size_t t = 0; t < triggerNames.size(); t++) {
		//LogWarning("triggerbit") << triggerNames[t] << " " <<  (*triggerBits)[t] << " " << (*triggerPrescale)[t];
		if ( TString( triggerNames[t] ).Contains( HLTtrigger ) ) {
			if ( ((*triggerPrescale)[t] == 1) || baselineTrigger ) triggerFired = (*triggerBits)[t];
		}
	}
	if ( HLTtrigger.Contains( "NOTRIGGER" ) ) triggerFired = 1;
	return triggerFired;
}	

inline bool checkORListOfTriggerBits( vector<string> triggerNames, Handle<vector<float>> triggerBits, Handle<vector<int>> triggerPrescale, vector<string> triggerPass, bool baselineTrigger ){

	bool ORTriggers = 0;
	if ( triggerNames.size() == triggerBits->size() ) {

		vector<bool> triggersFired;
		for (size_t t = 0; t < triggerPass.size(); t++) {
			bool triggerFired = checkTriggerBits( triggerNames, triggerBits, triggerPrescale, triggerPass[t], baselineTrigger );
			triggersFired.push_back( triggerFired );
			//LogWarning("test") << triggerPass[t] << " " << triggerFired;
			//if ( triggerFired ) LogWarning("test") << triggerPass[t] << " " << triggerFired;
		}
		
		ORTriggers = any_of(triggersFired.begin(), triggersFired.end(), [](bool v) { return v; }); 
	} else LogError("TriggerBits") << "triggerNameList has different size than triggerBit";

	return ORTriggers;
}

inline bool checkTriggerBitsMiniAOD( TriggerNames triggerNames, Handle<TriggerResults> triggerBits, Handle<PackedTriggerPrescales> triggerPrescales, TString HLTtrigger, bool baselineTrigger  ){

  	bool triggerFired = 0;
	for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
		//LogWarning("triggerbit") << triggerNames.triggerName(i) << " " <<  triggerBits->accept(i) << " " << triggerPrescales->getPrescaleForIndex(i);
		if (TString(triggerNames.triggerName(i)).Contains(HLTtrigger) && (triggerBits->accept(i))) {
			if ( (triggerPrescales->getPrescaleForIndex(i) == 1) || baselineTrigger ) triggerFired=1;
		}
	}

	return triggerFired;
}	

inline bool checkORListOfTriggerBitsMiniAOD( TriggerNames triggerNames, Handle<TriggerResults> triggerBits, Handle<PackedTriggerPrescales> triggerPrescales, vector<string>  triggerPass, bool baselineTrigger  ){

	vector<bool> triggersFired;
	for (size_t t = 0; t < triggerPass.size(); t++) {
		bool triggerFired = checkTriggerBitsMiniAOD( triggerNames, triggerBits, triggerPrescales, triggerPass[t], baselineTrigger );
		triggersFired.push_back( triggerFired );
		//LogWarning("test") << triggerPass[t] << " " << triggerFired;
		//if ( triggerFired ) LogWarning("test") << triggerPass[t] << " " << triggerFired;
	}
	
	bool ORTriggers = any_of(triggersFired.begin(), triggersFired.end(), [](bool v) { return v; }); 
	
	return ORTriggers;
}

inline double corrections( TLorentzVector rawJet, double jetArea, double Rho, int NPV, FactorizedJetCorrector* jetCorrector ){

	jetCorrector->setJetPt ( rawJet.Pt() );
	jetCorrector->setJetEta( rawJet.Eta() );
	jetCorrector->setJetPhi( rawJet.Phi() );
	jetCorrector->setJetE  ( rawJet.E() );
	jetCorrector->setJetA  ( jetArea );
	jetCorrector->setRho   ( Rho );
	jetCorrector->setNPV   ( NPV );
	double jecFactor= jetCorrector->getCorrection();

	return jecFactor;
}

inline double uncertainty( TLorentzVector rawJet, JetCorrectionUncertainty* jetUnc, bool getUnc ){

	jetUnc->setJetPt ( rawJet.Pt()  );
	jetUnc->setJetEta( rawJet.Eta() );
	double jetUncFactor = jetUnc->getUncertainty( getUnc );
	return jetUncFactor;
}

inline double getJER( double jetEta, int JERType ){

	double scaleNom = 1.0;
	double scaleUnc = 1.0;
	double eta = fabs(jetEta);
	if(eta>=0.0 && eta<0.8) { scaleNom = 1.061; scaleUnc = 0.023; }
	if(eta>=0.8 && eta<1.3) { scaleNom = 1.088; scaleUnc = 0.029; }
	if(eta>=1.3 && eta<1.9) { scaleNom = 1.106; scaleUnc = 0.030; }
	if(eta>=1.9 && eta<2.5) { scaleNom = 1.126; scaleUnc = 0.094; }
	if(eta>=2.5 && eta<3.0) { scaleNom = 1.343; scaleUnc = 0.123; }
	if(eta>=3.0 && eta<3.2) { scaleNom = 1.303; scaleUnc = 0.111; }
	if(eta>=3.2 && eta<5.0) { scaleNom = 1.320; scaleUnc = 0.286; }

	if ( JERType == 0 ) return scaleNom;
	else if ( JERType == 1 ) return (scaleNom + scaleUnc);
	else if ( JERType == -1 ) return (scaleNom - scaleUnc);
	else return 1.;
}

// all this section is based on https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L215-L281
inline void getWeights( Handle<LHEEventProduct> lheEvtInfo, int lha_pdf_id_, vector<float> & scaleWeights, vector<float> & pdfWeights, vector<float> & alphaWeights ){

	/* Renormalization/Factorization scale weights
	 * These are the first 9 weights for all generators
	 * mu_R and mu_F are varied independently (by a factor of 1, 2, 0.5) - check LHE header
	 * [0] is the default weight (no variation) - it has worse precision even
	 * --> I skip saving it (one can just use 1)
	 * --> Also do not save unphysical combinations as recommended
	 * (mu_R = 0.5, mu_F = 2 and mu_R = 2, mu_F = 0.5)
	 * Save only: 1,2,3,4,6,8
	 */
	double lheOrigWeight = lheEvtInfo->originalXWGTUP();
	vector <gen::WeightsInfo> weightsTemp = lheEvtInfo->weights();
	if ( weightsTemp.size()>=9) for (size_t i=0; i<9; ++i) if (i!=0&&i!=5&&i!=7) scaleWeights.push_back( weightsTemp[i].wgt/lheOrigWeight);

	/* PDF weights
	 * Usually a set of 100 weights (excluding default)
	 * Only default PDF variation is saved, but if needed
	 * likely others are available depending on the generator
	 * index of first weight varies, beware!
	 * Sometimes first weight is default=1 weight (has to skip!)
	 * Additional info: MC2Hessian conversion will soon be provided
	 */
	size_t first = 9;
	// Madgraph nf5 - have to skip 1 weight which is default
	if (lha_pdf_id_ == 263000) first = 10;
	// Madgraph nf4 uses NNPDF30_lo_as_0130_nf_4 (ID=263400)
	// Which is the second 101 pdf set, again has to skip first weight
	if (lha_pdf_id_ == 263400) first = 111;
	if (weightsTemp.size()>=first+100) for (size_t i=first; i<first+100; ++i) pdfWeights.push_back(weightsTemp[i].wgt/lheOrigWeight);

	/* Alpha_s weights (only for NLO!)
	 * A set of two weights for 
	 * alpha_s = 0.118 - 0.002 and
	 * alpha_s = 0.118 + 0.002
	 * is given --> scale result uncertainty by 0.75
	 */
	if ( weightsTemp.size()>=111 &&
	    ( (lha_pdf_id_ == 260000) || // Powheg 5nf
	     (lha_pdf_id_ == 260400) || // Powheg 4nf 
	     (lha_pdf_id_ == 292000) || // aMC@NLO 5nf
	     (lha_pdf_id_ == 292200)    // aMC@NLO 5nf
	     ) ) {
		alphaWeights.push_back(weightsTemp[109].wgt/lheOrigWeight);
		alphaWeights.push_back(weightsTemp[110].wgt/lheOrigWeight);
	}
}

//// Btagging scale factors
inline double btagSF( string csvFile, double jetPt, double jetEta, double flavour, string sysType, string measurementType ){

	double jetSF = 1;
        BTagCalibration calib("csv", csvFile );
	BTagCalibrationReader reader(BTagEntry::OP_LOOSE,	// operating point
			"central",				// central sys type
			{"up", "down"});			// other sys types

	if ( TMath::Abs( flavour ) == 5 ) { 
		reader.load( calib, BTagEntry::FLAV_B, measurementType);
		jetSF = reader.eval_auto_bounds( sysType, BTagEntry::FLAV_B, jetEta, jetPt );  

	} else if ( TMath::Abs( flavour ) == 4 ) { 
		reader.load( calib, BTagEntry::FLAV_C, measurementType);
		jetSF = reader.eval_auto_bounds( sysType, BTagEntry::FLAV_C, jetEta, jetPt ); 

	} else if ( TMath::Abs( flavour ) == 0 ) {
		reader.load( calib, BTagEntry::FLAV_UDSG, "incl");
		jetSF = reader.eval_auto_bounds( sysType, BTagEntry::FLAV_UDSG, jetEta, jetPt );  
	}

	return jetSF;
}

inline vector< myJet > pairing( vector< myJet > JETS, bool withDeltaR ) {

	vector< myJet > reorderedJETS;
	if ( JETS.size() > 3 ) {
		vector<double> tmpDijetR;
		bool invertOrder = false;
		double pair1234 = 9999, pair1324 = 9999, pair1423 = 9999;
		TLorentzVector j12, j34, j13, j24, j14, j23;

		if( withDeltaR ) {
			double dR12 = JETS[0].p4.DeltaR( JETS[1].p4 );
			double dR34 = JETS[2].p4.DeltaR( JETS[3].p4 );
			if ( dR12 < dR34 ) invertOrder = true;
			pair1234 = abs( dR12 - 0.8 )  + abs( dR34 - 0.8 );
		} else {
			j12 = JETS[0].p4 + JETS[1].p4;
			j34 = JETS[2].p4 + JETS[3].p4;
			if ( j12.Pt() < j34.Pt() ) invertOrder = true;
			pair1234 = TMath::Abs( j12.M() - j34.M() )/ (j12.M() + j34.M());
		}
		tmpDijetR.push_back( pair1234 );

		if( withDeltaR ) {
			double dR13 = JETS[0].p4.DeltaR( JETS[2].p4 );
			double dR24 = JETS[1].p4.DeltaR( JETS[3].p4 );
			if ( dR13 < dR24 ) invertOrder = true;
			pair1324 = abs( dR13 - 0.8 )  + abs( dR24 - 0.8 );
		} else {
			j13 = JETS[0].p4 + JETS[2].p4;
			j24 = JETS[1].p4 + JETS[3].p4;
			if ( j13.Pt() < j24.Pt() ) invertOrder = true;
			pair1324 = TMath::Abs( j13.M() - j24.M() )/ (j13.M() + j24.M());
		}
		tmpDijetR.push_back( pair1324 );


		if( withDeltaR ) {
			double dR14 = JETS[0].p4.DeltaR( JETS[3].p4 );
			double dR23 = JETS[1].p4.DeltaR( JETS[2].p4 );
			if ( dR14 < dR23 ) invertOrder = true;
			pair1423 = abs( dR14 - 0.8 )  + abs( dR23 - 0.8 );
		} else {
			j14 = JETS[0].p4 + JETS[3].p4;
			j23 = JETS[1].p4 + JETS[2].p4;
			if ( j14.Pt() < j23.Pt() ) invertOrder = true;
			pair1423 = TMath::Abs( j14.M() - j23.M() )/ (j14.M() + j23.M());
		}
		tmpDijetR.push_back( pair1423 );

		int tmpMin = min_element(tmpDijetR.begin(), tmpDijetR.end()) - tmpDijetR.begin();
		if ( invertOrder ) tmpMin = tmpMin + 10;
		/*
		 * Notation:
		 * Pair 12 34: 0 and 10 
		 * Pair 13 24: 1 and 11
		 * Pair 14 23: 2 and 12
		 */
		if ( tmpMin == 0 ) {
			reorderedJETS.push_back( JETS[0] );
			reorderedJETS.push_back( JETS[1] );
			reorderedJETS.push_back( JETS[2] );
			reorderedJETS.push_back( JETS[3] );
		} else if ( tmpMin == 10 ) {
			reorderedJETS.push_back( JETS[2] );
			reorderedJETS.push_back( JETS[3] );
			reorderedJETS.push_back( JETS[0] );
			reorderedJETS.push_back( JETS[1] );
		} else if ( tmpMin == 1 ) {
			reorderedJETS.push_back( JETS[0] );
			reorderedJETS.push_back( JETS[2] );
			reorderedJETS.push_back( JETS[1] );
			reorderedJETS.push_back( JETS[3] );
		} else if ( tmpMin == 11 ) {
			reorderedJETS.push_back( JETS[1] );
			reorderedJETS.push_back( JETS[3] );
			reorderedJETS.push_back( JETS[0] );
			reorderedJETS.push_back( JETS[2] );
		} else if ( tmpMin == 2 ) {
			reorderedJETS.push_back( JETS[0] );
			reorderedJETS.push_back( JETS[3] );
			reorderedJETS.push_back( JETS[1] );
			reorderedJETS.push_back( JETS[2] );
		} else if ( tmpMin == 12 ) {
			reorderedJETS.push_back( JETS[1] );
			reorderedJETS.push_back( JETS[2] );
			reorderedJETS.push_back( JETS[0] );
			reorderedJETS.push_back( JETS[3] );
		}

		/// storing the rest of the jets
		for (unsigned int ijet = 4; ijet < JETS.size(); ijet++) {
			reorderedJETS.push_back( JETS[ijet] );
		}

	} else LogWarning("JETS") << "Number of input jets is lower than 4."; 
	return reorderedJETS;

}

inline vector< myJet > pairingMinChi( vector< myJet > JETS, vector<float> * vectorChi2 ) {

	vector< myJet > reorderedJETS;
	vector< int > indices;
	double minChi2 = 99999999999999999;
	int bestInd[4] = { -1, -1, -1, -1};
	if( JETS.size() > 3 ) {
		for (unsigned int a = 0; a < JETS.size(); a++) {
			for (unsigned int b = 1; b < JETS.size(); b++) {
				for (unsigned int c = 2; c < JETS.size(); c++) {
					for (unsigned int d = 3; d < JETS.size(); d++) {
						if( (a!=b) && (a!=c) && (a!=d) && (b!=c) && (b!=d) && (c!=d) && (a<b) && (c<d) ) {
							//cout << a << " " << b << " " << c << " " << d << endl; 
							double dijet1Mass = ( JETS[a].p4 + JETS[b].p4 ).M();	
							double dijet2Mass = ( JETS[c].p4 + JETS[d].p4 ).M();	
							double aveMass = (dijet1Mass + dijet2Mass )/2 ;
							double sigma = 0.0724281 - (4.14372e-05 * aveMass);
							double tmpChi2 = TMath::Power( (dijet1Mass - dijet2Mass )/aveMass, 2 ) / TMath::Power( sigma, 2 );
							//LogWarning("calc chi2") << a << b << c << d << " " << dijet1Mass << " " << dijet2Mass << " " << tmpChi2 << " " << sigma  << " " << aveMass;
							vectorChi2->push_back( tmpChi2 );
							if( tmpChi2 < minChi2 ){
								minChi2 = tmpChi2;
								bestInd[0] = a;
								bestInd[1] = b;
								bestInd[2] = c;
								bestInd[3] = d;
								//LogWarning("min chi2") << a << b << c << d << " " << dijet1Mass << " " << dijet2Mass << " " << minChi2;
							}
						}
					}
				}
			}
		}
	}
	//LogWarning("index") << bestInd[0] << " " << bestInd[1] << " " << bestInd[2] << " " << bestInd[3];
	reorderedJETS.push_back( JETS[ bestInd[0] ] );
	reorderedJETS.push_back( JETS[ bestInd[1] ] );
	reorderedJETS.push_back( JETS[ bestInd[2] ] );
	reorderedJETS.push_back( JETS[ bestInd[3] ] );

	return reorderedJETS;
}
