// -*- C++ -*-
//
// Package:    Ntuples/Ntuples
// Class:      RUNONBoostedAnalysis
// 
/**\class RUNONBoostedAnalysis RUNONBoostedAnalysis.cc Ntuples/Ntuples/plugins/RUNONBoostedAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  alejandro gomez
//         Created:  Tue, 14 Oct 2014 23:13:13 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <TLorentzVector.h>
#include <TH2.h>
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

using namespace edm;
using namespace std;

//
// constants, enums and typedefs
//
typedef struct Jet_struc {
	TLorentzVector p4;
	TLorentzVector subjet0;
	TLorentzVector subjet1;
	double mass;
	double tau1;
	double tau2;
	double tau3;
	bool btagCSV;
} JETtype;
//
// class declaration
//
class RUNONBoostedAnalysis : public EDAnalyzer {
   public:
      explicit RUNONBoostedAnalysis(const ParameterSet&);
      ~RUNONBoostedAnalysis();

   private:
      virtual void beginJob() override;
      virtual void analyze(const Event&, const EventSetup&) override;
      virtual void endJob() override;
      virtual void clearVariables();

      //virtual void beginRun(Run const&, EventSetup const&) override;
      //virtual void endRun(Run const&, EventSetup const&) override;
      //virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
      //virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

      // ----------member data ---------------------------
      //GetterOfProducts< vector<float> > getterOfProducts_;
      //EDGetTokenT<TriggerResults> triggerBits_;
      Service<TFileService> fs_;
      TTree *RUNAtree;
      map< string, TH1D* > histos1D_;
      map< string, TH2D* > histos2D_;
      vector< string > cutLabels;
      map< string, double > cutmap;

      bool bjSample;
      double scale;
      double cutTrimmedMassvalue;
      double cutHTvalue;
      double cutAsymvalue;
      double cutCosThetavalue;
      double cutSubjetPtRatiovalue;
      double cutTau31value;
      double cutTau21value;

      vector<float> *jetsPt = new std::vector<float>();
      vector<float> *jetsEta = new std::vector<float>();
      vector<float> *jetsPhi = new std::vector<float>();
      vector<float> *jetsE = new std::vector<float>();
      vector<float> *jet1SubjetsPt = new std::vector<float>();
      vector<float> *jet1SubjetsEta = new std::vector<float>();
      vector<float> *jet1SubjetsPhi = new std::vector<float>();
      vector<float> *jet1SubjetsE = new std::vector<float>();
      vector<float> *jet2SubjetsPt = new std::vector<float>();
      vector<float> *jet2SubjetsEta = new std::vector<float>();
      vector<float> *jet2SubjetsPhi = new std::vector<float>();
      vector<float> *jet2SubjetsE = new std::vector<float>();

      EDGetTokenT<vector<float>> jetPt_;
      EDGetTokenT<vector<float>> jetEta_;
      EDGetTokenT<vector<float>> jetPhi_;
      EDGetTokenT<vector<float>> jetE_;
      EDGetTokenT<vector<float>> jetMass_;
      EDGetTokenT<vector<float>> jetTrimmedMass_;
      EDGetTokenT<vector<float>> jetTau1_;
      EDGetTokenT<vector<float>> jetTau2_;
      EDGetTokenT<vector<float>> jetTau3_;
      EDGetTokenT<vector<float>> jetNSubjets_;
      EDGetTokenT<vector<float>> jetSubjetIndex0_;
      EDGetTokenT<vector<float>> jetSubjetIndex1_;
      EDGetTokenT<vector<float>> jetSubjetIndex2_;
      EDGetTokenT<vector<float>> jetSubjetIndex3_;
      EDGetTokenT<vector<vector<int>>> jetKeys_;
      EDGetTokenT<vector<float>> jetCSV_;
      EDGetTokenT<vector<float>> jetCSVV1_;
      EDGetTokenT<int> NPV_;

      //Jet ID
      EDGetTokenT<vector<float>> jecFactor_;
      EDGetTokenT<vector<float>> neutralHadronEnergy_;
      EDGetTokenT<vector<float>> neutralEmEnergy_;
      EDGetTokenT<vector<float>> chargeEmEnergy_;
      EDGetTokenT<vector<float>> muonEnergy_; 

      // Subjets
      EDGetTokenT<vector<float>> subjetPt_;
      EDGetTokenT<vector<float>> subjetEta_;
      EDGetTokenT<vector<float>> subjetPhi_;
      EDGetTokenT<vector<float>> subjetE_;
      EDGetTokenT<vector<float>> subjetMass_;

};

//
// static data member definitions
//

//
// constructors and destructor
//
RUNONBoostedAnalysis::RUNONBoostedAnalysis(const ParameterSet& iConfig):
//	getterOfProducts_(ProcessMatch(*), this) {
//	triggerBits_(consumes<TriggerResults>(iConfig.getParameter<InputTag>("bits"))),
	jetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPt"))),
	jetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetEta"))),
	jetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPhi"))),
	jetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetE"))),
	jetMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetMass"))),
	jetTrimmedMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTrimmedMass"))),
	jetTau1_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTau1"))),
	jetTau2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTau2"))),
	jetTau3_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTau3"))),
	jetNSubjets_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetNSubjets"))),
	jetSubjetIndex0_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex0"))),
	jetSubjetIndex1_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex1"))),
	jetSubjetIndex2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex2"))),
	jetSubjetIndex3_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex3"))),
	jetKeys_(consumes<vector<vector<int>>>(iConfig.getParameter<InputTag>("jetKeys"))),
	jetCSV_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCSV"))),
	jetCSVV1_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCSVV1"))),
	NPV_(consumes<int>(iConfig.getParameter<InputTag>("NPV"))),
	//Jet ID,
	jecFactor_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jecFactor"))),
	neutralHadronEnergy_(consumes<vector<float>>(iConfig.getParameter<InputTag>("neutralHadronEnergy"))),
	neutralEmEnergy_(consumes<vector<float>>(iConfig.getParameter<InputTag>("neutralEmEnergy"))),
	chargeEmEnergy_(consumes<vector<float>>(iConfig.getParameter<InputTag>("chargeEmEnergy"))),
	muonEnergy_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonEnergy"))),
	// Subjets
	subjetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetPt"))),
	subjetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetEta"))),
	subjetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetPhi"))),
	subjetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetE"))),
	subjetMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetMass")))
{
	scale = iConfig.getParameter<double>("scale");
	bjSample = iConfig.getParameter<bool>("bjSample");
	cutTrimmedMassvalue = iConfig.getParameter<double>("cutTrimmedMassvalue");
	cutHTvalue = iConfig.getParameter<double>("cutHTvalue");
	cutAsymvalue = iConfig.getParameter<double>("cutAsymvalue");
	cutCosThetavalue = iConfig.getParameter<double>("cutCosThetavalue");
	cutSubjetPtRatiovalue = iConfig.getParameter<double>("cutSubjetPtRatiovalue");
	cutTau31value = iConfig.getParameter<double>("cutTau31value");
	cutTau21value = iConfig.getParameter<double>("cutTau21value");
}


RUNONBoostedAnalysis::~RUNONBoostedAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNONBoostedAnalysis::analyze(const Event& iEvent, const EventSetup& iSetup) {


	/*vector<Handle< vector<float> > > handles;
	getterOfProducts_.fillHandles(event, handles);
	*/

	Handle<vector<float> > jetPt;
	iEvent.getByToken(jetPt_, jetPt);

	Handle<vector<float> > jetEta;
	iEvent.getByToken(jetEta_, jetEta);

	Handle<vector<float> > jetPhi;
	iEvent.getByToken(jetPhi_, jetPhi);

	Handle<vector<float> > jetE;
	iEvent.getByToken(jetE_, jetE);

	Handle<vector<float> > jetMass;
	iEvent.getByToken(jetMass_, jetMass);

	Handle<vector<float> > jetTrimmedMass;
	iEvent.getByToken(jetTrimmedMass_, jetTrimmedMass);

	Handle<vector<float> > jetTau1;
	iEvent.getByToken(jetTau1_, jetTau1);

	Handle<vector<float> > jetTau2;
	iEvent.getByToken(jetTau2_, jetTau2);

	Handle<vector<float> > jetTau3;
	iEvent.getByToken(jetTau3_, jetTau3);

	Handle<vector<float> > jetNSubjets;
	iEvent.getByToken(jetNSubjets_, jetNSubjets);

	Handle<vector<float> > jetSubjetIndex0;
	iEvent.getByToken(jetSubjetIndex0_, jetSubjetIndex0);

	Handle<vector<float> > jetSubjetIndex1;
	iEvent.getByToken(jetSubjetIndex1_, jetSubjetIndex1);

	Handle<vector<float> > jetSubjetIndex2;
	iEvent.getByToken(jetSubjetIndex2_, jetSubjetIndex2);

	Handle<vector<float> > jetSubjetIndex3;
	iEvent.getByToken(jetSubjetIndex3_, jetSubjetIndex3);

	Handle<vector<vector<int> > > jetKeys;
	iEvent.getByToken(jetKeys_, jetKeys);

	Handle<vector<float> > jetCSV;
	iEvent.getByToken(jetCSV_, jetCSV);

	Handle<vector<float> > jetCSVV1;
	iEvent.getByToken(jetCSVV1_, jetCSVV1);

	Handle<int> NPV;
	iEvent.getByToken(NPV_, NPV);

	/// Jet ID
	Handle<vector<float> > jecFactor;
	iEvent.getByToken(jecFactor_, jecFactor);

	Handle<vector<float> > neutralHadronEnergy;
	iEvent.getByToken(neutralHadronEnergy_, neutralHadronEnergy);

	Handle<vector<float> > neutralEmEnergy;
	iEvent.getByToken(neutralEmEnergy_, neutralEmEnergy);

	Handle<vector<float> > chargeEmEnergy;
	iEvent.getByToken(chargeEmEnergy_, chargeEmEnergy);

	Handle<vector<float> > muonEnergy;
	iEvent.getByToken(muonEnergy_, muonEnergy);

	/// Subjets
	Handle<vector<float> > subjetPt;
	iEvent.getByToken(subjetPt_, subjetPt);

	Handle<vector<float> > subjetEta;
	iEvent.getByToken(subjetEta_, subjetEta);

	Handle<vector<float> > subjetPhi;
	iEvent.getByToken(subjetPhi_, subjetPhi);

	Handle<vector<float> > subjetE;
	iEvent.getByToken(subjetE_, subjetE);

	Handle<vector<float> > subjetMass;
	iEvent.getByToken(subjetMass_, subjetMass);


	cutmap["Processed"] += 1;

	int numPV = *NPV;
	vector< JETtype > JETS;
	vector< float > tmpTriggerMass;
	int numJets = 0;
	double HT = 0;
	double rawHT = 0;
	bool cutHT = 0;
	bool cutMass = 0;
	bool bTagCSV = 0;
	for (size_t i = 0; i < jetPt->size(); i++) {

		if( TMath::Abs( (*jetEta)[i] ) > 2.4 ) continue;

		tmpTriggerMass.push_back( (*jetTrimmedMass)[i] );

		rawHT += (*jetPt)[i];
		histos1D_[ "rawJetPt" ]->Fill( (*jetPt)[i], scale  );

		double jec = 1. / ( (*jecFactor)[i] * (*jetE)[i] );
		double nhf = (*neutralHadronEnergy)[i] * jec;
		double nEMf = (*neutralEmEnergy)[i] * jec;
		double cEMf = (*chargeEmEnergy)[i] * jec;
		double muf = (*muonEnergy)[i] * jec;
		//int npr = (*chargedHadronMultiplicity)[i] + (*neutralHadronMultiplicity)[i] ;  //// REMEMBER TO INCLUDE # of constituents

		bool idL = ( (nhf<0.99) && (nEMf<0.99) && (muf<0.8) && (cEMf<0.9) );

		//if( !idL ) LogWarning("jetID") << (*jetPt)[i] << " " << jec << " " << nhf << " " << nEMf << " " << muf << " " << cEMf;

		if( (*jetPt)[i] > 40  && idL ) { 
			//LogWarning("jetInfo") << i << " " << (*jetPt)[i] << " " << (*jetEta)[i] << " " << (*jetPhi)[i] << " " << (*jetMass)[i];

			HT += (*jetPt)[i];
			++numJets;

			TLorentzVector tmpJet;
			tmpJet.SetPtEtaPhiE( (*jetPt)[i], (*jetEta)[i], (*jetPhi)[i], (*jetE)[i] );

			/// Vector of zeros
			TLorentzVector tmpSubjet0, tmpSubjet1, tmpZeros;
			tmpZeros.SetPtEtaPhiE( 0, 0, 0, 0 );

			//LogWarning("jetSubjetIndex") << (*jetSubjetIndex0)[i] << " " <<  (*jetSubjetIndex1)[i] << " " << (*jetSubjetIndex2)[i] << " " << (*jetSubjetIndex3)[i];

			for (size_t j = 0; j < subjetPt->size(); j++) {
				if( j == (*jetSubjetIndex0)[i] ) {
					//LogWarning("subjets0") << j << " " << (*jetSubjetIndex0)[i] << " " <<  subjetPt->size() << " " << (*subjetPt)[j];
					tmpSubjet0.SetPtEtaPhiE( (*subjetPt)[j], (*subjetEta)[j], (*subjetPhi)[j], (*subjetE)[j] );
				} //else tmpSubjet0 = tmpZeros ; 
					
				if( j == (*jetSubjetIndex1)[i] ) {
					tmpSubjet1.SetPtEtaPhiE( (*subjetPt)[j], (*subjetEta)[j], (*subjetPhi)[j], (*subjetE)[j] );
				} //else tmpSubjet1 = tmpZeros ; 
			}

			//if ( (*jetCSV)[i] > 0.244 ) bTagCSV = 1; 	// CSVL
			if ( (*jetCSV)[i] > 0.679 ) bTagCSV = 1; 	// CSVM
			//if ( (*jetCSVV1)[i] > 0.405 ) bTagCSV = 1; 	// CSVV1L
			//if ( (*jetCSVV1)[i] > 0.783 ) bTagCSV = 1; 	// CSVV1M

			JETtype tmpJET;
			tmpJET.p4 = tmpJet;
			tmpJET.subjet0 = tmpSubjet0;
			tmpJET.subjet1 = tmpSubjet1;
			tmpJET.mass = (*jetMass)[i];
			tmpJET.tau1 = (*jetTau1)[i];
			tmpJET.tau2 = (*jetTau2)[i];
			tmpJET.tau3 = (*jetTau3)[i];
			tmpJET.btagCSV = bTagCSV;
			JETS.push_back( tmpJET );
	   
			histos1D_[ "jetPt" ]->Fill( (*jetPt)[i], scale  );
			histos1D_[ "jetEta" ]->Fill( (*jetEta)[i], scale  );
			histos1D_[ "jetMass" ]->Fill( (*jetMass)[i], scale  );
		}
	}

	//sort(JETS.begin(), JETS.end(), [](const JETtype &p1, const JETtype &p2) { TLorentzVector tmpP1, tmpP2; tmpP1 = p1.p4; tmpP2 = p2.p4;  return tmpP1.M() > tmpP2.M(); }); 
	histos1D_[ "jetNum" ]->Fill( numJets, scale );
	histos1D_[ "NPV" ]->Fill( numPV, scale );
	if ( HT > 0 ) histos1D_[ "HT" ]->Fill( HT, scale  );
	if ( rawHT > 0 ) histos1D_[ "rawHT" ]->Fill( rawHT, scale  );
	if ( HT > cutHTvalue ) cutHT = 1;

	sort(tmpTriggerMass.begin(), tmpTriggerMass.end(), [](const float p1, const float p2) { return p1 > p2; }); 
	if ( ( tmpTriggerMass.size()> 0 ) && ( tmpTriggerMass[0] > cutTrimmedMassvalue) ){
		cutMass = 1;
		histos1D_[ "jetTrimmedMass" ]->Fill( tmpTriggerMass[0], scale  );
	}
	//LogWarning("mass") << tmpTriggerMass[0] << " " << tmpTriggerMass[1] ;
	tmpTriggerMass.clear();
						
	clearVariables();

	if( cutMass && cutHT ){
				
		cutmap["Trigger"] += 1;

		if( ( numJets == 4 ) && ( HT > cutHTvalue ) ){
			cutmap["HT"] += 1;
			histos1D_[ "HT_cutHT" ]->Fill( HT, scale  );
			histos1D_[ "NPV_cutHT" ]->Fill( numPV, scale );
			histos1D_[ "jetNum_cutHT" ]->Fill( numJets, scale );
			histos1D_[ "jetPt_cutHT" ]->Fill( (*jetPt)[0], scale  );
			histos1D_[ "jetPt_cutHT" ]->Fill( (*jetPt)[1], scale  );
			histos1D_[ "jetEta_cutHT" ]->Fill( (*jetEta)[0], scale  );
			histos1D_[ "jetEta_cutHT" ]->Fill( (*jetEta)[1], scale  );
			histos1D_[ "jetMass_cutHT_Ptsort" ]->Fill( (*jetMass)[0], scale );
			histos1D_[ "jetMass_cutHT" ]->Fill( JETS[0].mass, scale );

		
			vector<double> tmpDijetR;
			double dR12 = JETS[0].p4.DeltaR( JETS[1].p4 );
			double dR34 = JETS[2].p4.DeltaR( JETS[3].p4 );
			double dijetR1234 = abs( dR12 - 1 )  + abs( dR34 - 1 );
			tmpDijetR.push_back( dijetR1234 );

			double dR13 = JETS[0].p4.DeltaR( JETS[2].p4 );
			double dR24 = JETS[1].p4.DeltaR( JETS[3].p4 );
			double dijetR1324 = abs( dR13 - 1 )  + abs( dR24 - 1 );
			tmpDijetR.push_back( dijetR1324 );

			double dR14 = JETS[0].p4.DeltaR( JETS[3].p4 );
			double dR23 = JETS[1].p4.DeltaR( JETS[2].p4 );
			double dijetR1423 = abs( dR14 - 1 )  + abs( dR23 - 1 );
			tmpDijetR.push_back( dijetR1423 );

			//LogWarning("test") << min_element(tmpDijetR.begin(), tmpDijetR.end()) - tmpDijetR.begin()  << " " << dijetR1234 << " " << dijetR1324 << " " << dijetR1423 ;
			int minDeltaR = min_element(tmpDijetR.begin(), tmpDijetR.end()) - tmpDijetR.begin();
			TLorentzVector j1, j2, j3, j4;

			if( minDeltaR == 0 ){
				j1 = JETS[0].p4;
				j2 = JETS[1].p4;
				j3 = JETS[2].p4;
				j4 = JETS[3].p4;
			} else if ( minDeltaR == 1 ) {
				j1 = JETS[0].p4;
				j2 = JETS[2].p4;
				j3 = JETS[1].p4;
				j4 = JETS[3].p4;
			} else if ( minDeltaR == 2 ) {
				j1 = JETS[0].p4;
				j2 = JETS[3].p4;
				j3 = JETS[1].p4;
				j4 = JETS[2].p4;
			}


			vector<double> dalitz1, dalitz2;
			double dalitzY1 = -9999;
			double dalitzY2 = -9999;
			double dalitzY3 = -9999;
			double dalitzY4 = -9999;
			double dalitzY5 = -9999;
			double dalitzY6 = -9999;
			double dalitzX1 = -9999; 
			double dalitzX2 = -9999; 
			double dalitzX3 = -9999; 
			double dalitzX4 = -9999; 
			double dalitzX5 = -9999; 
			double dalitzX6 = -9999; 


			
			double m1 = j1.M();
			double m2 = j2.M();
			//double m3 = JETS[2].p4.M();
			//double m4 = JETS[3].p4.M();

			double m12 = ( j1 + j2 ).M() ;
			double m34 = ( j3 + j4 ).M() ;
			double m134 = ( j1 + j3 + j4 ).M() ;
			double m123 = ( j1 + j2 + j3 ).M() ;
			double m124 = ( j1 + j2 + j4 ).M() ;
			double m234 = ( j2 + j3 + j4 ).M() ;
			double m1234 = ( j1 + j2 + j3 + j4 ).M() ;
			

			double tmptilde = pow( m1, 2 ) + pow( m2, 2) + pow( m34, 2 ) + pow( m1234, 2);
			double mtilde12 = pow( m12, 2 ) / tmptilde;
			double mtilde134 = pow( m134, 2 ) / tmptilde;
			double mtilde234 = pow( m234, 2 ) / tmptilde;
			//double tmpMtilde = mtilde12 + mtilde134 + mtilde234;
			//LogWarning("dalitz0") << tmpMtilde << " " << mtilde12 << " " << mtilde134 << " " <<  mtilde234;
			dalitz1.push_back( mtilde12 );
			dalitz1.push_back( mtilde134 );
			dalitz1.push_back( mtilde234 );
			sort( dalitz1.begin(), dalitz1.end(), [](const double &p1, const double &p2) { return p1 > p2; }); 
			//LogWarning("dalitz1") << dalitz1[0] << " " << dalitz1[1] << " " << dalitz1[2];
			histos1D_[ "mu1_cutHT" ]->Fill( dalitz1[0], scale );
			histos1D_[ "mu2_cutHT" ]->Fill( dalitz1[1], scale );
			histos1D_[ "mu3_cutHT" ]->Fill( dalitz1[2], scale );
			histos2D_[ "mu1234_cutHT" ]->Fill( dalitz1[0], dalitz1[2], scale );
			histos2D_[ "mu1234_cutHT" ]->Fill( dalitz1[1], dalitz1[2], scale );
			histos2D_[ "mu1234_cutHT" ]->Fill( dalitz1[0], dalitz1[1], scale );

			dalitzX1 = ( dalitz1[1] + ( 2 * dalitz1[0] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz1234_cutHT" ]->Fill( dalitzX1, dalitz1[1], scale );
			//LogWarning("X1") << dalitzX1 << " " << dalitz1[1] ;
			dalitzX2 = ( dalitz1[2] + ( 2 * dalitz1[0] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz1234_cutHT" ]->Fill( dalitzX2, dalitz1[2], scale );
			//LogWarning("X2") << dalitzX2 << " " << dalitz1[2] ;
			dalitzX3 = ( dalitz1[0] + ( 2 * dalitz1[1] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz1234_cutHT" ]->Fill( dalitzX3, dalitz1[0], scale );
			//LogWarning("X3") << dalitzX3 << " " << dalitz1[0] ;
			dalitzX4 = ( dalitz1[2] + ( 2 * dalitz1[1] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz1234_cutHT" ]->Fill( dalitzX4, dalitz1[2], scale );
			//LogWarning("X4") << dalitzX4 << " " << dalitz1[2] ;
			dalitzX5 = ( dalitz1[0] + ( 2 * dalitz1[2] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz1234_cutHT" ]->Fill( dalitzX5, dalitz1[0], scale );
			//LogWarning("X5") << dalitzX5 << " " << dalitz1[0] ;
			dalitzX6 = ( dalitz1[1] + ( 2 * dalitz1[2] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz1234_cutHT" ]->Fill( dalitzX6, dalitz1[1], scale );
			//LogWarning("X6") << dalitzX6 << " " << dalitz1[1] ;


			double mtilde34 = pow( m34, 2 ) / tmptilde;
			double mtilde123 = pow( m123, 2 ) / tmptilde;
			double mtilde124 = pow( m124, 2 ) / tmptilde;
			dalitz2.push_back( mtilde34 );
			dalitz2.push_back( mtilde123 );
			dalitz2.push_back( mtilde124 );
			sort( dalitz2.begin(), dalitz2.end(), [](const double &p1, const double &p2) { return p1 > p2; }); 
			histos1D_[ "mu4_cutHT" ]->Fill( dalitz2[0], scale );
			histos1D_[ "mu5_cutHT" ]->Fill( dalitz2[1], scale );
			histos1D_[ "mu6_cutHT" ]->Fill( dalitz2[2], scale );
			histos2D_[ "mu3412_cutHT" ]->Fill( dalitz2[0], dalitz2[2], scale );
			histos2D_[ "mu3412_cutHT" ]->Fill( dalitz2[1], dalitz2[2], scale );
			histos2D_[ "mu3412_cutHT" ]->Fill( dalitz2[0], dalitz2[1], scale );

			dalitzY1 = ( dalitz2[1] + ( 2 * dalitz2[0] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz3412_cutHT" ]->Fill( dalitzY1, dalitz2[1], scale );
			dalitzY2 = ( dalitz2[2] + ( 2 * dalitz2[0] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz3412_cutHT" ]->Fill( dalitzY2, dalitz2[2], scale );
			dalitzY3 = ( dalitz2[0] + ( 2 * dalitz2[1] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz3412_cutHT" ]->Fill( dalitzY3, dalitz2[0], scale );
			dalitzY4 = ( dalitz2[2] + ( 2 * dalitz2[1] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz3412_cutHT" ]->Fill( dalitzY4, dalitz2[2], scale );
			dalitzY5 = ( dalitz2[0] + ( 2 * dalitz2[2] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz3412_cutHT" ]->Fill( dalitzY5, dalitz2[0], scale );
			dalitzY6 = ( dalitz2[1] + ( 2 * dalitz2[2] ) ) / TMath::Sqrt(3);
			histos2D_[ "dalitz3412_cutHT" ]->Fill( dalitzY6, dalitz2[1], scale );



		}
	}
	JETS.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void RUNONBoostedAnalysis::beginJob() {

	RUNAtree = fs_->make< TTree >("RUNATree", "RUNATree"); 
	RUNAtree->Branch( "jetsPt", "vector<float>", &jetsPt);
	RUNAtree->Branch( "jetsEta", "vector<float>", &jetsEta);
	RUNAtree->Branch( "jetsPhi", "vector<float>", &jetsPhi);
	RUNAtree->Branch( "jetsE", "vector<float>", &jetsE);
	RUNAtree->Branch( "jet1SubjetsPt", "vector<float>", &jet1SubjetsPt);
	RUNAtree->Branch( "jet1SubjetsEta", "vector<float>", &jet1SubjetsEta);
	RUNAtree->Branch( "jet1SubjetsPhi", "vector<float>", &jet1SubjetsPhi);
	RUNAtree->Branch( "jet1SubjetsE", "vector<float>", &jet1SubjetsE);
	RUNAtree->Branch( "jet2SubjetsPt", "vector<float>", &jet2SubjetsPt);
	RUNAtree->Branch( "jet2SubjetsEta", "vector<float>", &jet2SubjetsEta);
	RUNAtree->Branch( "jet2SubjetsPhi", "vector<float>", &jet2SubjetsPhi);
	RUNAtree->Branch( "jet2SubjetsE", "vector<float>", &jet2SubjetsE);

	histos1D_[ "rawJetPt" ] = fs_->make< TH1D >( "rawJetPt", "rawJetPt", 100, 0., 1000. );
	histos1D_[ "rawJetPt" ]->Sumw2();
	histos1D_[ "rawHT" ] = fs_->make< TH1D >( "rawHT", "rawHT", 150, 0., 1500. );
	histos1D_[ "rawHT" ]->Sumw2();

	histos1D_[ "jetPt" ] = fs_->make< TH1D >( "jetPt", "jetPt", 100, 0., 1000. );
	histos1D_[ "jetPt" ]->Sumw2();
	histos1D_[ "jetEta" ] = fs_->make< TH1D >( "jetEta", "jetEta", 100, -5., 5. );
	histos1D_[ "jetEta" ]->Sumw2();
	histos1D_[ "jetNum" ] = fs_->make< TH1D >( "jetNum", "jetNum", 10, 0., 10. );
	histos1D_[ "jetNum" ]->Sumw2();
	histos1D_[ "jetMass" ] = fs_->make< TH1D >( "jetMass", "jetMass", 30, 0., 300. );
	histos1D_[ "jetMass" ]->Sumw2();
	histos1D_[ "jetTrimmedMass" ] = fs_->make< TH1D >( "jetTrimmedMass", "jetTrimmedMass", 30, 0., 300. );
	histos1D_[ "jetTrimmedMass" ]->Sumw2();
	histos1D_[ "HT" ] = fs_->make< TH1D >( "HT", "HT", 150, 0., 1500. );
	histos1D_[ "HT" ]->Sumw2();
	histos1D_[ "NPV" ] = fs_->make< TH1D >( "NPV", "NPV", 80, 0., 80. );
	histos1D_[ "NPV" ]->Sumw2();

	histos1D_[ "HT_cutHT" ] = fs_->make< TH1D >( "HT_cutHT", "HT_cutHT", 150, 0., 1500. );
	histos1D_[ "HT_cutHT" ]->Sumw2();
	histos1D_[ "NPV_cutHT" ] = fs_->make< TH1D >( "NPV_cutHT", "NPV_cutHT", 80, 0., 80. );
	histos1D_[ "NPV_cutHT" ]->Sumw2();
	histos1D_[ "jetPt_cutHT" ] = fs_->make< TH1D >( "jetPt_cutHT", "jetPt_cutHT", 100, 0., 1000. );
	histos1D_[ "jetPt_cutHT" ]->Sumw2();
	histos1D_[ "jetEta_cutHT" ] = fs_->make< TH1D >( "jetEta_cutHT", "jetEta_cutHT", 100, -5., 5. );
	histos1D_[ "jetEta_cutHT" ]->Sumw2();
	histos1D_[ "jetNum_cutHT" ] = fs_->make< TH1D >( "jetNum_cutHT", "jetNum_cutHT", 10, 0., 10. );
	histos1D_[ "jetNum_cutHT" ]->Sumw2();
	histos1D_[ "jetMass_cutHT_Ptsort" ] = fs_->make< TH1D >( "jetMass_cutHT_Ptsort", "jetMass_cutHT_Ptsort", 30, 0., 300. );
	histos1D_[ "jetMass_cutHT_Ptsort" ]->Sumw2();
	histos1D_[ "jetMass_cutHT" ] = fs_->make< TH1D >( "jetMass_cutHT", "jetMass_cutHT", 30, 0., 300. );
	histos1D_[ "jetMass_cutHT" ]->Sumw2();
	histos1D_[ "massAsymmetry_cutHT" ] = fs_->make< TH1D >( "massAsymmetry_cutHT", "massAsymmetry_cutHT", 20, 0., 1. );
	histos1D_[ "massAsymmetry_cutHT" ]->Sumw2();
	histos1D_[ "massAve_cutHT" ] = fs_->make< TH1D >( "massAve_cutHT", "massAve_cutHT", 30, 0., 300. );
	histos1D_[ "massAve_cutHT" ]->Sumw2();
	histos1D_[ "massAveLowPU_cutHT" ] = fs_->make< TH1D >( "massAveLowPU_cutHT", "massAveLowPU_cutHT", 30, 0., 300. );
	histos1D_[ "massAveLowPU_cutHT" ]->Sumw2();
	histos1D_[ "massAveMedPU_cutHT" ] = fs_->make< TH1D >( "massAveMedPU_cutHT", "massAveMedPU_cutHT", 30, 0., 300. );
	histos1D_[ "massAveMedPU_cutHT" ]->Sumw2();
	histos1D_[ "massAveHighPU_cutHT" ] = fs_->make< TH1D >( "massAveHighPU_cutHT", "massAveHighPU_cutHT", 30, 0., 300. );
	histos1D_[ "massAveHighPU_cutHT" ]->Sumw2();
	histos2D_[ "massAvevsJet1Mass_cutHT" ] = fs_->make< TH2D >( "massAvevsJet1Mass_cutHT", "massAvevsJet1Mass_cutHT", 30, 0., 300., 30, 0., 300. );
	histos2D_[ "massAvevsJet1Mass_cutHT" ]->Sumw2();
	histos2D_[ "massAvevsJet2Mass_cutHT" ] = fs_->make< TH2D >( "massAvevsJet2Mass_cutHT", "massAvevsJet2Mass_cutHT", 30, 0., 300., 30, 0., 300. );
	histos2D_[ "massAvevsJet2Mass_cutHT" ]->Sumw2();
	histos2D_[ "jet1vs2Mass_cutHT" ] = fs_->make< TH2D >( "jet1vs2Mass_cutHT", "jet1vs2Mass_cutHT", 30, 0., 300., 30, 0., 300. );
	histos2D_[ "jet1vs2Mass_cutHT" ]->Sumw2();
	histos1D_[ "cosThetaStar_cutHT" ] = fs_->make< TH1D >( "cosThetaStar_cutHT", "cosThetaStar_cutHT", 20, 0., 1. );
	histos1D_[ "cosThetaStar_cutHT" ]->Sumw2();
	histos1D_[ "jet1Tau1_cutHT" ] = fs_->make< TH1D >( "jet1Tau1_cutHT", "jet1Tau1_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Tau1_cutHT" ]->Sumw2();
	histos1D_[ "jet1Tau2_cutHT" ] = fs_->make< TH1D >( "jet1Tau2_cutHT", "jet1Tau2_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Tau2_cutHT" ]->Sumw2();
	histos1D_[ "jet1Tau3_cutHT" ] = fs_->make< TH1D >( "jet1Tau3_cutHT", "jet1Tau3_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Tau3_cutHT" ]->Sumw2();
	histos1D_[ "jet1Tau21_cutHT" ] = fs_->make< TH1D >( "jet1Tau21_cutHT", "jet1Tau21_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Tau21_cutHT" ]->Sumw2();
	histos1D_[ "jet1Tau31_cutHT" ] = fs_->make< TH1D >( "jet1Tau31_cutHT", "jet1Tau31_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Tau31_cutHT" ]->Sumw2();
	histos1D_[ "jet1Tau32_cutHT" ] = fs_->make< TH1D >( "jet1Tau32_cutHT", "jet1Tau32_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Tau32_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet1Pt_cutHT" ] = fs_->make< TH1D >( "jet1Subjet1Pt_cutHT", "jet1Subjet1Pt_cutHT", 100, 0., 1000. );
	histos1D_[ "jet1Subjet1Pt_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet2Pt_cutHT" ] = fs_->make< TH1D >( "jet1Subjet2Pt_cutHT", "jet1Subjet2Pt_cutHT", 100, 0., 1000. );
	histos1D_[ "jet1Subjet2Pt_cutHT" ]->Sumw2();
	histos1D_[ "jet1SubjetPtRatio_cutHT" ] = fs_->make< TH1D >( "jet1SubjetPtRatio_cutHT", "jet1SubjetPtRatio_cutHT", 20, 0, 1.);
	histos1D_[ "jet1SubjetPtRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet1SubjetMass21Ratio_cutHT" ] = fs_->make< TH1D >( "jet1SubjetMass21Ratio_cutHT", "jet1SubjetMass21Ratio_cutHT", 20, 0., 1. );
	histos1D_[ "jet1SubjetMass21Ratio_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet112MassRatio_cutHT" ] = fs_->make< TH1D >( "jet1Subjet112MassRatio_cutHT", "jet1Subjet112MassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Subjet112MassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet1JetMassRatio_cutHT" ] = fs_->make< TH1D >( "jet1Subjet1JetMassRatio_cutHT", "jet1Subjet1JetMassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Subjet1JetMassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet212MassRatio_cutHT" ] = fs_->make< TH1D >( "jet1Subjet212MassRatio_cutHT", "jet1Subjet212MassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Subjet212MassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet2JetMassRatio_cutHT" ] = fs_->make< TH1D >( "jet1Subjet2JetMassRatio_cutHT", "jet1Subjet2JetMassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet1Subjet2JetMassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet1Mass_cutHT" ] = fs_->make< TH1D >( "jet1Subjet1Mass_cutHT", "jet1Subjet1Mass_cutHT", 20, 0., 100. );
	histos1D_[ "jet1Subjet1Mass_cutHT" ]->Sumw2();
	histos1D_[ "jet1Subjet2Mass_cutHT" ] = fs_->make< TH1D >( "jet1Subjet2Mass_cutHT", "jet1Subjet2Mass_cutHT", 20, 0., 100. );
	histos1D_[ "jet1Subjet2Mass_cutHT" ]->Sumw2();
	histos2D_[ "jet1Subjet12Mass_cutHT" ] = fs_->make< TH2D >( "jet1Subjet12Mass_cutHT", "jet1Subjet12Mass_cutHT", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "jet1Subjet12Mass_cutHT" ]->Sumw2();
	histos2D_[ "jet1Subjet112vs212MassRatio_cutHT" ] = fs_->make< TH2D >( "jet1Subjet112vs212MassRatio_cutHT", "jet1Subjet112vs212MassRatio_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "jet1Subjet112vs212MassRatio_cutHT" ]->Sumw2();
	histos2D_[ "jet1Subjet1JetvsSubjet2JetMassRatio_cutHT" ] = fs_->make< TH2D >( "jet1Subjet1JetvsSubjet2JetMassRatio_cutHT", "jet1Subjet1JetvsSubjet2JetMassRatio_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "jet1Subjet1JetvsSubjet2JetMassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet1Pt_cutHT" ] = fs_->make< TH1D >( "jet2Subjet1Pt_cutHT", "jet2Subjet1Pt_cutHT", 100, 0., 1000. );
	histos1D_[ "jet2Subjet1Pt_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet2Pt_cutHT" ] = fs_->make< TH1D >( "jet2Subjet2Pt_cutHT", "jet2Subjet2Pt_cutHT", 100, 0., 1000. );
	histos1D_[ "jet2Subjet2Pt_cutHT" ]->Sumw2();
	histos1D_[ "jet2SubjetPtRatio_cutHT" ] = fs_->make< TH1D >( "jet2SubjetPtRatio_cutHT", "jet2SubjetPtRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet2SubjetPtRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet2SubjetMass21Ratio_cutHT" ] = fs_->make< TH1D >( "jet2SubjetMass21Ratio_cutHT", "jet2SubjetMass21Ratio_cutHT", 20, 0., 1. );
	histos1D_[ "jet2SubjetMass21Ratio_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet112MassRatio_cutHT" ] = fs_->make< TH1D >( "jet2Subjet112MassRatio_cutHT", "jet2Subjet112MassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet2Subjet112MassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet1JetMassRatio_cutHT" ] = fs_->make< TH1D >( "jet2Subjet1JetMassRatio_cutHT", "jet2Subjet1JetMassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet2Subjet1JetMassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet212MassRatio_cutHT" ] = fs_->make< TH1D >( "jet2Subjet212MassRatio_cutHT", "jet2Subjet212MassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet2Subjet212MassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet2JetMassRatio_cutHT" ] = fs_->make< TH1D >( "jet2Subjet2JetMassRatio_cutHT", "jet2Subjet2JetMassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "jet2Subjet2JetMassRatio_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet1Mass_cutHT" ] = fs_->make< TH1D >( "jet2Subjet1Mass_cutHT", "jet2Subjet1Mass_cutHT", 20, 0., 100.);
	histos1D_[ "jet2Subjet1Mass_cutHT" ]->Sumw2();
	histos1D_[ "jet2Subjet2Mass_cutHT" ] = fs_->make< TH1D >( "jet2Subjet2Mass_cutHT", "jet2Subjet2Mass_cutHT", 20, 0., 100. );
	histos1D_[ "jet2Subjet2Mass_cutHT" ]->Sumw2();
	histos2D_[ "jet2Subjet12Mass_cutHT" ] = fs_->make< TH2D >( "jet2Subjet12Mass_cutHT", "jet2Subjet12Mass_cutHT", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "jet2Subjet12Mass_cutHT" ]->Sumw2();
	histos2D_[ "jet2Subjet112vs212MassRatio_cutHT" ] = fs_->make< TH2D >( "jet2Subjet112vs212MassRatio_cutHT", "jet2Subjet112vs212MassRatio_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "jet2Subjet112vs212MassRatio_cutHT" ]->Sumw2();
	histos2D_[ "jet2Subjet1JetvsSubjet2JetMassRatio_cutHT" ] = fs_->make< TH2D >( "jet2Subjet1JetvsSubjet2JetMassRatio_cutHT", "jet2Subjet1JetvsSubjet2JetMassRatio_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "jet2Subjet1JetvsSubjet2JetMassRatio_cutHT" ]->Sumw2();

	histos1D_[ "subjetPtRatio_cutHT" ] = fs_->make< TH1D >( "subjetPtRatio_cutHT", "subjetPtRatio_cutHT", 20, 0., 1. );
	histos1D_[ "subjetPtRatio_cutHT" ]->Sumw2();
	histos1D_[ "subjetMass21Ratio_cutHT" ] = fs_->make< TH1D >( "subjetMass21Ratio_cutHT", "subjetMass21Ratio_cutHT", 20, 0., 1. );
	histos1D_[ "subjetMass21Ratio_cutHT" ]->Sumw2();
	histos1D_[ "subjet112MassRatio_cutHT" ] = fs_->make< TH1D >( "subjet112MassRatio_cutHT", "subjet112MassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "subjet112MassRatio_cutHT" ]->Sumw2();
	histos1D_[ "subjet212MassRatio_cutHT" ] = fs_->make< TH1D >( "subjet212MassRatio_cutHT", "subjet212MassRatio_cutHT", 20, 0., 1. );
	histos1D_[ "subjet212MassRatio_cutHT" ]->Sumw2();
	histos2D_[ "subjet12Mass_cutHT" ] = fs_->make< TH2D >( "subjet12Mass_cutHT", "subjet12Mass_cutHT", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "subjet12Mass_cutHT" ]->Sumw2();
	histos2D_[ "dijetCorr_cutHT" ] = fs_->make< TH2D >( "dijetCorr_cutHT", "dijetCorr_cutHT", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorr_cutHT" ]->Sumw2();
	histos2D_[ "dijetCorrPhi_cutHT" ] = fs_->make< TH2D >( "dijetCorrPhi_cutHT", "dijetCorrPhi_cutHT", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorrPhi_cutHT" ]->Sumw2();
	histos2D_[ "subjet112vs212MassRatio_cutHT" ] = fs_->make< TH2D >( "subjet112vs212MassRatio_cutHT", "subjet112vs212MassRatio_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet112vs212MassRatio_cutHT" ]->Sumw2();
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutHT" ] = fs_->make< TH2D >( "subjet1JetvsSubjet2JetMassRatio_cutHT", "subjet1JetvsSubjet2JetMassRatio_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutHT" ]->Sumw2();
	histos1D_[ "subjetPolAngle13412_cutHT" ] = fs_->make< TH1D >( "subjetPolAngle13412_cutHT", "subjetPolAngle13412_cutHT", 20, 0., 1. );
	histos1D_[ "subjetPolAngle13412_cutHT" ]->Sumw2();
	histos1D_[ "subjetPolAngle31234_cutHT" ] = fs_->make< TH1D >( "subjetPolAngle31234_cutHT", "subjetPolAngle31234_cutHT", 20, 0., 1. );
	histos1D_[ "subjetPolAngle31234_cutHT" ]->Sumw2();
	histos2D_[ "subjetPolAngle13412vs31234_cutHT" ] = fs_->make< TH2D >( "subjetPolAngle13412vs31234_cutHT", "subjetPolAngle13412vs31234_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjetPolAngle13412vs31234_cutHT" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle13412_cutHT" ] = fs_->make< TH1D >( "tmpSubjetPolAngle13412_cutHT", "tmpSubjetPolAngle13412_cutHT", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle13412_cutHT" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle31234_cutHT" ] = fs_->make< TH1D >( "tmpSubjetPolAngle31234_cutHT", "tmpSubjetPolAngle31234_cutHT", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle31234_cutHT" ]->Sumw2();
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutHT" ] = fs_->make< TH2D >( "tmpSubjetPolAngle13412vs31234_cutHT", "tmpSubjetPolAngle13412vs31234_cutHT", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutHT" ]->Sumw2();
	histos1D_[ "mu1_cutHT" ] = fs_->make< TH1D >( "mu1_cutHT", "mu1_cutHT", 150, 0., 1.5 );
	histos1D_[ "mu1_cutHT" ]->Sumw2();
	histos1D_[ "mu2_cutHT" ] = fs_->make< TH1D >( "mu2_cutHT", "mu2_cutHT", 150, 0., 1.5 );
	histos1D_[ "mu2_cutHT" ]->Sumw2();
	histos1D_[ "mu3_cutHT" ] = fs_->make< TH1D >( "mu3_cutHT", "mu3_cutHT", 150, 0., 1.5 );
	histos1D_[ "mu3_cutHT" ]->Sumw2();
	histos2D_[ "mu1234_cutHT" ] = fs_->make< TH2D >( "mu1234_cutHT", "mu1234_cutHT", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu1234_cutHT" ]->Sumw2();
	histos2D_[ "dalitz1234_1_cutHT" ] = fs_->make< TH2D >( "dalitz1234_1_cutHT", "dalitz1234_1_cutHT", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz1234_1_cutHT" ]->Sumw2();
	histos2D_[ "dalitz1234_2_cutHT" ] = fs_->make< TH2D >( "dalitz1234_2_cutHT", "dalitz1234_2_cutHT", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz1234_2_cutHT" ]->Sumw2();
	histos2D_[ "dalitz1234_3_cutHT" ] = fs_->make< TH2D >( "dalitz1234_3_cutHT", "dalitz1234_3_cutHT", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz1234_3_cutHT" ]->Sumw2();
	histos2D_[ "dalitz1234_cutHT" ] = fs_->make< TH2D >( "dalitz1234_cutHT", "dalitz1234_cutHT", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz1234_cutHT" ]->Sumw2();
	histos1D_[ "mu4_cutHT" ] = fs_->make< TH1D >( "mu4_cutHT", "mu4_cutHT", 150, 0., 1.5 );
	histos1D_[ "mu4_cutHT" ]->Sumw2();
	histos1D_[ "mu5_cutHT" ] = fs_->make< TH1D >( "mu5_cutHT", "mu5_cutHT", 150, 0., 1.5 );
	histos1D_[ "mu5_cutHT" ]->Sumw2();
	histos1D_[ "mu6_cutHT" ] = fs_->make< TH1D >( "mu6_cutHT", "mu6_cutHT", 150, 0., 1.5 );
	histos1D_[ "mu6_cutHT" ]->Sumw2();
	histos2D_[ "mu3412_cutHT" ] = fs_->make< TH2D >( "mu3412_cutHT", "mu3412_cutHT", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu3412_cutHT" ]->Sumw2();
	histos2D_[ "dalitz3412_cutHT" ] = fs_->make< TH2D >( "dalitz3412_cutHT", "dalitz3412_cutHT", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz3412_cutHT" ]->Sumw2();

	histos1D_[ "massAve_cutAsym" ] = fs_->make< TH1D >( "massAve_cutAsym", "massAve_cutAsym", 30, 0., 300. );
	histos1D_[ "massAve_cutAsym" ]->Sumw2();
	histos1D_[ "massAveLowPU_cutAsym" ] = fs_->make< TH1D >( "massAveLowPU_cutAsym", "massAveLowPU_cutAsym", 30, 0., 300. );
	histos1D_[ "massAveLowPU_cutAsym" ]->Sumw2();
	histos1D_[ "massAveMedPU_cutAsym" ] = fs_->make< TH1D >( "massAveMedPU_cutAsym", "massAveMedPU_cutAsym", 30, 0., 300. );
	histos1D_[ "massAveMedPU_cutAsym" ]->Sumw2();
	histos1D_[ "massAveHighPU_cutAsym" ] = fs_->make< TH1D >( "massAveHighPU_cutAsym", "massAveHighPU_cutAsym", 30, 0., 300. );
	histos1D_[ "massAveHighPU_cutAsym" ]->Sumw2();
	histos1D_[ "cosThetaStar_cutAsym" ] = fs_->make< TH1D >( "cosThetaStar_cutAsym", "cosThetaStar_cutAsym", 20, 0., 1. );
	histos1D_[ "cosThetaStar_cutAsym" ]->Sumw2();
	histos1D_[ "jet1Tau21_cutAsym" ] = fs_->make< TH1D >( "jet1Tau21_cutAsym", "jet1Tau21_cutAsym", 20, 0., 1. );
	histos1D_[ "jet1Tau21_cutAsym" ]->Sumw2();
	histos1D_[ "jet1Tau31_cutAsym" ] = fs_->make< TH1D >( "jet1Tau31_cutAsym", "jet1Tau31_cutAsym", 20, 0., 1. );
	histos1D_[ "jet1Tau31_cutAsym" ]->Sumw2();
	histos1D_[ "jet1Tau32_cutAsym" ] = fs_->make< TH1D >( "jet1Tau32_cutAsym", "jet1Tau32_cutAsym", 20, 0., 1. );
	histos1D_[ "jet1Tau32_cutAsym" ]->Sumw2();
	histos1D_[ "subjetPtRatio_cutAsym" ] = fs_->make< TH1D >( "subjetPtRatio_cutAsym", "subjetPtRatio_cutAsym", 20, 0., 1. );
	histos1D_[ "subjetPtRatio_cutAsym" ]->Sumw2();
	histos1D_[ "subjetMass21Ratio_cutAsym" ] = fs_->make< TH1D >( "subjetMass21Ratio_cutAsym", "subjetMass21Ratio_cutAsym", 20, 0., 1. );
	histos1D_[ "subjetMass21Ratio_cutAsym" ]->Sumw2();
	histos1D_[ "subjet112MassRatio_cutAsym" ] = fs_->make< TH1D >( "subjet112MassRatio_cutAsym", "subjet112MassRatio_cutAsym", 20, 0., 1. );
	histos1D_[ "subjet112MassRatio_cutAsym" ]->Sumw2();
	histos1D_[ "subjet212MassRatio_cutAsym" ] = fs_->make< TH1D >( "subjet212MassRatio_cutAsym", "subjet212MassRatio_cutAsym", 20, 0., 1. );
	histos1D_[ "subjet212MassRatio_cutAsym" ]->Sumw2();
	histos1D_[ "subjetPolAngle13412_cutAsym" ] = fs_->make< TH1D >( "subjetPolAngle13412_cutAsym", "subjetPolAngle13412_cutAsym", 20, 0., 1. );
	histos1D_[ "subjetPolAngle13412_cutAsym" ]->Sumw2();
	histos1D_[ "subjetPolAngle31234_cutAsym" ] = fs_->make< TH1D >( "subjetPolAngle31234_cutAsym", "subjetPolAngle31234_cutAsym", 20, 0., 1. );
	histos1D_[ "subjetPolAngle31234_cutAsym" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle13412_cutAsym" ] = fs_->make< TH1D >( "tmpSubjetPolAngle13412_cutAsym", "tmpSubjetPolAngle13412_cutAsym", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle13412_cutAsym" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle31234_cutAsym" ] = fs_->make< TH1D >( "tmpSubjetPolAngle31234_cutAsym", "tmpSubjetPolAngle31234_cutAsym", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle31234_cutAsym" ]->Sumw2();
	histos2D_[ "subjet12Mass_cutAsym" ] = fs_->make< TH2D >( "subjet12Mass_cutAsym", "subjet12Mass_cutAsym", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "subjet12Mass_cutAsym" ]->Sumw2();
	histos2D_[ "dijetCorr_cutAsym" ] = fs_->make< TH2D >( "dijetCorr_cutAsym", "dijetCorr_cutAsym", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorr_cutAsym" ]->Sumw2();
	histos2D_[ "dijetCorrPhi_cutAsym" ] = fs_->make< TH2D >( "dijetCorrPhi_cutAsym", "dijetCorrPhi_cutAsym", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorrPhi_cutAsym" ]->Sumw2();
	histos2D_[ "subjet112vs212MassRatio_cutAsym" ] = fs_->make< TH2D >( "subjet112vs212MassRatio_cutAsym", "subjet112vs212MassRatio_cutAsym", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet112vs212MassRatio_cutAsym" ]->Sumw2();
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutAsym" ] = fs_->make< TH2D >( "subjet1JetvsSubjet2JetMassRatio_cutAsym", "subjet1JetvsSubjet2JetMassRatio_cutAsym", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutAsym" ]->Sumw2();
	histos2D_[ "subjetPolAngle13412vs31234_cutAsym" ] = fs_->make< TH2D >( "subjetPolAngle13412vs31234_cutAsym", "subjetPolAngle13412vs31234_cutAsym", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjetPolAngle13412vs31234_cutAsym" ]->Sumw2();
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutAsym" ] = fs_->make< TH2D >( "tmpSubjetPolAngle13412vs31234_cutAsym", "tmpSubjetPolAngle13412vs31234_cutAsym", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutAsym" ]->Sumw2();
	histos2D_[ "mu1234_cutAsym" ] = fs_->make< TH2D >( "mu1234_cutAsym", "mu1234_cutAsym", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu1234_cutAsym" ]->Sumw2();
	histos2D_[ "mu3412_cutAsym" ] = fs_->make< TH2D >( "mu3412_cutAsym", "mu3412_cutAsym", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu3412_cutAsym" ]->Sumw2();
	histos2D_[ "dalitz1234_cutAsym" ] = fs_->make< TH2D >( "dalitz1234_cutAsym", "dalitz1234_cutAsym", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz1234_cutAsym" ]->Sumw2();
	histos2D_[ "dalitz3412_cutAsym" ] = fs_->make< TH2D >( "dalitz3412_cutAsym", "dalitz3412_cutAsym", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz3412_cutAsym" ]->Sumw2();

	histos1D_[ "massAve_cutCosTheta" ] = fs_->make< TH1D >( "massAve_cutCosTheta", "massAve_cutCosTheta", 30, 0., 300. );
	histos1D_[ "massAve_cutCosTheta" ]->Sumw2();
	histos1D_[ "massAveLowPU_cutCosTheta" ] = fs_->make< TH1D >( "massAveLowPU_cutCosTheta", "massAveLowPU_cutCosTheta", 30, 0., 300. );
	histos1D_[ "massAveLowPU_cutCosTheta" ]->Sumw2();
	histos1D_[ "massAveMedPU_cutCosTheta" ] = fs_->make< TH1D >( "massAveMedPU_cutCosTheta", "massAveMedPU_cutCosTheta", 30, 0., 300. );
	histos1D_[ "massAveMedPU_cutCosTheta" ]->Sumw2();
	histos1D_[ "massAveHighPU_cutCosTheta" ] = fs_->make< TH1D >( "massAveHighPU_cutCosTheta", "massAveHighPU_cutCosTheta", 30, 0., 300. );
	histos1D_[ "massAveHighPU_cutCosTheta" ]->Sumw2();
	histos1D_[ "jet1Tau21_cutCosTheta" ] = fs_->make< TH1D >( "jet1Tau21_cutCosTheta", "jet1Tau21_cutCosTheta", 20, 0., 1. );
	histos1D_[ "jet1Tau21_cutCosTheta" ]->Sumw2();
	histos1D_[ "jet1Tau31_cutCosTheta" ] = fs_->make< TH1D >( "jet1Tau31_cutCosTheta", "jet1Tau31_cutCosTheta", 20, 0., 1. );
	histos1D_[ "jet1Tau31_cutCosTheta" ]->Sumw2();
	histos1D_[ "jet1Tau32_cutCosTheta" ] = fs_->make< TH1D >( "jet1Tau32_cutCosTheta", "jet1Tau32_cutCosTheta", 20, 0., 1. );
	histos1D_[ "jet1Tau32_cutCosTheta" ]->Sumw2();
	histos1D_[ "subjetPtRatio_cutCosTheta" ] = fs_->make< TH1D >( "subjetPtRatio_cutCosTheta", "subjetPtRatio_cutCosTheta", 20, 0., 1. );
	histos1D_[ "subjetPtRatio_cutCosTheta" ]->Sumw2();
	histos1D_[ "subjetMass21Ratio_cutCosTheta" ] = fs_->make< TH1D >( "subjetMass21Ratio_cutCosTheta", "subjetMass21Ratio_cutCosTheta", 20, 0., 1. );
	histos1D_[ "subjetMass21Ratio_cutCosTheta" ]->Sumw2();
	histos1D_[ "subjet112MassRatio_cutCosTheta" ] = fs_->make< TH1D >( "subjet112MassRatio_cutCosTheta", "subjet112MassRatio_cutCosTheta", 20, 0., 1. );
	histos1D_[ "subjet112MassRatio_cutCosTheta" ]->Sumw2();
	histos1D_[ "subjet212MassRatio_cutCosTheta" ] = fs_->make< TH1D >( "subjet212MassRatio_cutCosTheta", "subjet212MassRatio_cutCosTheta", 20, 0., 1. );
	histos1D_[ "subjet212MassRatio_cutCosTheta" ]->Sumw2();
	histos1D_[ "subjetPolAngle13412_cutCosTheta" ] = fs_->make< TH1D >( "subjetPolAngle13412_cutCosTheta", "subjetPolAngle13412_cutCosTheta", 20, 0., 1. );
	histos1D_[ "subjetPolAngle13412_cutCosTheta" ]->Sumw2();
	histos1D_[ "subjetPolAngle31234_cutCosTheta" ] = fs_->make< TH1D >( "subjetPolAngle31234_cutCosTheta", "subjetPolAngle31234_cutCosTheta", 20, 0., 1. );
	histos1D_[ "subjetPolAngle31234_cutCosTheta" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle13412_cutCosTheta" ] = fs_->make< TH1D >( "tmpSubjetPolAngle13412_cutCosTheta", "tmpSubjetPolAngle13412_cutCosTheta", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle13412_cutCosTheta" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle31234_cutCosTheta" ] = fs_->make< TH1D >( "tmpSubjetPolAngle31234_cutCosTheta", "tmpSubjetPolAngle31234_cutCosTheta", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle31234_cutCosTheta" ]->Sumw2();
	histos2D_[ "subjet12Mass_cutCosTheta" ] = fs_->make< TH2D >( "subjet12Mass_cutCosTheta", "subjet12Mass_cutCosTheta", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "subjet12Mass_cutCosTheta" ]->Sumw2();
	histos2D_[ "dijetCorr_cutCosTheta" ] = fs_->make< TH2D >( "dijetCorr_cutCosTheta", "dijetCorr_cutCosTheta", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorr_cutCosTheta" ]->Sumw2();
	histos2D_[ "dijetCorrPhi_cutCosTheta" ] = fs_->make< TH2D >( "dijetCorrPhi_cutCosTheta", "dijetCorrPhi_cutCosTheta", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorrPhi_cutCosTheta" ]->Sumw2();
	histos2D_[ "subjet112vs212MassRatio_cutCosTheta" ] = fs_->make< TH2D >( "subjet112vs212MassRatio_cutCosTheta", "subjet112vs212MassRatio_cutCosTheta", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet112vs212MassRatio_cutCosTheta" ]->Sumw2();
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutCosTheta" ] = fs_->make< TH2D >( "subjet1JetvsSubjet2JetMassRatio_cutCosTheta", "subjet1JetvsSubjet2JetMassRatio_cutCosTheta", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutCosTheta" ]->Sumw2();
	histos2D_[ "subjetPolAngle13412vs31234_cutCosTheta" ] = fs_->make< TH2D >( "subjetPolAngle13412vs31234_cutCosTheta", "subjetPolAngle13412vs31234_cutCosTheta", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjetPolAngle13412vs31234_cutCosTheta" ]->Sumw2();
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutCosTheta" ] = fs_->make< TH2D >( "tmpSubjetPolAngle13412vs31234_cutCosTheta", "tmpSubjetPolAngle13412vs31234_cutCosTheta", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutCosTheta" ]->Sumw2();
	histos2D_[ "mu1234_cutCosTheta" ] = fs_->make< TH2D >( "mu1234_cutCosTheta", "mu1234_cutCosTheta", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu1234_cutCosTheta" ]->Sumw2();
	histos2D_[ "mu3412_cutCosTheta" ] = fs_->make< TH2D >( "mu3412_cutCosTheta", "mu3412_cutCosTheta", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu3412_cutCosTheta" ]->Sumw2();
	histos2D_[ "dalitz1234_cutCosTheta" ] = fs_->make< TH2D >( "dalitz1234_cutCosTheta", "dalitz1234_cutCosTheta", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz1234_cutCosTheta" ]->Sumw2();
	histos2D_[ "dalitz3412_cutCosTheta" ] = fs_->make< TH2D >( "dalitz3412_cutCosTheta", "dalitz3412_cutCosTheta", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz3412_cutCosTheta" ]->Sumw2();

	histos1D_[ "massAve_cutSubjetPtRatio" ] = fs_->make< TH1D >( "massAve_cutSubjetPtRatio", "massAve_cutSubjetPtRatio", 30, 0., 300. );
	histos1D_[ "massAve_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "massAveLowPU_cutSubjetPtRatio" ] = fs_->make< TH1D >( "massAveLowPU_cutSubjetPtRatio", "massAveLowPU_cutSubjetPtRatio", 30, 0., 300. );
	histos1D_[ "massAveLowPU_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "massAveMedPU_cutSubjetPtRatio" ] = fs_->make< TH1D >( "massAveMedPU_cutSubjetPtRatio", "massAveMedPU_cutSubjetPtRatio", 30, 0., 300. );
	histos1D_[ "massAveMedPU_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "massAveHighPU_cutSubjetPtRatio" ] = fs_->make< TH1D >( "massAveHighPU_cutSubjetPtRatio", "massAveHighPU_cutSubjetPtRatio", 30, 0., 300. );
	histos1D_[ "massAveHighPU_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "jet1Tau21_cutSubjetPtRatio" ] = fs_->make< TH1D >( "jet1Tau21_cutSubjetPtRatio", "jet1Tau21_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "jet1Tau21_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "jet1Tau31_cutSubjetPtRatio" ] = fs_->make< TH1D >( "jet1Tau31_cutSubjetPtRatio", "jet1Tau31_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "jet1Tau31_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "jet1Tau32_cutSubjetPtRatio" ] = fs_->make< TH1D >( "jet1Tau32_cutSubjetPtRatio", "jet1Tau32_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "jet1Tau32_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "subjetMass21Ratio_cutSubjetPtRatio" ] = fs_->make< TH1D >( "subjetMass21Ratio_cutSubjetPtRatio", "subjetMass21Ratio_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "subjetMass21Ratio_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "subjet112MassRatio_cutSubjetPtRatio" ] = fs_->make< TH1D >( "subjet112MassRatio_cutSubjetPtRatio", "subjet112MassRatio_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "subjet112MassRatio_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "subjet212MassRatio_cutSubjetPtRatio" ] = fs_->make< TH1D >( "subjet212MassRatio_cutSubjetPtRatio", "subjet212MassRatio_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "subjet212MassRatio_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "subjetPolAngle13412_cutSubjetPtRatio" ] = fs_->make< TH1D >( "subjetPolAngle13412_cutSubjetPtRatio", "subjetPolAngle13412_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "subjetPolAngle13412_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "subjetPolAngle31234_cutSubjetPtRatio" ] = fs_->make< TH1D >( "subjetPolAngle31234_cutSubjetPtRatio", "subjetPolAngle31234_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "subjetPolAngle31234_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle13412_cutSubjetPtRatio" ] = fs_->make< TH1D >( "tmpSubjetPolAngle13412_cutSubjetPtRatio", "tmpSubjetPolAngle13412_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle13412_cutSubjetPtRatio" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle31234_cutSubjetPtRatio" ] = fs_->make< TH1D >( "tmpSubjetPolAngle31234_cutSubjetPtRatio", "tmpSubjetPolAngle31234_cutSubjetPtRatio", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle31234_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "subjet12Mass_cutSubjetPtRatio" ] = fs_->make< TH2D >( "subjet12Mass_cutSubjetPtRatio", "subjet12Mass_cutSubjetPtRatio", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "subjet12Mass_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "dijetCorr_cutSubjetPtRatio" ] = fs_->make< TH2D >( "dijetCorr_cutSubjetPtRatio", "dijetCorr_cutSubjetPtRatio", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorr_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "dijetCorrPhi_cutSubjetPtRatio" ] = fs_->make< TH2D >( "dijetCorrPhi_cutSubjetPtRatio", "dijetCorrPhi_cutSubjetPtRatio", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorrPhi_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "subjet112vs212MassRatio_cutSubjetPtRatio" ] = fs_->make< TH2D >( "subjet112vs212MassRatio_cutSubjetPtRatio", "subjet112vs212MassRatio_cutSubjetPtRatio", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet112vs212MassRatio_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutSubjetPtRatio" ] = fs_->make< TH2D >( "subjet1JetvsSubjet2JetMassRatio_cutSubjetPtRatio", "subjet1JetvsSubjet2JetMassRatio_cutSubjetPtRatio", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "subjetPolAngle13412vs31234_cutSubjetPtRatio" ] = fs_->make< TH2D >( "subjetPolAngle13412vs31234_cutSubjetPtRatio", "subjetPolAngle13412vs31234_cutSubjetPtRatio", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjetPolAngle13412vs31234_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutSubjetPtRatio" ] = fs_->make< TH2D >( "tmpSubjetPolAngle13412vs31234_cutSubjetPtRatio", "tmpSubjetPolAngle13412vs31234_cutSubjetPtRatio", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "mu1234_cutSubjetPtRatio" ] = fs_->make< TH2D >( "mu1234_cutSubjetPtRatio", "mu1234_cutSubjetPtRatio", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu1234_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "mu3412_cutSubjetPtRatio" ] = fs_->make< TH2D >( "mu3412_cutSubjetPtRatio", "mu3412_cutSubjetPtRatio", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "mu3412_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "dalitz1234_cutSubjetPtRatio" ] = fs_->make< TH2D >( "dalitz1234_cutSubjetPtRatio", "dalitz1234_cutSubjetPtRatio", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz1234_cutSubjetPtRatio" ]->Sumw2();
	histos2D_[ "dalitz3412_cutSubjetPtRatio" ] = fs_->make< TH2D >( "dalitz3412_cutSubjetPtRatio", "dalitz3412_cutSubjetPtRatio", 150, 0., 1.5, 150, 0., 1.5 );
	histos2D_[ "dalitz3412_cutSubjetPtRatio" ]->Sumw2();

	histos1D_[ "massAve_cutTau31" ] = fs_->make< TH1D >( "massAve_cutTau31", "massAve_cutTau31", 30, 0., 300. );
	histos1D_[ "massAve_cutTau31" ]->Sumw2();
	histos1D_[ "subjetMass21Ratio_cutTau31" ] = fs_->make< TH1D >( "subjetMass21Ratio_cutTau31", "subjetMass21Ratio_cutTau31", 20, 0., 1. );
	histos1D_[ "subjetMass21Ratio_cutTau31" ]->Sumw2();
	histos1D_[ "subjet112MassRatio_cutTau31" ] = fs_->make< TH1D >( "subjet112MassRatio_cutTau31", "subjet112MassRatio_cutTau31", 20, 0., 1. );
	histos1D_[ "subjet112MassRatio_cutTau31" ]->Sumw2();
	histos1D_[ "subjet212MassRatio_cutTau31" ] = fs_->make< TH1D >( "subjet212MassRatio_cutTau31", "subjet212MassRatio_cutTau31", 20, 0., 1. );
	histos1D_[ "subjet212MassRatio_cutTau31" ]->Sumw2();
	histos1D_[ "subjetPolAngle13412_cutTau31" ] = fs_->make< TH1D >( "subjetPolAngle13412_cutTau31", "subjetPolAngle13412_cutTau31", 20, 0., 1. );
	histos1D_[ "subjetPolAngle13412_cutTau31" ]->Sumw2();
	histos1D_[ "subjetPolAngle31234_cutTau31" ] = fs_->make< TH1D >( "subjetPolAngle31234_cutTau31", "subjetPolAngle31234_cutTau31", 20, 0., 1. );
	histos1D_[ "subjetPolAngle31234_cutTau31" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle13412_cutTau31" ] = fs_->make< TH1D >( "tmpSubjetPolAngle13412_cutTau31", "tmpSubjetPolAngle13412_cutTau31", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle13412_cutTau31" ]->Sumw2();
	histos1D_[ "tmpSubjetPolAngle31234_cutTau31" ] = fs_->make< TH1D >( "tmpSubjetPolAngle31234_cutTau31", "tmpSubjetPolAngle31234_cutTau31", 20, 0., 1. );
	histos1D_[ "tmpSubjetPolAngle31234_cutTau31" ]->Sumw2();
	histos2D_[ "subjet12Mass_cutTau31" ] = fs_->make< TH2D >( "subjet12Mass_cutTau31", "subjet12Mass_cutTau31", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "subjet12Mass_cutTau31" ]->Sumw2();
	histos2D_[ "dijetCorr_cutTau31" ] = fs_->make< TH2D >( "dijetCorr_cutTau31", "dijetCorr_cutTau31", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorr_cutTau31" ]->Sumw2();
	histos2D_[ "dijetCorrPhi_cutTau31" ] = fs_->make< TH2D >( "dijetCorrPhi_cutTau31", "dijetCorrPhi_cutTau31", 14, -3.5, 3.5, 14, -3.5, 3.5 );
	histos2D_[ "dijetCorrPhi_cutTau31" ]->Sumw2();
	histos2D_[ "subjet112vs212MassRatio_cutTau31" ] = fs_->make< TH2D >( "subjet112vs212MassRatio_cutTau31", "subjet112vs212MassRatio_cutTau31", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet112vs212MassRatio_cutTau31" ]->Sumw2();
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutTau31" ] = fs_->make< TH2D >( "subjet1JetvsSubjet2JetMassRatio_cutTau31", "subjet1JetvsSubjet2JetMassRatio_cutTau31", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjet1JetvsSubjet2JetMassRatio_cutTau31" ]->Sumw2();
	histos2D_[ "subjetPolAngle13412vs31234_cutTau31" ] = fs_->make< TH2D >( "subjetPolAngle13412vs31234_cutTau31", "subjetPolAngle13412vs31234_cutTau31", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "subjetPolAngle13412vs31234_cutTau31" ]->Sumw2();
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutTau31" ] = fs_->make< TH2D >( "tmpSubjetPolAngle13412vs31234_cutTau31", "tmpSubjetPolAngle13412vs31234_cutTau31", 20, 0., 1., 20, 0., 1. );
	histos2D_[ "tmpSubjetPolAngle13412vs31234_cutTau31" ]->Sumw2();

	histos1D_[ "massAve_cutTau21" ] = fs_->make< TH1D >( "massAve_cutTau21", "massAve_cutTau21", 30, 0., 300. );
	histos1D_[ "massAve_cutTau21" ]->Sumw2();
	histos2D_[ "subjet12Mass_cutTau21" ] = fs_->make< TH2D >( "subjet12Mass_cutTau21", "subjet12Mass_cutTau21", 20, 0., 100., 20, 0., 100. );
	histos2D_[ "subjet12Mass_cutTau21" ]->Sumw2();

	histos1D_[ "massAve_cutBtagAfterSubjetPtRatio" ] = fs_->make< TH1D >( "massAve_cutBtagAfterSubjetPtRatio", "massAve_cutBtagAfterSubjetPtRatio", 30, 0., 300. );
	histos1D_[ "massAve_cutBtagAfterSubjetPtRatio" ]->Sumw2();
	histos1D_[ "massAve_cutBtagAfterTau31" ] = fs_->make< TH1D >( "massAve_cutBtagAfterTau31", "massAve_cutBtagAfterTau31", 30, 0., 300. );
	histos1D_[ "massAve_cutBtagAfterTau31" ]->Sumw2();
	histos1D_[ "massAve_cutBtagAfterTau21" ] = fs_->make< TH1D >( "massAve_cutBtagAfterTau21", "massAve_cutBtagAfterTau21", 30, 0., 300. );
	histos1D_[ "massAve_cutBtagAfterTau21" ]->Sumw2();

	cutLabels.push_back("Processed");
	cutLabels.push_back("Trigger");
	cutLabels.push_back("HT");
	cutLabels.push_back("Asymmetry");
	cutLabels.push_back("CosTheta");
	cutLabels.push_back("SubjetPtRatio");
	cutLabels.push_back("btagAfterSubjetPtRatio");
	cutLabels.push_back("Tau31");
	cutLabels.push_back("btagAfterTau31");
	cutLabels.push_back("Tau21");
	cutLabels.push_back("btagAfterTau21");
	histos1D_[ "hcutflow" ] = fs_->make< TH1D >("cutflow","cut flow", cutLabels.size(), 0.5, cutLabels.size() +0.5 );
	histos1D_[ "hcutflow" ]->Sumw2();
	histos1D_[ "hcutflowSimple" ] = fs_->make< TH1D >("cutflowSimple","simple cut flow", cutLabels.size(), 0.5, cutLabels.size() +0.5 );
	histos1D_[ "hcutflowSimple" ]->Sumw2();
	for( const string &ivec : cutLabels ) cutmap[ ivec ] = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void RUNONBoostedAnalysis::endJob() {

	int ibin = 1;
	for( const string &ivec : cutLabels ) {
		histos1D_["hcutflow"]->SetBinContent( ibin, cutmap[ ivec ] * scale );
		histos1D_["hcutflow"]->GetXaxis()->SetBinLabel( ibin, ivec.c_str() );
		histos1D_["hcutflowSimple"]->SetBinContent( ibin, cutmap[ ivec ] );
		histos1D_["hcutflowSimple"]->GetXaxis()->SetBinLabel( ibin, ivec.c_str() );
		ibin++;
	}

}

void RUNONBoostedAnalysis::clearVariables() {

	jetsPt->clear();
	jetsEta->clear();
	jetsPhi->clear();
	jetsE->clear();
	jet1SubjetsPt->clear();
	jet1SubjetsEta->clear();
	jet1SubjetsPhi->clear();
	jet1SubjetsE->clear();
	jet2SubjetsPt->clear();
	jet2SubjetsEta->clear();
	jet2SubjetsPhi->clear();
	jet2SubjetsE->clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(RUNONBoostedAnalysis);
