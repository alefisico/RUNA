// -*- C++ -*-
//
// Package:    RUNA/RUNTriggerEfficiency
// Class:      RUNBoostedTriggerEfficiency
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

#include "RUNA/RUNAnalysis/interface/CommonVariablesStructure.h"

using namespace edm;
using namespace std;

//
// constants, enums and typedefs
//

//
// class declaration
//
class RUNBoostedTriggerEfficiency : public EDAnalyzer {
	public:
		explicit RUNBoostedTriggerEfficiency(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~RUNBoostedTriggerEfficiency();

	private:
		virtual void beginJob() override;
		virtual void analyze(const Event&, const EventSetup&) override;
		virtual void endJob() override;

		virtual void beginRun(Run const&, EventSetup const&) override;
		virtual void endRun(Run const&, EventSetup const&) override;
		//virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
		//virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

		// ----------member data ---------------------------
		Service<TFileService> fs_;
		TTree *RUNAtree;
		map< string, TH1D* > histos1D_;
		map< string, TH2D* > histos2D_;
		vector< string > cutLabels;

		double cutAK8jetPt;
		double cutAK8HT;
		double cutAK8jet1Pt;
		double cutAK8jet2Pt;
		double cutAK8jet1Mass;
		string jecVersion;
		string PUMethod;
		TString baseTrigger;
		vector<string> triggerPass, triggerNamesList;
		vector<JetCorrectorParameters> jetPar;
		FactorizedJetCorrector * jetJECAK8;
		vector<JetCorrectorParameters> massPar;
		FactorizedJetCorrector * massJECAK8;

		ULong64_t event = 0;
		unsigned int lumi = 0, run=0;

		EDGetTokenT<vector<float>> jetPt_;
		EDGetTokenT<vector<float>> jetEta_;
		EDGetTokenT<vector<float>> jetPhi_;
		EDGetTokenT<vector<float>> jetE_;
		EDGetTokenT<vector<float>> jetTrimmedMass_;
		EDGetTokenT<vector<float>> jetPrunedMass_;
		EDGetTokenT<vector<float>> jetSoftDropMass_;
		EDGetTokenT<vector<float>> jetCSVv2_;
		EDGetTokenT<vector<float>> jetCSVv2V1_;
		EDGetTokenT<vector<float>> jetArea_;
		EDGetTokenT<int> NPV_;
		EDGetTokenT<vector<float>> rho_;
		EDGetTokenT<unsigned int> lumi_;
		EDGetTokenT<unsigned int> run_;
		EDGetTokenT<ULong64_t> event_;

		// Trigger
		EDGetTokenT<vector<int>> triggerPrescale_;
		EDGetTokenT<vector<float>> triggerBit_;
		EDGetTokenT<vector<string>> triggerName_;

		//Jet ID
		EDGetTokenT<vector<float>> jecFactor_;
		EDGetTokenT<vector<float>> neutralHadronEnergyFrac_;
		EDGetTokenT<vector<float>> neutralEmEnergyFrac_;
		EDGetTokenT<vector<float>> chargedHadronEnergyFrac_;
		EDGetTokenT<vector<float>> chargedEmEnergyFrac_;
		EDGetTokenT<vector<float>> neutralMultiplicity_;
		EDGetTokenT<vector<float>> chargedMultiplicity_;
		EDGetTokenT<vector<float>> muonEnergy_; 

};

//
// static data member definitions
//

//
// constructors and destructor
//
RUNBoostedTriggerEfficiency::RUNBoostedTriggerEfficiency(const ParameterSet& iConfig):
	jetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPt"))),
	jetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetEta"))),
	jetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPhi"))),
	jetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetE"))),
	jetTrimmedMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTrimmedMass"))),
	jetPrunedMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPrunedMass"))),
	jetSoftDropMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSoftDropMass"))),
	jetCSVv2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCSVv2"))),
	jetArea_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetArea"))),
	NPV_(consumes<int>(iConfig.getParameter<InputTag>("NPV"))),
	rho_(consumes<vector<float>>(iConfig.getParameter<InputTag>("rho"))),
	lumi_(consumes<unsigned int>(iConfig.getParameter<InputTag>("Lumi"))),
	run_(consumes<unsigned int>(iConfig.getParameter<InputTag>("Run"))),
	event_(consumes<ULong64_t>(iConfig.getParameter<InputTag>("Event"))),
	// Trigger
	triggerPrescale_(consumes<vector<int>>(iConfig.getParameter<InputTag>("triggerPrescale"))),
	triggerBit_(consumes<vector<float>>(iConfig.getParameter<InputTag>("triggerBit"))),
	triggerName_(consumes<vector<string>,InRun>(iConfig.getParameter<InputTag>("triggerName"))),
	//Jet ID,
	jecFactor_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jecFactor"))),
	neutralHadronEnergyFrac_(consumes<vector<float>>(iConfig.getParameter<InputTag>("neutralHadronEnergyFrac"))),
	neutralEmEnergyFrac_(consumes<vector<float>>(iConfig.getParameter<InputTag>("neutralEmEnergyFrac"))),
	chargedHadronEnergyFrac_(consumes<vector<float>>(iConfig.getParameter<InputTag>("chargedHadronEnergyFrac"))),
	chargedEmEnergyFrac_(consumes<vector<float>>(iConfig.getParameter<InputTag>("chargedEmEnergyFrac"))),
	neutralMultiplicity_(consumes<vector<float>>(iConfig.getParameter<InputTag>("neutralMultiplicity"))),
	chargedMultiplicity_(consumes<vector<float>>(iConfig.getParameter<InputTag>("chargedMultiplicity"))),
	muonEnergy_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonEnergy")))
{
	cutAK8jetPt 	= iConfig.getParameter<double>("cutAK8jetPt");
	cutAK8HT 	= iConfig.getParameter<double>("cutAK8HT");
	cutAK8jet1Pt 	= iConfig.getParameter<double>("cutAK8jet1Pt");
	cutAK8jet2Pt 	= iConfig.getParameter<double>("cutAK8jet2Pt");
	cutAK8jet1Mass 	= iConfig.getParameter<double>("cutAK8jet1Mass");
	baseTrigger 	= iConfig.getParameter<string>("baseTrigger");
	triggerPass 	= iConfig.getParameter<vector<string>>("triggerPass");
	jecVersion 	= iConfig.getParameter<string>("jecVersion");
	PUMethod 	= iConfig.getParameter<string>("PUMethod");

	/////// JECs
	string prefix;
	prefix = jecVersion + "_DATA_";

	// all jet
	vector<string> jecAK8PayloadNames_;
	jecAK8PayloadNames_.push_back(prefix + "L1FastJet_AK8PF"+PUMethod+".txt");
	jecAK8PayloadNames_.push_back(prefix + "L2Relative_AK8PF"+PUMethod+".txt");
	jecAK8PayloadNames_.push_back(prefix + "L3Absolute_AK8PF"+PUMethod+".txt");
	jecAK8PayloadNames_.push_back(prefix + "L2L3Residual_AK8PF"+PUMethod+".txt");

	for ( vector<string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		jetPar.push_back(pars);
	}
	jetJECAK8 = new FactorizedJetCorrector(jetPar);

	// jet mass
	vector<string> massjecAK8PayloadNames_;
	massjecAK8PayloadNames_.push_back(prefix + "L2Relative_AK8PF"+PUMethod+".txt");
	massjecAK8PayloadNames_.push_back(prefix + "L3Absolute_AK8PF"+PUMethod+".txt");
	massjecAK8PayloadNames_.push_back(prefix + "L2L3Residual_AK8PF"+PUMethod+".txt");

	for ( vector<string>::const_iterator payloadBegin = massjecAK8PayloadNames_.begin(), payloadEnd = massjecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters massPars(*ipayload);
		massPar.push_back(massPars);
	}
	massJECAK8 = new FactorizedJetCorrector(massPar);
}


RUNBoostedTriggerEfficiency::~RUNBoostedTriggerEfficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNBoostedTriggerEfficiency::analyze(const Event& iEvent, const EventSetup& iSetup) {

	Handle<vector<float> > jetPt;
	iEvent.getByToken(jetPt_, jetPt);

	Handle<vector<float> > jetEta;
	iEvent.getByToken(jetEta_, jetEta);

	Handle<vector<float> > jetPhi;
	iEvent.getByToken(jetPhi_, jetPhi);

	Handle<vector<float> > jetE;
	iEvent.getByToken(jetE_, jetE);

	Handle<vector<float> > jetTrimmedMass;
	iEvent.getByToken(jetTrimmedMass_, jetTrimmedMass);

	Handle<vector<float> > jetPrunedMass;
	iEvent.getByToken(jetPrunedMass_, jetPrunedMass);

	Handle<vector<float> > jetSoftDropMass;
	iEvent.getByToken(jetSoftDropMass_, jetSoftDropMass);

	Handle<vector<float> > jetCSVv2;
	iEvent.getByToken(jetCSVv2_, jetCSVv2);

	Handle<vector<float> > jetArea;
	iEvent.getByToken(jetArea_, jetArea);

	Handle<int> NPV;
	iEvent.getByToken(NPV_, NPV);

	Handle<vector<float>> rho;
	iEvent.getByToken(rho_, rho);

	Handle<unsigned int> Lumi;
	iEvent.getByToken(lumi_, Lumi);

	Handle<unsigned int> Run;
	iEvent.getByToken(run_, Run);

	Handle<ULong64_t> ievent;
	iEvent.getByToken(event_, ievent);

	/// Trigger
	Handle<vector<int> > triggerPrescale;
	iEvent.getByToken(triggerPrescale_, triggerPrescale);

	Handle<vector<float> > triggerBit;
	iEvent.getByToken(triggerBit_, triggerBit);

	/// Jet ID
	Handle<vector<float> > jecFactor;
	iEvent.getByToken(jecFactor_, jecFactor);

	Handle<vector<float> > neutralHadronEnergyFrac;
	iEvent.getByToken(neutralHadronEnergyFrac_, neutralHadronEnergyFrac);

	Handle<vector<float> > neutralEmEnergyFrac;
	iEvent.getByToken(neutralEmEnergyFrac_, neutralEmEnergyFrac);

	Handle<vector<float> > chargedHadronEnergyFrac;
	iEvent.getByToken(chargedHadronEnergyFrac_, chargedHadronEnergyFrac);

	Handle<vector<float> > chargedEmEnergyFrac;
	iEvent.getByToken(chargedEmEnergyFrac_, chargedEmEnergyFrac);

	Handle<vector<float> > neutralMultiplicity;
	iEvent.getByToken(neutralMultiplicity_, neutralMultiplicity);

	Handle<vector<float> > chargedMultiplicity;
	iEvent.getByToken(chargedMultiplicity_, chargedMultiplicity);

	Handle<vector<float> > muonEnergy;
	iEvent.getByToken(muonEnergy_, muonEnergy);

	bool basedTriggerFired = checkTriggerBits( triggerNamesList, triggerBit, triggerPrescale, baseTrigger, true  );
	bool ORTriggers = checkORListOfTriggerBits( triggerNamesList, triggerBit, triggerPrescale, triggerPass, false );

	/// Applying kinematic, trigger and jet ID
	vector< myJet > JETS;
	float HT = 0;

	for (size_t i = 0; i < jetPt->size(); i++) {

		if( TMath::Abs( (*jetEta)[i] ) > 2.4 ) continue;
		string typeOfJetID = "looseJetID";	// check trigger with looser jet id
		bool jetId = jetID( (*jetEta)[i], (*jetE)[i], (*jecFactor)[i], (*neutralHadronEnergyFrac)[i], (*neutralEmEnergyFrac)[i], (*chargedHadronEnergyFrac)[i], (*muonEnergy)[i], (*chargedEmEnergyFrac)[i], (*chargedMultiplicity)[i], (*neutralMultiplicity)[i], typeOfJetID ); 

		TLorentzVector tmpJet, rawJet, corrJet, genJet, smearJet;
		tmpJet.SetPtEtaPhiE( (*jetPt)[i], (*jetEta)[i], (*jetPhi)[i], (*jetE)[i] );
		rawJet = tmpJet* (*jecFactor)[i] ;

		double JEC = corrections( rawJet, (*jetArea)[i], (*rho)[i] ,*NPV, jetJECAK8); 
		corrJet = rawJet* ( JEC ); 

		//if( ( (*jetPt)[i] > cutAK8jetPt ) && idL ) { 
		if( ( corrJet.Pt() > cutAK8jetPt ) && jetId ) { 

			//HT += (*jetPt)[i];
			HT += corrJet.Pt();

			//TLorentzVector tmpJet;
			//tmpJet.SetPtEtaPhiE( (*jetPt)[i], (*jetEta)[i], (*jetPhi)[i], (*jetE)[i] );
			double massJEC = corrections( rawJet, (*jetArea)[i], (*rho)[i] ,*NPV, massJECAK8); 
			double corrTrimmedMass = (*jetTrimmedMass)[i] * massJEC ;
			double corrPrunedMass = (*jetPrunedMass)[i] * massJEC ;
			double corrSoftDropMass = (*jetSoftDropMass)[i] * massJEC ;


			/// Calculating RHO
			double rhoMass = TMath::Log( ( corrJet.M() * corrJet.M() )/( corrJet.Pt() * corrJet.Pt() ) );
			double rhoPrunedMass = TMath::Log( ( corrPrunedMass * corrPrunedMass )/( corrJet.Pt() * corrJet.Pt() ) );

			myJet tmpJET;
			tmpJET.p4 = corrJet; //tmpJet;
			tmpJET.trimmedMass = corrTrimmedMass;
			tmpJET.prunedMass = corrPrunedMass;
			tmpJET.softDropMass = corrSoftDropMass;
			tmpJET.btagCSVv2 = (*jetCSVv2)[i];
			tmpJET.tau1 = rhoMass;		// Just for convenience, rho stored in tau1
			tmpJET.tau2 = rhoPrunedMass;
			JETS.push_back( tmpJET );
		}
	}


	if ( JETS.size() > 1 ) {

		// Mass average 
		float prunedMassAve = massAverage( JETS[0].prunedMass, JETS[1].prunedMass );
		float softDropMassAve = massAverage( JETS[0].softDropMass, JETS[1].softDropMass );
		//////////////////////////////////////////////////////////////////////////

		
		histos2D_[ "jet1PrunedMassHT_cutDijet_noTrigger" ]->Fill( JETS[0].prunedMass, HT );
		histos2D_[ "jet1SoftDropMassHT_cutDijet_noTrigger" ]->Fill( JETS[0].softDropMass, HT );

		if ( basedTriggerFired || ORTriggers ) {
			if ( basedTriggerFired && ORTriggers ) {
				histos2D_[ "jet1PrunedMassHT_cutDijet_triggerOneAndTwo" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerOneAndTwo" ]->Fill( JETS[0].softDropMass, HT );
			} else if ( basedTriggerFired ) {
				histos2D_[ "jet1PrunedMassHT_cutDijet_triggerOne" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerOne" ]->Fill( JETS[0].softDropMass, HT );
			} else if ( ORTriggers ) {
				histos2D_[ "jet1PrunedMassHT_cutDijet_triggerTwo" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerTwo" ]->Fill( JETS[0].softDropMass, HT );
			}
		} else {
			histos2D_[ "jet1PrunedMassHT_cutDijet_noTriggerOneOrTwo" ]->Fill( JETS[0].prunedMass, HT );
			histos2D_[ "jet1SoftDropMassHT_cutDijet_noTriggerOneOrTwo" ]->Fill( JETS[0].softDropMass, HT );
		}


		if ( basedTriggerFired ) {
			histos1D_[ "HTDenom_cutDijet" ]->Fill( HT  );
			histos1D_[ "prunedMassAveDenom_cutDijet" ]->Fill( prunedMassAve );
			histos1D_[ "softDropMassAveDenom_cutDijet" ]->Fill( softDropMassAve );
			histos1D_[ "jet1PrunedMassDenom_cutDijet" ]->Fill( JETS[0].prunedMass  );
			histos1D_[ "jet1SoftDropMassDenom_cutDijet" ]->Fill( JETS[0].softDropMass  );
			histos1D_[ "jet1PtDenom_cutDijet" ]->Fill( JETS[0].p4.Pt()   );
			histos1D_[ "jet1RhoDenom_cutDijet" ]->Fill( JETS[0].tau1 );	// tau1 is rho!!!!!
			histos1D_[ "jet1RhoPrunedDenom_cutDijet" ]->Fill( JETS[0].tau2 );
			histos1D_[ "jet2PrunedMassDenom_cutDijet" ]->Fill( JETS[1].prunedMass  );
			histos1D_[ "jet2SoftDropMassDenom_cutDijet" ]->Fill( JETS[1].softDropMass  );
			histos1D_[ "jet2PtDenom_cutDijet" ]->Fill( JETS[1].p4.Pt() );
			histos1D_[ "jet2RhoDenom_cutDijet" ]->Fill( JETS[1].tau1 );
			histos1D_[ "jet2RhoPrunedDenom_cutDijet" ]->Fill( JETS[1].tau2 );

			histos2D_[ "jet1PrunedMassHTDenom_cutDijet" ]->Fill( JETS[0].prunedMass, HT );
			histos2D_[ "jet1PrunedMassjet1PtDenom_cutDijet" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
			histos2D_[ "jet1SoftDropMassHTDenom_cutDijet" ]->Fill( JETS[0].softDropMass, HT );
			histos2D_[ "jet1SoftDropMassjet1PtDenom_cutDijet" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

			histos2D_[ "jet2PrunedMassHTDenom_cutDijet" ]->Fill( JETS[1].prunedMass, HT );
			histos2D_[ "jet2PrunedMassjet2PtDenom_cutDijet" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
			histos2D_[ "jet2SoftDropMassHTDenom_cutDijet" ]->Fill( JETS[1].softDropMass, HT );
			histos2D_[ "jet2SoftDropMassjet2PtDenom_cutDijet" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

			histos2D_[ "prunedMassAveHTDenom_cutDijet" ]->Fill( prunedMassAve, HT );
			histos2D_[ "prunedMassAvejet1PtDenom_cutDijet" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
			histos2D_[ "prunedMassAvejet2PtDenom_cutDijet" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
			histos2D_[ "softDropMassAveHTDenom_cutDijet" ]->Fill( softDropMassAve, HT );
			histos2D_[ "softDropMassAvejet1PtDenom_cutDijet" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
			histos2D_[ "softDropMassAvejet2PtDenom_cutDijet" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );

			if ( ORTriggers ){
				histos1D_[ "HTPassing_cutDijet" ]->Fill( HT  );
				histos1D_[ "prunedMassAvePassing_cutDijet" ]->Fill( prunedMassAve );
				histos1D_[ "softDropMassAvePassing_cutDijet" ]->Fill( softDropMassAve );
				histos1D_[ "jet1PrunedMassPassing_cutDijet" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassPassing_cutDijet" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtPassing_cutDijet" ]->Fill( JETS[0].p4.Pt()   );
				histos1D_[ "jet1RhoPassing_cutDijet" ]->Fill( JETS[0].tau1 );	// tau1 is rho!!!!!
				histos1D_[ "jet1RhoPrunedPassing_cutDijet" ]->Fill( JETS[0].tau2 );
				histos1D_[ "jet2PrunedMassPassing_cutDijet" ]->Fill( JETS[1].prunedMass  );
				histos1D_[ "jet2SoftDropMassPassing_cutDijet" ]->Fill( JETS[1].softDropMass  );
				histos1D_[ "jet2PtPassing_cutDijet" ]->Fill( JETS[1].p4.Pt() );
				histos1D_[ "jet2RhoPassing_cutDijet" ]->Fill( JETS[1].tau1 );
				histos1D_[ "jet2RhoPrunedPassing_cutDijet" ]->Fill( JETS[1].tau2 );

				histos2D_[ "jet1PrunedMassHTPassing_cutDijet" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1PrunedMassjet1PtPassing_cutDijet" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassHTPassing_cutDijet" ]->Fill( JETS[0].softDropMass, HT );
				histos2D_[ "jet1SoftDropMassjet1PtPassing_cutDijet" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

				histos2D_[ "jet2PrunedMassHTPassing_cutDijet" ]->Fill( JETS[1].prunedMass, HT );
				histos2D_[ "jet2PrunedMassjet2PtPassing_cutDijet" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
				histos2D_[ "jet2SoftDropMassHTPassing_cutDijet" ]->Fill( JETS[1].softDropMass, HT );
				histos2D_[ "jet2SoftDropMassjet2PtPassing_cutDijet" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

				histos2D_[ "prunedMassAveHTPassing_cutDijet" ]->Fill( prunedMassAve, HT );
				histos2D_[ "prunedMassAvejet1PtPassing_cutDijet" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
				histos2D_[ "prunedMassAvejet2PtPassing_cutDijet" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
				histos2D_[ "softDropMassAveHTPassing_cutDijet" ]->Fill( softDropMassAve, HT );
				histos2D_[ "softDropMassAvejet1PtPassing_cutDijet" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
				histos2D_[ "softDropMassAvejet2PtPassing_cutDijet" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );
			}

			if ( HT > cutAK8HT ) {
				histos1D_[ "HTDenom_cutHT" ]->Fill( HT  );
				histos1D_[ "prunedMassAveDenom_cutHT" ]->Fill( prunedMassAve );
				histos1D_[ "softDropMassAveDenom_cutHT" ]->Fill( softDropMassAve );
				histos1D_[ "jet1PrunedMassDenom_cutHT" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassDenom_cutHT" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtDenom_cutHT" ]->Fill( JETS[0].p4.Pt()   );
				histos1D_[ "jet1RhoDenom_cutHT" ]->Fill( JETS[0].tau1 );	// tau1 is rho!!!!!
				histos1D_[ "jet1RhoPrunedDenom_cutHT" ]->Fill( JETS[0].tau2 );
				histos1D_[ "jet2PrunedMassDenom_cutHT" ]->Fill( JETS[1].prunedMass  );
				histos1D_[ "jet2SoftDropMassDenom_cutHT" ]->Fill( JETS[1].softDropMass  );
				histos1D_[ "jet2PtDenom_cutHT" ]->Fill( JETS[1].p4.Pt() );
				histos1D_[ "jet2RhoDenom_cutHT" ]->Fill( JETS[1].tau1 );
				histos1D_[ "jet2RhoPrunedDenom_cutHT" ]->Fill( JETS[1].tau2 );

				histos2D_[ "jet1PrunedMassHTDenom_cutHT" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1PrunedMassjet1PtDenom_cutHT" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassHTDenom_cutHT" ]->Fill( JETS[0].softDropMass, HT );
				histos2D_[ "jet1SoftDropMassjet1PtDenom_cutHT" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

				histos2D_[ "jet2PrunedMassHTDenom_cutHT" ]->Fill( JETS[1].prunedMass, HT );
				histos2D_[ "jet2PrunedMassjet2PtDenom_cutHT" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
				histos2D_[ "jet2SoftDropMassHTDenom_cutHT" ]->Fill( JETS[1].softDropMass, HT );
				histos2D_[ "jet2SoftDropMassjet2PtDenom_cutHT" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

				histos2D_[ "prunedMassAveHTDenom_cutHT" ]->Fill( prunedMassAve, HT );
				histos2D_[ "prunedMassAvejet1PtDenom_cutHT" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
				histos2D_[ "prunedMassAvejet2PtDenom_cutHT" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
				histos2D_[ "softDropMassAveHTDenom_cutHT" ]->Fill( softDropMassAve, HT );
				histos2D_[ "softDropMassAvejet1PtDenom_cutHT" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
				histos2D_[ "softDropMassAvejet2PtDenom_cutHT" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );

				if ( ORTriggers ){
					histos1D_[ "HTPassing_cutHT" ]->Fill( HT  );
					histos1D_[ "prunedMassAvePassing_cutHT" ]->Fill( prunedMassAve );
					histos1D_[ "softDropMassAvePassing_cutHT" ]->Fill( softDropMassAve );
					histos1D_[ "jet1PrunedMassPassing_cutHT" ]->Fill( JETS[0].prunedMass  );
					histos1D_[ "jet1SoftDropMassPassing_cutHT" ]->Fill( JETS[0].softDropMass  );
					histos1D_[ "jet1PtPassing_cutHT" ]->Fill( JETS[0].p4.Pt()   );
					histos1D_[ "jet1RhoPassing_cutHT" ]->Fill( JETS[0].tau1 );	// tau1 is rho!!!!!
					histos1D_[ "jet1RhoPrunedPassing_cutHT" ]->Fill( JETS[0].tau2 );
					histos1D_[ "jet2PrunedMassPassing_cutHT" ]->Fill( JETS[1].prunedMass  );
					histos1D_[ "jet2SoftDropMassPassing_cutHT" ]->Fill( JETS[1].softDropMass  );
					histos1D_[ "jet2PtPassing_cutHT" ]->Fill( JETS[1].p4.Pt() );
					histos1D_[ "jet2RhoPassing_cutHT" ]->Fill( JETS[1].tau1 );
					histos1D_[ "jet2RhoPrunedPassing_cutHT" ]->Fill( JETS[1].tau2 );

					histos2D_[ "jet1PrunedMassHTPassing_cutHT" ]->Fill( JETS[0].prunedMass, HT );
					histos2D_[ "jet1PrunedMassjet1PtPassing_cutHT" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
					histos2D_[ "jet1SoftDropMassHTPassing_cutHT" ]->Fill( JETS[0].softDropMass, HT );
					histos2D_[ "jet1SoftDropMassjet1PtPassing_cutHT" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

					histos2D_[ "jet2PrunedMassHTPassing_cutHT" ]->Fill( JETS[1].prunedMass, HT );
					histos2D_[ "jet2PrunedMassjet2PtPassing_cutHT" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
					histos2D_[ "jet2SoftDropMassHTPassing_cutHT" ]->Fill( JETS[1].softDropMass, HT );
					histos2D_[ "jet2SoftDropMassjet2PtPassing_cutHT" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

					histos2D_[ "prunedMassAveHTPassing_cutHT" ]->Fill( prunedMassAve, HT );
					histos2D_[ "prunedMassAvejet1PtPassing_cutHT" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
					histos2D_[ "prunedMassAvejet2PtPassing_cutHT" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
					histos2D_[ "softDropMassAveHTPassing_cutHT" ]->Fill( softDropMassAve, HT );
					histos2D_[ "softDropMassAvejet1PtPassing_cutHT" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
					histos2D_[ "softDropMassAvejet2PtPassing_cutHT" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );
				}
			}


			if ( JETS[0].prunedMass > cutAK8jet1Mass ) {
				histos1D_[ "HTDenom_cutjet1PrunedMass" ]->Fill( HT  );
				histos1D_[ "prunedMassAveDenom_cutjet1PrunedMass" ]->Fill( prunedMassAve );
				histos1D_[ "softDropMassAveDenom_cutjet1PrunedMass" ]->Fill( softDropMassAve );
				histos1D_[ "jet1PrunedMassDenom_cutjet1PrunedMass" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassDenom_cutjet1PrunedMass" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtDenom_cutjet1PrunedMass" ]->Fill( JETS[0].p4.Pt()   );
				histos1D_[ "jet2PrunedMassDenom_cutjet1PrunedMass" ]->Fill( JETS[1].prunedMass  );
				histos1D_[ "jet2SoftDropMassDenom_cutjet1PrunedMass" ]->Fill( JETS[1].softDropMass  );
				histos1D_[ "jet2PtDenom_cutjet1PrunedMass" ]->Fill( JETS[1].p4.Pt() );

				histos2D_[ "jet1PrunedMassHTDenom_cutjet1PrunedMass" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1PrunedMassjet1PtDenom_cutjet1PrunedMass" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassHTDenom_cutjet1PrunedMass" ]->Fill( JETS[0].softDropMass, HT );
				histos2D_[ "jet1SoftDropMassjet1PtDenom_cutjet1PrunedMass" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

				histos2D_[ "jet2PrunedMassHTDenom_cutjet1PrunedMass" ]->Fill( JETS[1].prunedMass, HT );
				histos2D_[ "jet2PrunedMassjet2PtDenom_cutjet1PrunedMass" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
				histos2D_[ "jet2SoftDropMassHTDenom_cutjet1PrunedMass" ]->Fill( JETS[1].softDropMass, HT );
				histos2D_[ "jet2SoftDropMassjet2PtDenom_cutjet1PrunedMass" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

				histos2D_[ "prunedMassAveHTDenom_cutjet1PrunedMass" ]->Fill( prunedMassAve, HT );
				histos2D_[ "prunedMassAvejet1PtDenom_cutjet1PrunedMass" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
				histos2D_[ "prunedMassAvejet2PtDenom_cutjet1PrunedMass" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
				histos2D_[ "softDropMassAveHTDenom_cutjet1PrunedMass" ]->Fill( softDropMassAve, HT );
				histos2D_[ "softDropMassAvejet1PtDenom_cutjet1PrunedMass" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
				histos2D_[ "softDropMassAvejet2PtDenom_cutjet1PrunedMass" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );

				if ( ORTriggers ){
					histos1D_[ "HTPassing_cutjet1PrunedMass" ]->Fill( HT  );
					histos1D_[ "prunedMassAvePassing_cutjet1PrunedMass" ]->Fill( prunedMassAve );
					histos1D_[ "softDropMassAvePassing_cutjet1PrunedMass" ]->Fill( softDropMassAve );
					histos1D_[ "jet1PrunedMassPassing_cutjet1PrunedMass" ]->Fill( JETS[0].prunedMass  );
					histos1D_[ "jet1SoftDropMassPassing_cutjet1PrunedMass" ]->Fill( JETS[0].softDropMass  );
					histos1D_[ "jet1PtPassing_cutjet1PrunedMass" ]->Fill( JETS[0].p4.Pt()   );
					histos1D_[ "jet2PrunedMassPassing_cutjet1PrunedMass" ]->Fill( JETS[1].prunedMass  );
					histos1D_[ "jet2SoftDropMassPassing_cutjet1PrunedMass" ]->Fill( JETS[1].softDropMass  );
					histos1D_[ "jet2PtPassing_cutjet1PrunedMass" ]->Fill( JETS[1].p4.Pt() );

					histos2D_[ "jet1PrunedMassHTPassing_cutjet1PrunedMass" ]->Fill( JETS[0].prunedMass, HT );
					histos2D_[ "jet1PrunedMassjet1PtPassing_cutjet1PrunedMass" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
					histos2D_[ "jet1SoftDropMassHTPassing_cutjet1PrunedMass" ]->Fill( JETS[0].softDropMass, HT );
					histos2D_[ "jet1SoftDropMassjet1PtPassing_cutjet1PrunedMass" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

					histos2D_[ "jet2PrunedMassHTPassing_cutjet1PrunedMass" ]->Fill( JETS[1].prunedMass, HT );
					histos2D_[ "jet2PrunedMassjet2PtPassing_cutjet1PrunedMass" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
					histos2D_[ "jet2SoftDropMassHTPassing_cutjet1PrunedMass" ]->Fill( JETS[1].softDropMass, HT );
					histos2D_[ "jet2SoftDropMassjet2PtPassing_cutjet1PrunedMass" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

					histos2D_[ "prunedMassAveHTPassing_cutjet1PrunedMass" ]->Fill( prunedMassAve, HT );
					histos2D_[ "prunedMassAvejet1PtPassing_cutjet1PrunedMass" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
					histos2D_[ "prunedMassAvejet2PtPassing_cutjet1PrunedMass" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
					histos2D_[ "softDropMassAveHTPassing_cutjet1PrunedMass" ]->Fill( softDropMassAve, HT );
					histos2D_[ "softDropMassAvejet1PtPassing_cutjet1PrunedMass" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
					histos2D_[ "softDropMassAvejet2PtPassing_cutjet1PrunedMass" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );
				}
			}

			if ( JETS[0].softDropMass > cutAK8jet1Mass ) {
				histos1D_[ "HTDenom_cutjet1SoftDropMass" ]->Fill( HT  );
				histos1D_[ "prunedMassAveDenom_cutjet1SoftDropMass" ]->Fill( prunedMassAve );
				histos1D_[ "softDropMassAveDenom_cutjet1SoftDropMass" ]->Fill( softDropMassAve );
				histos1D_[ "jet1PrunedMassDenom_cutjet1SoftDropMass" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassDenom_cutjet1SoftDropMass" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtDenom_cutjet1SoftDropMass" ]->Fill( JETS[0].p4.Pt()   );
				histos1D_[ "jet2PrunedMassDenom_cutjet1SoftDropMass" ]->Fill( JETS[1].prunedMass  );
				histos1D_[ "jet2SoftDropMassDenom_cutjet1SoftDropMass" ]->Fill( JETS[1].softDropMass  );
				histos1D_[ "jet2PtDenom_cutjet1SoftDropMass" ]->Fill( JETS[1].p4.Pt() );

				histos2D_[ "jet1PrunedMassHTDenom_cutjet1SoftDropMass" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1PrunedMassjet1PtDenom_cutjet1SoftDropMass" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassHTDenom_cutjet1SoftDropMass" ]->Fill( JETS[0].softDropMass, HT );
				histos2D_[ "jet1SoftDropMassjet1PtDenom_cutjet1SoftDropMass" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

				histos2D_[ "jet2PrunedMassHTDenom_cutjet1SoftDropMass" ]->Fill( JETS[1].prunedMass, HT );
				histos2D_[ "jet2PrunedMassjet2PtDenom_cutjet1SoftDropMass" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
				histos2D_[ "jet2SoftDropMassHTDenom_cutjet1SoftDropMass" ]->Fill( JETS[1].softDropMass, HT );
				histos2D_[ "jet2SoftDropMassjet2PtDenom_cutjet1SoftDropMass" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

				histos2D_[ "prunedMassAveHTDenom_cutjet1SoftDropMass" ]->Fill( prunedMassAve, HT );
				histos2D_[ "prunedMassAvejet1PtDenom_cutjet1SoftDropMass" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
				histos2D_[ "prunedMassAvejet2PtDenom_cutjet1SoftDropMass" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
				histos2D_[ "softDropMassAveHTDenom_cutjet1SoftDropMass" ]->Fill( softDropMassAve, HT );
				histos2D_[ "softDropMassAvejet1PtDenom_cutjet1SoftDropMass" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
				histos2D_[ "softDropMassAvejet2PtDenom_cutjet1SoftDropMass" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );

				if ( ORTriggers ){
					histos1D_[ "HTPassing_cutjet1SoftDropMass" ]->Fill( HT  );
					histos1D_[ "prunedMassAvePassing_cutjet1SoftDropMass" ]->Fill( prunedMassAve );
					histos1D_[ "softDropMassAvePassing_cutjet1SoftDropMass" ]->Fill( softDropMassAve );
					histos1D_[ "jet1PrunedMassPassing_cutjet1SoftDropMass" ]->Fill( JETS[0].prunedMass  );
					histos1D_[ "jet1SoftDropMassPassing_cutjet1SoftDropMass" ]->Fill( JETS[0].softDropMass  );
					histos1D_[ "jet1PtPassing_cutjet1SoftDropMass" ]->Fill( JETS[0].p4.Pt()   );
					histos1D_[ "jet2PrunedMassPassing_cutjet1SoftDropMass" ]->Fill( JETS[1].prunedMass  );
					histos1D_[ "jet2SoftDropMassPassing_cutjet1SoftDropMass" ]->Fill( JETS[1].softDropMass  );
					histos1D_[ "jet2PtPassing_cutjet1SoftDropMass" ]->Fill( JETS[1].p4.Pt() );

					histos2D_[ "jet1PrunedMassHTPassing_cutjet1SoftDropMass" ]->Fill( JETS[0].prunedMass, HT );
					histos2D_[ "jet1PrunedMassjet1PtPassing_cutjet1SoftDropMass" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
					histos2D_[ "jet1SoftDropMassHTPassing_cutjet1SoftDropMass" ]->Fill( JETS[0].softDropMass, HT );
					histos2D_[ "jet1SoftDropMassjet1PtPassing_cutjet1SoftDropMass" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

					histos2D_[ "jet2PrunedMassHTPassing_cutjet1SoftDropMass" ]->Fill( JETS[1].prunedMass, HT );
					histos2D_[ "jet2PrunedMassjet2PtPassing_cutjet1SoftDropMass" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
					histos2D_[ "jet2SoftDropMassHTPassing_cutjet1SoftDropMass" ]->Fill( JETS[1].softDropMass, HT );
					histos2D_[ "jet2SoftDropMassjet2PtPassing_cutjet1SoftDropMass" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

					histos2D_[ "prunedMassAveHTPassing_cutjet1SoftDropMass" ]->Fill( prunedMassAve, HT );
					histos2D_[ "prunedMassAvejet1PtPassing_cutjet1SoftDropMass" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
					histos2D_[ "prunedMassAvejet2PtPassing_cutjet1SoftDropMass" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
					histos2D_[ "softDropMassAveHTPassing_cutjet1SoftDropMass" ]->Fill( softDropMassAve, HT );
					histos2D_[ "softDropMassAvejet1PtPassing_cutjet1SoftDropMass" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
					histos2D_[ "softDropMassAvejet2PtPassing_cutjet1SoftDropMass" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );
				}
			}

			if ( ( JETS[0].p4.Pt() > cutAK8jet1Pt ) && ( JETS[1].p4.Pt() > cutAK8jet2Pt ) ) {
				histos1D_[ "HTDenom_cut2jetsPt" ]->Fill( HT  );
				histos1D_[ "prunedMassAveDenom_cut2jetsPt" ]->Fill( prunedMassAve );
				histos1D_[ "softDropMassAveDenom_cut2jetsPt" ]->Fill( softDropMassAve );
				histos1D_[ "jet1PrunedMassDenom_cut2jetsPt" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassDenom_cut2jetsPt" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtDenom_cut2jetsPt" ]->Fill( JETS[0].p4.Pt()   );
				histos1D_[ "jet2PrunedMassDenom_cut2jetsPt" ]->Fill( JETS[1].prunedMass  );
				histos1D_[ "jet2SoftDropMassDenom_cut2jetsPt" ]->Fill( JETS[1].softDropMass  );
				histos1D_[ "jet2PtDenom_cut2jetsPt" ]->Fill( JETS[1].p4.Pt() );

				histos2D_[ "jet1PrunedMassHTDenom_cut2jetsPt" ]->Fill( JETS[0].prunedMass, HT );
				histos2D_[ "jet1PrunedMassjet1PtDenom_cut2jetsPt" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassHTDenom_cut2jetsPt" ]->Fill( JETS[0].softDropMass, HT );
				histos2D_[ "jet1SoftDropMassjet1PtDenom_cut2jetsPt" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

				histos2D_[ "jet2PrunedMassHTDenom_cut2jetsPt" ]->Fill( JETS[1].prunedMass, HT );
				histos2D_[ "jet2PrunedMassjet2PtDenom_cut2jetsPt" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
				histos2D_[ "jet2SoftDropMassHTDenom_cut2jetsPt" ]->Fill( JETS[1].softDropMass, HT );
				histos2D_[ "jet2SoftDropMassjet2PtDenom_cut2jetsPt" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

				histos2D_[ "prunedMassAveHTDenom_cut2jetsPt" ]->Fill( prunedMassAve, HT );
				histos2D_[ "prunedMassAvejet1PtDenom_cut2jetsPt" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
				histos2D_[ "prunedMassAvejet2PtDenom_cut2jetsPt" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
				histos2D_[ "softDropMassAveHTDenom_cut2jetsPt" ]->Fill( softDropMassAve, HT );
				histos2D_[ "softDropMassAvejet1PtDenom_cut2jetsPt" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
				histos2D_[ "softDropMassAvejet2PtDenom_cut2jetsPt" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );

				if ( ORTriggers ){
					histos1D_[ "HTPassing_cut2jetsPt" ]->Fill( HT  );
					histos1D_[ "prunedMassAvePassing_cut2jetsPt" ]->Fill( prunedMassAve );
					histos1D_[ "softDropMassAvePassing_cut2jetsPt" ]->Fill( softDropMassAve );
					histos1D_[ "jet1PrunedMassPassing_cut2jetsPt" ]->Fill( JETS[0].prunedMass  );
					histos1D_[ "jet1SoftDropMassPassing_cut2jetsPt" ]->Fill( JETS[0].softDropMass  );
					histos1D_[ "jet1PtPassing_cut2jetsPt" ]->Fill( JETS[0].p4.Pt()   );
					histos1D_[ "jet2PrunedMassPassing_cut2jetsPt" ]->Fill( JETS[1].prunedMass  );
					histos1D_[ "jet2SoftDropMassPassing_cut2jetsPt" ]->Fill( JETS[1].softDropMass  );
					histos1D_[ "jet2PtPassing_cut2jetsPt" ]->Fill( JETS[1].p4.Pt() );

					histos2D_[ "jet1PrunedMassHTPassing_cut2jetsPt" ]->Fill( JETS[0].prunedMass, HT );
					histos2D_[ "jet1PrunedMassjet1PtPassing_cut2jetsPt" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
					histos2D_[ "jet1SoftDropMassHTPassing_cut2jetsPt" ]->Fill( JETS[0].softDropMass, HT );
					histos2D_[ "jet1SoftDropMassjet1PtPassing_cut2jetsPt" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

					histos2D_[ "jet2PrunedMassHTPassing_cut2jetsPt" ]->Fill( JETS[1].prunedMass, HT );
					histos2D_[ "jet2PrunedMassjet2PtPassing_cut2jetsPt" ]->Fill( JETS[1].prunedMass, JETS[1].p4.Pt() );
					histos2D_[ "jet2SoftDropMassHTPassing_cut2jetsPt" ]->Fill( JETS[1].softDropMass, HT );
					histos2D_[ "jet2SoftDropMassjet2PtPassing_cut2jetsPt" ]->Fill( JETS[1].softDropMass, JETS[1].p4.Pt() );

					histos2D_[ "prunedMassAveHTPassing_cut2jetsPt" ]->Fill( prunedMassAve, HT );
					histos2D_[ "prunedMassAvejet1PtPassing_cut2jetsPt" ]->Fill( prunedMassAve, JETS[0].p4.Pt() );
					histos2D_[ "prunedMassAvejet2PtPassing_cut2jetsPt" ]->Fill( prunedMassAve, JETS[1].p4.Pt() );
					histos2D_[ "softDropMassAveHTPassing_cut2jetsPt" ]->Fill( softDropMassAve, HT );
					histos2D_[ "softDropMassAvejet1PtPassing_cut2jetsPt" ]->Fill( softDropMassAve, JETS[0].p4.Pt() );
					histos2D_[ "softDropMassAvejet2PtPassing_cut2jetsPt" ]->Fill( softDropMassAve, JETS[1].p4.Pt() );
				}
			}
		}
	}
	JETS.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void RUNBoostedTriggerEfficiency::beginJob() {

	//////// test plots
	histos2D_[ "jet1PrunedMassHT_cutDijet_noTrigger" ] = fs_->make< TH2D >( "jet1PrunedMassHT_cutDijet_noTrigger", "jet1PrunedMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHT_cutDijet_noTrigger" ]->Sumw2();
	histos2D_[ "jet1PrunedMassHT_cutDijet_triggerOne" ] = fs_->make< TH2D >( "jet1PrunedMassHT_cutDijet_triggerOne", "jet1PrunedMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHT_cutDijet_triggerOne" ]->Sumw2();
	histos2D_[ "jet1PrunedMassHT_cutDijet_triggerTwo" ] = fs_->make< TH2D >( "jet1PrunedMassHT_cutDijet_triggerTwo", "jet1PrunedMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHT_cutDijet_triggerTwo" ]->Sumw2();
	histos2D_[ "jet1PrunedMassHT_cutDijet_triggerOneAndTwo" ] = fs_->make< TH2D >( "jet1PrunedMassHT_cutDijet_triggerOneAndTwo", "jet1PrunedMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHT_cutDijet_triggerOneAndTwo" ]->Sumw2();
	histos2D_[ "jet1PrunedMassHT_cutDijet_noTriggerOneOrTwo" ] = fs_->make< TH2D >( "jet1PrunedMassHT_cutDijet_noTriggerOneOrTwo", "jet1PrunedMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHT_cutDijet_noTriggerOneOrTwo" ]->Sumw2();


	histos2D_[ "jet1SoftDropMassHT_cutDijet_noTrigger" ] = fs_->make< TH2D >( "jet1SoftDropMassHT_cutDijet_noTrigger", "jet1SoftDropMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHT_cutDijet_noTrigger" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerOne" ] = fs_->make< TH2D >( "jet1SoftDropMassHT_cutDijet_triggerOne", "jet1SoftDropMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerOne" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassHT_cutDijet_triggerTwo", "jet1SoftDropMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerTwo" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerOneAndTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassHT_cutDijet_triggerOneAndTwo", "jet1SoftDropMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHT_cutDijet_triggerOneAndTwo" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHT_cutDijet_noTriggerOneOrTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassHT_cutDijet_noTriggerOneOrTwo", "jet1SoftDropMassHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHT_cutDijet_noTriggerOneOrTwo" ]->Sumw2();
	///////////////////////////////////////////////////////////////////

	/////// Denom cutDijet	
	histos1D_[ "HTDenom_cutDijet" ] = fs_->make< TH1D >( "HTDenom_cutDijet", "HTDenom_cutDijet", 500, 0., 5000. );
	histos1D_[ "HTDenom_cutDijet" ]->Sumw2();
	histos1D_[ "prunedMassAveDenom_cutDijet" ] = fs_->make< TH1D >( "prunedMassAveDenom_cutDijet", "prunedMassAveDenom_cutDijet", 60, 0., 600. );
	histos1D_[ "prunedMassAveDenom_cutDijet" ]->Sumw2();
	histos1D_[ "softDropMassAveDenom_cutDijet" ] = fs_->make< TH1D >( "softDropMassAveDenom_cutDijet", "softDropMassAveDenom_cutDijet", 60, 0., 600. );
	histos1D_[ "softDropMassAveDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet1PrunedMassDenom_cutDijet" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cutDijet", "jet1PrunedMassDenom_cutDijet", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cutDijet" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutDijet", "jet1SoftDropMassDenom_cutDijet", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cutDijet" ] = fs_->make< TH1D >( "jet1PtDenom_cutDijet", "jet1PtDenom_cutDijet", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet1RhoDenom_cutDijet" ] = fs_->make< TH1D >( "jet1RhoDenom_cutDijet", "jet1RhoDenom_cutDijet", 40, -10., 10. );
	histos1D_[ "jet1RhoDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet1RhoPrunedDenom_cutDijet" ] = fs_->make< TH1D >( "jet1RhoPrunedDenom_cutDijet", "jet1RhoPrunedDenom_cutDijet", 40, -10., 10. );
	histos1D_[ "jet1RhoPrunedDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet2PrunedMassDenom_cutDijet" ] = fs_->make< TH1D >( "jet2PrunedMassDenom_cutDijet", "jet2PrunedMassDenom_cutDijet", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassDenom_cutDijet" ] = fs_->make< TH1D >( "jet2SoftDropMassDenom_cutDijet", "jet2SoftDropMassDenom_cutDijet", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet2PtDenom_cutDijet" ] = fs_->make< TH1D >( "jet2PtDenom_cutDijet", "jet2PtDenom_cutDijet", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet2RhoDenom_cutDijet" ] = fs_->make< TH1D >( "jet2RhoDenom_cutDijet", "jet2RhoDenom_cutDijet", 40, -10., 10. );
	histos1D_[ "jet2RhoDenom_cutDijet" ]->Sumw2();
	histos1D_[ "jet2RhoPrunedDenom_cutDijet" ] = fs_->make< TH1D >( "jet2RhoPrunedDenom_cutDijet", "jet2RhoPrunedDenom_cutDijet", 40, -10., 10. );
	histos1D_[ "jet2RhoPrunedDenom_cutDijet" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTDenom_cutDijet" ] = fs_->make< TH2D >( "jet1PrunedMassHTDenom_cutDijet", "jet1PrunedMassHTDenom_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTDenom_cutDijet" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutDijet" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cutDijet", "jet1PrunedMassjet1PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutDijet" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTDenom_cutDijet" ] = fs_->make< TH2D >( "jet1SoftDropMassHTDenom_cutDijet", "jet1SoftDropMassHTDenom_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTDenom_cutDijet" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutDijet" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutDijet", "jet1SoftDropMassjet1PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutDijet" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTDenom_cutDijet" ] = fs_->make< TH2D >( "jet2PrunedMassHTDenom_cutDijet", "jet2PrunedMassHTDenom_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTDenom_cutDijet" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutDijet" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtDenom_cutDijet", "jet2PrunedMassjet2PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutDijet" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTDenom_cutDijet" ] = fs_->make< TH2D >( "jet2SoftDropMassHTDenom_cutDijet", "jet2SoftDropMassHTDenom_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTDenom_cutDijet" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutDijet" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtDenom_cutDijet", "jet2SoftDropMassjet2PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutDijet" ]->Sumw2();

	histos2D_[ "prunedMassAveHTDenom_cutDijet" ] = fs_->make< TH2D >( "prunedMassAveHTDenom_cutDijet", "prunedMassAveHTDenom_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTDenom_cutDijet" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtDenom_cutDijet" ] = fs_->make< TH2D >( "prunedMassAvejet1PtDenom_cutDijet", "prunedMassAvejet1PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtDenom_cutDijet" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtDenom_cutDijet" ] = fs_->make< TH2D >( "prunedMassAvejet2PtDenom_cutDijet", "prunedMassAvejet2PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtDenom_cutDijet" ]->Sumw2();
	histos2D_[ "softDropMassAveHTDenom_cutDijet" ] = fs_->make< TH2D >( "softDropMassAveHTDenom_cutDijet", "softDropMassAveHTDenom_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTDenom_cutDijet" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtDenom_cutDijet" ] = fs_->make< TH2D >( "softDropMassAvejet1PtDenom_cutDijet", "softDropMassAvejet1PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtDenom_cutDijet" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtDenom_cutDijet" ] = fs_->make< TH2D >( "softDropMassAvejet2PtDenom_cutDijet", "softDropMassAvejet2PtDenom_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtDenom_cutDijet" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cutDijet	
	histos1D_[ "HTPassing_cutDijet" ] = fs_->make< TH1D >( "HTPassing_cutDijet", "HTPassing_cutDijet", 500, 0., 5000. );
	histos1D_[ "HTPassing_cutDijet" ]->Sumw2();
	histos1D_[ "prunedMassAvePassing_cutDijet" ] = fs_->make< TH1D >( "prunedMassAvePassing_cutDijet", "prunedMassAvePassing_cutDijet", 60, 0., 600. );
	histos1D_[ "prunedMassAvePassing_cutDijet" ]->Sumw2();
	histos1D_[ "softDropMassAvePassing_cutDijet" ] = fs_->make< TH1D >( "softDropMassAvePassing_cutDijet", "softDropMassAvePassing_cutDijet", 60, 0., 600. );
	histos1D_[ "softDropMassAvePassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet1PrunedMassPassing_cutDijet" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cutDijet", "jet1PrunedMassPassing_cutDijet", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cutDijet" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutDijet", "jet1SoftDropMassPassing_cutDijet", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutDijet" ] = fs_->make< TH1D >( "jet1PtPassing_cutDijet", "jet1PtPassing_cutDijet", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet1RhoPassing_cutDijet" ] = fs_->make< TH1D >( "jet1RhoPassing_cutDijet", "jet1RhoPassing_cutDijet", 40, -10., 10. );
	histos1D_[ "jet1RhoPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet1RhoPrunedPassing_cutDijet" ] = fs_->make< TH1D >( "jet1RhoPrunedPassing_cutDijet", "jet1RhoPrunedPassing_cutDijet", 40, -10., 10. );
	histos1D_[ "jet1RhoPrunedPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet2PrunedMassPassing_cutDijet" ] = fs_->make< TH1D >( "jet2PrunedMassPassing_cutDijet", "jet2PrunedMassPassing_cutDijet", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassPassing_cutDijet" ] = fs_->make< TH1D >( "jet2SoftDropMassPassing_cutDijet", "jet2SoftDropMassPassing_cutDijet", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cutDijet" ] = fs_->make< TH1D >( "jet2PtPassing_cutDijet", "jet2PtPassing_cutDijet", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet2RhoPassing_cutDijet" ] = fs_->make< TH1D >( "jet2RhoPassing_cutDijet", "jet2RhoPassing_cutDijet", 40, -10., 10. );
	histos1D_[ "jet2RhoPassing_cutDijet" ]->Sumw2();
	histos1D_[ "jet2RhoPrunedPassing_cutDijet" ] = fs_->make< TH1D >( "jet2RhoPrunedPassing_cutDijet", "jet2RhoPrunedPassing_cutDijet", 40, -10., 10. );
	histos1D_[ "jet2RhoPrunedPassing_cutDijet" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTPassing_cutDijet" ] = fs_->make< TH2D >( "jet1PrunedMassHTPassing_cutDijet", "jet1PrunedMassHTPassing_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTPassing_cutDijet" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutDijet" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cutDijet", "jet1PrunedMassjet1PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutDijet" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTPassing_cutDijet" ] = fs_->make< TH2D >( "jet1SoftDropMassHTPassing_cutDijet", "jet1SoftDropMassHTPassing_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTPassing_cutDijet" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutDijet" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutDijet", "jet1SoftDropMassjet1PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutDijet" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTPassing_cutDijet" ] = fs_->make< TH2D >( "jet2PrunedMassHTPassing_cutDijet", "jet2PrunedMassHTPassing_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTPassing_cutDijet" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutDijet" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtPassing_cutDijet", "jet2PrunedMassjet2PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutDijet" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTPassing_cutDijet" ] = fs_->make< TH2D >( "jet2SoftDropMassHTPassing_cutDijet", "jet2SoftDropMassHTPassing_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTPassing_cutDijet" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutDijet" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtPassing_cutDijet", "jet2SoftDropMassjet2PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutDijet" ]->Sumw2();

	histos2D_[ "prunedMassAveHTPassing_cutDijet" ] = fs_->make< TH2D >( "prunedMassAveHTPassing_cutDijet", "prunedMassAveHTPassing_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTPassing_cutDijet" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtPassing_cutDijet" ] = fs_->make< TH2D >( "prunedMassAvejet1PtPassing_cutDijet", "prunedMassAvejet1PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtPassing_cutDijet" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtPassing_cutDijet" ] = fs_->make< TH2D >( "prunedMassAvejet2PtPassing_cutDijet", "prunedMassAvejet2PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtPassing_cutDijet" ]->Sumw2();
	histos2D_[ "softDropMassAveHTPassing_cutDijet" ] = fs_->make< TH2D >( "softDropMassAveHTPassing_cutDijet", "softDropMassAveHTPassing_cutDijet", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTPassing_cutDijet" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtPassing_cutDijet" ] = fs_->make< TH2D >( "softDropMassAvejet1PtPassing_cutDijet", "softDropMassAvejet1PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtPassing_cutDijet" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtPassing_cutDijet" ] = fs_->make< TH2D >( "softDropMassAvejet2PtPassing_cutDijet", "softDropMassAvejet2PtPassing_cutDijet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtPassing_cutDijet" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Denom cutHT	
	histos1D_[ "HTDenom_cutHT" ] = fs_->make< TH1D >( "HTDenom_cutHT", "HTDenom_cutHT", 500, 0., 5000. );
	histos1D_[ "HTDenom_cutHT" ]->Sumw2();
	histos1D_[ "prunedMassAveDenom_cutHT" ] = fs_->make< TH1D >( "prunedMassAveDenom_cutHT", "prunedMassAveDenom_cutHT", 60, 0., 600. );
	histos1D_[ "prunedMassAveDenom_cutHT" ]->Sumw2();
	histos1D_[ "softDropMassAveDenom_cutHT" ] = fs_->make< TH1D >( "softDropMassAveDenom_cutHT", "softDropMassAveDenom_cutHT", 60, 0., 600. );
	histos1D_[ "softDropMassAveDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet1PrunedMassDenom_cutHT" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cutHT", "jet1PrunedMassDenom_cutHT", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cutHT" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutHT", "jet1SoftDropMassDenom_cutHT", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cutHT" ] = fs_->make< TH1D >( "jet1PtDenom_cutHT", "jet1PtDenom_cutHT", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet1RhoDenom_cutHT" ] = fs_->make< TH1D >( "jet1RhoDenom_cutHT", "jet1RhoDenom_cutHT", 40, -10., 10. );
	histos1D_[ "jet1RhoDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet1RhoPrunedDenom_cutHT" ] = fs_->make< TH1D >( "jet1RhoPrunedDenom_cutHT", "jet1RhoPrunedDenom_cutHT", 40, -10., 10. );
	histos1D_[ "jet1RhoPrunedDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet2PrunedMassDenom_cutHT" ] = fs_->make< TH1D >( "jet2PrunedMassDenom_cutHT", "jet2PrunedMassDenom_cutHT", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassDenom_cutHT" ] = fs_->make< TH1D >( "jet2SoftDropMassDenom_cutHT", "jet2SoftDropMassDenom_cutHT", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet2PtDenom_cutHT" ] = fs_->make< TH1D >( "jet2PtDenom_cutHT", "jet2PtDenom_cutHT", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet2RhoDenom_cutHT" ] = fs_->make< TH1D >( "jet2RhoDenom_cutHT", "jet2RhoDenom_cutHT", 40, -10., 10. );
	histos1D_[ "jet2RhoDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet2RhoPrunedDenom_cutHT" ] = fs_->make< TH1D >( "jet2RhoPrunedDenom_cutHT", "jet2RhoPrunedDenom_cutHT", 40, -10., 10. );
	histos1D_[ "jet2RhoPrunedDenom_cutHT" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTDenom_cutHT" ] = fs_->make< TH2D >( "jet1PrunedMassHTDenom_cutHT", "jet1PrunedMassHTDenom_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTDenom_cutHT" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutHT" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cutHT", "jet1PrunedMassjet1PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutHT" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTDenom_cutHT" ] = fs_->make< TH2D >( "jet1SoftDropMassHTDenom_cutHT", "jet1SoftDropMassHTDenom_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTDenom_cutHT" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutHT" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutHT", "jet1SoftDropMassjet1PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutHT" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTDenom_cutHT" ] = fs_->make< TH2D >( "jet2PrunedMassHTDenom_cutHT", "jet2PrunedMassHTDenom_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTDenom_cutHT" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutHT" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtDenom_cutHT", "jet2PrunedMassjet2PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutHT" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTDenom_cutHT" ] = fs_->make< TH2D >( "jet2SoftDropMassHTDenom_cutHT", "jet2SoftDropMassHTDenom_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTDenom_cutHT" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutHT" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtDenom_cutHT", "jet2SoftDropMassjet2PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutHT" ]->Sumw2();

	histos2D_[ "prunedMassAveHTDenom_cutHT" ] = fs_->make< TH2D >( "prunedMassAveHTDenom_cutHT", "prunedMassAveHTDenom_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTDenom_cutHT" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtDenom_cutHT" ] = fs_->make< TH2D >( "prunedMassAvejet1PtDenom_cutHT", "prunedMassAvejet1PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtDenom_cutHT" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtDenom_cutHT" ] = fs_->make< TH2D >( "prunedMassAvejet2PtDenom_cutHT", "prunedMassAvejet2PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtDenom_cutHT" ]->Sumw2();
	histos2D_[ "softDropMassAveHTDenom_cutHT" ] = fs_->make< TH2D >( "softDropMassAveHTDenom_cutHT", "softDropMassAveHTDenom_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTDenom_cutHT" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtDenom_cutHT" ] = fs_->make< TH2D >( "softDropMassAvejet1PtDenom_cutHT", "softDropMassAvejet1PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtDenom_cutHT" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtDenom_cutHT" ] = fs_->make< TH2D >( "softDropMassAvejet2PtDenom_cutHT", "softDropMassAvejet2PtDenom_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtDenom_cutHT" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cutHT	
	histos1D_[ "HTPassing_cutHT" ] = fs_->make< TH1D >( "HTPassing_cutHT", "HTPassing_cutHT", 500, 0., 5000. );
	histos1D_[ "HTPassing_cutHT" ]->Sumw2();
	histos1D_[ "prunedMassAvePassing_cutHT" ] = fs_->make< TH1D >( "prunedMassAvePassing_cutHT", "prunedMassAvePassing_cutHT", 60, 0., 600. );
	histos1D_[ "prunedMassAvePassing_cutHT" ]->Sumw2();
	histos1D_[ "softDropMassAvePassing_cutHT" ] = fs_->make< TH1D >( "softDropMassAvePassing_cutHT", "softDropMassAvePassing_cutHT", 60, 0., 600. );
	histos1D_[ "softDropMassAvePassing_cutHT" ]->Sumw2();
	histos1D_[ "jet1PrunedMassPassing_cutHT" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cutHT", "jet1PrunedMassPassing_cutHT", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cutHT" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutHT", "jet1SoftDropMassPassing_cutHT", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutHT" ] = fs_->make< TH1D >( "jet1PtPassing_cutHT", "jet1PtPassing_cutHT", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet1RhoPassing_cutHT" ] = fs_->make< TH1D >( "jet1RhoPassing_cutHT", "jet1RhoPassing_cutHT", 40, -10., 10. );
	histos1D_[ "jet1RhoPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet1RhoPrunedPassing_cutHT" ] = fs_->make< TH1D >( "jet1RhoPrunedPassing_cutHT", "jet1RhoPrunedPassing_cutHT", 40, -10., 10. );
	histos1D_[ "jet1RhoPrunedPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet2PrunedMassPassing_cutHT" ] = fs_->make< TH1D >( "jet2PrunedMassPassing_cutHT", "jet2PrunedMassPassing_cutHT", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassPassing_cutHT" ] = fs_->make< TH1D >( "jet2SoftDropMassPassing_cutHT", "jet2SoftDropMassPassing_cutHT", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cutHT" ] = fs_->make< TH1D >( "jet2PtPassing_cutHT", "jet2PtPassing_cutHT", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet2RhoPassing_cutHT" ] = fs_->make< TH1D >( "jet2RhoPassing_cutHT", "jet2RhoPassing_cutHT", 40, -10., 10. );
	histos1D_[ "jet2RhoPassing_cutHT" ]->Sumw2();
	histos1D_[ "jet2RhoPrunedPassing_cutHT" ] = fs_->make< TH1D >( "jet2RhoPrunedPassing_cutHT", "jet2RhoPrunedPassing_cutHT", 40, -10., 10. );
	histos1D_[ "jet2RhoPrunedPassing_cutHT" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTPassing_cutHT" ] = fs_->make< TH2D >( "jet1PrunedMassHTPassing_cutHT", "jet1PrunedMassHTPassing_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTPassing_cutHT" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutHT" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cutHT", "jet1PrunedMassjet1PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutHT" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTPassing_cutHT" ] = fs_->make< TH2D >( "jet1SoftDropMassHTPassing_cutHT", "jet1SoftDropMassHTPassing_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTPassing_cutHT" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutHT" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutHT", "jet1SoftDropMassjet1PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutHT" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTPassing_cutHT" ] = fs_->make< TH2D >( "jet2PrunedMassHTPassing_cutHT", "jet2PrunedMassHTPassing_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTPassing_cutHT" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutHT" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtPassing_cutHT", "jet2PrunedMassjet2PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutHT" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTPassing_cutHT" ] = fs_->make< TH2D >( "jet2SoftDropMassHTPassing_cutHT", "jet2SoftDropMassHTPassing_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTPassing_cutHT" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutHT" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtPassing_cutHT", "jet2SoftDropMassjet2PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutHT" ]->Sumw2();

	histos2D_[ "prunedMassAveHTPassing_cutHT" ] = fs_->make< TH2D >( "prunedMassAveHTPassing_cutHT", "prunedMassAveHTPassing_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTPassing_cutHT" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtPassing_cutHT" ] = fs_->make< TH2D >( "prunedMassAvejet1PtPassing_cutHT", "prunedMassAvejet1PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtPassing_cutHT" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtPassing_cutHT" ] = fs_->make< TH2D >( "prunedMassAvejet2PtPassing_cutHT", "prunedMassAvejet2PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtPassing_cutHT" ]->Sumw2();
	histos2D_[ "softDropMassAveHTPassing_cutHT" ] = fs_->make< TH2D >( "softDropMassAveHTPassing_cutHT", "softDropMassAveHTPassing_cutHT", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTPassing_cutHT" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtPassing_cutHT" ] = fs_->make< TH2D >( "softDropMassAvejet1PtPassing_cutHT", "softDropMassAvejet1PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtPassing_cutHT" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtPassing_cutHT" ] = fs_->make< TH2D >( "softDropMassAvejet2PtPassing_cutHT", "softDropMassAvejet2PtPassing_cutHT", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtPassing_cutHT" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Denom cutjet1PrunedMass	
	histos1D_[ "HTDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "HTDenom_cutjet1PrunedMass", "HTDenom_cutjet1PrunedMass", 500, 0., 5000. );
	histos1D_[ "HTDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "prunedMassAveDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "prunedMassAveDenom_cutjet1PrunedMass", "prunedMassAveDenom_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "prunedMassAveDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "softDropMassAveDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "softDropMassAveDenom_cutjet1PrunedMass", "softDropMassAveDenom_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "softDropMassAveDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet1PrunedMassDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cutjet1PrunedMass", "jet1PrunedMassDenom_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutjet1PrunedMass", "jet1SoftDropMassDenom_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet1PtDenom_cutjet1PrunedMass", "jet1PtDenom_cutjet1PrunedMass", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet2PrunedMassDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet2PrunedMassDenom_cutjet1PrunedMass", "jet2PrunedMassDenom_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet2SoftDropMassDenom_cutjet1PrunedMass", "jet2SoftDropMassDenom_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassDenom_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet2PtDenom_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet2PtDenom_cutjet1PrunedMass", "jet2PtDenom_cutjet1PrunedMass", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cutjet1PrunedMass" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1PrunedMassHTDenom_cutjet1PrunedMass", "jet1PrunedMassHTDenom_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cutjet1PrunedMass", "jet1PrunedMassjet1PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1SoftDropMassHTDenom_cutjet1PrunedMass", "jet1SoftDropMassHTDenom_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutjet1PrunedMass", "jet1SoftDropMassjet1PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutjet1PrunedMass" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2PrunedMassHTDenom_cutjet1PrunedMass", "jet2PrunedMassHTDenom_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtDenom_cutjet1PrunedMass", "jet2PrunedMassjet2PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2SoftDropMassHTDenom_cutjet1PrunedMass", "jet2SoftDropMassHTDenom_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtDenom_cutjet1PrunedMass", "jet2SoftDropMassjet2PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutjet1PrunedMass" ]->Sumw2();

	histos2D_[ "prunedMassAveHTDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "prunedMassAveHTDenom_cutjet1PrunedMass", "prunedMassAveHTDenom_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "prunedMassAvejet1PtDenom_cutjet1PrunedMass", "prunedMassAvejet1PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "prunedMassAvejet2PtDenom_cutjet1PrunedMass", "prunedMassAvejet2PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "softDropMassAveHTDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "softDropMassAveHTDenom_cutjet1PrunedMass", "softDropMassAveHTDenom_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "softDropMassAvejet1PtDenom_cutjet1PrunedMass", "softDropMassAvejet1PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtDenom_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtDenom_cutjet1PrunedMass" ] = fs_->make< TH2D >( "softDropMassAvejet2PtDenom_cutjet1PrunedMass", "softDropMassAvejet2PtDenom_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtDenom_cutjet1PrunedMass" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cutjet1PrunedMass	
	histos1D_[ "HTPassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "HTPassing_cutjet1PrunedMass", "HTPassing_cutjet1PrunedMass", 500, 0., 5000. );
	histos1D_[ "HTPassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "prunedMassAvePassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "prunedMassAvePassing_cutjet1PrunedMass", "prunedMassAvePassing_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "prunedMassAvePassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "softDropMassAvePassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "softDropMassAvePassing_cutjet1PrunedMass", "softDropMassAvePassing_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "softDropMassAvePassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet1PrunedMassPassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cutjet1PrunedMass", "jet1PrunedMassPassing_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutjet1PrunedMass", "jet1SoftDropMassPassing_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet1PtPassing_cutjet1PrunedMass", "jet1PtPassing_cutjet1PrunedMass", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet2PrunedMassPassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet2PrunedMassPassing_cutjet1PrunedMass", "jet2PrunedMassPassing_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassPassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassPassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet2SoftDropMassPassing_cutjet1PrunedMass", "jet2SoftDropMassPassing_cutjet1PrunedMass", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassPassing_cutjet1PrunedMass" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cutjet1PrunedMass" ] = fs_->make< TH1D >( "jet2PtPassing_cutjet1PrunedMass", "jet2PtPassing_cutjet1PrunedMass", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cutjet1PrunedMass" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1PrunedMassHTPassing_cutjet1PrunedMass", "jet1PrunedMassHTPassing_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cutjet1PrunedMass", "jet1PrunedMassjet1PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1SoftDropMassHTPassing_cutjet1PrunedMass", "jet1SoftDropMassHTPassing_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutjet1PrunedMass", "jet1SoftDropMassjet1PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutjet1PrunedMass" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2PrunedMassHTPassing_cutjet1PrunedMass", "jet2PrunedMassHTPassing_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtPassing_cutjet1PrunedMass", "jet2PrunedMassjet2PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2SoftDropMassHTPassing_cutjet1PrunedMass", "jet2SoftDropMassHTPassing_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtPassing_cutjet1PrunedMass", "jet2SoftDropMassjet2PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutjet1PrunedMass" ]->Sumw2();

	histos2D_[ "prunedMassAveHTPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "prunedMassAveHTPassing_cutjet1PrunedMass", "prunedMassAveHTPassing_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "prunedMassAvejet1PtPassing_cutjet1PrunedMass", "prunedMassAvejet1PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "prunedMassAvejet2PtPassing_cutjet1PrunedMass", "prunedMassAvejet2PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "softDropMassAveHTPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "softDropMassAveHTPassing_cutjet1PrunedMass", "softDropMassAveHTPassing_cutjet1PrunedMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "softDropMassAvejet1PtPassing_cutjet1PrunedMass", "softDropMassAvejet1PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtPassing_cutjet1PrunedMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtPassing_cutjet1PrunedMass" ] = fs_->make< TH2D >( "softDropMassAvejet2PtPassing_cutjet1PrunedMass", "softDropMassAvejet2PtPassing_cutjet1PrunedMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtPassing_cutjet1PrunedMass" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Denom cutjet1SoftDropMass	
	histos1D_[ "HTDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "HTDenom_cutjet1SoftDropMass", "HTDenom_cutjet1SoftDropMass", 500, 0., 5000. );
	histos1D_[ "HTDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "prunedMassAveDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "prunedMassAveDenom_cutjet1SoftDropMass", "prunedMassAveDenom_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "prunedMassAveDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "softDropMassAveDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "softDropMassAveDenom_cutjet1SoftDropMass", "softDropMassAveDenom_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "softDropMassAveDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet1PrunedMassDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cutjet1SoftDropMass", "jet1PrunedMassDenom_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutjet1SoftDropMass", "jet1SoftDropMassDenom_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet1PtDenom_cutjet1SoftDropMass", "jet1PtDenom_cutjet1SoftDropMass", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet2PrunedMassDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet2PrunedMassDenom_cutjet1SoftDropMass", "jet2PrunedMassDenom_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet2SoftDropMassDenom_cutjet1SoftDropMass", "jet2SoftDropMassDenom_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet2PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet2PtDenom_cutjet1SoftDropMass", "jet2PtDenom_cutjet1SoftDropMass", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cutjet1SoftDropMass" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1PrunedMassHTDenom_cutjet1SoftDropMass", "jet1PrunedMassHTDenom_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cutjet1SoftDropMass", "jet1PrunedMassjet1PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1SoftDropMassHTDenom_cutjet1SoftDropMass", "jet1SoftDropMassHTDenom_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutjet1SoftDropMass", "jet1SoftDropMassjet1PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutjet1SoftDropMass" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2PrunedMassHTDenom_cutjet1SoftDropMass", "jet2PrunedMassHTDenom_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtDenom_cutjet1SoftDropMass", "jet2PrunedMassjet2PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2SoftDropMassHTDenom_cutjet1SoftDropMass", "jet2SoftDropMassHTDenom_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtDenom_cutjet1SoftDropMass", "jet2SoftDropMassjet2PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cutjet1SoftDropMass" ]->Sumw2();

	histos2D_[ "prunedMassAveHTDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "prunedMassAveHTDenom_cutjet1SoftDropMass", "prunedMassAveHTDenom_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "prunedMassAvejet1PtDenom_cutjet1SoftDropMass", "prunedMassAvejet1PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "prunedMassAvejet2PtDenom_cutjet1SoftDropMass", "prunedMassAvejet2PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "softDropMassAveHTDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "softDropMassAveHTDenom_cutjet1SoftDropMass", "softDropMassAveHTDenom_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "softDropMassAvejet1PtDenom_cutjet1SoftDropMass", "softDropMassAvejet1PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtDenom_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtDenom_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "softDropMassAvejet2PtDenom_cutjet1SoftDropMass", "softDropMassAvejet2PtDenom_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtDenom_cutjet1SoftDropMass" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cutjet1SoftDropMass	
	histos1D_[ "HTPassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "HTPassing_cutjet1SoftDropMass", "HTPassing_cutjet1SoftDropMass", 500, 0., 5000. );
	histos1D_[ "HTPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "prunedMassAvePassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "prunedMassAvePassing_cutjet1SoftDropMass", "prunedMassAvePassing_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "prunedMassAvePassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "softDropMassAvePassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "softDropMassAvePassing_cutjet1SoftDropMass", "softDropMassAvePassing_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "softDropMassAvePassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet1PrunedMassPassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cutjet1SoftDropMass", "jet1PrunedMassPassing_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutjet1SoftDropMass", "jet1SoftDropMassPassing_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet1PtPassing_cutjet1SoftDropMass", "jet1PtPassing_cutjet1SoftDropMass", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet2PrunedMassPassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet2PrunedMassPassing_cutjet1SoftDropMass", "jet2PrunedMassPassing_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassPassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet2SoftDropMassPassing_cutjet1SoftDropMass", "jet2SoftDropMassPassing_cutjet1SoftDropMass", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH1D >( "jet2PtPassing_cutjet1SoftDropMass", "jet2PtPassing_cutjet1SoftDropMass", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cutjet1SoftDropMass" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1PrunedMassHTPassing_cutjet1SoftDropMass", "jet1PrunedMassHTPassing_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cutjet1SoftDropMass", "jet1PrunedMassjet1PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1SoftDropMassHTPassing_cutjet1SoftDropMass", "jet1SoftDropMassHTPassing_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutjet1SoftDropMass", "jet1SoftDropMassjet1PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutjet1SoftDropMass" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2PrunedMassHTPassing_cutjet1SoftDropMass", "jet2PrunedMassHTPassing_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtPassing_cutjet1SoftDropMass", "jet2PrunedMassjet2PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2SoftDropMassHTPassing_cutjet1SoftDropMass", "jet2SoftDropMassHTPassing_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtPassing_cutjet1SoftDropMass", "jet2SoftDropMassjet2PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cutjet1SoftDropMass" ]->Sumw2();

	histos2D_[ "prunedMassAveHTPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "prunedMassAveHTPassing_cutjet1SoftDropMass", "prunedMassAveHTPassing_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "prunedMassAvejet1PtPassing_cutjet1SoftDropMass", "prunedMassAvejet1PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "prunedMassAvejet2PtPassing_cutjet1SoftDropMass", "prunedMassAvejet2PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "softDropMassAveHTPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "softDropMassAveHTPassing_cutjet1SoftDropMass", "softDropMassAveHTPassing_cutjet1SoftDropMass", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "softDropMassAvejet1PtPassing_cutjet1SoftDropMass", "softDropMassAvejet1PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtPassing_cutjet1SoftDropMass" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtPassing_cutjet1SoftDropMass" ] = fs_->make< TH2D >( "softDropMassAvejet2PtPassing_cutjet1SoftDropMass", "softDropMassAvejet2PtPassing_cutjet1SoftDropMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtPassing_cutjet1SoftDropMass" ]->Sumw2();
	///////////////////////////////////////////////////////////////////



	/////// Denom cut2jetsPt	
	histos1D_[ "HTDenom_cut2jetsPt" ] = fs_->make< TH1D >( "HTDenom_cut2jetsPt", "HTDenom_cut2jetsPt", 500, 0., 5000. );
	histos1D_[ "HTDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "prunedMassAveDenom_cut2jetsPt" ] = fs_->make< TH1D >( "prunedMassAveDenom_cut2jetsPt", "prunedMassAveDenom_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "prunedMassAveDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "softDropMassAveDenom_cut2jetsPt" ] = fs_->make< TH1D >( "softDropMassAveDenom_cut2jetsPt", "softDropMassAveDenom_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "softDropMassAveDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet1PrunedMassDenom_cut2jetsPt" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cut2jetsPt", "jet1PrunedMassDenom_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cut2jetsPt" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cut2jetsPt", "jet1SoftDropMassDenom_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cut2jetsPt" ] = fs_->make< TH1D >( "jet1PtDenom_cut2jetsPt", "jet1PtDenom_cut2jetsPt", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet2PrunedMassDenom_cut2jetsPt" ] = fs_->make< TH1D >( "jet2PrunedMassDenom_cut2jetsPt", "jet2PrunedMassDenom_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassDenom_cut2jetsPt" ] = fs_->make< TH1D >( "jet2SoftDropMassDenom_cut2jetsPt", "jet2SoftDropMassDenom_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassDenom_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet2PtDenom_cut2jetsPt" ] = fs_->make< TH1D >( "jet2PtDenom_cut2jetsPt", "jet2PtDenom_cut2jetsPt", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cut2jetsPt" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet1PrunedMassHTDenom_cut2jetsPt", "jet1PrunedMassHTDenom_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cut2jetsPt", "jet1PrunedMassjet1PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet1SoftDropMassHTDenom_cut2jetsPt", "jet1SoftDropMassHTDenom_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cut2jetsPt", "jet1SoftDropMassjet1PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cut2jetsPt" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet2PrunedMassHTDenom_cut2jetsPt", "jet2PrunedMassHTDenom_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtDenom_cut2jetsPt", "jet2PrunedMassjet2PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet2SoftDropMassHTDenom_cut2jetsPt", "jet2SoftDropMassHTDenom_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtDenom_cut2jetsPt", "jet2SoftDropMassjet2PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtDenom_cut2jetsPt" ]->Sumw2();

	histos2D_[ "prunedMassAveHTDenom_cut2jetsPt" ] = fs_->make< TH2D >( "prunedMassAveHTDenom_cut2jetsPt", "prunedMassAveHTDenom_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "prunedMassAvejet1PtDenom_cut2jetsPt", "prunedMassAvejet1PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "prunedMassAvejet2PtDenom_cut2jetsPt", "prunedMassAvejet2PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "softDropMassAveHTDenom_cut2jetsPt" ] = fs_->make< TH2D >( "softDropMassAveHTDenom_cut2jetsPt", "softDropMassAveHTDenom_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "softDropMassAvejet1PtDenom_cut2jetsPt", "softDropMassAvejet1PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtDenom_cut2jetsPt" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtDenom_cut2jetsPt" ] = fs_->make< TH2D >( "softDropMassAvejet2PtDenom_cut2jetsPt", "softDropMassAvejet2PtDenom_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtDenom_cut2jetsPt" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cut2jetsPt	
	histos1D_[ "HTPassing_cut2jetsPt" ] = fs_->make< TH1D >( "HTPassing_cut2jetsPt", "HTPassing_cut2jetsPt", 500, 0., 5000. );
	histos1D_[ "HTPassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "prunedMassAvePassing_cut2jetsPt" ] = fs_->make< TH1D >( "prunedMassAvePassing_cut2jetsPt", "prunedMassAvePassing_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "prunedMassAvePassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "softDropMassAvePassing_cut2jetsPt" ] = fs_->make< TH1D >( "softDropMassAvePassing_cut2jetsPt", "softDropMassAvePassing_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "softDropMassAvePassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet1PrunedMassPassing_cut2jetsPt" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cut2jetsPt", "jet1PrunedMassPassing_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cut2jetsPt" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cut2jetsPt", "jet1SoftDropMassPassing_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cut2jetsPt" ] = fs_->make< TH1D >( "jet1PtPassing_cut2jetsPt", "jet1PtPassing_cut2jetsPt", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet2PrunedMassPassing_cut2jetsPt" ] = fs_->make< TH1D >( "jet2PrunedMassPassing_cut2jetsPt", "jet2PrunedMassPassing_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet2PrunedMassPassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet2SoftDropMassPassing_cut2jetsPt" ] = fs_->make< TH1D >( "jet2SoftDropMassPassing_cut2jetsPt", "jet2SoftDropMassPassing_cut2jetsPt", 60, 0., 600. );
	histos1D_[ "jet2SoftDropMassPassing_cut2jetsPt" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cut2jetsPt" ] = fs_->make< TH1D >( "jet2PtPassing_cut2jetsPt", "jet2PtPassing_cut2jetsPt", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cut2jetsPt" ]->Sumw2();

	histos2D_[ "jet1PrunedMassHTPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet1PrunedMassHTPassing_cut2jetsPt", "jet1PrunedMassHTPassing_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1PrunedMassHTPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet1PrunedMassjet1PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cut2jetsPt", "jet1PrunedMassjet1PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassHTPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet1SoftDropMassHTPassing_cut2jetsPt", "jet1SoftDropMassHTPassing_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet1SoftDropMassHTPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cut2jetsPt", "jet1SoftDropMassjet1PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cut2jetsPt" ]->Sumw2();

	histos2D_[ "jet2PrunedMassHTPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet2PrunedMassHTPassing_cut2jetsPt", "jet2PrunedMassHTPassing_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2PrunedMassHTPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet2PrunedMassjet2PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet2PrunedMassjet2PtPassing_cut2jetsPt", "jet2PrunedMassjet2PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2PrunedMassjet2PtPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassHTPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet2SoftDropMassHTPassing_cut2jetsPt", "jet2SoftDropMassHTPassing_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "jet2SoftDropMassHTPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "jet2SoftDropMassjet2PtPassing_cut2jetsPt", "jet2SoftDropMassjet2PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet2SoftDropMassjet2PtPassing_cut2jetsPt" ]->Sumw2();

	histos2D_[ "prunedMassAveHTPassing_cut2jetsPt" ] = fs_->make< TH2D >( "prunedMassAveHTPassing_cut2jetsPt", "prunedMassAveHTPassing_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "prunedMassAveHTPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "prunedMassAvejet1PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "prunedMassAvejet1PtPassing_cut2jetsPt", "prunedMassAvejet1PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet1PtPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "prunedMassAvejet2PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "prunedMassAvejet2PtPassing_cut2jetsPt", "prunedMassAvejet2PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "prunedMassAvejet2PtPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "softDropMassAveHTPassing_cut2jetsPt" ] = fs_->make< TH2D >( "softDropMassAveHTPassing_cut2jetsPt", "softDropMassAveHTPassing_cut2jetsPt", 60, 0., 600., 500, 0., 5000.);
	histos2D_[ "softDropMassAveHTPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "softDropMassAvejet1PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "softDropMassAvejet1PtPassing_cut2jetsPt", "softDropMassAvejet1PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet1PtPassing_cut2jetsPt" ]->Sumw2();
	histos2D_[ "softDropMassAvejet2PtPassing_cut2jetsPt" ] = fs_->make< TH2D >( "softDropMassAvejet2PtPassing_cut2jetsPt", "softDropMassAvejet2PtPassing_cut2jetsPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "softDropMassAvejet2PtPassing_cut2jetsPt" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


}

// ------------ method called once each job just after ending the event loop  ------------
void RUNBoostedTriggerEfficiency::endJob() {

}

void RUNBoostedTriggerEfficiency::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<double>("cutAK8jetPt", 150.);
	desc.add<double>("cutAK8HT", 900.);
	desc.add<double>("cutAK8jet1Pt", 500.);
	desc.add<double>("cutAK8jet2Pt", 450.);
	desc.add<double>("cutAK8jet1Mass", 60.);
	desc.add<string>("baseTrigger", "HLT_PFHT475");
	desc.add<string>("jecVersion", "supportFiles/Fall15_25nsV2");
	desc.add<string>("PUMethod", "chs");
	vector<string> HLTPass;
	HLTPass.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");
	desc.add<vector<string>>("triggerPass",	HLTPass);

	desc.add<InputTag>("Lumi", 	InputTag("eventInfo:evtInfoLumiBlock"));
	desc.add<InputTag>("Run", 	InputTag("eventInfo:evtInfoRunNumber"));
	desc.add<InputTag>("Event", 	InputTag("eventInfo:evtInfoEventNumber"));
	desc.add<InputTag>("NPV", 	InputTag("vertexInfo:npv"));
	desc.add<InputTag>("rho", 	InputTag("vertexInfo:rho"));
	desc.add<InputTag>("jetPt", 	InputTag("jetsAK8CHS:jetAK8CHSPt"));
	desc.add<InputTag>("jetEta", 	InputTag("jetsAK8CHS:jetAK8CHSEta"));
	desc.add<InputTag>("jetPhi", 	InputTag("jetsAK8CHS:jetAK8CHSPhi"));
	desc.add<InputTag>("jetE", 	InputTag("jetsAK8CHS:jetAK8CHSE"));
	desc.add<InputTag>("jetArea", 	InputTag("jetsAK8CHS:jetAK8CHSjetArea"));
	desc.add<InputTag>("jetPrunedMass", 	InputTag("jetsAK8CHS:jetAK8CHSprunedMassCHS"));
	desc.add<InputTag>("jetSoftDropMass", 	InputTag("jetsAK8CHS:jetAK8CHSsoftDropMassCHS"));
	desc.add<InputTag>("jetTrimmedMass", 	InputTag("jetsAK8CHS:jetAK8CHStrimmedMassCHS"));
	desc.add<InputTag>("jetCSVv2", 	InputTag("jetsAK8CHS:jetAK8CHSCSVv2"));
	// JetID
	desc.add<InputTag>("jecFactor", 		InputTag("jetsAK8CHS:jetAK8CHSjecFactor0"));
	desc.add<InputTag>("neutralHadronEnergyFrac", 	InputTag("jetsAK8CHS:jetAK8CHSneutralHadronEnergyFrac"));
	desc.add<InputTag>("neutralEmEnergyFrac", 	InputTag("jetsAK8CHS:jetAK8CHSneutralEmEnergyFrac"));
	desc.add<InputTag>("chargedEmEnergyFrac", 	InputTag("jetsAK8CHS:jetAK8CHSchargedEmEnergyFrac"));
	desc.add<InputTag>("muonEnergy", 		InputTag("jetsAK8CHS:jetAK8CHSMuonEnergy"));
	desc.add<InputTag>("chargedHadronEnergyFrac", 	InputTag("jetsAK8CHS:jetAK8CHSchargedHadronEnergyFrac"));
	desc.add<InputTag>("neutralMultiplicity",	InputTag("jetsAK8CHS:jetAK8CHSneutralMultiplicity"));
	desc.add<InputTag>("chargedMultiplicity", 	InputTag("jetsAK8CHS:jetAK8CHSchargedMultiplicity"));
	// Trigger
	desc.add<InputTag>("triggerPrescale",		InputTag("TriggerUserData:triggerPrescaleTree"));
	desc.add<InputTag>("triggerBit",		InputTag("TriggerUserData:triggerBitTree"));
	desc.add<InputTag>("triggerName",		InputTag("TriggerUserData:triggerNameTree"));
	descriptions.addDefault(desc);
}
      
void RUNBoostedTriggerEfficiency::beginRun(const Run& iRun, const EventSetup& iSetup){

	/// Getting the names of the triggers from Run
	Handle<vector<string> > triggerName;
	iRun.getByToken(triggerName_, triggerName);
	LogWarning("TriggerNames") << "List of triggers found:";
	for (size_t q = 0; q < triggerName->size(); q++) {
		triggerNamesList.push_back( (*triggerName)[q] );
		cout << (*triggerName)[q] << endl; 
	}
	if ( triggerNamesList.size() == 0 ) LogError("TriggerNames") << "No triggers found.";
		
}

void RUNBoostedTriggerEfficiency::endRun(const Run& iRun, const EventSetup& iSetup){
	triggerNamesList.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(RUNBoostedTriggerEfficiency);
