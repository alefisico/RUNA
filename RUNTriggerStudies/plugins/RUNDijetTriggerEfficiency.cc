// -*- C++ -*-
//
// Package:    RUNA/RUNTriggerEfficiency
// Class:      RUNDijetTriggerEfficiency
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
class RUNDijetTriggerEfficiency : public EDAnalyzer {
	public:
		explicit RUNDijetTriggerEfficiency(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~RUNDijetTriggerEfficiency();

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
		double cutAK8jet1Pt;
		double cutAK8jet1Mass;
		TString baseTrigger;
		string jecVersion;
		string PUMethod;
		vector<string> triggerPass, triggerNamesList;
		vector<JetCorrectorParameters> jetPar;
		FactorizedJetCorrector * jetJECAK8;
		vector<JetCorrectorParameters> massPar;
		FactorizedJetCorrector * massJECAK8;

		EDGetTokenT<vector<float>> jetPt_;
		EDGetTokenT<vector<float>> jetEta_;
		EDGetTokenT<vector<float>> jetPhi_;
		EDGetTokenT<vector<float>> jetE_;
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

		/// Muons
		EDGetTokenT<vector<float>> muonPt_;
		EDGetTokenT<vector<float>> muonEta_;
		EDGetTokenT<vector<float>> muonPhi_;
		EDGetTokenT<vector<float>> muonE_;
		EDGetTokenT<vector<float>> muonIsLoose_;
		EDGetTokenT<vector<float>> muonIsGlobal_;
};

//
// static data member definitions
//

//
// constructors and destructor
//
RUNDijetTriggerEfficiency::RUNDijetTriggerEfficiency(const ParameterSet& iConfig):
	jetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPt"))),
	jetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetEta"))),
	jetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPhi"))),
	jetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetE"))),
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
	muonEnergy_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonEnergy"))),
	// Muons
	muonPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonPt"))), 
	muonEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonEta"))), 
	muonPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonPhi"))), 
	muonE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonE"))), 
	muonIsLoose_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonIsLoose"))), 
	muonIsGlobal_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonIsGlobal")))
{
	cutAK8jetPt = iConfig.getParameter<double>("cutAK8jetPt");
	cutAK8jet1Pt = iConfig.getParameter<double>("cutAK8jet1Pt");
	cutAK8jet1Mass = iConfig.getParameter<double>("cutAK8jet1Mass");
	baseTrigger = iConfig.getParameter<string>("baseTrigger");
	triggerPass = iConfig.getParameter<vector<string>>("triggerPass");
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


RUNDijetTriggerEfficiency::~RUNDijetTriggerEfficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNDijetTriggerEfficiency::analyze(const Event& iEvent, const EventSetup& iSetup) {

	Handle<vector<float> > jetPt;
	iEvent.getByToken(jetPt_, jetPt);

	Handle<vector<float> > jetEta;
	iEvent.getByToken(jetEta_, jetEta);

	Handle<vector<float> > jetPhi;
	iEvent.getByToken(jetPhi_, jetPhi);

	Handle<vector<float> > jetE;
	iEvent.getByToken(jetE_, jetE);

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

	/// Muons 
	Handle<vector<float> > muonPt;
	iEvent.getByToken(muonPt_, muonPt);

	Handle<vector<float> > muonEta;
	iEvent.getByToken(muonEta_, muonEta);

	Handle<vector<float> > muonPhi;
	iEvent.getByToken(muonPhi_, muonPhi);

	Handle<vector<float> > muonE;
	iEvent.getByToken(muonE_, muonE);

	Handle<vector<float> > muonIsLoose;
	iEvent.getByToken(muonIsLoose_, muonIsLoose);

	Handle<vector<float> > muonIsGlobal;
	iEvent.getByToken(muonIsGlobal_, muonIsGlobal);


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

		if( ( corrJet.Pt() > cutAK8jetPt ) && jetId ) { 

			HT += corrJet.Pt();
			double massJEC = corrections( rawJet, (*jetArea)[i], (*rho)[i] ,*NPV, massJECAK8); 
			double corrSoftDropMass = (*jetSoftDropMass)[i] * massJEC ;
			double corrPrunedMass = (*jetPrunedMass)[i] * massJEC ;

			myJet tmpJET;
			tmpJET.p4 = corrJet;
			tmpJET.prunedMass = corrPrunedMass;
			tmpJET.softDropMass = corrSoftDropMass;
			JETS.push_back( tmpJET );
		}
	}

	////// Muon veto
	vector< TLorentzVector > MUONS;
	if ( muonPt->size() > 0 ) {
		for (size_t m = 0; m < muonPt->size(); m++) {
			if( ( (*muonIsLoose)[m] ) && ( (*muonPt)[m] > 10 ) && ( TMath::Abs((*muonEta)[m]) < 2.5 ) ) {
				//LogWarning("muon") << (*muonPt)[m] << " " << (*muonEta)[m] << " " << (*muonIsLoose)[m] << " " << (*muonIsGlobal)[m] ;
				//muonsPt->push_back( (*muonPt)[m] );
				TLorentzVector tmpMuon;
				tmpMuon.SetPtEtaPhiE( (*muonPt)[m], (*muonEta)[m], (*muonPhi)[m], (*muonE)[m] );
				MUONS.push_back( tmpMuon );
			}
		}
	}


	bool isoMuon = true;
	if ( JETS.size() > 0 ) {

		histos2D_[ "jet1PrunedMassPt_cutJet_noTrigger" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
		histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

		if ( basedTriggerFired || ORTriggers ) {
			if ( basedTriggerFired && ORTriggers ) {
				histos2D_[ "jet1PrunedMassPt_cutJet_triggerOneAndTwo" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOneAndTwo" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );
			} else if ( basedTriggerFired ) {
				histos2D_[ "jet1PrunedMassPt_cutJet_triggerOne" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOne" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );
			} else if ( ORTriggers ) {
				histos2D_[ "jet1PrunedMassPt_cutJet_triggerTwo" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassPt_cutJet_triggerTwo" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );
			}
		} else {
			histos2D_[ "jet1PrunedMassPt_cutJet_noTriggerOneOrTwo" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
			histos2D_[ "jet1SoftDropMassPt_cutJet_noTriggerOneOrTwo" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );
		}

		if ( MUONS.size() > 0 ) {
			double jetMuonDeltaR = JETS[0].p4.DeltaR( MUONS[0] );
			if ( jetMuonDeltaR < 0.8 ) isoMuon = false;
			//LogWarning( "muons" ) << MUONS[0].Pt() << " " << jetMuonDeltaR << " " << isoMuon << " " << basedTriggerFired;
		}


		if ( basedTriggerFired && isoMuon ) {
			//LogWarning( "passed" ) << basedTriggerFired << " " << isoMuon ;
			histos1D_[ "jet1PrunedMassDenom_cutJet" ]->Fill( JETS[0].prunedMass  );
			histos1D_[ "jet1SoftDropMassDenom_cutJet" ]->Fill( JETS[0].softDropMass  );
			histos1D_[ "jet1PtDenom_cutJet" ]->Fill( JETS[0].p4.Pt()   );

			histos2D_[ "jet1PrunedMassjet1PtDenom_cutJet" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
			histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJet" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

			if ( ORTriggers ){
				histos1D_[ "jet1PrunedMassPassing_cutJet" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassPassing_cutJet" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtPassing_cutJet" ]->Fill( JETS[0].p4.Pt()   );

				histos2D_[ "jet1PrunedMassjet1PtPassing_cutJet" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJet" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

			}

			if ( JETS[0].softDropMass > cutAK8jet1Mass ) {
				histos1D_[ "jet1PrunedMassDenom_cutJetMass" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassDenom_cutJetMass" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtDenom_cutJetMass" ]->Fill( JETS[0].p4.Pt()   );

				histos2D_[ "jet1PrunedMassjet1PtDenom_cutJetMass" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetMass" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

				if ( ORTriggers ){
					histos1D_[ "jet1PrunedMassPassing_cutJetMass" ]->Fill( JETS[0].prunedMass  );
					histos1D_[ "jet1SoftDropMassPassing_cutJetMass" ]->Fill( JETS[0].softDropMass  );
					histos1D_[ "jet1PtPassing_cutJetMass" ]->Fill( JETS[0].p4.Pt()   );

					histos2D_[ "jet1PrunedMassjet1PtPassing_cutJetMass" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
					histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetMass" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );
				}
			}

			if ( JETS[0].p4.Pt() > cutAK8jet1Pt ) {
				histos1D_[ "jet1PrunedMassDenom_cutJetPt" ]->Fill( JETS[0].prunedMass  );
				histos1D_[ "jet1SoftDropMassDenom_cutJetPt" ]->Fill( JETS[0].softDropMass  );
				histos1D_[ "jet1PtDenom_cutJetPt" ]->Fill( JETS[0].p4.Pt()   );

				histos2D_[ "jet1PrunedMassjet1PtDenom_cutJetPt" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
				histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetPt" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );

				if ( ORTriggers ){
					histos1D_[ "jet1PrunedMassPassing_cutJetPt" ]->Fill( JETS[0].prunedMass  );
					histos1D_[ "jet1SoftDropMassPassing_cutJetPt" ]->Fill( JETS[0].softDropMass  );
					histos1D_[ "jet1PtPassing_cutJetPt" ]->Fill( JETS[0].p4.Pt()   );

					histos2D_[ "jet1PrunedMassjet1PtPassing_cutJetPt" ]->Fill( JETS[0].prunedMass, JETS[0].p4.Pt() );
					histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetPt" ]->Fill( JETS[0].softDropMass, JETS[0].p4.Pt() );
				}
			}

		}
	}
	JETS.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void RUNDijetTriggerEfficiency::beginJob() {

	//////// test plots
	histos2D_[ "jet1PrunedMassPt_cutJet_noTrigger" ] = fs_->make< TH2D >( "jet1PrunedMassPt_cutJet_noTrigger", "jet1PrunedMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassPt_cutJet_noTrigger" ]->Sumw2();
	histos2D_[ "jet1PrunedMassPt_cutJet_noTrigger" ] = fs_->make< TH2D >( "jet1PrunedMassPt_cutJet_noTrigger", "jet1PrunedMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassPt_cutJet_noTrigger" ]->Sumw2();
	histos2D_[ "jet1PrunedMassPt_cutJet_triggerOne" ] = fs_->make< TH2D >( "jet1PrunedMassPt_cutJet_triggerOne", "jet1PrunedMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassPt_cutJet_triggerOne" ]->Sumw2();
	histos2D_[ "jet1PrunedMassPt_cutJet_triggerTwo" ] = fs_->make< TH2D >( "jet1PrunedMassPt_cutJet_triggerTwo", "jet1PrunedMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassPt_cutJet_triggerTwo" ]->Sumw2();
	histos2D_[ "jet1PrunedMassPt_cutJet_triggerOneAndTwo" ] = fs_->make< TH2D >( "jet1PrunedMassPt_cutJet_triggerOneAndTwo", "jet1PrunedMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassPt_cutJet_triggerOneAndTwo" ]->Sumw2();
	histos2D_[ "jet1PrunedMassPt_cutJet_noTriggerOneOrTwo" ] = fs_->make< TH2D >( "jet1PrunedMassPt_cutJet_noTriggerOneOrTwo", "jet1PrunedMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassPt_cutJet_noTriggerOneOrTwo" ]->Sumw2();


	histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_noTrigger", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_noTrigger", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOne" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_triggerOne", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOne" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_triggerTwo", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerTwo" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOneAndTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_triggerOneAndTwo", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOneAndTwo" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassPt_cutJet_noTriggerOneOrTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_noTriggerOneOrTwo", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_noTriggerOneOrTwo" ]->Sumw2();
	///////////////////////////////////////////////////////////////////

	/////// Denom cutJet	
	histos1D_[ "jet1PrunedMassDenom_cutJet" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cutJet", "jet1PrunedMassDenom_cutJet", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cutJet" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cutJet" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutJet", "jet1SoftDropMassDenom_cutJet", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cutJet" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cutJet" ] = fs_->make< TH1D >( "jet1PtDenom_cutJet", "jet1PtDenom_cutJet", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutJet" ]->Sumw2();

	histos2D_[ "jet1PrunedMassjet1PtDenom_cutJet" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cutJet", "jet1PrunedMassjet1PtDenom_cutJet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutJet" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJet" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutJet", "jet1SoftDropMassjet1PtDenom_cutJet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJet" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cutJet	
	histos1D_[ "jet1PrunedMassPassing_cutJet" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cutJet", "jet1PrunedMassPassing_cutJet", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cutJet" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cutJet" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutJet", "jet1SoftDropMassPassing_cutJet", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cutJet" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutJet" ] = fs_->make< TH1D >( "jet1PtPassing_cutJet", "jet1PtPassing_cutJet", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutJet" ]->Sumw2();

	histos2D_[ "jet1PrunedMassjet1PtPassing_cutJet" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cutJet", "jet1PrunedMassjet1PtPassing_cutJet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutJet" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJet" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutJet", "jet1SoftDropMassjet1PtPassing_cutJet", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJet" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Denom cutJetMass	
	histos1D_[ "jet1PrunedMassDenom_cutJetMass" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cutJetMass", "jet1PrunedMassDenom_cutJetMass", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cutJetMass" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cutJetMass" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutJetMass", "jet1SoftDropMassDenom_cutJetMass", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cutJetMass" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cutJetMass" ] = fs_->make< TH1D >( "jet1PtDenom_cutJetMass", "jet1PtDenom_cutJetMass", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutJetMass" ]->Sumw2();

	histos2D_[ "jet1PrunedMassjet1PtDenom_cutJetMass" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cutJetMass", "jet1PrunedMassjet1PtDenom_cutJetMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutJetMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutJetMass", "jet1SoftDropMassjet1PtDenom_cutJetMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetMass" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cutJetMass	
	histos1D_[ "jet1PrunedMassPassing_cutJetMass" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cutJetMass", "jet1PrunedMassPassing_cutJetMass", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cutJetMass" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cutJetMass" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutJetMass", "jet1SoftDropMassPassing_cutJetMass", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cutJetMass" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutJetMass" ] = fs_->make< TH1D >( "jet1PtPassing_cutJetMass", "jet1PtPassing_cutJetMass", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutJetMass" ]->Sumw2();

	histos2D_[ "jet1PrunedMassjet1PtPassing_cutJetMass" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cutJetMass", "jet1PrunedMassjet1PtPassing_cutJetMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutJetMass" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutJetMass", "jet1SoftDropMassjet1PtPassing_cutJetMass", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetMass" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Denom cutJetPt	
	histos1D_[ "jet1PrunedMassDenom_cutJetPt" ] = fs_->make< TH1D >( "jet1PrunedMassDenom_cutJetPt", "jet1PrunedMassDenom_cutJetPt", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassDenom_cutJetPt" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassDenom_cutJetPt" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutJetPt", "jet1SoftDropMassDenom_cutJetPt", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassDenom_cutJetPt" ]->Sumw2();
	histos1D_[ "jet1PtDenom_cutJetPt" ] = fs_->make< TH1D >( "jet1PtDenom_cutJetPt", "jet1PtDenom_cutJetPt", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutJetPt" ]->Sumw2();

	histos2D_[ "jet1PrunedMassjet1PtDenom_cutJetPt" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtDenom_cutJetPt", "jet1PrunedMassjet1PtDenom_cutJetPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtDenom_cutJetPt" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetPt" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutJetPt", "jet1SoftDropMassjet1PtDenom_cutJetPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetPt" ]->Sumw2();
	///////////////////////////////////////////////////////////////////


	/////// Passing cutJetPt	
	histos1D_[ "jet1PrunedMassPassing_cutJetPt" ] = fs_->make< TH1D >( "jet1PrunedMassPassing_cutJetPt", "jet1PrunedMassPassing_cutJetPt", 60, 0., 600. );
	histos1D_[ "jet1PrunedMassPassing_cutJetPt" ]->Sumw2();
	histos1D_[ "jet1SoftDropMassPassing_cutJetPt" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutJetPt", "jet1SoftDropMassPassing_cutJetPt", 60, 0., 600. );
	histos1D_[ "jet1SoftDropMassPassing_cutJetPt" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutJetPt" ] = fs_->make< TH1D >( "jet1PtPassing_cutJetPt", "jet1PtPassing_cutJetPt", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutJetPt" ]->Sumw2();

	histos2D_[ "jet1PrunedMassjet1PtPassing_cutJetPt" ] = fs_->make< TH2D >( "jet1PrunedMassjet1PtPassing_cutJetPt", "jet1PrunedMassjet1PtPassing_cutJetPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1PrunedMassjet1PtPassing_cutJetPt" ]->Sumw2();
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetPt" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutJetPt", "jet1SoftDropMassjet1PtPassing_cutJetPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetPt" ]->Sumw2();
	///////////////////////////////////////////////////////////////////

}

// ------------ method called once each job just after ending the event loop  ------------
void RUNDijetTriggerEfficiency::endJob() {

}

void RUNDijetTriggerEfficiency::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<double>("cutAK8jetPt", 150.);
	desc.add<double>("cutAK8jet1Pt", 500.);
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
	desc.add<InputTag>("rho", 	InputTag("vertexInfo:rho"));
	desc.add<InputTag>("NPV", 	InputTag("vertexInfo:npv"));
	desc.add<InputTag>("jetArea", 	InputTag("jetsAK8Puppi:jetAK8PuppijetArea"));
	desc.add<InputTag>("jetPt", 	InputTag("jetsAK8Puppi:jetAK8PuppiPt"));
	desc.add<InputTag>("jetEta", 	InputTag("jetsAK8Puppi:jetAK8PuppiEta"));
	desc.add<InputTag>("jetPhi", 	InputTag("jetsAK8Puppi:jetAK8PuppiPhi"));
	desc.add<InputTag>("jetE", 	InputTag("jetsAK8Puppi:jetAK8PuppiE"));
	desc.add<InputTag>("jetPrunedMass", 	InputTag("jetsAK8Puppi:jetAK8PuppiprunedMass"));
	desc.add<InputTag>("jetSoftDropMass", 	InputTag("jetsAK8Puppi:jetAK8PuppisoftDropMass"));
	desc.add<InputTag>("jetCSVv2", 	InputTag("jetsAK8Puppi:jetAK8PuppiCSVv2"));
	// JetID
	desc.add<InputTag>("jecFactor", 		InputTag("jetsAK8Puppi:jetAK8PuppijecFactor0"));
	desc.add<InputTag>("neutralHadronEnergyFrac", 	InputTag("jetsAK8Puppi:jetAK8PuppineutralHadronEnergyFrac"));
	desc.add<InputTag>("neutralEmEnergyFrac", 	InputTag("jetsAK8Puppi:jetAK8PuppineutralEmEnergyFrac"));
	desc.add<InputTag>("chargedEmEnergyFrac", 	InputTag("jetsAK8Puppi:jetAK8PuppichargedEmEnergyFrac"));
	desc.add<InputTag>("muonEnergy", 		InputTag("jetsAK8Puppi:jetAK8PuppiMuonEnergy"));
	desc.add<InputTag>("chargedHadronEnergyFrac", 	InputTag("jetsAK8Puppi:jetAK8PuppichargedHadronEnergyFrac"));
	desc.add<InputTag>("neutralMultiplicity",	InputTag("jetsAK8Puppi:jetAK8PuppineutralMultiplicity"));
	desc.add<InputTag>("chargedMultiplicity", 	InputTag("jetsAK8Puppi:jetAK8PuppichargedMultiplicity"));
	// Trigger
	desc.add<InputTag>("triggerPrescale",		InputTag("TriggerUserData:triggerPrescaleTree"));
	desc.add<InputTag>("triggerBit",		InputTag("TriggerUserData:triggerBitTree"));
	desc.add<InputTag>("triggerName",		InputTag("TriggerUserData:triggerNameTree"));
	// Muons
	desc.add<InputTag>("muonPt", 		InputTag("muons:muPt"));
	desc.add<InputTag>("muonEta", 		InputTag("muons:muEta"));
	desc.add<InputTag>("muonPhi", 		InputTag("muons:muPhi"));
	desc.add<InputTag>("muonE", 		InputTag("muons:muE"));
	desc.add<InputTag>("muonIsLoose", 		InputTag("muons:muIsLooseMuon"));
	desc.add<InputTag>("muonIsGlobal", 		InputTag("muons:muIsGlobalMuon"));
	descriptions.addDefault(desc);
}
      
void RUNDijetTriggerEfficiency::beginRun(const Run& iRun, const EventSetup& iSetup){

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

void RUNDijetTriggerEfficiency::endRun(const Run& iRun, const EventSetup& iSetup){
	triggerNamesList.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(RUNDijetTriggerEfficiency);
