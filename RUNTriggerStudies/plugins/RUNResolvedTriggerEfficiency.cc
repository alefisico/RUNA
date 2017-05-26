// -*- C++ -*-
//
// Package:    RUNA/RUNTriggerEfficiency
// Class:      RUNResolvedTriggerEfficiency
// Original Author:  Alejandro Gomez Espinosa
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
class RUNResolvedTriggerEfficiency : public EDAnalyzer {
	public:
		explicit RUNResolvedTriggerEfficiency(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~RUNResolvedTriggerEfficiency();

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

		TString baseTrigger;
		double cutAK4jetPt;
		double cutAK4jet4Pt;
		double cutAK4HT;
		vector<string> triggerPass, triggerNamesList;

		ULong64_t event = 0;
		int numJets = 0, numPV = 0;
		unsigned int lumi = 0, run=0;

		EDGetTokenT<vector<float>> jetPt_;
		EDGetTokenT<vector<float>> jetEta_;
		EDGetTokenT<vector<float>> jetPhi_;
		EDGetTokenT<vector<float>> jetE_;
		EDGetTokenT<vector<float>> jetCSV_;
		EDGetTokenT<vector<float>> jetCSVV1_;
		EDGetTokenT<int> NPV_;
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
RUNResolvedTriggerEfficiency::RUNResolvedTriggerEfficiency(const ParameterSet& iConfig):
	jetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPt"))),
	jetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetEta"))),
	jetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPhi"))),
	jetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetE"))),
	jetCSV_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCSV"))),
	jetCSVV1_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCSVV1"))),
	NPV_(consumes<int>(iConfig.getParameter<InputTag>("NPV"))),
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
	baseTrigger = iConfig.getParameter<string>("baseTrigger");
	cutAK4jetPt = iConfig.getParameter<double>("cutAK4jetPt");
	cutAK4jet4Pt = iConfig.getParameter<double>("cutAK4jet4Pt");
	cutAK4HT = iConfig.getParameter<double>("cutAK4HT");
	triggerPass = iConfig.getParameter<vector<string>>("triggerPass");
}


RUNResolvedTriggerEfficiency::~RUNResolvedTriggerEfficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNResolvedTriggerEfficiency::analyze(const Event& iEvent, const EventSetup& iSetup) {

	Handle<vector<float> > jetPt;
	iEvent.getByToken(jetPt_, jetPt);

	Handle<vector<float> > jetEta;
	iEvent.getByToken(jetEta_, jetEta);

	Handle<vector<float> > jetPhi;
	iEvent.getByToken(jetPhi_, jetPhi);

	Handle<vector<float> > jetE;
	iEvent.getByToken(jetE_, jetE);

	Handle<vector<float> > jetCSV;
	iEvent.getByToken(jetCSV_, jetCSV);

	Handle<vector<float> > jetCSVV1;
	iEvent.getByToken(jetCSVV1_, jetCSVV1);

	Handle<int> NPV;
	iEvent.getByToken(NPV_, NPV);

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
	vector< TLorentzVector > JETS;
	vector< myJet > myJETS;
	float HT = 0;
	for (size_t i = 0; i < jetPt->size(); i++) {

		if( TMath::Abs( (*jetEta)[i] ) > 2.4 ) continue;

		string typeOfJetID = "looseJetID";	// check trigger with looser jet id
		bool idL = jetID( (*jetEta)[i], (*jetE)[i], (*jecFactor)[i], (*neutralHadronEnergyFrac)[i], (*neutralEmEnergyFrac)[i], (*chargedHadronEnergyFrac)[i], (*muonEnergy)[i], (*chargedEmEnergyFrac)[i], (*chargedMultiplicity)[i], (*neutralMultiplicity)[i], typeOfJetID ); 

		if( ( (*jetPt)[i] > cutAK4jetPt ) && idL ) { 
			//LogWarning("jetInfo") << i << " " << (*jetPt)[i] << " " << (*jetEta)[i] << " " << (*jetPhi)[i] << " " << (*jetMass)[i];
			HT += (*jetPt)[i];
			TLorentzVector tmpJet;
			tmpJet.SetPtEtaPhiE( (*jetPt)[i], (*jetEta)[i], (*jetPhi)[i], (*jetE)[i] );
			JETS.push_back( tmpJet );

			myJet tmpMyJet;
			tmpMyJet.p4 = tmpJet;
			myJETS.push_back( tmpMyJet );

		}
	}


	if ( JETS.size() > 3 ) {

		if ( basedTriggerFired ) {

			vector< myJet > myPairedJets;
			myPairedJets = pairing( myJETS, true );
			double mass1 = ( myPairedJets[0].p4 + myPairedJets[1].p4 ).M();
			double mass2 = ( myPairedJets[2].p4 + myPairedJets[3].p4 ).M();
			double massAve = ( mass1 + mass2 ) / 2;

			histos1D_[ "jet1PtDenom_cut4Jet" ]->Fill( JETS[0].Pt() );
			histos1D_[ "jet2PtDenom_cut4Jet" ]->Fill( JETS[1].Pt() );
			histos1D_[ "jet3PtDenom_cut4Jet" ]->Fill( JETS[2].Pt() );
			histos1D_[ "jet4PtDenom_cut4Jet" ]->Fill( JETS[3].Pt() );
			histos1D_[ "HTDenom_cut4Jet" ]->Fill( HT  );
			histos2D_[ "jet4PtHTDenom_cut4Jet" ]->Fill( JETS[3].Pt(), HT );
			histos1D_[ "massAveDenom_cut4Jet" ]->Fill( massAve );

			if ( ORTriggers ){
				histos1D_[ "jet1PtPassing_cut4Jet" ]->Fill( JETS[0].Pt() );
				histos1D_[ "jet2PtPassing_cut4Jet" ]->Fill( JETS[1].Pt() );
				histos1D_[ "jet3PtPassing_cut4Jet" ]->Fill( JETS[2].Pt() );
				histos1D_[ "jet4PtPassing_cut4Jet" ]->Fill( JETS[3].Pt() );
				histos1D_[ "HTPassing_cut4Jet" ]->Fill( HT  );
				histos2D_[ "jet4PtHTPassing_cut4Jet" ]->Fill( JETS[3].Pt(), HT );
				histos1D_[ "massAvePassing_cut4Jet" ]->Fill( massAve );
			}

			if ( JETS[3].Pt() > cutAK4jet4Pt ) {
				histos1D_[ "jet1PtDenom_cutJet4Pt" ]->Fill( JETS[0].Pt() );
				histos1D_[ "jet2PtDenom_cutJet4Pt" ]->Fill( JETS[1].Pt() );
				histos1D_[ "jet3PtDenom_cutJet4Pt" ]->Fill( JETS[2].Pt() );
				histos1D_[ "jet4PtDenom_cutJet4Pt" ]->Fill( JETS[3].Pt() );
				histos1D_[ "HTDenom_cutJet4Pt" ]->Fill( HT  );
				histos2D_[ "jet4PtHTDenom_cutJet4Pt" ]->Fill( JETS[3].Pt(), HT );

				if ( ORTriggers ){
					histos1D_[ "jet1PtPassing_cutJet4Pt" ]->Fill( JETS[0].Pt() );
					histos1D_[ "jet2PtPassing_cutJet4Pt" ]->Fill( JETS[1].Pt() );
					histos1D_[ "jet3PtPassing_cutJet4Pt" ]->Fill( JETS[2].Pt() );
					histos1D_[ "jet4PtPassing_cutJet4Pt" ]->Fill( JETS[3].Pt() );
					histos1D_[ "HTPassing_cutJet4Pt" ]->Fill( HT  );
					histos2D_[ "jet4PtHTPassing_cutJet4Pt" ]->Fill( JETS[3].Pt(), HT );
				}
			}

			if ( HT > cutAK4HT ) {
				histos1D_[ "jet1PtDenom_cutHT" ]->Fill( JETS[0].Pt() );
				histos1D_[ "jet2PtDenom_cutHT" ]->Fill( JETS[1].Pt() );
				histos1D_[ "jet3PtDenom_cutHT" ]->Fill( JETS[2].Pt() );
				histos1D_[ "jet4PtDenom_cutHT" ]->Fill( JETS[3].Pt() );
				histos1D_[ "HTDenom_cutHT" ]->Fill( HT  );
				histos2D_[ "jet4PtHTDenom_cutHT" ]->Fill( JETS[3].Pt(), HT );
				histos1D_[ "massAveDenom_cutHT" ]->Fill( massAve );

				if ( ORTriggers ){
					histos1D_[ "jet1PtPassing_cutHT" ]->Fill( JETS[0].Pt() );
					histos1D_[ "jet2PtPassing_cutHT" ]->Fill( JETS[1].Pt() );
					histos1D_[ "jet3PtPassing_cutHT" ]->Fill( JETS[2].Pt() );
					histos1D_[ "jet4PtPassing_cutHT" ]->Fill( JETS[3].Pt() );
					histos1D_[ "HTPassing_cutHT" ]->Fill( HT  );
					histos2D_[ "jet4PtHTPassing_cutHT" ]->Fill( JETS[3].Pt(), HT );
					histos1D_[ "massAvePassing_cutHT" ]->Fill( massAve );
				}

			}
		}
	}
	JETS.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void RUNResolvedTriggerEfficiency::beginJob() {

	histos1D_[ "HTDenom_cut4Jet" ] = fs_->make< TH1D >( "HTDenom_cut4Jet", "HTDenom_cut4Jet", 500, 0., 5000. );
	histos1D_[ "HTDenom_cut4Jet" ]->Sumw2();
	histos1D_[ "HTPassing_cut4Jet" ] = fs_->make< TH1D >( "HTPassing_cut4Jet", "HTPassing_cut4Jet", 500, 0., 5000. );
	histos1D_[ "HTPassing_cut4Jet" ]->Sumw2();

	histos1D_[ "jet1PtDenom_cut4Jet" ] = fs_->make< TH1D >( "jet1PtDenom_cut4Jet", "jet1PtDenom_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cut4Jet" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cut4Jet" ] = fs_->make< TH1D >( "jet1PtPassing_cut4Jet", "jet1PtPassing_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cut4Jet" ]->Sumw2();

	histos1D_[ "jet2PtDenom_cut4Jet" ] = fs_->make< TH1D >( "jet2PtDenom_cut4Jet", "jet2PtDenom_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cut4Jet" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cut4Jet" ] = fs_->make< TH1D >( "jet2PtPassing_cut4Jet", "jet2PtPassing_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cut4Jet" ]->Sumw2();

	histos1D_[ "jet3PtDenom_cut4Jet" ] = fs_->make< TH1D >( "jet3PtDenom_cut4Jet", "jet3PtDenom_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet3PtDenom_cut4Jet" ]->Sumw2();
	histos1D_[ "jet3PtPassing_cut4Jet" ] = fs_->make< TH1D >( "jet3PtPassing_cut4Jet", "jet3PtPassing_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet3PtPassing_cut4Jet" ]->Sumw2();

	histos1D_[ "jet4PtDenom_cut4Jet" ] = fs_->make< TH1D >( "jet4PtDenom_cut4Jet", "jet4PtDenom_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet4PtDenom_cut4Jet" ]->Sumw2();
	histos1D_[ "jet4PtPassing_cut4Jet" ] = fs_->make< TH1D >( "jet4PtPassing_cut4Jet", "jet4PtPassing_cut4Jet", 150, 0., 1500. );
	histos1D_[ "jet4PtPassing_cut4Jet" ]->Sumw2();

	histos1D_[ "massAveDenom_cut4Jet" ] = fs_->make< TH1D >( "massAveDenom_cut4Jet", "massAveDenom_cut4Jet", 1500, 0., 1500. );
	histos1D_[ "massAveDenom_cut4Jet" ]->Sumw2();
	histos1D_[ "massAvePassing_cut4Jet" ] = fs_->make< TH1D >( "massAvePassing_cut4Jet", "massAvePassing_cut4Jet", 1500, 0., 1500. );
	histos1D_[ "massAvePassing_cut4Jet" ]->Sumw2();

	histos2D_[ "jet4PtHTDenom_cut4Jet" ] = fs_->make< TH2D >( "jet4PtHTDenom_cut4Jet", "HT vs 4th Leading Jet Pt", 150, 0., 1500., 500, 0., 5000.);
	histos2D_[ "jet4PtHTDenom_cut4Jet" ]->Sumw2();

	histos2D_[ "jet4PtHTPassing_cut4Jet" ] = fs_->make< TH2D >( "jet4PtHTPassing_cut4Jet", "HT vs 4th Leading Jet Pt", 150, 0., 1500., 500, 0., 5000.);
	histos2D_[ "jet4PtHTPassing_cut4Jet" ]->Sumw2();

	histos1D_[ "HTDenom_cutJet4Pt" ] = fs_->make< TH1D >( "HTDenom_cutJet4Pt", "HTDenom_cutJet4Pt", 500, 0., 5000. );
	histos1D_[ "HTDenom_cutJet4Pt" ]->Sumw2();
	histos1D_[ "HTPassing_cutJet4Pt" ] = fs_->make< TH1D >( "HTPassing_cutJet4Pt", "HTPassing_cutJet4Pt", 500, 0., 5000. );
	histos1D_[ "HTPassing_cutJet4Pt" ]->Sumw2();

	histos1D_[ "jet1PtDenom_cutJet4Pt" ] = fs_->make< TH1D >( "jet1PtDenom_cutJet4Pt", "jet1PtDenom_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutJet4Pt" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutJet4Pt" ] = fs_->make< TH1D >( "jet1PtPassing_cutJet4Pt", "jet1PtPassing_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutJet4Pt" ]->Sumw2();

	histos1D_[ "jet2PtDenom_cutJet4Pt" ] = fs_->make< TH1D >( "jet2PtDenom_cutJet4Pt", "jet2PtDenom_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cutJet4Pt" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cutJet4Pt" ] = fs_->make< TH1D >( "jet2PtPassing_cutJet4Pt", "jet2PtPassing_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cutJet4Pt" ]->Sumw2();

	histos1D_[ "jet3PtDenom_cutJet4Pt" ] = fs_->make< TH1D >( "jet3PtDenom_cutJet4Pt", "jet3PtDenom_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet3PtDenom_cutJet4Pt" ]->Sumw2();
	histos1D_[ "jet3PtPassing_cutJet4Pt" ] = fs_->make< TH1D >( "jet3PtPassing_cutJet4Pt", "jet3PtPassing_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet3PtPassing_cutJet4Pt" ]->Sumw2();

	histos1D_[ "jet4PtDenom_cutJet4Pt" ] = fs_->make< TH1D >( "jet4PtDenom_cutJet4Pt", "jet4PtDenom_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet4PtDenom_cutJet4Pt" ]->Sumw2();
	histos1D_[ "jet4PtPassing_cutJet4Pt" ] = fs_->make< TH1D >( "jet4PtPassing_cutJet4Pt", "jet4PtPassing_cutJet4Pt", 150, 0., 1500. );
	histos1D_[ "jet4PtPassing_cutJet4Pt" ]->Sumw2();

	histos2D_[ "jet4PtHTDenom_cutJet4Pt" ] = fs_->make< TH2D >( "jet4PtHTDenom_cutJet4Pt", "HT vs 4th Leading Jet Pt", 150, 0., 1500., 500, 0., 5000.);
	histos2D_[ "jet4PtHTDenom_cutJet4Pt" ]->Sumw2();

	histos2D_[ "jet4PtHTPassing_cutJet4Pt" ] = fs_->make< TH2D >( "jet4PtHTPassing_cutJet4Pt", "HT vs 4th Leading Jet Pt", 150, 0., 1500., 500, 0., 5000.);
	histos2D_[ "jet4PtHTPassing_cutJet4Pt" ]->Sumw2();


	histos1D_[ "HTDenom_cutHT" ] = fs_->make< TH1D >( "HTDenom_cutHT", "HTDenom_cutHT", 500, 0., 5000. );
	histos1D_[ "HTDenom_cutHT" ]->Sumw2();
	histos1D_[ "HTPassing_cutHT" ] = fs_->make< TH1D >( "HTPassing_cutHT", "HTPassing_cutHT", 500, 0., 5000. );
	histos1D_[ "HTPassing_cutHT" ]->Sumw2();

	histos1D_[ "jet1PtDenom_cutHT" ] = fs_->make< TH1D >( "jet1PtDenom_cutHT", "jet1PtDenom_cutHT", 150, 0., 1500. );
	histos1D_[ "jet1PtDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet1PtPassing_cutHT" ] = fs_->make< TH1D >( "jet1PtPassing_cutHT", "jet1PtPassing_cutHT", 150, 0., 1500. );
	histos1D_[ "jet1PtPassing_cutHT" ]->Sumw2();

	histos1D_[ "jet2PtDenom_cutHT" ] = fs_->make< TH1D >( "jet2PtDenom_cutHT", "jet2PtDenom_cutHT", 150, 0., 1500. );
	histos1D_[ "jet2PtDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet2PtPassing_cutHT" ] = fs_->make< TH1D >( "jet2PtPassing_cutHT", "jet2PtPassing_cutHT", 150, 0., 1500. );
	histos1D_[ "jet2PtPassing_cutHT" ]->Sumw2();

	histos1D_[ "jet3PtDenom_cutHT" ] = fs_->make< TH1D >( "jet3PtDenom_cutHT", "jet3PtDenom_cutHT", 150, 0., 1500. );
	histos1D_[ "jet3PtDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet3PtPassing_cutHT" ] = fs_->make< TH1D >( "jet3PtPassing_cutHT", "jet3PtPassing_cutHT", 150, 0., 1500. );
	histos1D_[ "jet3PtPassing_cutHT" ]->Sumw2();

	histos1D_[ "jet4PtDenom_cutHT" ] = fs_->make< TH1D >( "jet4PtDenom_cutHT", "jet4PtDenom_cutHT", 150, 0., 1500. );
	histos1D_[ "jet4PtDenom_cutHT" ]->Sumw2();
	histos1D_[ "jet4PtPassing_cutHT" ] = fs_->make< TH1D >( "jet4PtPassing_cutHT", "jet4PtPassing_cutHT", 150, 0., 1500. );
	histos1D_[ "jet4PtPassing_cutHT" ]->Sumw2();

	histos2D_[ "jet4PtHTDenom_cutHT" ] = fs_->make< TH2D >( "jet4PtHTDenom_cutHT", "jet4PtHTDenom_cutHT", 150, 0., 1500., 500, 0., 5000.);
	histos2D_[ "jet4PtHTDenom_cutHT" ]->Sumw2();

	histos2D_[ "jet4PtHTPassing_cutHT" ] = fs_->make< TH2D >( "jet4PtHTPassing_cutHT", "jet4PtHTPassing_cutHT", 150, 0., 1500., 500, 0., 5000.);
	histos2D_[ "jet4PtHTPassing_cutHT" ]->Sumw2();

	histos1D_[ "massAveDenom_cutHT" ] = fs_->make< TH1D >( "massAveDenom_cutHT", "massAveDenom_cutHT", 1500, 0., 1500. );
	histos1D_[ "massAveDenom_cutHT" ]->Sumw2();
	histos1D_[ "massAvePassing_cutHT" ] = fs_->make< TH1D >( "massAvePassing_cutHT", "massAvePassing_cutHT", 1500, 0., 1500. );
	histos1D_[ "massAvePassing_cutHT" ]->Sumw2();


}

// ------------ method called once each job just after ending the event loop  ------------
void RUNResolvedTriggerEfficiency::endJob() {

}

void RUNResolvedTriggerEfficiency::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<double>("cutAK4jetPt", 50);
	desc.add<double>("cutAK4jet4Pt", 80);
	desc.add<double>("cutAK4HT", 800);
	desc.add<string>("baseTrigger", "HLT_PFHT475");
	vector<string> HLTPass;
	HLTPass.push_back("HLT_PFHT800");
	desc.add<vector<string>>("triggerPass",	HLTPass);

	desc.add<InputTag>("Lumi", 	InputTag("eventInfo:evtInfoLumiBlock"));
	desc.add<InputTag>("Run", 	InputTag("eventInfo:evtInfoRunNumber"));
	desc.add<InputTag>("Event", 	InputTag("eventInfo:evtInfoEventNumber"));
	desc.add<InputTag>("NPV", 	InputTag("eventUserData:npv"));
	desc.add<InputTag>("jetPt", 	InputTag("jetsAK4CHS:jetAK4CHSPt"));
	desc.add<InputTag>("jetEta", 	InputTag("jetsAK4CHS:jetAK4CHSEta"));
	desc.add<InputTag>("jetPhi", 	InputTag("jetsAK4CHS:jetAK4CHSPhi"));
	desc.add<InputTag>("jetE", 	InputTag("jetsAK4CHS:jetAK4CHSE"));
	desc.add<InputTag>("jetCSV", 	InputTag("jetsAK4CHS:jetAK4CHSCSV"));
	desc.add<InputTag>("jetCSVV1", 	InputTag("jetsAK4CHS:jetAK4CHSCSVV1"));
	// JetID
	desc.add<InputTag>("jecFactor", 		InputTag("jetsAK4CHS:jetAK4CHSjecFactor0"));
	desc.add<InputTag>("neutralHadronEnergyFrac", 	InputTag("jetsAK4CHS:jetAK4CHSneutralHadronEnergyFrac"));
	desc.add<InputTag>("neutralEmEnergyFrac", 		InputTag("jetsAK4CHS:jetAK4CHSneutralEmEnergyFrac"));
	desc.add<InputTag>("chargedEmEnergyFrac", 		InputTag("jetsAK4CHS:jetAK4CHSchargedEmEnergyFrac"));
	desc.add<InputTag>("muonEnergy", 		InputTag("jetsAK4CHS:jetAK4CHSMuonEnergy"));
	desc.add<InputTag>("chargedHadronEnergyFrac", 	InputTag("jetsAK4CHS:jetAK4CHSchargedHadronEnergyFrac"));
	desc.add<InputTag>("neutralMultiplicity",	InputTag("jetsAK4CHS:jetAK4CHSneutralMultiplicity"));
	desc.add<InputTag>("chargedMultiplicity", 	InputTag("jetsAK4CHS:jetAK4CHSchargedMultiplicity"));
	// Trigger
	desc.add<InputTag>("triggerPrescale",		InputTag("TriggerUserData:triggerPrescaleTree"));
	desc.add<InputTag>("triggerBit",		InputTag("TriggerUserData:triggerBitTree"));
	desc.add<InputTag>("triggerName",		InputTag("TriggerUserData:triggerNameTree"));
	descriptions.addDefault(desc);
}
      
void RUNResolvedTriggerEfficiency::beginRun(const Run& iRun, const EventSetup& iSetup){

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
void RUNResolvedTriggerEfficiency::endRun(const Run& iRun, const EventSetup& iSetup){
	triggerNamesList.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(RUNResolvedTriggerEfficiency);
