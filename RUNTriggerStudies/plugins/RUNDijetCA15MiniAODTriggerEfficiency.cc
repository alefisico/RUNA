// -*- C++ -*-
//
// Package:    RUNA/RUNTriggerEfficiency
// Class:      RUNDijetCA15MiniAODTriggerEfficiency
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
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "RUNA/RUNAnalysis/interface/CommonVariablesStructure.h"

using namespace edm;
using namespace std;

//
// constants, enums and typedefs
//

//
// class declaration
//
class RUNDijetCA15MiniAODTriggerEfficiency : public EDAnalyzer {
	public:
		explicit RUNDijetCA15MiniAODTriggerEfficiency(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~RUNDijetCA15MiniAODTriggerEfficiency();

	private:
		virtual void beginJob() override;
		virtual void analyze(const Event&, const EventSetup&) override;
		virtual void endJob() override;

		//virtual void beginRun(Run const&, EventSetup const&) override;
		//virtual void endRun(Run const&, EventSetup const&) override;
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
		vector<string> triggerPass, triggerNamesList;

		edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
		edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
		edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;
};

//
// static data member definitions
//

//
// constructors and destructor
//
RUNDijetCA15MiniAODTriggerEfficiency::RUNDijetCA15MiniAODTriggerEfficiency(const ParameterSet& iConfig):
	triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
	triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
	triggerEvent_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("hltTrigger"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("recoJets")))
{
	cutAK8jetPt = iConfig.getParameter<double>("cutAK8jetPt");
	cutAK8jet1Pt = iConfig.getParameter<double>("cutAK8jet1Pt");
	cutAK8jet1Mass = iConfig.getParameter<double>("cutAK8jet1Mass");
	baseTrigger = iConfig.getParameter<string>("baseTrigger");
	triggerPass = iConfig.getParameter<vector<string>>("triggerPass");
}


RUNDijetCA15MiniAODTriggerEfficiency::~RUNDijetCA15MiniAODTriggerEfficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNDijetCA15MiniAODTriggerEfficiency::analyze(const Event& iEvent, const EventSetup& iSetup) {

	edm::Handle<edm::TriggerResults> triggerBits;
	edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
	edm::Handle<trigger::TriggerEvent> trigEvent; 
	edm::Handle<pat::JetCollection> jets;

	iEvent.getByToken(triggerBits_, triggerBits);
	iEvent.getByToken(triggerPrescales_, triggerPrescales);
	iEvent.getByToken(triggerEvent_,trigEvent);
	iEvent.getByToken(jetToken_, jets);

	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  	bool ORTriggers = checkORListOfTriggerBitsMiniAOD( names, triggerBits, triggerPrescales, triggerPass, false );
  	bool basedTriggerFired = checkTriggerBitsMiniAOD( names, triggerBits, triggerPrescales, baseTrigger, true );

	/// Applying kinematic, trigger and jet ID
	pat::JetCollection JETS;
	float HT = 0, jet1SoftDropMass = -999;

	for (const pat::Jet &jet : *jets) {

		if ( TMath::Abs( jet.eta() ) > 2.4 ) continue;
		string typeOfJetID = "looseJetID";	// check trigger with looser jet id
		bool idL = jetID( jet.eta(), jet.energy(), jet.jecFactor(0), jet.neutralHadronEnergyFraction(), jet.neutralEmEnergyFraction(), jet.chargedHadronEnergyFraction(), jet.muonEnergy(), jet.chargedEmEnergyFraction(), jet.chargedMultiplicity(), jet.neutralMultiplicity(), typeOfJetID ); 

		
		///// Puppi jets
		if ( jet.pt() > cutAK8jetPt && idL ){
			HT += jet.pt();
			JETS.push_back( jet );
			jet1SoftDropMass = jet.userFloat("ca15PFJetsPuppiSoftDropMass");
		}
	}

	if ( JETS.size() > 0 ) {

		histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ]->Fill( jet1SoftDropMass, JETS[0].pt() );

		if ( basedTriggerFired || ORTriggers ) {
			if ( basedTriggerFired && ORTriggers ) {
				histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOneAndTwo" ]->Fill( jet1SoftDropMass, JETS[0].pt() );
			} else if ( basedTriggerFired ) {
				histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOne" ]->Fill( jet1SoftDropMass, JETS[0].pt() );
			} else if ( ORTriggers ) {
				histos2D_[ "jet1SoftDropMassPt_cutJet_triggerTwo" ]->Fill( jet1SoftDropMass, JETS[0].pt() );
			}
		} else {
			histos2D_[ "jet1SoftDropMassPt_cutJet_noTriggerOneOrTwo" ]->Fill( jet1SoftDropMass, JETS[0].pt() );
		}


		if ( basedTriggerFired ) {
			histos1D_[ "jet1SoftDropMassDenom_cutJet" ]->Fill( jet1SoftDropMass  );
			histos1D_[ "jet1PtDenom_cutJet" ]->Fill( JETS[0].pt()   );

			histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJet" ]->Fill( jet1SoftDropMass, JETS[0].pt() );

			if ( ORTriggers ){
				histos1D_[ "jet1SoftDropMassPassing_cutJet" ]->Fill( jet1SoftDropMass  );
				histos1D_[ "jet1PtPassing_cutJet" ]->Fill( JETS[0].pt()   );

				histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJet" ]->Fill( jet1SoftDropMass, JETS[0].pt() );

			}

			if ( jet1SoftDropMass > cutAK8jet1Mass ) {
				histos1D_[ "jet1SoftDropMassDenom_cutJetMass" ]->Fill( jet1SoftDropMass  );
				histos1D_[ "jet1PtDenom_cutJetMass" ]->Fill( JETS[0].pt()   );

				histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetMass" ]->Fill( jet1SoftDropMass, JETS[0].pt() );

				if ( ORTriggers ){
					histos1D_[ "jet1SoftDropMassPassing_cutJetMass" ]->Fill( jet1SoftDropMass  );
					histos1D_[ "jet1PtPassing_cutJetMass" ]->Fill( JETS[0].pt()   );

					histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetMass" ]->Fill( jet1SoftDropMass, JETS[0].pt() );
				}
			}

			if ( JETS[0].pt() > cutAK8jet1Pt ) {
				histos1D_[ "jet1SoftDropMassDenom_cutJetPt" ]->Fill( jet1SoftDropMass  );
				histos1D_[ "jet1PtDenom_cutJetPt" ]->Fill( JETS[0].pt()   );

				histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetPt" ]->Fill( jet1SoftDropMass, JETS[0].pt() );

				if ( ORTriggers ){
					histos1D_[ "jet1SoftDropMassPassing_cutJetPt" ]->Fill( jet1SoftDropMass  );
					histos1D_[ "jet1PtPassing_cutJetPt" ]->Fill( JETS[0].pt()   );

					histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetPt" ]->Fill( jet1SoftDropMass, JETS[0].pt() );
				}
			}

		}
	}



}


// ------------ method called once each job just before starting event loop  ------------
void RUNDijetCA15MiniAODTriggerEfficiency::beginJob() {

	//////// test plots


	histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_noTrigger", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_noTrigger" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_noTrigger", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOne" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_triggerOne", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_triggerTwo", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_triggerOneAndTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_triggerOneAndTwo", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	histos2D_[ "jet1SoftDropMassPt_cutJet_noTriggerOneOrTwo" ] = fs_->make< TH2D >( "jet1SoftDropMassPt_cutJet_noTriggerOneOrTwo", "jet1SoftDropMassPt", 60, 0., 600., 150, 0., 1500.);
	///////////////////////////////////////////////////////////////////

	/////// Denom cutJet	
	histos1D_[ "jet1SoftDropMassDenom_cutJet" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutJet", "jet1SoftDropMassDenom_cutJet", 60, 0., 600. );
	histos1D_[ "jet1PtDenom_cutJet" ] = fs_->make< TH1D >( "jet1PtDenom_cutJet", "jet1PtDenom_cutJet", 150, 0., 1500. );

	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJet" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutJet", "jet1SoftDropMassjet1PtDenom_cutJet", 60, 0., 600., 150, 0., 1500.);
	///////////////////////////////////////////////////////////////////


	/////// Passing cutJet	
	histos1D_[ "jet1SoftDropMassPassing_cutJet" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutJet", "jet1SoftDropMassPassing_cutJet", 60, 0., 600. );
	histos1D_[ "jet1PtPassing_cutJet" ] = fs_->make< TH1D >( "jet1PtPassing_cutJet", "jet1PtPassing_cutJet", 150, 0., 1500. );

	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJet" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutJet", "jet1SoftDropMassjet1PtPassing_cutJet", 60, 0., 600., 150, 0., 1500.);
	///////////////////////////////////////////////////////////////////


	/////// Denom cutJetMass	
	histos1D_[ "jet1SoftDropMassDenom_cutJetMass" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutJetMass", "jet1SoftDropMassDenom_cutJetMass", 60, 0., 600. );
	histos1D_[ "jet1PtDenom_cutJetMass" ] = fs_->make< TH1D >( "jet1PtDenom_cutJetMass", "jet1PtDenom_cutJetMass", 150, 0., 1500. );

	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutJetMass", "jet1SoftDropMassjet1PtDenom_cutJetMass", 60, 0., 600., 150, 0., 1500.);
	///////////////////////////////////////////////////////////////////


	/////// Passing cutJetMass	
	histos1D_[ "jet1SoftDropMassPassing_cutJetMass" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutJetMass", "jet1SoftDropMassPassing_cutJetMass", 60, 0., 600. );
	histos1D_[ "jet1PtPassing_cutJetMass" ] = fs_->make< TH1D >( "jet1PtPassing_cutJetMass", "jet1PtPassing_cutJetMass", 150, 0., 1500. );

	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetMass" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutJetMass", "jet1SoftDropMassjet1PtPassing_cutJetMass", 60, 0., 600., 150, 0., 1500.);
	///////////////////////////////////////////////////////////////////


	/////// Denom cutJetPt	
	histos1D_[ "jet1SoftDropMassDenom_cutJetPt" ] = fs_->make< TH1D >( "jet1SoftDropMassDenom_cutJetPt", "jet1SoftDropMassDenom_cutJetPt", 60, 0., 600. );
	histos1D_[ "jet1PtDenom_cutJetPt" ] = fs_->make< TH1D >( "jet1PtDenom_cutJetPt", "jet1PtDenom_cutJetPt", 150, 0., 1500. );

	histos2D_[ "jet1SoftDropMassjet1PtDenom_cutJetPt" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtDenom_cutJetPt", "jet1SoftDropMassjet1PtDenom_cutJetPt", 60, 0., 600., 150, 0., 1500.);
	///////////////////////////////////////////////////////////////////


	/////// Passing cutJetPt	
	histos1D_[ "jet1SoftDropMassPassing_cutJetPt" ] = fs_->make< TH1D >( "jet1SoftDropMassPassing_cutJetPt", "jet1SoftDropMassPassing_cutJetPt", 60, 0., 600. );
	histos1D_[ "jet1PtPassing_cutJetPt" ] = fs_->make< TH1D >( "jet1PtPassing_cutJetPt", "jet1PtPassing_cutJetPt", 150, 0., 1500. );

	histos2D_[ "jet1SoftDropMassjet1PtPassing_cutJetPt" ] = fs_->make< TH2D >( "jet1SoftDropMassjet1PtPassing_cutJetPt", "jet1SoftDropMassjet1PtPassing_cutJetPt", 60, 0., 600., 150, 0., 1500.);
	///////////////////////////////////////////////////////////////////
	
	///// Sumw2 all the histos
	for( auto const& histo : histos1D_ ) histos1D_[ histo.first ]->Sumw2();
	for( auto const& histo : histos2D_ ) histos2D_[ histo.first ]->Sumw2();


}

// ------------ method called once each job just after ending the event loop  ------------
void RUNDijetCA15MiniAODTriggerEfficiency::endJob() {

}

void RUNDijetCA15MiniAODTriggerEfficiency::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<double>("cutAK8jetPt", 150.);
	desc.add<double>("cutAK8jet1Pt", 500.);
	desc.add<double>("cutAK8jet1Mass", 60.);
	desc.add<string>("baseTrigger", "HLT_PFHT475");
	vector<string> HLTPass;
	HLTPass.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");
	desc.add<vector<string>>("triggerPass",	HLTPass);
	desc.add<InputTag>("bits", 	InputTag("TriggerResults", "", "HLT"));
	desc.add<InputTag>("prescales", 	InputTag("patTrigger"));
	desc.add<InputTag>("hltTrigger", 	InputTag("hltTriggerSummaryAOD","","HLT"));
	desc.add<InputTag>("recoJets", 	InputTag("slimmedJetsAK8"));

	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RUNDijetCA15MiniAODTriggerEfficiency);
