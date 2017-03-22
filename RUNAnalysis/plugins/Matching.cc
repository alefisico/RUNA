// -*- C++ -*-
//
// Package:    RUNA/RUNtuples
// Class:      Matching
// 
/**\class Matching Matching.cc RUNA/RUNtuples/plugins/Matching.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  alejandro gomez
//         Created:  Fri, 31 Oct 2014 16:15:59 GMT
//
//


// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TH1D.h"
#include <TLorentzVector.h>
//
// class declaration
//
using namespace edm;
using namespace reco;
using namespace std;


class Matching : public edm::EDAnalyzer {
	public:
		explicit Matching(const edm::ParameterSet&);
		~Matching();

		//static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<pat::JetCollection> AK4jets_;
		edm::EDGetTokenT<pat::JetCollection> AK8jets_;
		edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
		int particle1, particle2, particle3;

		int oneBoostedP1 = 0;
		int oneBoostedNoResolved = 0;
		int oneResolvedP1 = 0;
		int boostedAndResolvedP1 = 0;
		int twoBoostedP1 = 0;
		int fourResolvedP1 = 0;
		int fourResolvedwithtwoBoostedP1 = 0;
		int noBoostedNoResolved = 0;
		int none = 0;
		std::map< std::string, TH1D* > histos1D_;
		//std::map< std::size_t , reco::Candidate > decayStory_;
};

//
// constants, enums and typedefs
//
/*typedef struct {
	bool pass;
	double deltaR;
	int indexJet;
	pat::Jet matchJet;
} matched;*/

typedef struct {
	pat::JetCollection AK8matchedJets;
	pat::JetCollection AK4matchedJets;
	reco::CandidateCollection daughters;
	TLorentzVector genPartP4;
	int genPartId;
} fullParentInfo;

//
// static data member definitions
//

//
// constructors and destructor
//
Matching::Matching(const edm::ParameterSet& iConfig):
    AK4jets_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("AK4jets"))),
    AK8jets_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("AK8jets"))),
    genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles")))
{
    particle1 = iConfig.getParameter<int>("particle1");
    particle2 = iConfig.getParameter<int>("particle2");
    particle3 = iConfig.getParameter<int>("particle3");

}


Matching::~Matching()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
//Check recursively if any ancestor of particle is the given one
bool isAncestor(reco::Candidate & ancestor, const reco::Candidate * particle) {

	//particle is already the ancestor
	//edm::LogWarning("testing") << ancestor->pdgId() << " " << ancestor->pt() << " " << particle->pdgId() << " " << particle->pt();
        if( ( ancestor.pdgId() == particle->pdgId() ) && ( ancestor.mass() == particle->mass() ) ) {  
		//edm::LogWarning("is ancestor") << ancestor->pdgId() << " " << ancestor->status() << " " << particle->mother()->pdgId() << " " << particle->pdgId();
		if ( ( ancestor.pt() == particle->pt() ) && ( particle->status() == 22 ) ) return true;
	}

	//otherwise loop on mothers, if any and return true if the ancestor is found
	if( particle->mother() != nullptr ) {
		if( isAncestor( ancestor, particle->mother()) ) return true;
	}
	//if we did not return yet, then particle and ancestor are not relatives
	return false;
}

pat::Jet checkDeltaR(reco::Candidate & p1, Handle<pat::JetCollection> jets, double minDeltaR, TH1D * allDeltaR, TH1D * histoMinDeltaR){

	pat::Jet matchedJet;
	double deltaR = 99999;

	for( unsigned int j=0; j<jets->size(); j++ ) {
		const pat::Jet & p2 = (*jets)[j];
		//double tmpdeltaR2 = reco::deltaR2( p1.rapidity(), p1.phi(), p2.rapidity(), p2.phi() );
		double tmpdeltaR = reco::deltaR2( p1.eta(), p1.phi(), p2.eta(), p2.phi() );
		//TLorentzVector tmp1, tmp2;
		//tmp1.SetPtEtaPhiE( p1.pt(), p1.eta(), p1.phi(), p1.energy() );
		//tmp2.SetPtEtaPhiE( p2.pt(), p2.eta(), p2.phi(), p2.energy() );
		//double tmpdeltaR = tmp1.DeltaR( tmp2 );
		//double tmpdeltaR3 = TMath::Sqrt( TMath::Power( (p1.eta()-p2.eta()), 2) + TMath::Power( (p1.phi()-p2.phi()), 2) );
		allDeltaR->Fill( tmpdeltaR );
		//edm::LogWarning("calc deltaR") << j << " "  << tmpdeltaR << " " << p2.pt();
		if( tmpdeltaR < deltaR ) {
			deltaR = tmpdeltaR;
			if( deltaR < minDeltaR ) { 
				matchedJet = p2;
				//edm::LogWarning("final deltaR") << j << " "  << tmpdeltaR << " " << matchedJet.pt();
			}
		}
	}
	//edm::LogWarning("deltaR") << deltaR; // << " " << ind ;
	if( minDeltaR < 999 ) histoMinDeltaR->Fill( minDeltaR );
	return matchedJet; 
}

reco::CandidateCollection checkDaughters( reco::Candidate & p1, reco::CandidateCollection finalParticlesCollection ){

	reco::CandidateCollection daughters;
	for( auto & fp : finalParticlesCollection ) {
		const reco::Candidate * finalMother = fp.mother();
		if( isAncestor( p1, finalMother ) ) daughters.push_back( fp ); // LogWarning("Particle found 1") << jp1->pdgId() << " " << fp.pdgId(); }
	}
	return daughters;
}

// ------------ method called for each event  ------------
void
Matching::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

	edm::Handle<pat::JetCollection> AK4jets;
	iEvent.getByToken(AK4jets_, AK4jets);
	
	edm::Handle<pat::JetCollection> AK8jets;
	iEvent.getByToken(AK8jets_, AK8jets);
	
	edm::Handle<reco::GenParticleCollection> genParticles;
	iEvent.getByToken(genParticles_, genParticles);


	reco::CandidateCollection p1Collection, p2Collection, p3Collection, finalParticlesCollection;

	for( size_t i = 0; i < genParticles->size(); i++ ) {

		const reco::Candidate &p = ( *genParticles )[i];

		if( ( TMath::Abs( p.pdgId() ) == particle1 ) && p.status() == 22 ) { 
			p1Collection.push_back( p );
			histos1D_[ "p1FathersPdgId" ]->Fill( TMath::Abs( p.pdgId() ) );
			//LogWarning("mother") << p.pdgId();
		}
		if( ( TMath::Abs( p.pdgId() ) == particle2 ) && p.status() == 22 ) { 
			p2Collection.push_back( p );
			histos1D_[ "p2FathersPdgId" ]->Fill( TMath::Abs( p.pdgId() ) );
		}
		if( ( TMath::Abs( p.pdgId() ) == particle3 ) && p.status() == 22 ) {
			p3Collection.push_back( p );
			histos1D_[ "p3FathersPdgId" ]->Fill( TMath::Abs( p.pdgId() ) );
		}

		bool parton = ( ( TMath::Abs( p.pdgId() ) < 6 ) || ( p.pdgId() == 21 )  );
		if( p.status() == 23 && parton ) { 
			finalParticlesCollection.push_back( p );
			//LogWarning("daughter") << p.pdgId();
		}

	}

	vector< fullParentInfo > parents1; 
	for( auto & part : p1Collection ) {

		reco::CandidateCollection daughtersCollection = checkDaughters( part, finalParticlesCollection ); 
		
		pat::JetCollection ak8JetsMatched, ak4JetsMatched;
		//double tmpAK8JetPt = 999, tmpAK4JetPt = 999;
		//int tmpNumDauAK8 = 0; 
		for( auto & dau : daughtersCollection ) {
			//LogWarning( "daughters") << "Parent" << part.pdgId() << " daughter " << dau.pdgId();
			//LogWarning("AK8 check") << " ";
			pat::Jet tmpAK8Jet = checkDeltaR( dau, AK8jets, 0.6, histos1D_[ "p1AK8DeltaR" ], histos1D_[ "minP1AK8DeltaR" ] );
			if( tmpAK8Jet.pt() > 0 ) ak8JetsMatched.push_back( tmpAK8Jet );
			/*if( tmpAK8Jet.pt() > 0 ) { 
				tmpNumDauAK8+=1;
				if( ( tmpNumDauAK8 > 1 ) && ( tmpAK8JetPt == tmpAK8Jet.pt() ) ) ak8JetsMatched.push_back( tmpAK8Jet );
				tmpAK8JetPt = tmpAK8Jet.pt();
			}*/

			//LogWarning("AK4 check") << " ";
			pat::Jet tmpAK4Jet = checkDeltaR( dau, AK4jets, 0.3, histos1D_[ "p1AK4DeltaR" ], histos1D_[ "minP1AK4DeltaR" ] );
			if( tmpAK4Jet.pt() > 0 ) ak4JetsMatched.push_back( tmpAK4Jet );
			/*if( ( tmpAK4Jet.pt() > 0 ) && ( tmpAK4JetPt != tmpAK4Jet.pt() ) ) {
				//LogWarning("testAK4") << tmpAK4JetPt << " " << tmpAK4Jet.pt();
				ak4JetsMatched.push_back( tmpAK4Jet );
				tmpAK4JetPt = tmpAK4Jet.pt();
			}*/
		}

		//for( auto & ak8J : ak8JetsMatched ) LogWarning("matched AK8") << ak8J.pt();
		//for( auto & ak4J : ak4JetsMatched ) LogWarning("matched AK4") << ak4J.pt();
		TLorentzVector tmpParentP4;
		tmpParentP4.SetPtEtaPhiE( part.pt(), part.eta(), part.phi(), part.energy() );
		
		fullParentInfo tmpParent;
		tmpParent.genPartP4 = tmpParentP4;
		tmpParent.genPartId = part.pdgId();
		tmpParent.AK8matchedJets = ak8JetsMatched;
		tmpParent.AK4matchedJets = ak4JetsMatched;
		tmpParent.daughters = daughtersCollection;
		parents1.push_back( tmpParent );
	}
	
	double numBoosted = 0;
	double numResolved = 0;
	for( auto & parent : parents1 ){
		//LogWarning("parent") << parent.genPartId << " size ak8jets " << parent.AK8matchedJets.size() << " size ak4 jets" << parent.AK4matchedJets.size();
		if ( parent.AK8matchedJets.size() == 2 ) {
			if ( parent.AK8matchedJets[0].pt() ==  parent.AK8matchedJets[1].pt() ) { 
				numBoosted += 1;
				//LogWarning("matched ak8") << parent.AK8matchedJets[0].pt() << " " << parent.AK8matchedJets[1].pt();
			}
		}

		if ( parent.AK4matchedJets.size() == 2 ) {
			if( parent.AK4matchedJets[0].pt() !=  parent.AK4matchedJets[1].pt() ) { 
				numResolved += parent.AK4matchedJets.size();
				//LogWarning("matched ak4") << parent.AK4matchedJets[0].pt() << " " << parent.AK4matchedJets[1].pt();
			}
		}

		for( auto & dau : parent.daughters ){
			histos1D_[ "p1DaughtersPdgId" ]->Fill( dau.pdgId() );
		}
	}
	//LogWarning("count") << numBoosted << " " << numResolved;

	if ( ( numBoosted==2 ) && ( numResolved == 4 ) ) {
		
		boostedAndResolvedP1+=1;
		for( auto & parent : parents1 ){
			histos1D_[ "jet1ak8Mass_TwoBoostedFourResolved" ]->Fill( parent.AK8matchedJets[0].mass() );
			histos1D_[ "jet1ak8MPt_TwoBoostedFourResolved" ]->Fill( parent.AK8matchedJets[0].pt() );

			TLorentzVector tmpJ1, tmpJ2, tmpDijet;
			tmpJ1.SetPtEtaPhiE( parent.AK4matchedJets[0].pt(), parent.AK4matchedJets[0].eta(), parent.AK4matchedJets[0].phi(), parent.AK4matchedJets[0].energy());
			tmpJ2.SetPtEtaPhiE( parent.AK4matchedJets[1].pt(), parent.AK4matchedJets[1].eta(), parent.AK4matchedJets[1].phi(), parent.AK4matchedJets[1].energy());
			tmpDijet = tmpJ1 + tmpJ2;
			histos1D_[ "jetSysak4Mass_TwoBoostedFourResolved" ]->Fill( tmpDijet.M() );
			histos1D_[ "jetSysak4MPt_TwoBoostedFourResolved" ]->Fill( tmpDijet.Pt() );
		}

	} else if ( numBoosted == 2 ) {
		
		twoBoostedP1+=1;
		for( auto & parent : parents1 ){
			
			//LogWarning("size of ak8 matched") << parent.genPartId << " " << parent.AK8matchedJets.size() << " " << parent.daughters.size();
			/*double deltaR = 999;
			for( auto & dau : parent.daughters ){ 
				deltaR = reco::deltaR2( parent.AK8matchedJets[0].eta(), parent.AK8matchedJets[0].phi(), dau.eta(), dau.phi() );
				//LogWarning("ak8 matched deltaR") << deltaR << " " <<  parent.AK8matchedJets[0].pt();
			}*/
			//LogWarning("ak8 matched") << parent.AK8matchedJets[0].mass();
			histos1D_[ "jet1ak8Mass_TwoBoosted" ]->Fill( parent.AK8matchedJets[0].mass() );
			histos1D_[ "jet1ak8MPt_TwoBoosted" ]->Fill( parent.AK8matchedJets[0].pt() );

		}

	} else if ( numResolved == 4 ) {
		
		fourResolvedP1+=1;

		for( auto & parent : parents1 ){
			TLorentzVector tmpJ1, tmpJ2, tmpDijet;
			tmpJ1.SetPtEtaPhiE( parent.AK4matchedJets[0].pt(), parent.AK4matchedJets[0].eta(), parent.AK4matchedJets[0].phi(), parent.AK4matchedJets[0].energy());
			tmpJ2.SetPtEtaPhiE( parent.AK4matchedJets[1].pt(), parent.AK4matchedJets[1].eta(), parent.AK4matchedJets[1].phi(), parent.AK4matchedJets[1].energy());
			tmpDijet = tmpJ1 + tmpJ2;
			histos1D_[ "jetSysak4Mass_FourResolved" ]->Fill( tmpDijet.M() );
			histos1D_[ "jetSysak4MPt_FourResolved" ]->Fill( tmpDijet.Pt() );
		}
	
	} else if ( ( numBoosted==1 ) && ( numResolved == 0 ) ) {
		oneBoostedNoResolved+=1;
		LogWarning("One boosted No Resolved") << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event();
		histos1D_[ "numBoostedJets_OneBoostedNoResolved" ]->Fill( AK8jets->size() );
		if( AK8jets->size() > 0 ) histos1D_[ "jet1ak8Pt_OneBoostedNoResolved" ]->Fill( (*AK8jets)[0].pt() );
		if( AK8jets->size() > 1 ) histos1D_[ "jet2ak8Pt_OneBoostedNoResolved" ]->Fill( (*AK8jets)[1].pt() );
		histos1D_[ "numResolvedJets_OneBoostedNoResolved" ]->Fill( AK4jets->size() );
		if( AK4jets->size() > 0 ) histos1D_[ "jet1ak4Pt_OneBoostedNoResolved" ]->Fill( (*AK4jets)[0].pt() );
		if( AK4jets->size() > 1 ) histos1D_[ "jet2ak4Pt_OneBoostedNoResolved" ]->Fill( (*AK4jets)[1].pt() );
		if( AK4jets->size() > 2 ) histos1D_[ "jet3ak4Pt_OneBoostedNoResolved" ]->Fill( (*AK4jets)[2].pt() );
		if( AK4jets->size() > 3 ) histos1D_[ "jet4ak4Pt_OneBoostedNoResolved" ]->Fill( (*AK4jets)[3].pt() );

	} else if ( ( numBoosted==1 ) && ( numResolved == 2 ) ) {
		oneBoostedP1+=1;
		LogWarning("One boosted One Resolved") << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event();
		histos1D_[ "numBoostedJets_OneBoostedTwoResolved" ]->Fill( AK8jets->size() );
		if( AK8jets->size() > 0 ) histos1D_[ "jet1ak8Pt_OneBoostedTwoResolved" ]->Fill( (*AK8jets)[0].pt() );
		if( AK8jets->size() > 1 ) histos1D_[ "jet2ak8Pt_OneBoostedTwoResolved" ]->Fill( (*AK8jets)[1].pt() );
		histos1D_[ "numResolvedJets_OneBoostedTwoResolved" ]->Fill( AK4jets->size() );
		if( AK4jets->size() > 0 ) histos1D_[ "jet1ak4Pt_OneBoostedTwoResolved" ]->Fill( (*AK4jets)[0].pt() );
		if( AK4jets->size() > 1 ) histos1D_[ "jet2ak4Pt_OneBoostedTwoResolved" ]->Fill( (*AK4jets)[1].pt() );
		if( AK4jets->size() > 2 ) histos1D_[ "jet3ak4Pt_OneBoostedTwoResolved" ]->Fill( (*AK4jets)[2].pt() );
		if( AK4jets->size() > 3 ) histos1D_[ "jet4ak4Pt_OneBoostedTwoResolved" ]->Fill( (*AK4jets)[3].pt() );

	} else if ( ( numBoosted==0 ) && ( numResolved == 2 ) ) { 
		oneResolvedP1+=1;
		LogWarning("No boosted One Resolved") << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event();
		histos1D_[ "numBoostedJets_noBoostedTwoResolved" ]->Fill( AK8jets->size() );
		if( AK8jets->size() > 0 ) histos1D_[ "jet1ak8Pt_noBoostedTwoResolved" ]->Fill( (*AK8jets)[0].pt() );
		if( AK8jets->size() > 1 ) histos1D_[ "jet2ak8Pt_noBoostedTwoResolved" ]->Fill( (*AK8jets)[1].pt() );
		histos1D_[ "numResolvedJets_noBoostedTwoResolved" ]->Fill( AK4jets->size() );
		if( AK4jets->size() > 0 ) histos1D_[ "jet1ak4Pt_noBoostedTwoResolved" ]->Fill( (*AK4jets)[0].pt() );
		if( AK4jets->size() > 1 ) histos1D_[ "jet2ak4Pt_noBoostedTwoResolved" ]->Fill( (*AK4jets)[1].pt() );
		if( AK4jets->size() > 2 ) histos1D_[ "jet3ak4Pt_noBoostedTwoResolved" ]->Fill( (*AK4jets)[2].pt() );
		if( AK4jets->size() > 3 ) histos1D_[ "jet4ak4Pt_noBoostedTwoResolved" ]->Fill( (*AK4jets)[3].pt() );

	} else if ( ( numBoosted==0 ) && ( numResolved == 0 ) ) { 
		noBoostedNoResolved+=1;
		LogWarning("No boosted no Resolved") << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event();
		histos1D_[ "numBoostedJets_noBoostedNoResolved" ]->Fill( AK8jets->size() );
		if( AK8jets->size() > 0 ) histos1D_[ "jet1ak8Pt_noBoostedNoResolved" ]->Fill( (*AK8jets)[0].pt() );
		if( AK8jets->size() > 1 ) histos1D_[ "jet2ak8Pt_noBoostedNoResolved" ]->Fill( (*AK8jets)[1].pt() );
		histos1D_[ "numResolvedJets_noBoostedNoResolved" ]->Fill( AK4jets->size() );
		if( AK4jets->size() > 0 ) histos1D_[ "jet1ak4Pt_noBoostedNoResolved" ]->Fill( (*AK4jets)[0].pt() );
		if( AK4jets->size() > 1 ) histos1D_[ "jet2ak4Pt_noBoostedNoResolved" ]->Fill( (*AK4jets)[1].pt() );
		if( AK4jets->size() > 2 ) histos1D_[ "jet3ak4Pt_noBoostedNoResolved" ]->Fill( (*AK4jets)[2].pt() );
		if( AK4jets->size() > 3 ) histos1D_[ "jet4ak4Pt_noBoostedNoResolved" ]->Fill( (*AK4jets)[3].pt() );
	} else {
		none+=1;
		LogWarning("None") << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event();
	}


}


// ------------ method called once each job just before starting event loop  ------------
void Matching::beginJob() {

	edm::Service< TFileService > fileService;

	histos1D_[ "p1FathersPdgId" ] = fileService->make< TH1D >( "p1FathersPdgId", "p1FathersPdgId", 50, 1000000, 1000050 );
	histos1D_[ "p1FathersPdgId" ]->SetXTitle( "abs(pdfId) p1Collection Fathers" );
	histos1D_[ "p2FathersPdgId" ] = fileService->make< TH1D >( "p2FathersPdgId", "p2FathersPdgId", 50, 1000000, 1000050 );
	histos1D_[ "p2FathersPdgId" ]->SetXTitle( "abs(pdfId) p2Collection Fathers" );
	histos1D_[ "p3FathersPdgId" ] = fileService->make< TH1D >( "p3FathersPdgId", "p3FathersPdgId", 50, 1000000, 1000050 );
	histos1D_[ "p3FathersPdgId" ]->SetXTitle( "abs(pdfId) p3Collection Fathers" );

	histos1D_[ "p1DaughtersPdgId" ] = fileService->make< TH1D >( "p1DaughtersPdgId", "p1DaughtersPdgId", 61, -30.5, 30.5 );
	histos1D_[ "p1DaughtersPdgId" ]->SetXTitle( "p1Collection daughters pdgId" );
	histos1D_[ "p2DaughtersPdgId" ] = fileService->make< TH1D >( "p2DaughtersPdgId", "p2DaughtersPdgId", 60, -30, 30 );
	histos1D_[ "p2DaughtersPdgId" ]->SetXTitle( "p2Collection daughters pdgId" );
	histos1D_[ "p3DaughtersPdgId" ] = fileService->make< TH1D >( "p3DaughtersPdgId", "p3DaughtersPdgId", 60, -30, 30 );
	histos1D_[ "p3DaughtersPdgId" ]->SetXTitle( "p3Collection daughters pdgId" );

	histos1D_[ "p1AK8DeltaR" ] = fileService->make< TH1D >( "p1AK8DeltaR", "p1AK8DeltaR", 150, 0., 1.5 );
	histos1D_[ "p1AK8DeltaR" ]->SetXTitle( "#Delta R( jet, parton)" );
	histos1D_[ "minP1AK8DeltaR" ] = fileService->make< TH1D >( "minP1AK8DeltaR", "minP1AK8DeltaR", 50, 0., 5. );
	histos1D_[ "minP1AK8DeltaR" ]->SetXTitle( "min #Delta R( jet, parton)" );
	histos1D_[ "p1AK4DeltaR" ] = fileService->make< TH1D >( "p1AK4DeltaR", "p1AK4DeltaR", 150, 0., 1.5 );
	histos1D_[ "p1AK4DeltaR" ]->SetXTitle( "#Delta R( jet, parton)" );
	histos1D_[ "minP1AK4DeltaR" ] = fileService->make< TH1D >( "minP1AK4DeltaR", "minP1AK4DeltaR", 50, 0., 5. );
	histos1D_[ "minP1AK4DeltaR" ]->SetXTitle( "min #Delta R( jet, parton)" );

	histos1D_[ "jet1ak8Mass_TwoBoostedFourResolved" ] = fileService->make< TH1D >( "jet1ak8Mass_TwoBoostedFourResolved", "jet1ak8Mass_TwoBoostedFourResolved", 1200, 0., 1200. );
	histos1D_[ "jet1ak8Mass_TwoBoostedFourResolved" ]->SetXTitle( "2 boosted matched Mass Jet [GeV]" );
	histos1D_[ "jet1ak8MPt_TwoBoostedFourResolved" ] = fileService->make< TH1D >( "jet1ak8MPt_TwoBoostedFourResolved", "jet1ak8MPt_TwoBoostedFourResolved", 1200, 0., 1200. );
	histos1D_[ "jet1ak8MPt_TwoBoostedFourResolved" ]->SetXTitle( "2 boosted matched Pt Jet [GeV]" );
	histos1D_[ "jetSysak4Mass_TwoBoostedFourResolved" ] = fileService->make< TH1D >( "jetSysak4Mass_TwoBoostedFourResolved", "jetSysak4Mass_TwoBoostedFourResolved", 1200, 0., 1200. );
	histos1D_[ "jetSysak4Mass_TwoBoostedFourResolved" ]->SetXTitle( "2 boosted matched Mass DiJetSys [GeV]" );
	histos1D_[ "jetSysak4MPt_TwoBoostedFourResolved" ] = fileService->make< TH1D >( "jetSysak4MPt_TwoBoostedFourResolved", "jetSysak4MPt_TwoBoostedFourResolved", 1200, 0., 1200. );
	histos1D_[ "jetSysak4MPt_TwoBoostedFourResolved" ]->SetXTitle( "2 boosted matched Pt DiJetSys [GeV]" );

	histos1D_[ "jet1ak8Mass_TwoBoosted" ] = fileService->make< TH1D >( "jet1ak8Mass_TwoBoosted", "jet1ak8Mass_TwoBoosted", 1200, 0., 1200. );
	histos1D_[ "jet1ak8Mass_TwoBoosted" ]->SetXTitle( "2 boosted matched Mass Jet [GeV]" );
	histos1D_[ "jet1ak8MPt_TwoBoosted" ] = fileService->make< TH1D >( "jet1ak8MPt_TwoBoosted", "jet1ak8MPt_TwoBoosted", 1200, 0., 1200. );
	histos1D_[ "jet1ak8MPt_TwoBoosted" ]->SetXTitle( "2 boosted matched Pt Jet [GeV]" );
	histos1D_[ "jetSysak4Mass_FourResolved" ] = fileService->make< TH1D >( "jetSysak4Mass_FourResolved", "jetSysak4Mass_FourResolved", 1200, 0., 1200. );
	histos1D_[ "jetSysak4Mass_FourResolved" ]->SetXTitle( "2 boosted matched Mass DiJetSys [GeV]" );
	histos1D_[ "jetSysak4MPt_FourResolved" ] = fileService->make< TH1D >( "jetSysak4MPt_FourResolved", "jetSysak4MPt_FourResolved", 1200, 0., 1200. );
	histos1D_[ "jetSysak4MPt_FourResolved" ]->SetXTitle( "2 boosted matched Pt DiJetSys [GeV]" );

	histos1D_[ "numBoostedJets_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "numBoostedJets_OneBoostedTwoResolved", "numBoostedJets_OneBoostedTwoResolved", 10, 0., 10 );
	histos1D_[ "numBoostedJets_OneBoostedTwoResolved" ]->SetXTitle( "Number of Boosted Jets" );
	histos1D_[ "jet1ak8Pt_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "jet1ak8Pt_OneBoostedTwoResolved", "jet1ak8Pt_OneBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak8Pt_OneBoostedTwoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak8Pt_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "jet2ak8Pt_OneBoostedTwoResolved", "jet2ak8Pt_OneBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak8Pt_OneBoostedTwoResolved" ]->SetXTitle( "2 Leading Jet Pt" );

	histos1D_[ "numResolvedJets_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "numResolvedJets_OneBoostedTwoResolved", "numResolvedJets_OneBoostedTwoResolved", 10, 0., 10 );
	histos1D_[ "numResolvedJets_OneBoostedTwoResolved" ]->SetXTitle( "Number of Resolved Jets" );
	histos1D_[ "jet1ak4Pt_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "jet1ak4Pt_OneBoostedTwoResolved", "jet1ak4Pt_OneBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak4Pt_OneBoostedTwoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak4Pt_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "jet2ak4Pt_OneBoostedTwoResolved", "jet2ak4Pt_OneBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak4Pt_OneBoostedTwoResolved" ]->SetXTitle( "2 Leading Jet Pt" );
	histos1D_[ "jet3ak4Pt_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "jet3ak4Pt_OneBoostedTwoResolved", "jet3ak4Pt_OneBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet3ak4Pt_OneBoostedTwoResolved" ]->SetXTitle( "3 Leading Jet Pt" );
	histos1D_[ "jet4ak4Pt_OneBoostedTwoResolved" ] = fileService->make< TH1D >( "jet4ak4Pt_OneBoostedTwoResolved", "jet4ak4Pt_OneBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet4ak4Pt_OneBoostedTwoResolved" ]->SetXTitle( "4 Leading Jet Pt" );

	histos1D_[ "numBoostedJets_OneBoostedNoResolved" ] = fileService->make< TH1D >( "numBoostedJets_OneBoostedNoResolved", "numBoostedJets_OneBoostedNoResolved", 10, 0., 10 );
	histos1D_[ "numBoostedJets_OneBoostedNoResolved" ]->SetXTitle( "Number of Boosted Jets" );
	histos1D_[ "jet1ak8Pt_OneBoostedNoResolved" ] = fileService->make< TH1D >( "jet1ak8Pt_OneBoostedNoResolved", "jet1ak8Pt_OneBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak8Pt_OneBoostedNoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak8Pt_OneBoostedNoResolved" ] = fileService->make< TH1D >( "jet2ak8Pt_OneBoostedNoResolved", "jet2ak8Pt_OneBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak8Pt_OneBoostedNoResolved" ]->SetXTitle( "2 Leading Jet Pt" );

	histos1D_[ "numResolvedJets_OneBoostedNoResolved" ] = fileService->make< TH1D >( "numResolvedJets_OneBoostedNoResolved", "numResolvedJets_OneBoostedNoResolved", 10, 0., 10 );
	histos1D_[ "numResolvedJets_OneBoostedNoResolved" ]->SetXTitle( "Number of Resolved Jets" );
	histos1D_[ "jet1ak4Pt_OneBoostedNoResolved" ] = fileService->make< TH1D >( "jet1ak4Pt_OneBoostedNoResolved", "jet1ak4Pt_OneBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak4Pt_OneBoostedNoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak4Pt_OneBoostedNoResolved" ] = fileService->make< TH1D >( "jet2ak4Pt_OneBoostedNoResolved", "jet2ak4Pt_OneBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak4Pt_OneBoostedNoResolved" ]->SetXTitle( "2 Leading Jet Pt" );
	histos1D_[ "jet3ak4Pt_OneBoostedNoResolved" ] = fileService->make< TH1D >( "jet3ak4Pt_OneBoostedNoResolved", "jet3ak4Pt_OneBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet3ak4Pt_OneBoostedNoResolved" ]->SetXTitle( "3 Leading Jet Pt" );
	histos1D_[ "jet4ak4Pt_OneBoostedNoResolved" ] = fileService->make< TH1D >( "jet4ak4Pt_OneBoostedNoResolved", "jet4ak4Pt_OneBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet4ak4Pt_OneBoostedNoResolved" ]->SetXTitle( "4 Leading Jet Pt" );

	histos1D_[ "numBoostedJets_noBoostedNoResolved" ] = fileService->make< TH1D >( "numBoostedJets_noBoostedNoResolved", "numBoostedJets_noBoostedNoResolved", 10, 0., 10 );
	histos1D_[ "numBoostedJets_noBoostedNoResolved" ]->SetXTitle( "Number of Boosted Jets" );
	histos1D_[ "jet1ak8Pt_noBoostedNoResolved" ] = fileService->make< TH1D >( "jet1ak8Pt_noBoostedNoResolved", "jet1ak8Pt_noBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak8Pt_noBoostedNoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak8Pt_noBoostedNoResolved" ] = fileService->make< TH1D >( "jet2ak8Pt_noBoostedNoResolved", "jet2ak8Pt_noBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak8Pt_noBoostedNoResolved" ]->SetXTitle( "2 Leading Jet Pt" );

	histos1D_[ "numResolvedJets_noBoostedNoResolved" ] = fileService->make< TH1D >( "numResolvedJets_noBoostedNoResolved", "numResolvedJets_noBoostedNoResolved", 10, 0., 10 );
	histos1D_[ "numResolvedJets_noBoostedNoResolved" ]->SetXTitle( "Number of Resolved Jets" );
	histos1D_[ "jet1ak4Pt_noBoostedNoResolved" ] = fileService->make< TH1D >( "jet1ak4Pt_noBoostedNoResolved", "jet1ak4Pt_noBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak4Pt_noBoostedNoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak4Pt_noBoostedNoResolved" ] = fileService->make< TH1D >( "jet2ak4Pt_noBoostedNoResolved", "jet2ak4Pt_noBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak4Pt_noBoostedNoResolved" ]->SetXTitle( "2 Leading Jet Pt" );
	histos1D_[ "jet3ak4Pt_noBoostedNoResolved" ] = fileService->make< TH1D >( "jet3ak4Pt_noBoostedNoResolved", "jet3ak4Pt_noBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet3ak4Pt_noBoostedNoResolved" ]->SetXTitle( "3 Leading Jet Pt" );
	histos1D_[ "jet4ak4Pt_noBoostedNoResolved" ] = fileService->make< TH1D >( "jet4ak4Pt_noBoostedNoResolved", "jet4ak4Pt_noBoostedNoResolved", 100, 0., 1000 );
	histos1D_[ "jet4ak4Pt_noBoostedNoResolved" ]->SetXTitle( "4 Leading Jet Pt" );

	histos1D_[ "numBoostedJets_noBoostedTwoResolved" ] = fileService->make< TH1D >( "numBoostedJets_noBoostedTwoResolved", "numBoostedJets_noBoostedTwoResolved", 10, 0., 10 );
	histos1D_[ "numBoostedJets_noBoostedTwoResolved" ]->SetXTitle( "Number of Boosted Jets" );
	histos1D_[ "jet1ak8Pt_noBoostedTwoResolved" ] = fileService->make< TH1D >( "jet1ak8Pt_noBoostedTwoResolved", "jet1ak8Pt_noBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak8Pt_noBoostedTwoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak8Pt_noBoostedTwoResolved" ] = fileService->make< TH1D >( "jet2ak8Pt_noBoostedTwoResolved", "jet2ak8Pt_noBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak8Pt_noBoostedTwoResolved" ]->SetXTitle( "2 Leading Jet Pt" );

	histos1D_[ "numResolvedJets_noBoostedTwoResolved" ] = fileService->make< TH1D >( "numResolvedJets_noBoostedTwoResolved", "numResolvedJets_noBoostedTwoResolved", 10, 0., 10 );
	histos1D_[ "numResolvedJets_noBoostedTwoResolved" ]->SetXTitle( "Number of Resolved Jets" );
	histos1D_[ "jet1ak4Pt_noBoostedTwoResolved" ] = fileService->make< TH1D >( "jet1ak4Pt_noBoostedTwoResolved", "jet1ak4Pt_noBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet1ak4Pt_noBoostedTwoResolved" ]->SetXTitle( "Leading Jet Pt" );
	histos1D_[ "jet2ak4Pt_noBoostedTwoResolved" ] = fileService->make< TH1D >( "jet2ak4Pt_noBoostedTwoResolved", "jet2ak4Pt_noBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet2ak4Pt_noBoostedTwoResolved" ]->SetXTitle( "2 Leading Jet Pt" );
	histos1D_[ "jet3ak4Pt_noBoostedTwoResolved" ] = fileService->make< TH1D >( "jet3ak4Pt_noBoostedTwoResolved", "jet3ak4Pt_noBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet3ak4Pt_noBoostedTwoResolved" ]->SetXTitle( "3 Leading Jet Pt" );
	histos1D_[ "jet4ak4Pt_noBoostedTwoResolved" ] = fileService->make< TH1D >( "jet4ak4Pt_noBoostedTwoResolved", "jet4ak4Pt_noBoostedTwoResolved", 100, 0., 1000 );
	histos1D_[ "jet4ak4Pt_noBoostedTwoResolved" ]->SetXTitle( "4 Leading Jet Pt" );
}

// ------------ method called once each job just after ending the event loop  ------------
void Matching::endJob() {
	double total = 0;
	total = ( boostedAndResolvedP1 + twoBoostedP1 + fourResolvedP1 + noBoostedNoResolved + oneBoostedP1 + oneBoostedNoResolved + oneResolvedP1 + none );
	LogWarning( "Summary" ) << "boosted AND resolved " << boostedAndResolvedP1/total
		<< "\ntwo boosted " << twoBoostedP1/total 
		<< "\nfour resolved " << fourResolvedP1/total 
		<< "\nno Boosted No Resolved " << noBoostedNoResolved/total 
		<< "\none boosted two resolved " << oneBoostedP1/total 
		<< "\none boosted no resolved " << oneBoostedNoResolved/total 
		<< "\ntwo resolved " << oneResolvedP1/total 
		<< "\nnone " << none/total  ;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Matching::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Matching::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Matching::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Matching::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
Matching::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}*/

//define this as a plug-in
DEFINE_FWK_MODULE(Matching);
