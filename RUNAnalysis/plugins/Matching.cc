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
		int twoBoostedP1 = 0;
		int fourResolvedP1 = 0;
		int away = 0;
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

/*matched checkDeltaR(reco::Candidate & p1, Handle<pat::JetCollection> jets, double minDeltaR, TH1D * histo){

	std::vector<double> deltaRVec;
	pat::Jet matchedJet;
	double deltaR = 99999;
	//int ind = -1;

	//edm::LogWarning("genParticle ")  << p1.pdgId() << " " << jets->size();

	for( unsigned int j=0; j<jets->size(); j++ ) {
		const pat::Jet & p2 = (*jets)[j];
	//for( auto & p2 : jets ) {
		double tmpdeltaR2 = reco::deltaR2( p1.rapidity(), p1.phi(), p2.rapidity(), p2.phi() );
		double tmpdeltaR4 = reco::deltaR2( p1.eta(), p1.phi(), p2.eta(), p2.phi() );
		TLorentzVector tmp1, tmp2;
		tmp1.SetPtEtaPhiE( p1.pt(), p1.eta(), p1.phi(), p1.energy() );
		tmp2.SetPtEtaPhiE( p2.pt(), p2.eta(), p2.phi(), p2.energy() );
		double tmpdeltaR = tmp1.DeltaR( tmp2 );
		double tmpdeltaR3 = TMath::Sqrt( TMath::Power( (p1.eta()-p2.eta()), 2) + TMath::Power( (p1.phi()-p2.phi()), 2) );
		histo->Fill( tmpdeltaR );
		//edm::LogWarning("calc deltaR") << " " << tmpdeltaR << " " << tmpdeltaR2 << " " << tmpdeltaR3 << " " << tmpdeltaR4;
		if( tmpdeltaR < deltaR ) {
			deltaR = tmpdeltaR;
			//ind = j;
			matchedJet = p2;
		}
		//dummy+=1;
	}
	//edm::LogWarning("deltaR") << deltaR; // << " " << ind ;
	matched results;
	if( deltaR < minDeltaR) {
		//edm::LogWarning("pass") << deltaR << " " << ind ;
		results.pass = true;
		results.deltaR = deltaR;
		//results.indexJet = ind;
		results.matchJet = matchedJet;
	} else {
		results.pass = false;
	}

	return results;
}*/

pat::Jet checkDeltaR(reco::Candidate & p1, Handle<pat::JetCollection> jets, double minDeltaR, TH1D * histo){

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
		histo->Fill( tmpdeltaR );
		//edm::LogWarning("calc deltaR") << " " << tmpdeltaR << " " << tmpdeltaR2 << " " << tmpdeltaR3 << " " << tmpdeltaR4;
		if( tmpdeltaR < deltaR ) {
			deltaR = tmpdeltaR;
			if( deltaR < minDeltaR ) matchedJet = p2;
		}
	}
	//edm::LogWarning("deltaR") << deltaR; // << " " << ind ;
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

	//reco::CandidateCollection  daughtersParticle1, daughtersParticle2, daughtersParticle3;
	vector< fullParentInfo > parents1; 
	for( auto & part : p1Collection ) {

		reco::CandidateCollection daughtersCollection = checkDaughters( part, finalParticlesCollection ); 
		
		pat::JetCollection ak8JetsMatched, ak4JetsMatched;
		double tmpJetPt = -999;
		for( auto & dau : daughtersCollection ) {
			//LogWarning( "daughters") << part.pdgId() << " " << dau.pdgId();
			pat::Jet tmpAK8Jet = checkDeltaR( dau, AK8jets, 0.8, histos1D_[ "p1DaughtersAK8DeltaR" ] );
			if( ( tmpAK8Jet.pt() > 0 ) && ( tmpJetPt != tmpAK8Jet.pt() ) ) {
				LogWarning("test") << tmpJetPt << " " << tmpAK8Jet.pt();
				ak8JetsMatched.push_back( tmpAK8Jet );
				tmpJetPt = tmpAK8Jet.pt();
			}

			pat::Jet tmpAK4Jet = checkDeltaR( dau, AK4jets, 0.4, histos1D_[ "p1DaughtersAK4DeltaR" ] );
			ak4JetsMatched.push_back( tmpAK4Jet );
		}

		for( auto & ak8J : ak8JetsMatched ) LogWarning("matched AK8") << ak8J.pt();
		for( auto & ak4J : ak4JetsMatched ) LogWarning("matched AK4") << ak4J.pt();
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
	
	/*for( auto & parent : parents1 ){
		LogWarning("test") << parent.genPartId << " " << parent.AK8matchedJets.size() << " " << parent.daughters.size();
	}*/

	/*
	LogWarning("number of daughters 1") << daughtersParticle1.size() << " " << daughtersParticle1[0].size();
	LogWarning("number of daughters 2") << daughtersParticle2.size() << " " << daughtersParticle2[0].size();
	LogWarning("number of daughters 3") << daughtersParticle3.size() << " " << daughtersParticle3[0].size();
	*/
	
	/*/ Particle 1
	map< int, matched > AK8Particle1, AK4Particle1;
	for ( auto set : daughtersParticle1 ) {
		int numAK8matched = 0;
		for ( reco::Candidate & p : set.second ) {
			//LogWarning("test") << set.first << " " << p.pdgId();
			histos1D_[ "p1DaughtersPdgId" ]->Fill( p.pdgId() );
	       		matched infoMatchedAK8;
			//LogWarning("START AK8") << AK8jets->size();
			infoMatchedAK8 = checkDeltaR( p, AK8jets, 0.8, histos1D_[ "p1DaughtersAK8DeltaR" ] );
			histos1D_[ "p1AK8DeltaR" ]->Fill( infoMatchedAK8.deltaR  );
			bool passParAK8 = infoMatchedAK8.pass;
			if ( passParAK8 ) {
				AK8Particle1[ p.pdgId() ] = infoMatchedAK8;
				numAK8matched+=1;
				//LogWarning("pdgId") << set.first << " "  <<  p.pdgId();
			}
		}

		if ( numAK8matched < 2 ) {
			for ( reco::Candidate & p : set.second ) {
				matched infoMatchedAK4;
				//LogWarning("START AK4") << AK4jets->size();
				infoMatchedAK4 = checkDeltaR( p, AK4jets, 0.4, histos1D_[ "p1DaughtersAK4DeltaR" ]);
				bool passParAK4 = infoMatchedAK4.pass;
				if ( passParAK4 ) {
					AK4Particle1[ p.pdgId() ] = infoMatchedAK4;
					histos1D_[ "p1AK4DeltaR" ]->Fill( infoMatchedAK4.deltaR  );
				}
			}
		}
	}
	//LogWarning("test") << AK8Particle1.size() << " " << AK4Particle1.size(); 
	
	if (( AK8jets->size() < 2 ) && ( AK4jets->size() < 2 ) ){
		away+=1;
	} else if ( AK8Particle1.size() == 4 ) {
		twoBoostedP1+=1;
	} else if ( AK4Particle1.size() == 4 ) {
		fourResolvedP1+=1;
	} else if ( ( AK8Particle1.size() == 2 ) && ( AK4Particle1.size() == 2 ) ){
		oneBoostedP1+=1;
	} else none+=1;

	*/

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

	histos1D_[ "p1DaughtersAK8DeltaR" ] = fileService->make< TH1D >( "p1DaughtersAK8DeltaR", "p1DaughtersAK8DeltaR", 150, 0., 1.5 );
	histos1D_[ "p1DaughtersAK8DeltaR" ]->SetXTitle( "#Delta R( parton, parton)" );
	histos1D_[ "p1AK8DeltaR" ] = fileService->make< TH1D >( "p1AK8DeltaR", "p1AK8DeltaR", 50, 0., 5. );
	histos1D_[ "p1AK8DeltaR" ]->SetXTitle( "min #Delta R( jet, parton)" );
	histos1D_[ "p1DaughtersAK4DeltaR" ] = fileService->make< TH1D >( "p1DaughtersAK4DeltaR", "p1DaughtersAK4DeltaR", 150, 0., 1.5 );
	histos1D_[ "p1DaughtersAK4DeltaR" ]->SetXTitle( "#Delta R( parton, parton)" );
	histos1D_[ "p1AK4DeltaR" ] = fileService->make< TH1D >( "p1AK4DeltaR", "p1AK4DeltaR", 50, 0., 5. );
	histos1D_[ "p1AK4DeltaR" ]->SetXTitle( "min #Delta R( jet, parton)" );

	histos1D_[ "p2DaughtersDeltaR" ] = fileService->make< TH1D >( "p2DaughtersDeltaR", "p2DaughtersDeltaR", 150, 0., 1.5 );
	histos1D_[ "p2DaughtersDeltaR" ]->SetXTitle( "#Delta R( parton, parton)" );
	histos1D_[ "p1JetMass" ] = fileService->make< TH1D >( "p1JetMass", "p1JetMass", 120, 0., 1200. );
	histos1D_[ "p1JetMass" ]->SetXTitle( "Mass Jet matched [GeV]" );

	histos1D_[ "p2DeltaR" ] = fileService->make< TH1D >( "p2DeltaR", "p2DeltaR", 50, 0., 5. );
	histos1D_[ "p2DeltaR" ]->SetXTitle( "min #Delta R( jet, parton)" );
	histos1D_[ "p2JetMass" ] = fileService->make< TH1D >( "p2JetMass", "p2JetMass", 120, 0., 1200. );
	histos1D_[ "p2JetMass" ]->SetXTitle( "Mass Jet matched [GeV]" );

	histos1D_[ "p3DeltaR" ] = fileService->make< TH1D >( "p3DeltaR", "p3DeltaR", 50, 0., 5. );
	histos1D_[ "p3DeltaR" ]->SetXTitle( "min #Delta R( jet, parton)" );
	histos1D_[ "p3JetMass" ] = fileService->make< TH1D >( "p3JetMass", "p3JetMass", 120, 0., 1200. );
	histos1D_[ "p3JetMass" ]->SetXTitle( "Mass Jet matched [GeV]" );

}

// ------------ method called once each job just after ending the event loop  ------------
void Matching::endJob() {
	LogWarning( "Summary" ) << "two boosted " << twoBoostedP1 << " four resolved " << fourResolvedP1 << "  away " << away << " one boosted " << oneBoostedP1 << " none " << none  ;
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
