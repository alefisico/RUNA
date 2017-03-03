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
typedef struct {
	bool pass;
	double deltaR;
	int indexJet;
	pat::Jet matchJet;
} matched;

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

matched checkDeltaR(reco::Candidate & p1, Handle<pat::JetCollection> jets, double minDeltaR){

	std::vector<double> deltaRVec;
	pat::Jet matchedJet;
	double deltaR = 99999;
	int ind = -1;

	//edm::LogWarning("genParticle ")  << p1.pdgId() << " " << jets->size();

	for( unsigned int j=0; j<jets->size(); j++ ) {
		const pat::Jet & p2 = (*jets)[j];
		double tmpdeltaR = reco::deltaR2( p1.rapidity(), p1.phi(), p2.rapidity(), p2.phi() );
		//edm::LogWarning("calc deltaR") << tmpdeltaR << " " << j;
		if( tmpdeltaR < deltaR ) {
			deltaR = tmpdeltaR;
			ind = j;
			matchedJet = p2;
		}
		//dummy+=1;
	}
	//edm::LogWarning("deltaR") << deltaR << " " << ind ;
	matched results;
	if( deltaR < minDeltaR) {
		results.pass = true;
		results.deltaR = deltaR;
		results.indexJet = ind;
		results.matchJet = matchedJet;
	} else {
		results.pass = false;
	}

	return results;
}

std::map< int, reco::CandidateCollection > checkDaughters( reco::CandidateCollection pCollection, reco::CandidateCollection finalParticlesCollection ){

	std::map< int, reco::CandidateCollection >  daughtersParticle;
	int dummy = 0;
	for( auto & p : pCollection ) {
		reco::CandidateCollection tmp1;
		for( auto & fp : finalParticlesCollection ) {
			const reco::Candidate * finalMother = fp.mother();
			if( isAncestor( p, finalMother ) ) tmp1.push_back( fp ); // LogWarning("Particle found 1") << jp1->pdgId() << " " << fp.pdgId(); }
		}
		daughtersParticle[ dummy ] = tmp1;
		dummy+=1;
	}
	return daughtersParticle;
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
			//LogWarning("mother") << p.pdgId();
		}
		if( ( TMath::Abs( p.pdgId() ) == particle2 ) && p.status() == 22 ) p2Collection.push_back( p );
		if( ( TMath::Abs( p.pdgId() ) == particle3 ) && p.status() == 22 ) p3Collection.push_back( p );

		bool parton = ( ( TMath::Abs( p.pdgId() ) < 6 ) || ( p.pdgId() == 21 )  );
		if( p.status() == 23 && parton ) { 
			finalParticlesCollection.push_back( p );
			//LogWarning("daughter") << p.pdgId();
		}

	}

	std::map< int, reco::CandidateCollection  > daughtersParticle1, daughtersParticle2, daughtersParticle3;
	daughtersParticle1 = checkDaughters( p1Collection, finalParticlesCollection ); 
	daughtersParticle2 = checkDaughters( p2Collection, finalParticlesCollection ); 
	daughtersParticle3 = checkDaughters( p3Collection, finalParticlesCollection ); 

	/*
	LogWarning("number of daughters 1") << daughtersParticle1.size() << " " << daughtersParticle1[0].size();
	LogWarning("number of daughters 2") << daughtersParticle2.size() << " " << daughtersParticle2[0].size();
	LogWarning("number of daughters 3") << daughtersParticle3.size() << " " << daughtersParticle3[0].size();
	*/
	
	// Particle 1
	map< int, matched > AK8Particle1, AK4Particle1;
	for ( auto set : daughtersParticle1 ) {
		int numAK8matched = 0;
		for ( reco::Candidate & p : set.second ) {
			//LogWarning("test") << set.first << " " << p.pdgId();
	       		matched infoMatchedAK8;
			//LogWarning("START AK8") << AK8jets->size();
			infoMatchedAK8 = checkDeltaR( p, AK8jets, 0.8 );
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
				infoMatchedAK4 = checkDeltaR( p, AK4jets, 0.4 );
				bool passParAK4 = infoMatchedAK4.pass;
				if ( passParAK4 ) AK4Particle1[ p.pdgId() ] = infoMatchedAK4;
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



}


// ------------ method called once each job just before starting event loop  ------------
void Matching::beginJob() {

	edm::Service< TFileService > fileService;

	histos1D_[ "p1DaughtersDeltaR" ] = fileService->make< TH1D >( "p1DaughtersDeltaR", "p1DaughtersDeltaR", 150, 0., 1.5 );
	histos1D_[ "p1DaughtersDeltaR" ]->SetXTitle( "#Delta R( parton, parton)" );
	histos1D_[ "p2DaughtersDeltaR" ] = fileService->make< TH1D >( "p2DaughtersDeltaR", "p2DaughtersDeltaR", 150, 0., 1.5 );
	histos1D_[ "p2DaughtersDeltaR" ]->SetXTitle( "#Delta R( parton, parton)" );
	histos1D_[ "p1DeltaR" ] = fileService->make< TH1D >( "p1DeltaR", "p1DeltaR", 50, 0., 5. );
	histos1D_[ "p1DeltaR" ]->SetXTitle( "min #Delta R( jet, parton)" );
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
