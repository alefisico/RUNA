// -*- C++ -*-
//
// Package:    RUNA/RUNAnalysis
// Class:      RUNResolvedResolutionCalc
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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "RUNA/RUNAnalysis/interface/CommonVariablesStructure.h"
#include "RUNA/RUNAnalysis/interface/PUReweighter.h"


using namespace edm;
using namespace std;
using namespace reco;

//
// constants, enums and typedefs
//

//
// class declaration
//
class RUNResolvedResolutionCalc : public EDAnalyzer {
	public:
		explicit RUNResolvedResolutionCalc(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~RUNResolvedResolutionCalc();

	private:
		virtual void beginJob() override;
		virtual void analyze(const Event&, const EventSetup&) override;
		virtual void endJob() override;
		virtual void clearVariables();
		//virtual void beginRun(const Run&, const EventSetup&) override;
		//virtual void endRun(Run const&, EventSetup const&) override;

		//virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
		//virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

		// ----------member data ---------------------------
		PUReweighter PUWeight_;
		//int lhaPdfId;

		Service<TFileService> fs_;
		TTree *RUNAtree;
		map< string, TH1D* > histos1D_;
		map< string, TH2D* > histos2D_;
		vector< string > cutLabels;
		map< string, double > cutmap;

		string dataPUFile;
		string jecVersion;
		double cutAK4jetPt;
		double cutAK4HT;
		double cutAK4MassAsym;
		double cutDelta;
		double cutDeltaEtaDijetSyst;

		vector<JetCorrectorParameters> jetPar;
		FactorizedJetCorrector * jetJEC;
		JetCorrectionUncertainty *jetCorrUnc;

		vector<float> *jetsPt = new std::vector<float>();
		vector<float> *jetsEta = new std::vector<float>();
		vector<float> *jetsPhi = new std::vector<float>();
		vector<float> *jetsE = new std::vector<float>();
		vector<float> *jetsQGL = new std::vector<float>();
		vector<float> *jetsCSVv2 = new std::vector<float>();
		vector<float> *jetsCSVv2SF = new std::vector<float>();
		vector<float> *jetsCMVAv2 = new std::vector<float>();
		vector<float> *muonsPt = new std::vector<float>();
		vector<float> *elesPt  = new std::vector<float>();
		ULong64_t event = 0;
		int numJets = 0, numPV = 0, numEle = 0, numMuon = 0;
		unsigned int lumi = 0, run=0;
		float HT = 0, mass1 = -999, mass2 = -999, massAve = -999, MET = -999,
		      jet1Pt = -9999, jet2Pt = -9999, jet3Pt = -9999, jet4Pt = -9999,
		      delta1 = -999, delta2 = -999, massAsym = -999, eta1 = -999, eta2 = -999, deltaEta = -999, 
		      deltaR = -999, cosThetaStar1 = -999, cosThetaStar2 = -999,
		      puWeight = -999, genWeight = -999, lumiWeight = -999, pdfWeight = -999 ;
		bool btagCSVv2Pair12 = 0, btagCSVv2Pair34 = 0;
		vector<float> scaleWeights, pdfWeights, alphaWeights;

		/// Jets
		EDGetTokenT<vector<float>> jetPt_;
		EDGetTokenT<vector<float>> jetEta_;
		EDGetTokenT<vector<float>> jetPhi_;
		EDGetTokenT<vector<float>> jetE_;
		EDGetTokenT<vector<float>> jetQGL_;
		EDGetTokenT<vector<float>> jetCSVv2_;
		EDGetTokenT<vector<float>> jetCMVAv2_;
		EDGetTokenT<vector<float>> jetArea_;
		EDGetTokenT<vector<float>> jetGenPt_;
		EDGetTokenT<vector<float>> jetGenEta_;
		EDGetTokenT<vector<float>> jetGenPhi_;
		EDGetTokenT<vector<float>> jetGenE_;
		EDGetTokenT<vector<float>> jetHadronFlavour_;

		/// Event variables
		EDGetTokenT<int> NPV_;
		EDGetTokenT<vector<float>> metPt_;
		EDGetTokenT<int> trueNInt_;
		EDGetTokenT<vector<int>> bunchCross_;
		EDGetTokenT<vector<float>> rho_;
		EDGetTokenT<vector<int>> puNumInt_;
		EDGetTokenT<unsigned int> lumi_;
		EDGetTokenT<unsigned int> run_;
		EDGetTokenT<ULong64_t> event_;
		EDGetTokenT<GenEventInfoProduct> generator_;
		EDGetTokenT<GenParticleCollection> genParticles_;

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
RUNResolvedResolutionCalc::RUNResolvedResolutionCalc(const ParameterSet& iConfig):
	jetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPt"))),
	jetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetEta"))),
	jetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPhi"))),
	jetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetE"))),
	jetQGL_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetQGL"))),
	jetCSVv2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCSVv2"))),
	jetCMVAv2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCMVAv2"))),
	jetArea_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetArea"))),
	jetGenPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenPt"))),
	jetGenEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenEta"))),
	jetGenPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenPhi"))),
	jetGenE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenE"))),
	jetHadronFlavour_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetHadronFlavour"))),
	NPV_(consumes<int>(iConfig.getParameter<InputTag>("NPV"))),
	metPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("metPt"))),
	trueNInt_(consumes<int>(iConfig.getParameter<InputTag>("trueNInt"))),
	bunchCross_(consumes<vector<int>>(iConfig.getParameter<InputTag>("bunchCross"))),
	rho_(consumes<vector<float>>(iConfig.getParameter<InputTag>("rho"))),
	puNumInt_(consumes<vector<int>>(iConfig.getParameter<InputTag>("puNumInt"))),
	lumi_(consumes<unsigned int>(iConfig.getParameter<InputTag>("Lumi"))),
	run_(consumes<unsigned int>(iConfig.getParameter<InputTag>("Run"))),
	event_(consumes<ULong64_t>(iConfig.getParameter<InputTag>("Event"))),
	generator_(consumes<GenEventInfoProduct>(iConfig.getParameter<InputTag>("generator"))),
	genParticles_(consumes<GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles"))),
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
	cutAK4jetPt 	= iConfig.getParameter<double>("cutAK4jetPt");
	cutAK4HT	= iConfig.getParameter<double>("cutAK4HT");
	cutAK4MassAsym	= iConfig.getParameter<double>("cutAK4MassAsym");
	cutDelta        = iConfig.getParameter<double>("cutDelta");
	cutDeltaEtaDijetSyst	= iConfig.getParameter<double>("cutDeltaEtaDijetSyst");
	dataPUFile 	= iConfig.getParameter<string>("dataPUFile");
	jecVersion 	= iConfig.getParameter<string>("jecVersion");

	/////// JECs
	string tmpPrefix = jecVersion;
	string prefix;
	prefix = tmpPrefix + "_MC_";

	// all jet
	vector<string> jecPayloadNames_;
	jecPayloadNames_.push_back(prefix + "L1FastJet_AK4PFchs.txt");
	jecPayloadNames_.push_back(prefix + "L2Relative_AK4PFchs.txt");
	jecPayloadNames_.push_back(prefix + "L3Absolute_AK4PFchs.txt");

	for ( vector<string>::const_iterator payloadBegin = jecPayloadNames_.begin(), payloadEnd = jecPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		jetPar.push_back(pars);
	}
	jetJEC = new FactorizedJetCorrector(jetPar);

	// jec uncertainty
	JetCorrectorParameters jecUncParam( prefix + "Uncertainty_AK4PFchs.txt");
	jetCorrUnc  = new JetCorrectionUncertainty( jecUncParam);

}


RUNResolvedResolutionCalc::~RUNResolvedResolutionCalc()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNResolvedResolutionCalc::analyze(const Event& iEvent, const EventSetup& iSetup) {


	Handle<vector<float> > jetPt;
	iEvent.getByToken(jetPt_, jetPt);

	Handle<vector<float> > jetEta;
	iEvent.getByToken(jetEta_, jetEta);

	Handle<vector<float> > jetPhi;
	iEvent.getByToken(jetPhi_, jetPhi);

	Handle<vector<float> > jetE;
	iEvent.getByToken(jetE_, jetE);

	Handle<vector<float> > jetQGL;
	iEvent.getByToken(jetQGL_, jetQGL);

	Handle<vector<float> > jetCSVv2;
	iEvent.getByToken(jetCSVv2_, jetCSVv2);

	Handle<vector<float> > jetCMVAv2;
	iEvent.getByToken(jetCMVAv2_, jetCMVAv2);

	Handle<vector<float> > jetArea;
	iEvent.getByToken(jetArea_, jetArea);

	Handle<vector<float> > jetGenPt;
	iEvent.getByToken(jetGenPt_, jetGenPt);

	Handle<vector<float> > jetGenEta;
	iEvent.getByToken(jetGenEta_, jetGenEta);

	Handle<vector<float> > jetGenPhi;
	iEvent.getByToken(jetGenPhi_, jetGenPhi);

	Handle<vector<float> > jetGenE;
	iEvent.getByToken(jetGenE_, jetGenE);

	Handle<vector<float> > jetHadronFlavour;
	iEvent.getByToken(jetHadronFlavour_, jetHadronFlavour);

	Handle<int> NPV;
	iEvent.getByToken(NPV_, NPV);

	Handle<int> trueNInt;
	iEvent.getByToken(trueNInt_, trueNInt);

	Handle<vector<int>> bunchCross;
	iEvent.getByToken(bunchCross_, bunchCross);

	Handle<vector<float>> rho;
	iEvent.getByToken(rho_, rho);

	Handle<vector<int>> puNumInt;
	iEvent.getByToken(puNumInt_, puNumInt);

	Handle<unsigned int> Lumi;
	iEvent.getByToken(lumi_, Lumi);

	Handle<unsigned int> Run;
	iEvent.getByToken(run_, Run);

	Handle<ULong64_t> ievent;
	iEvent.getByToken(event_, ievent);

	Handle<vector<float> > metPt;
	iEvent.getByToken(metPt_, metPt);

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


	Handle< GenParticleCollection > genParticles;
	iEvent.getByToken( genParticles_, genParticles );

	vector<TLorentzVector> dauStop1, dauStop2;
	for( size_t i = 0; i < genParticles->size(); i++ ) {

		const reco::Candidate &p = ( *genParticles )[i];
		const reco::Candidate * pMother = p.mother();
		if (pMother){
			if( ( pMother->pdgId() == 1000002 ) && ( p.status() == 23 ) ) { 
				//LogWarning("daughter") << p.pdgId();
				TLorentzVector tmpStopDau;
				tmpStopDau.SetPtEtaPhiE( p.pt(), p.eta(), p.phi(), p.energy() );
				dauStop1.push_back( tmpStopDau );
			}
			
			if( ( pMother->pdgId() == -1000002 ) && ( p.status() == 23 ) ) { 
				//LogWarning("daughter") << p.pdgId();
				TLorentzVector tmpStopDau;
				tmpStopDau.SetPtEtaPhiE( p.pt(), p.eta(), p.phi(), p.energy() );
				dauStop2.push_back( tmpStopDau );
			}
		}

	}


	// PU Reweight
	puWeight = PUWeight_.getPUWeight( *trueNInt, *bunchCross );
	//puWeight = 1;
	histos1D_[ "PUWeight" ]->Fill( puWeight );
	double totalWeight = puWeight;
	///////////////////////////////////////////////////*/
	
	cutmap["Processed"] += 1;

	int numPV = *NPV;
	vector< myJet > JETS;
	int numberJets = 0;
	double rawHT = 0;
	HT = 0;

	/////// Preselect jets
	for (size_t i = 0; i < jetPt->size(); i++) {

		if( TMath::Abs( (*jetEta)[i] ) > 2.4 ) continue;

		rawHT += (*jetPt)[i];
		histos1D_[ "rawJetPt" ]->Fill( (*jetPt)[i] , totalWeight );

		string typeOfJetID = "tightLepVetoJetID";
		bool jetId = jetID( (*jetEta)[i], (*jetE)[i], (*jecFactor)[i], (*neutralHadronEnergyFrac)[i], (*neutralEmEnergyFrac)[i], (*chargedHadronEnergyFrac)[i], (*muonEnergy)[i], (*chargedEmEnergyFrac)[i], (*chargedMultiplicity)[i], (*neutralMultiplicity)[i],  typeOfJetID ); 

		TLorentzVector tmpJet, rawJet, corrJet, genJet, smearJet;
		tmpJet.SetPtEtaPhiE( (*jetPt)[i], (*jetEta)[i], (*jetPhi)[i], (*jetE)[i] );
		rawJet = tmpJet* (*jecFactor)[i] ;

		double JEC = corrections( rawJet, (*jetArea)[i], (*rho)[i] ,*NPV, jetJEC); 
		corrJet = rawJet* JEC;

		if( ( corrJet.Pt() > cutAK4jetPt ) && jetId ) { 

			HT += corrJet.Pt();
			++numberJets;
			//if ( (*jetCSVv2)[i] > 0.244 ) bTagCSVv2 = 1; 	// CSVv2L
			//if ( (*jetCSVv2)[i] > 0.679 ) bTagCSVv2 = 1; 	// CSVv2M
			//if ( (*jetCSVv2V1)[i] > 0.405 ) bTagCSVv2 = 1; 	// CSVv2V1L
			//if ( (*jetCSVv2V1)[i] > 0.783 ) bTagCSVv2 = 1; 	// CSVv2V1M
			double jec = 1. / (*jecFactor)[i];
			histos1D_[ "jetPt" ]->Fill( corrJet.Pt() , totalWeight );
			histos1D_[ "jetEta" ]->Fill( corrJet.Eta() , totalWeight );
			histos1D_[ "neutralHadronEnergyFrac" ]->Fill( (*neutralHadronEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "neutralEmEnergyFrac" ]->Fill( (*neutralEmEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "chargedHadronEnergyFrac" ]->Fill( (*chargedHadronEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "chargedEmEnergyFrac" ]->Fill( (*chargedEmEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "numConst" ]->Fill( (*chargedMultiplicity)[i] + (*neutralMultiplicity)[i], totalWeight );
			histos1D_[ "chargedMultiplicity" ]->Fill( (*chargedMultiplicity)[i] * jec, totalWeight );


			myJet tmpJET;
			tmpJET.p4 = corrJet;
			tmpJET.btagCSVv2 = (*jetCSVv2)[i];
			tmpJET.btagCMVAv2 = (*jetCMVAv2)[i];
			tmpJET.qgl = (*jetQGL)[i];
			tmpJET.nhf = (*neutralHadronEnergyFrac)[i] * jec;
			tmpJET.nEMf = (*neutralEmEnergyFrac)[i] * jec;
			tmpJET.chf = (*chargedHadronEnergyFrac)[i] * jec;
			tmpJET.cEMf = (*chargedEmEnergyFrac)[i] * jec;
			tmpJET.numConst = (*chargedMultiplicity)[i] + (*neutralMultiplicity)[i];
			tmpJET.chm = (*chargedMultiplicity)[i] * jec;
			tmpJET.hadronFlavour = (*jetHadronFlavour)[i];
			JETS.push_back( tmpJET );
		}
	}
	///////////////////////////////////////////////////*/

	
	sort(JETS.begin(), JETS.end(), [](const myJet &p1, const myJet &p2) { return p1.p4.Pt() > p2.p4.Pt(); });  /// after corrections jets are not 100% sorted in Pt

	numJets = numberJets;
	histos1D_[ "jetNum" ]->Fill( numJets, totalWeight );
	histos1D_[ "NPV" ]->Fill( numPV, totalWeight );
	histos1D_[ "NPV_NOPUWeight" ]->Fill( numPV );
	if ( HT > 0 ) histos1D_[ "HT" ]->Fill( HT , totalWeight );
	if ( rawHT > 0 ) histos1D_[ "rawHT" ]->Fill( rawHT , totalWeight );
	MET = (*metPt)[0];

	clearVariables();
	
	if ( ( dauStop1.size() == 2 ) && ( dauStop2.size() == 2 ) ) { 

		if( numJets > 3 ) { 
		
			cutmap["4Jets"] += 1;
			histos1D_[ "HT_cut4Jets" ]->Fill( HT, totalWeight );
			histos1D_[ "jetNum_cut4Jets" ]->Fill( numJets, totalWeight );
			histos1D_[ "jet1Pt_cut4Jets" ]->Fill( JETS[0].p4.Pt(), totalWeight );
			histos1D_[ "jet2Pt_cut4Jets" ]->Fill( JETS[1].p4.Pt(), totalWeight );
			histos1D_[ "jet3Pt_cut4Jets" ]->Fill( JETS[2].p4.Pt(), totalWeight );
			histos1D_[ "jet4Pt_cut4Jets" ]->Fill( JETS[3].p4.Pt(), totalWeight );
			histos1D_[ "MET_cut4Jets" ]->Fill( MET, totalWeight );
			histos1D_[ "METHT_cut4Jets" ]->Fill( MET/HT, totalWeight );

			if( ( HT > cutAK4HT ) && ( JETS[3].p4.Pt() > cutAK4jetPt ) ){

				myJet j1, j2, j3, j4;
				double tmpMinDR1 = 99999, tmpMinDR2 = 99999, tmpMinDR3 = 99999, tmpMinDR4 = 99999;
				for( auto & jet : JETS ){
					double tmpDeltaRjetDau00 = jet.p4.DeltaR( dauStop1[0] );
					if( tmpDeltaRjetDau00 < tmpMinDR1 ) {
						tmpMinDR1 = tmpDeltaRjetDau00;
						//LogWarning("test deltaR") << tmpMinDR1;
						if ( tmpMinDR1 < 0.3 ) {
							j1 = jet;
							//LogWarning("test j1 mindeltaR") << tmpMinDR1 << " " << j1.p4.Pt();
						}
					}

					double tmpDeltaRjetDau01 = jet.p4.DeltaR( dauStop1[1] );
					if( tmpDeltaRjetDau01 < tmpMinDR2 ) {
						tmpMinDR2 = tmpDeltaRjetDau01;
						if ( tmpMinDR2 < 0.3 ) {
							j2 = jet;
							//LogWarning("test j2 mindeltaR") << tmpMinDR2 << " " << j2.p4.Pt();
						}
					}

					double tmpDeltaRjetDau10 = jet.p4.DeltaR( dauStop2[0] );
					if( tmpDeltaRjetDau10 < tmpMinDR3 ) {
						tmpMinDR3 = tmpDeltaRjetDau10;
						if ( tmpMinDR3 < 0.3 ) j3 = jet;
					}

					double tmpDeltaRjetDau11 = jet.p4.DeltaR( dauStop2[1] );
					if( tmpDeltaRjetDau11 < tmpMinDR4 ) {
						tmpMinDR4 = tmpDeltaRjetDau11;
						if ( tmpMinDR4 < 0.3 ) j4 = jet;
					}
				}	
				//LogWarning("test all jet") << j1.p4.Pt() << " " << j2.p4.Pt() << " " << j3.p4.Pt() << " " << j4.p4.Pt() ;

				if( ( j1.p4.Pt() > 0 ) && ( j2.p4.Pt() > 0 ) && ( j3.p4.Pt() > 0 ) && ( j4.p4.Pt() > 0 ) ) {

					if( (j1.p4.Pt() != j2.p4.Pt()) && (j1.p4.Pt() != j3.p4.Pt()) && (j1.p4.Pt() != j4.p4.Pt()) 
							&& (j2.p4.Pt() != j3.p4.Pt()) && (j2.p4.Pt() != j4.p4.Pt())
							&& (j3.p4.Pt() != j4.p4.Pt()) ) {

						mass1 = ( j1.p4 + j2.p4 ).M();
						mass2 = ( j3.p4 + j4.p4 ).M();
						massAve = ( mass1 + mass2 ) / 2;
						delta1 = ( j1.p4.Pt() + j2.p4.Pt() ) - massAve;
						delta2 = ( j3.p4.Pt() + j4.p4.Pt() ) - massAve;
						massAsym = TMath::Abs( mass1 - mass2 ) / (mass1 + mass2) ;
						eta1 = ( j1.p4 + j2.p4 ).Eta();
						eta2 = ( j3.p4 + j4.p4 ).Eta();
						deltaEta = TMath::Abs( eta1 - eta2 );
						deltaR = abs( ( j1.p4.DeltaR( j2.p4 ) - 0.8 )  + abs( ( j3.p4.DeltaR( j4.p4 ) - 0.8 ) ) );
						cosThetaStar1 = cosThetaStar( j1.p4, j2.p4 );
						cosThetaStar2 = cosThetaStar( j3.p4, j4.p4 );
						btagCSVv2Pair12 = ( (j1.btagCSVv2 > 0.8) || (j2.btagCSVv2 > 0.8) );
						btagCSVv2Pair34 = ( (j3.btagCSVv2 > 0.8) || (j4.btagCSVv2 > 0.8) );
					
						cutmap["BestPair"] += 1;
						histos1D_[ "HT_cutBestPair" ]->Fill( HT, totalWeight );
						histos1D_[ "MET_cutBestPair" ]->Fill( MET, totalWeight );
						histos1D_[ "METHT_cutBestPair" ]->Fill( MET/HT, totalWeight );
						histos1D_[ "jetNum_cutBestPair" ]->Fill( numJets, totalWeight );
						histos1D_[ "jet1Pt_cutBestPair" ]->Fill( JETS[0].p4.Pt(), totalWeight );
						histos1D_[ "jet2Pt_cutBestPair" ]->Fill( JETS[1].p4.Pt(), totalWeight );
						histos1D_[ "jet3Pt_cutBestPair" ]->Fill( JETS[2].p4.Pt(), totalWeight );
						histos1D_[ "jet4Pt_cutBestPair" ]->Fill( JETS[3].p4.Pt(), totalWeight );
						histos1D_[ "jet1QGL_cutBestPair" ]->Fill( JETS[0].qgl, totalWeight );
						histos1D_[ "jet2QGL_cutBestPair" ]->Fill( JETS[1].qgl, totalWeight );
						histos1D_[ "jet3QGL_cutBestPair" ]->Fill( JETS[2].qgl, totalWeight );
						histos1D_[ "jet4QGL_cutBestPair" ]->Fill( JETS[3].qgl, totalWeight );
						histos1D_[ "NPV_cutBestPair" ]->Fill( numPV, totalWeight );

						histos1D_[ "massAve_cutBestPair" ]->Fill( massAve, totalWeight );
						histos1D_[ "massAsym_cutBestPair" ]->Fill( massAsym, totalWeight );
						histos1D_[ "deltaEta_cutBestPair" ]->Fill( deltaEta, totalWeight );
						histos1D_[ "deltaR_cutBestPair" ]->Fill( deltaR , totalWeight );
						histos1D_[ "cosThetaStar1_cutBestPair" ]->Fill( cosThetaStar1 , totalWeight );
						histos1D_[ "cosThetaStar2_cutBestPair" ]->Fill( cosThetaStar2 , totalWeight );
						histos2D_[ "deltavsMassAve_cutBestPair" ]->Fill( massAve, delta1 , totalWeight );
						histos2D_[ "deltavsMassAve_cutBestPair" ]->Fill( massAve, delta2 , totalWeight );
						histos2D_[ "dijetsEta_cutBestPair" ]->Fill( eta1, eta2, totalWeight );

						if ( massAsym < cutAK4MassAsym ) { 

							cutmap["MassAsym"] += 1;
							histos1D_[ "HT_cutMassAsym" ]->Fill( HT, totalWeight );
							histos1D_[ "jetNum_cutMassAsym" ]->Fill( numJets, totalWeight );
							histos1D_[ "jet1Pt_cutMassAsym" ]->Fill( JETS[0].p4.Pt(), totalWeight );
							histos1D_[ "jet2Pt_cutMassAsym" ]->Fill( JETS[1].p4.Pt(), totalWeight );
							histos1D_[ "jet3Pt_cutMassAsym" ]->Fill( JETS[2].p4.Pt(), totalWeight );
							histos1D_[ "jet4Pt_cutMassAsym" ]->Fill( JETS[3].p4.Pt(), totalWeight );
							histos1D_[ "massAve_cutMassAsym" ]->Fill( massAve, totalWeight );
							histos1D_[ "massAsym_cutMassAsym" ]->Fill( massAsym, totalWeight );
							histos1D_[ "deltaEta_cutMassAsym" ]->Fill( deltaEta, totalWeight );
							histos2D_[ "deltavsMassAve_cutMassAsym" ]->Fill( massAve, delta1 , totalWeight );
							histos2D_[ "deltavsMassAve_cutMassAsym" ]->Fill( massAve, delta2 , totalWeight );
							histos2D_[ "dijetsEta_cutMassAsym" ]->Fill( eta1, eta2, totalWeight );
						
							if ( deltaEta <  cutDeltaEtaDijetSyst ) {
								cutmap["DeltaEtaDijetSyst"] += 1;
								histos1D_[ "HT_cutDeltaEtaDijetSyst" ]->Fill( HT, totalWeight );
								histos1D_[ "jetNum_cutDeltaEtaDijetSyst" ]->Fill( numJets, totalWeight );
								histos1D_[ "jet1Pt_cutDeltaEtaDijetSyst" ]->Fill( JETS[0].p4.Pt(), totalWeight );
								histos1D_[ "jet2Pt_cutDeltaEtaDijetSyst" ]->Fill( JETS[1].p4.Pt(), totalWeight );
								histos1D_[ "jet3Pt_cutDeltaEtaDijetSyst" ]->Fill( JETS[2].p4.Pt(), totalWeight );
								histos1D_[ "jet4Pt_cutDeltaEtaDijetSyst" ]->Fill( JETS[3].p4.Pt(), totalWeight );
								histos1D_[ "massAve_cutDeltaEtaDijetSyst" ]->Fill( massAve, totalWeight );
								histos1D_[ "massAsym_cutDeltaEtaDijetSyst" ]->Fill( massAsym, totalWeight );
								histos1D_[ "deltaEta_cutDeltaEtaDijetSyst" ]->Fill( deltaEta, totalWeight );
								histos2D_[ "deltavsMassAve_cutDeltaEtaDijetSyst" ]->Fill( massAve, delta1 , totalWeight );
								histos2D_[ "deltavsMassAve_cutDeltaEtaDijetSyst" ]->Fill( massAve, delta2 , totalWeight );
								histos2D_[ "dijetsEta_cutDeltaEtaDijetSyst" ]->Fill( eta1, eta2, totalWeight );
								
								if ( ( delta1 > cutDelta ) && ( delta2  > cutDelta ) ) {
									cutmap["Delta"] += 1;
									histos1D_[ "HT_cutDelta" ]->Fill( HT, totalWeight );
									histos1D_[ "jetNum_cutDelta" ]->Fill( numJets, totalWeight );
									histos1D_[ "jet1Pt_cutDelta" ]->Fill( JETS[0].p4.Pt(), totalWeight );
									histos1D_[ "jet2Pt_cutDelta" ]->Fill( JETS[1].p4.Pt(), totalWeight );
									histos1D_[ "jet3Pt_cutDelta" ]->Fill( JETS[2].p4.Pt(), totalWeight );
									histos1D_[ "jet4Pt_cutDelta" ]->Fill( JETS[3].p4.Pt(), totalWeight );
									histos1D_[ "massAve_cutDelta" ]->Fill( massAve, totalWeight );
									histos1D_[ "massAsym_cutDelta" ]->Fill( massAsym, totalWeight );
									histos1D_[ "deltaEta_cutDelta" ]->Fill( deltaEta, totalWeight );
									histos2D_[ "deltavsMassAve_cutDelta" ]->Fill( massAve, delta1 , totalWeight );
									histos2D_[ "deltavsMassAve_cutDelta" ]->Fill( massAve, delta2 , totalWeight );
									histos2D_[ "dijetsEta_cutDelta" ]->Fill( eta1, eta2, totalWeight );

									if ( btagCSVv2Pair12 && btagCSVv2Pair34 ) { 

										cutmap["Btag"] += 1;
										histos1D_[ "HT_cutBtag" ]->Fill( HT, totalWeight );
										histos1D_[ "jetNum_cutBtag" ]->Fill( numJets, totalWeight );
										histos1D_[ "jet1Pt_cutBtag" ]->Fill( JETS[0].p4.Pt(), totalWeight );
										histos1D_[ "jet2Pt_cutBtag" ]->Fill( JETS[1].p4.Pt(), totalWeight );
										histos1D_[ "jet3Pt_cutBtag" ]->Fill( JETS[2].p4.Pt(), totalWeight );
										histos1D_[ "jet4Pt_cutBtag" ]->Fill( JETS[3].p4.Pt(), totalWeight );
										histos1D_[ "massAve_cutBtag" ]->Fill( massAve, totalWeight );
										histos1D_[ "massAsym_cutBtag" ]->Fill( massAsym, totalWeight );
										histos1D_[ "deltaEta_cutBtag" ]->Fill( deltaEta, totalWeight );
										histos2D_[ "deltavsMassAve_cutBtag" ]->Fill( massAve, delta1 , totalWeight );
										histos2D_[ "deltavsMassAve_cutBtag" ]->Fill( massAve, delta2 , totalWeight );
										histos2D_[ "dijetsEta_cutBtag" ]->Fill( eta1, eta2, totalWeight );
									}
								}
							}
						}
					}
				}
			}
		}
	}
	JETS.clear();
}


// ------------ method called once each job just before starting event loop  ------------
void RUNResolvedResolutionCalc::beginJob() {

	// Calculate PUWeight
	PUWeight_.generateWeights( dataPUFile );

	histos1D_[ "rawJetPt" ] = fs_->make< TH1D >( "rawJetPt", "rawJetPt", 100, 0., 1000. );
	histos1D_[ "rawJetPt" ]->Sumw2();
	histos1D_[ "rawHT" ] = fs_->make< TH1D >( "rawHT", "rawHT", 300, 0., 3000. );
	histos1D_[ "rawHT" ]->Sumw2();

	histos1D_[ "jetPt" ] = fs_->make< TH1D >( "jetPt", "jetPt", 100, 0., 1000. );
	histos1D_[ "jetPt" ]->Sumw2();
	histos1D_[ "jetEta" ] = fs_->make< TH1D >( "jetEta", "jetEta", 100, -5., 5. );
	histos1D_[ "jetEta" ]->Sumw2();
	histos1D_[ "jetNum" ] = fs_->make< TH1D >( "jetNum", "jetNum", 10, 0., 10. );
	histos1D_[ "jetNum" ]->Sumw2();
	histos1D_[ "HT" ] = fs_->make< TH1D >( "HT", "HT", 300, 0., 3000. );
	histos1D_[ "HT" ]->Sumw2();
	histos1D_[ "NPV" ] = fs_->make< TH1D >( "NPV", "NPV", 80, 0., 80. );
	histos1D_[ "NPV" ]->Sumw2();
	histos1D_[ "NPV_NOPUWeight" ] = fs_->make< TH1D >( "NPV_NOPUWeight", "NPV_NOPUWeight", 80, 0., 80. );
	histos1D_[ "NPV_NOPUWeight" ]->Sumw2();
	histos1D_[ "PUWeight" ] = fs_->make< TH1D >( "PUWeight", "PUWeight", 50, 0., 5. );
	histos1D_[ "PUWeight" ]->Sumw2();
	histos1D_[ "neutralHadronEnergyFrac" ] = fs_->make< TH1D >( "neutralHadronEnergyFrac", "neutralHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "neutralHadronEnergyFrac" ]->Sumw2();
	histos1D_[ "neutralEmEnergyFrac" ] = fs_->make< TH1D >( "neutralEmEnergyFrac", "neutralEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "neutralEmEnergyFrac" ]->Sumw2();
	histos1D_[ "chargedHadronEnergyFrac" ] = fs_->make< TH1D >( "chargedHadronEnergyFrac", "chargedHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedHadronEnergyFrac" ]->Sumw2();
	histos1D_[ "chargedEmEnergyFrac" ] = fs_->make< TH1D >( "chargedEmEnergyFrac", "chargedEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedEmEnergyFrac" ]->Sumw2();
	histos1D_[ "chargedMultiplicity" ] = fs_->make< TH1D >( "chargedMultiplicity", "chargedMultiplicity", 50, 0., 1. );
	histos1D_[ "chargedMultiplicity" ]->Sumw2();
	histos1D_[ "numConst" ] = fs_->make< TH1D >( "numConst", "numConst", 100, 0., 100. );
	histos1D_[ "numConst" ]->Sumw2();



	histos1D_[ "jet1Pt_cut4Jets" ] = fs_->make< TH1D >( "jet1Pt_cut4Jets", "jet1Pt_cut4Jets", 100, 0., 1000. );
	histos1D_[ "jet1Pt_cut4Jets" ]->Sumw2();
	histos1D_[ "jet2Pt_cut4Jets" ] = fs_->make< TH1D >( "jet2Pt_cut4Jets", "jet2Pt_cut4Jets", 100, 0., 1000. );
	histos1D_[ "jet2Pt_cut4Jets" ]->Sumw2();
	histos1D_[ "jet3Pt_cut4Jets" ] = fs_->make< TH1D >( "jet3Pt_cut4Jets", "jet3Pt_cut4Jets", 100, 0., 1000. );
	histos1D_[ "jet3Pt_cut4Jets" ]->Sumw2();
	histos1D_[ "jet4Pt_cut4Jets" ] = fs_->make< TH1D >( "jet4Pt_cut4Jets", "jet4Pt_cut4Jets", 100, 0., 1000. );
	histos1D_[ "jet4Pt_cut4Jets" ]->Sumw2();
	histos1D_[ "jetNum_cut4Jets" ] = fs_->make< TH1D >( "jetNum_cut4Jets", "jetNum_cut4Jets", 10, 0., 10. );
	histos1D_[ "jetNum_cut4Jets" ]->Sumw2();
	histos1D_[ "HT_cut4Jets" ] = fs_->make< TH1D >( "HT_cut4Jets", "HT_cut4Jets", 300, 0., 3000. );
	histos1D_[ "HT_cut4Jets" ]->Sumw2();
	histos1D_[ "MET_cut4Jets" ] = fs_->make< TH1D >( "MET_cut4Jets", "MET_cut4Jets", 20, 0., 200. );
	histos1D_[ "MET_cut4Jets" ]->Sumw2();
	histos1D_[ "METHT_cut4Jets" ] = fs_->make< TH1D >( "METHT_cut4Jets", "METHT_cut4Jets", 50, 0., 1. );
	histos1D_[ "METHT_cut4Jets" ]->Sumw2();

	histos1D_[ "jet1Pt_cutBestPair" ] = fs_->make< TH1D >( "jet1Pt_cutBestPair", "jet1Pt_cutBestPair", 100, 0., 1000. );
	histos1D_[ "jet1Pt_cutBestPair" ]->Sumw2();
	histos1D_[ "jet2Pt_cutBestPair" ] = fs_->make< TH1D >( "jet2Pt_cutBestPair", "jet2Pt_cutBestPair", 100, 0., 1000. );
	histos1D_[ "jet2Pt_cutBestPair" ]->Sumw2();
	histos1D_[ "jet3Pt_cutBestPair" ] = fs_->make< TH1D >( "jet3Pt_cutBestPair", "jet3Pt_cutBestPair", 100, 0., 1000. );
	histos1D_[ "jet3Pt_cutBestPair" ]->Sumw2();
	histos1D_[ "jet4Pt_cutBestPair" ] = fs_->make< TH1D >( "jet4Pt_cutBestPair", "jet4Pt_cutBestPair", 100, 0., 1000. );
	histos1D_[ "jet4Pt_cutBestPair" ]->Sumw2();
	histos1D_[ "jetNum_cutBestPair" ] = fs_->make< TH1D >( "jetNum_cutBestPair", "jetNum_cutBestPair", 10, 0., 10. );
	histos1D_[ "jetNum_cutBestPair" ]->Sumw2();
	histos1D_[ "HT_cutBestPair" ] = fs_->make< TH1D >( "HT_cutBestPair", "HT_cutBestPair", 300, 0., 3000. );
	histos1D_[ "HT_cutBestPair" ]->Sumw2();
	histos1D_[ "MET_cutBestPair" ] = fs_->make< TH1D >( "MET_cutBestPair", "MET_cutBestPair", 20, 0., 200. );
	histos1D_[ "MET_cutBestPair" ]->Sumw2();
	histos1D_[ "METHT_cutBestPair" ] = fs_->make< TH1D >( "METHT_cutBestPair", "METHT_cutBestPair", 50, 0., 1. );
	histos1D_[ "METHT_cutBestPair" ]->Sumw2();
	histos1D_[ "NPV_cutBestPair" ] = fs_->make< TH1D >( "NPV_cutBestPair", "NPV_cutBestPair", 80, 0., 80. );
	histos1D_[ "NPV_cutBestPair" ]->Sumw2();
	histos1D_[ "neutralHadronEnergyFrac_cutBestPair" ] = fs_->make< TH1D >( "neutralHadronEnergyFrac_cutBestPair", "neutralHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "neutralHadronEnergyFrac_cutBestPair" ]->Sumw2();
	histos1D_[ "neutralEmEnergyFrac_cutBestPair" ] = fs_->make< TH1D >( "neutralEmEnergyFrac_cutBestPair", "neutralEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "neutralEmEnergyFrac_cutBestPair" ]->Sumw2();
	histos1D_[ "chargedHadronEnergyFrac_cutBestPair" ] = fs_->make< TH1D >( "chargedHadronEnergyFrac_cutBestPair", "chargedHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedHadronEnergyFrac_cutBestPair" ]->Sumw2();
	histos1D_[ "chargedEmEnergyFrac_cutBestPair" ] = fs_->make< TH1D >( "chargedEmEnergyFrac_cutBestPair", "chargedEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedEmEnergyFrac_cutBestPair" ]->Sumw2();
	histos1D_[ "chargedMultiplicity_cutBestPair" ] = fs_->make< TH1D >( "chargedMultiplicity_cutBestPair", "chargedMultiplicity", 50, 0., 1. );
	histos1D_[ "chargedMultiplicity_cutBestPair" ]->Sumw2();
	histos1D_[ "numConst_cutBestPair" ] = fs_->make< TH1D >( "numConst_cutBestPair", "numConst", 100, 0., 100. );
	histos1D_[ "numConst_cutBestPair" ]->Sumw2();
	histos1D_[ "jet1QGL_cutBestPair" ] = fs_->make< TH1D >( "jet1QGL_cutBestPair", "jet1QGL_cutBestPair", 10, 0., 1. );
	histos1D_[ "jet1QGL_cutBestPair" ]->Sumw2();
	histos1D_[ "jet2QGL_cutBestPair" ] = fs_->make< TH1D >( "jet2QGL_cutBestPair", "jet2QGL_cutBestPair", 10, 0., 1. );
	histos1D_[ "jet2QGL_cutBestPair" ]->Sumw2();
	histos1D_[ "jet3QGL_cutBestPair" ] = fs_->make< TH1D >( "jet3QGL_cutBestPair", "jet3QGL_cutBestPair", 10, 0., 1. );
	histos1D_[ "jet3QGL_cutBestPair" ]->Sumw2();
	histos1D_[ "jet4QGL_cutBestPair" ] = fs_->make< TH1D >( "jet4QGL_cutBestPair", "jet4QGL_cutBestPair", 10, 0., 1. );
	histos1D_[ "jet4QGL_cutBestPair" ]->Sumw2();
	histos1D_[ "massAve_cutBestPair" ] = fs_->make< TH1D >( "massAve_cutBestPair", "massAve_cutBestPair", 2000, 0., 2000.);
	histos1D_[ "massAve_cutBestPair" ]->Sumw2();
	histos1D_[ "massAsym_cutBestPair" ] = fs_->make< TH1D >( "massAsym_cutBestPair", "massAsym_cutBestPair", 50, 0., 2. );
	histos1D_[ "massAsym_cutBestPair" ]->Sumw2();
	histos1D_[ "deltaEta_cutBestPair" ] = fs_->make< TH1D >( "deltaEta_cutBestPair", "deltaEta_cutBestPair", 50, 0., 10. );
	histos1D_[ "deltaEta_cutBestPair" ]->Sumw2();
	histos1D_[ "deltaR_cutBestPair" ] = fs_->make< TH1D >( "deltaR_cutBestPair", "deltaR_cutBestPair", 50, 0., 5. );
	histos1D_[ "deltaR_cutBestPair" ]->Sumw2();
	histos1D_[ "cosThetaStar1_cutBestPair" ] = fs_->make< TH1D >( "cosThetaStar1_cutBestPair", "cosThetaStar1_cutBestPair", 10, 0., 1. );
	histos1D_[ "cosThetaStar1_cutBestPair" ]->Sumw2();
	histos1D_[ "cosThetaStar2_cutBestPair" ] = fs_->make< TH1D >( "cosThetaStar2_cutBestPair", "cosThetaStar2_cutBestPair", 10, 0., 1. );
	histos1D_[ "cosThetaStar2_cutBestPair" ]->Sumw2();
	histos2D_[ "deltavsMassAve_cutBestPair" ] = fs_->make< TH2D >( "deltavsMassAve_cutBestPair", "deltavsMassAve_cutBestPair", 200, 0., 2000.,  300, -500., 1000. );
	histos2D_[ "deltavsMassAve_cutBestPair" ]->Sumw2();
	histos2D_[ "dijetsEta_cutBestPair" ] = fs_->make< TH2D >( "dijetsEta_cutBestPair", "dijetsEta_cutBestPair", 48, -3., 3., 48, -3., 3. );
	histos2D_[ "dijetsEta_cutBestPair" ]->Sumw2();

	histos1D_[ "jet1Pt_cutMassAsym" ] = fs_->make< TH1D >( "jet1Pt_cutMassAsym", "jet1Pt_cutMassAsym", 100, 0., 1000. );
	histos1D_[ "jet1Pt_cutMassAsym" ]->Sumw2();
	histos1D_[ "jet2Pt_cutMassAsym" ] = fs_->make< TH1D >( "jet2Pt_cutMassAsym", "jet2Pt_cutMassAsym", 100, 0., 1000. );
	histos1D_[ "jet2Pt_cutMassAsym" ]->Sumw2();
	histos1D_[ "jet3Pt_cutMassAsym" ] = fs_->make< TH1D >( "jet3Pt_cutMassAsym", "jet3Pt_cutMassAsym", 100, 0., 1000. );
	histos1D_[ "jet3Pt_cutMassAsym" ]->Sumw2();
	histos1D_[ "jet4Pt_cutMassAsym" ] = fs_->make< TH1D >( "jet4Pt_cutMassAsym", "jet4Pt_cutMassAsym", 100, 0., 1000. );
	histos1D_[ "jet4Pt_cutMassAsym" ]->Sumw2();
	histos1D_[ "jetNum_cutMassAsym" ] = fs_->make< TH1D >( "jetNum_cutMassAsym", "jetNum_cutMassAsym", 10, 0., 10. );
	histos1D_[ "jetNum_cutMassAsym" ]->Sumw2();
	histos1D_[ "HT_cutMassAsym" ] = fs_->make< TH1D >( "HT_cutMassAsym", "HT_cutMassAsym", 300, 0., 3000. );
	histos1D_[ "HT_cutMassAsym" ]->Sumw2();
	histos1D_[ "massAve_cutMassAsym" ] = fs_->make< TH1D >( "massAve_cutMassAsym", "massAve_cutMassAsym", 2000, 0., 2000.);
	histos1D_[ "massAve_cutMassAsym" ]->Sumw2();
	histos1D_[ "massAsym_cutMassAsym" ] = fs_->make< TH1D >( "massAsym_cutMassAsym", "massAsym_cutMassAsym", 50, 0., 2. );
	histos1D_[ "massAsym_cutMassAsym" ]->Sumw2();
	histos1D_[ "deltaEta_cutMassAsym" ] = fs_->make< TH1D >( "deltaEta_cutMassAsym", "deltaEta_cutMassAsym", 50, 0., 10. );
	histos1D_[ "deltaEta_cutMassAsym" ]->Sumw2();
	histos2D_[ "deltavsMassAve_cutMassAsym" ] = fs_->make< TH2D >( "deltavsMassAve_cutMassAsym", "deltavsMassAve_cutMassAsym", 200, 0., 2000.,  300, -500., 1000. );
	histos2D_[ "deltavsMassAve_cutMassAsym" ]->Sumw2();
	histos2D_[ "dijetsEta_cutMassAsym" ] = fs_->make< TH2D >( "dijetsEta_cutMassAsym", "dijetsEta_cutMassAsym", 48, -3., 3., 48, -3., 3. );
	histos2D_[ "dijetsEta_cutMassAsym" ]->Sumw2();

	histos1D_[ "jet1Pt_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "jet1Pt_cutDeltaEtaDijetSyst", "jet1Pt_cutDeltaEtaDijetSyst", 100, 0., 1000. );
	histos1D_[ "jet1Pt_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "jet2Pt_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "jet2Pt_cutDeltaEtaDijetSyst", "jet2Pt_cutDeltaEtaDijetSyst", 100, 0., 1000. );
	histos1D_[ "jet2Pt_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "jet3Pt_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "jet3Pt_cutDeltaEtaDijetSyst", "jet3Pt_cutDeltaEtaDijetSyst", 100, 0., 1000. );
	histos1D_[ "jet3Pt_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "jet4Pt_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "jet4Pt_cutDeltaEtaDijetSyst", "jet4Pt_cutDeltaEtaDijetSyst", 100, 0., 1000. );
	histos1D_[ "jet4Pt_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "jetNum_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "jetNum_cutDeltaEtaDijetSyst", "jetNum_cutDeltaEtaDijetSyst", 10, 0., 10. );
	histos1D_[ "jetNum_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "HT_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "HT_cutDeltaEtaDijetSyst", "HT_cutDeltaEtaDijetSyst", 300, 0., 3000. );
	histos1D_[ "HT_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "massAve_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "massAve_cutDeltaEtaDijetSyst", "massAve_cutDeltaEtaDijetSyst", 2000, 0., 2000.);
	histos1D_[ "massAve_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "massAsym_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "massAsym_cutDeltaEtaDijetSyst", "massAsym_cutDeltaEtaDijetSyst", 50, 0., 2. );
	histos1D_[ "massAsym_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "deltaEta_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "deltaEta_cutDeltaEtaDijetSyst", "deltaEta_cutDeltaEtaDijetSyst", 50, 0., 10. );
	histos1D_[ "deltaEta_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos1D_[ "minDeltaEtaDijetSystR_cutDeltaEtaDijetSyst" ] = fs_->make< TH1D >( "minDeltaEtaDijetSystR_cutDeltaEtaDijetSyst", "minDeltaEtaDijetSystR_cutDeltaEtaDijetSyst", 50, 0., 5. );
	histos1D_[ "minDeltaEtaDijetSystR_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos2D_[ "deltavsMassAve_cutDeltaEtaDijetSyst" ] = fs_->make< TH2D >( "deltavsMassAve_cutDeltaEtaDijetSyst", "deltavsMassAve_cutDeltaEtaDijetSyst", 200, 0., 2000.,  300, -500., 1000. );
	histos2D_[ "deltavsMassAve_cutDeltaEtaDijetSyst" ]->Sumw2();
	histos2D_[ "dijetsEta_cutDeltaEtaDijetSyst" ] = fs_->make< TH2D >( "dijetsEta_cutDeltaEtaDijetSyst", "dijetsEta_cutDeltaEtaDijetSyst", 48, -3., 3., 48, -3., 3. );
	histos2D_[ "dijetsEta_cutDeltaEtaDijetSyst" ]->Sumw2();
	
	
	histos1D_[ "jet1Pt_cutDelta" ] = fs_->make< TH1D >( "jet1Pt_cutDelta", "jet1Pt_cutDelta", 100, 0., 1000. );
	histos1D_[ "jet1Pt_cutDelta" ]->Sumw2();
	histos1D_[ "jet2Pt_cutDelta" ] = fs_->make< TH1D >( "jet2Pt_cutDelta", "jet2Pt_cutDelta", 100, 0., 1000. );
	histos1D_[ "jet2Pt_cutDelta" ]->Sumw2();
	histos1D_[ "jet3Pt_cutDelta" ] = fs_->make< TH1D >( "jet3Pt_cutDelta", "jet3Pt_cutDelta", 100, 0., 1000. );
	histos1D_[ "jet3Pt_cutDelta" ]->Sumw2();
	histos1D_[ "jet4Pt_cutDelta" ] = fs_->make< TH1D >( "jet4Pt_cutDelta", "jet4Pt_cutDelta", 100, 0., 1000. );
	histos1D_[ "jet4Pt_cutDelta" ]->Sumw2();
	histos1D_[ "jetNum_cutDelta" ] = fs_->make< TH1D >( "jetNum_cutDelta", "jetNum_cutDelta", 10, 0., 10. );
	histos1D_[ "jetNum_cutDelta" ]->Sumw2();
	histos1D_[ "HT_cutDelta" ] = fs_->make< TH1D >( "HT_cutDelta", "HT_cutDelta", 300, 0., 3000. );
	histos1D_[ "HT_cutDelta" ]->Sumw2();
	histos1D_[ "massAve_cutDelta" ] = fs_->make< TH1D >( "massAve_cutDelta", "massAve_cutDelta", 2000, 0., 2000.);
	histos1D_[ "massAve_cutDelta" ]->Sumw2();
	histos1D_[ "massAsym_cutDelta" ] = fs_->make< TH1D >( "massAsym_cutDelta", "massAsym_cutDelta", 50, 0., 2. );
	histos1D_[ "massAsym_cutDelta" ]->Sumw2();
	histos1D_[ "deltaEta_cutDelta" ] = fs_->make< TH1D >( "deltaEta_cutDelta", "deltaEta_cutDelta", 50, 0., 10. );
	histos1D_[ "deltaEta_cutDelta" ]->Sumw2();
	histos2D_[ "deltavsMassAve_cutDelta" ] = fs_->make< TH2D >( "deltavsMassAve_cutDelta", "deltavsMassAve_cutDelta", 200, 0., 2000.,  300, -500., 1000. );
	histos2D_[ "deltavsMassAve_cutDelta" ]->Sumw2();
	histos2D_[ "dijetsEta_cutDelta" ] = fs_->make< TH2D >( "dijetsEta_cutDelta", "dijetsEta_cutDelta", 48, -3., 3., 48, -3., 3. );
	histos2D_[ "dijetsEta_cutDelta" ]->Sumw2();

	histos1D_[ "jet1Pt_cutBtag" ] = fs_->make< TH1D >( "jet1Pt_cutBtag", "jet1Pt_cutBtag", 100, 0., 1000. );
	histos1D_[ "jet1Pt_cutBtag" ]->Sumw2();
	histos1D_[ "jet2Pt_cutBtag" ] = fs_->make< TH1D >( "jet2Pt_cutBtag", "jet2Pt_cutBtag", 100, 0., 1000. );
	histos1D_[ "jet2Pt_cutBtag" ]->Sumw2();
	histos1D_[ "jet3Pt_cutBtag" ] = fs_->make< TH1D >( "jet3Pt_cutBtag", "jet3Pt_cutBtag", 100, 0., 1000. );
	histos1D_[ "jet3Pt_cutBtag" ]->Sumw2();
	histos1D_[ "jet4Pt_cutBtag" ] = fs_->make< TH1D >( "jet4Pt_cutBtag", "jet4Pt_cutBtag", 100, 0., 1000. );
	histos1D_[ "jet4Pt_cutBtag" ]->Sumw2();
	histos1D_[ "jetNum_cutBtag" ] = fs_->make< TH1D >( "jetNum_cutBtag", "jetNum_cutBtag", 10, 0., 10. );
	histos1D_[ "jetNum_cutBtag" ]->Sumw2();
	histos1D_[ "HT_cutBtag" ] = fs_->make< TH1D >( "HT_cutBtag", "HT_cutBtag", 300, 0., 3000. );
	histos1D_[ "HT_cutBtag" ]->Sumw2();
	histos1D_[ "massAve_cutBtag" ] = fs_->make< TH1D >( "massAve_cutBtag", "massAve_cutBtag", 2000, 0., 2000.);
	histos1D_[ "massAve_cutBtag" ]->Sumw2();
	histos1D_[ "massAsym_cutBtag" ] = fs_->make< TH1D >( "massAsym_cutBtag", "massAsym_cutBtag", 50, 0., 2. );
	histos1D_[ "massAsym_cutBtag" ]->Sumw2();
	histos1D_[ "deltaEta_cutBtag" ] = fs_->make< TH1D >( "deltaEta_cutBtag", "deltaEta_cutBtag", 50, 0., 10. );
	histos1D_[ "deltaEta_cutBtag" ]->Sumw2();
	histos2D_[ "deltavsMassAve_cutBtag" ] = fs_->make< TH2D >( "deltavsMassAve_cutBtag", "deltavsMassAve_cutBtag", 200, 0., 2000.,  300, -500., 1000. );
	histos2D_[ "deltavsMassAve_cutBtag" ]->Sumw2();
	histos2D_[ "dijetsEta_cutBtag" ] = fs_->make< TH2D >( "dijetsEta_cutBtag", "dijetsEta_cutBtag", 48, -3., 3., 48, -3., 3. );
	histos2D_[ "dijetsEta_cutBtag" ]->Sumw2();


	cutLabels.push_back("Processed");
	cutLabels.push_back("4Jets");
	cutLabels.push_back("BestPair");
	cutLabels.push_back("MassAsym");
	cutLabels.push_back("Btag");
	cutLabels.push_back("DeltaEtaDijetSyst");
	cutLabels.push_back("Delta");
	histos1D_[ "hcutflow" ] = fs_->make< TH1D >("cutflow","cut flow", cutLabels.size(), 0.5, cutLabels.size() +0.5 );
	histos1D_[ "hcutflow" ]->Sumw2();
	for( const string &ivec : cutLabels ) cutmap[ ivec ] = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void RUNResolvedResolutionCalc::endJob() {

	int ibin = 1;
	for( const string &ivec : cutLabels ) {
		histos1D_["hcutflow"]->SetBinContent( ibin, cutmap[ ivec ] );
		histos1D_["hcutflow"]->GetXaxis()->SetBinLabel( ibin, ivec.c_str() );
		ibin++;
	}

}

void RUNResolvedResolutionCalc::clearVariables() {

	jetsPt->clear();
	jetsEta->clear();
	jetsPhi->clear();
	jetsE->clear();
	jetsQGL->clear();
	jetsCSVv2->clear();
	jetsCSVv2SF->clear();
	jetsCMVAv2->clear();

}

void RUNResolvedResolutionCalc::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<double>("cutAK4jetPt", 80);
	desc.add<double>("cutAK4HT", 800);
	desc.add<double>("cutAK4MassAsym", 0.20);
	desc.add<double>("cutDelta", 180);
	desc.add<double>("cutDeltaEtaDijetSyst", .75);
	desc.add<string>("dataPUFile", "supportFiles/PileupData2015D_JSON_10-23-2015.root");
	desc.add<string>("jecVersion", "supportFiles/Summer15_25nsV6");

	desc.add<InputTag>("Lumi", 	InputTag("eventInfo:evtInfoLumiBlock"));
	desc.add<InputTag>("Run", 	InputTag("eventInfo:evtInfoRunNumber"));
	desc.add<InputTag>("Event", 	InputTag("eventInfo:evtInfoEventNumber"));
	desc.add<InputTag>("generator", 	InputTag("generator"));
	desc.add<InputTag>("genParticles", 	InputTag("filteredPrunedGenParticles"));
	desc.add<InputTag>("bunchCross", 	InputTag("eventUserData:puBX"));
	desc.add<InputTag>("rho", 	InputTag("vertexInfo:rho"));
	desc.add<InputTag>("puNumInt", 	InputTag("eventUserData:puNInt"));
	desc.add<InputTag>("trueNInt", 	InputTag("eventUserData:puNtrueInt"));
	desc.add<InputTag>("jetPt", 	InputTag("jetsAK4CHS:jetAK4CHSPt"));
	desc.add<InputTag>("jetEta", 	InputTag("jetsAK4CHS:jetAK4CHSEta"));
	desc.add<InputTag>("jetPhi", 	InputTag("jetsAK4CHS:jetAK4CHSPhi"));
	desc.add<InputTag>("jetE", 	InputTag("jetsAK4CHS:jetAK4CHSE"));
	desc.add<InputTag>("jetQGL", 	InputTag("jetsAK4CHS:jetAK4CHSQGL"));
	desc.add<InputTag>("jetCSVv2", 	InputTag("jetsAK4CHS:jetAK4CHSCSVv2"));
	desc.add<InputTag>("jetCMVAv2", 	InputTag("jetsAK4CHS:jetAK4CHSCMVAv2"));
	desc.add<InputTag>("jetArea", 	InputTag("jetsAK4CHS:jetAK4CHSjetArea"));
	desc.add<InputTag>("jetGenPt", 	InputTag("jetsAK4CHS:jetAK4CHSGenJetPt"));
	desc.add<InputTag>("jetGenEta", 	InputTag("jetsAK4CHS:jetAK4CHSGenJetEta"));
	desc.add<InputTag>("jetGenPhi", 	InputTag("jetsAK4CHS:jetAK4CHSGenJetPhi"));
	desc.add<InputTag>("jetGenE", 	InputTag("jetsAK4CHS:jetAK4CHSGenJetE"));
	desc.add<InputTag>("jetHadronFlavour", 	InputTag("jetsAK4CHS:jetAK4CHSHadronFlavour"));
	desc.add<InputTag>("NPV", 	InputTag("vertexInfo:npv"));
	desc.add<InputTag>("metPt", 	InputTag("metFull:metFullPt"));
	// JetID
	desc.add<InputTag>("jecFactor", 		InputTag("jetsAK4CHS:jetAK4CHSjecFactor0"));
	desc.add<InputTag>("neutralHadronEnergyFrac", 	InputTag("jetsAK4CHS:jetAK4CHSneutralHadronEnergyFrac"));
	desc.add<InputTag>("neutralEmEnergyFrac", 	InputTag("jetsAK4CHS:jetAK4CHSneutralEmEnergyFrac"));
	desc.add<InputTag>("chargedEmEnergyFrac", 	InputTag("jetsAK4CHS:jetAK4CHSchargedEmEnergyFrac"));
	desc.add<InputTag>("muonEnergy", 		InputTag("jetsAK4CHS:jetAK4CHSMuonEnergy"));
	desc.add<InputTag>("chargedHadronEnergyFrac", 	InputTag("jetsAK4CHS:jetAK4CHSchargedHadronEnergyFrac"));
	desc.add<InputTag>("neutralMultiplicity",	InputTag("jetsAK4CHS:jetAK4CHSneutralMultiplicity"));
	desc.add<InputTag>("chargedMultiplicity", 	InputTag("jetsAK4CHS:jetAK4CHSchargedMultiplicity"));
	// Muons
	desc.add<InputTag>("muonPt", 		InputTag("muons:muPt"));
	desc.add<InputTag>("muonEta", 		InputTag("muons:muEta"));
	desc.add<InputTag>("muonIsLoose", 		InputTag("muons:muIsLooseMuon"));
	desc.add<InputTag>("muonIsGlobal", 		InputTag("muons:muIsGlobalMuon"));
	// Electrons
	desc.add<InputTag>("elePt", 		InputTag("electrons:elPt"));
	desc.add<InputTag>("eleEta", 		InputTag("electrons:elEta"));
	desc.add<InputTag>("eleLoose", 		InputTag("electrons:elvidLoose"));
	descriptions.addDefault(desc);
}
      
/*void RUNResolvedResolutionCalc::beginRun(const Run& iRun, const EventSetup& iSetup){
}*/
      
/*void RUNResolvedResolutionCalc::endRun(const Run& iRun, const EventSetup& iSetup){
}*/

//define this as a plug-in
DEFINE_FWK_MODULE(RUNResolvedResolutionCalc);
