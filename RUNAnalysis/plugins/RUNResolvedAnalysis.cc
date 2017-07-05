// -*- C++ -*-
//
// Package:    RUNA/RUNAnalysis
// Class:      RUNResolvedAnalysis
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

//
// constants, enums and typedefs
//

//
// class declaration
//
class RUNResolvedAnalysis : public EDAnalyzer {
	public:
		explicit RUNResolvedAnalysis(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~RUNResolvedAnalysis();

	private:
		virtual void beginJob() override;
		virtual void analyze(const Event&, const EventSetup&) override;
		virtual void endJob() override;
		virtual void clearVariables();
		virtual void beginRun(const Run&, const EventSetup&) override;
		virtual void endRun(Run const&, EventSetup const&) override;

		//virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
		//virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

		// ----------member data ---------------------------
		PUReweighter PUWeight_;
		int lhaPdfId;

		Service<TFileService> fs_;
		TTree *RUNAtree;
		map< string, TH1D* > histos1D_;
		map< string, TH2D* > histos2D_;
		vector< string > cutLabels;
		map< string, double > cutmap;

		bool isData;
		bool LHEcont;
		bool mkTree;
		bool massPairing;
		string dataPUFile;
		string jecVersion;
		string btagCSVFile;
		TString systematics;
		double scale;
		double cutAK4jetPt;
		double cutAK4HT;
		double cutAK4MassAsym;
		double cutDelta;
		double cutDeltaEtaDijetSyst;

		vector<string> triggerPass, triggerNamesList;
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
		EDGetTokenT<LHEEventProduct> extLHEProducer_;

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
		EDGetTokenT<vector<float>> muonIsLoose_;
		EDGetTokenT<vector<float>> muonIsGlobal_;

		/// Electrons
		EDGetTokenT<vector<float>> elePt_;
		EDGetTokenT<vector<float>> eleEta_;
		EDGetTokenT<vector<float>> eleLoose_;

};

//
// static data member definitions
//

//
// constructors and destructor
//
RUNResolvedAnalysis::RUNResolvedAnalysis(const ParameterSet& iConfig):
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
	extLHEProducer_(consumes<LHEEventProduct>(iConfig.getParameter<InputTag>("extLHEProducer"))),
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
	muonIsLoose_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonIsLoose"))), 
	muonIsGlobal_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonIsGlobal"))),
	//Electrons
	elePt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("elePt"))), 
	eleEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("eleEta"))), 
	eleLoose_(consumes<vector<float>>(iConfig.getParameter<InputTag>("eleLoose"))) 
{
	consumes<LHERunInfoProduct,edm::InRun> (edm::InputTag("externalLHEProducer"));
	scale 		= iConfig.getParameter<double>("scale");
	cutAK4jetPt 	= iConfig.getParameter<double>("cutAK4jetPt");
	cutAK4HT	= iConfig.getParameter<double>("cutAK4HT");
	cutAK4MassAsym	= iConfig.getParameter<double>("cutAK4MassAsym");
	cutDelta        = iConfig.getParameter<double>("cutDelta");
	cutDeltaEtaDijetSyst	= iConfig.getParameter<double>("cutDeltaEtaDijetSyst");
	triggerPass 	= iConfig.getParameter<vector<string>>("triggerPass");
	isData 		= iConfig.getParameter<bool>("isData");
	LHEcont		= iConfig.getParameter<bool>("LHEcont");
	mkTree 		= iConfig.getParameter<bool>("mkTree");
	massPairing	= iConfig.getParameter<bool>("massPairing");
	dataPUFile 	= iConfig.getParameter<string>("dataPUFile");
	jecVersion 	= iConfig.getParameter<string>("jecVersion");
	btagCSVFile 	= iConfig.getParameter<string>("btagCSVFile");
	systematics 	= iConfig.getParameter<string>("systematics");

	/////// JECs
	string tmpPrefix = jecVersion;
	string prefix;
	if (isData) prefix = tmpPrefix + "_DATA_";
	else prefix = tmpPrefix + "_MC_";

	// all jet
	vector<string> jecPayloadNames_;
	jecPayloadNames_.push_back(prefix + "L1FastJet_AK4PFchs.txt");
	jecPayloadNames_.push_back(prefix + "L2Relative_AK4PFchs.txt");
	jecPayloadNames_.push_back(prefix + "L3Absolute_AK4PFchs.txt");
	if (isData) jecPayloadNames_.push_back(prefix + "L2L3Residual_AK4PFchs.txt");

	for ( vector<string>::const_iterator payloadBegin = jecPayloadNames_.begin(), payloadEnd = jecPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		jetPar.push_back(pars);
	}
	jetJEC = new FactorizedJetCorrector(jetPar);

	// jec uncertainty
	JetCorrectorParameters jecUncParam( prefix + "Uncertainty_AK4PFchs.txt");
	jetCorrUnc  = new JetCorrectionUncertainty( jecUncParam);

}


RUNResolvedAnalysis::~RUNResolvedAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNResolvedAnalysis::analyze(const Event& iEvent, const EventSetup& iSetup) {


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

	Handle<vector<float> > muonIsLoose;
	iEvent.getByToken(muonIsLoose_, muonIsLoose);

	Handle<vector<float> > muonIsGlobal;
	iEvent.getByToken(muonIsGlobal_, muonIsGlobal);

	/// Electrons
	Handle<vector<float> > elePt;
	iEvent.getByToken(elePt_, elePt);

	Handle<vector<float> > eleEta;
	iEvent.getByToken(eleEta_, eleEta);

	Handle<vector<float> > eleLoose;
	iEvent.getByToken(eleLoose_, eleLoose);

	////// Muon veto
	if ( muonPt->size() > 0 ) {
		for (size_t m = 0; m < muonPt->size(); m++) {
			if( ( (*muonIsLoose)[m] ) && ( (*muonPt)[m] > 10 ) && ( TMath::Abs((*muonEta)[m]) < 2.5 ) ) {
				//LogWarning("muon") << (*muonPt)[m] << " " << (*muonEta)[m] << " " << (*muonIsLoose)[m] << " " << (*muonIsGlobal)[m] ;
				muonsPt->push_back( (*muonPt)[m] );
				numMuon++;
			}
		}
	}
	///////////////////////////////////////////////////*/
	
	////// Electron veto
	if ( elePt->size() > 0 ) {
		for (size_t m = 0; m < elePt->size(); m++) {
			if( ( (*eleLoose)[m] ) && ( (*elePt)[m] > 10 ) && ( TMath::Abs((*eleEta)[m]) < 2.5 ) ){
				//LogWarning("ele") << (*elePt)[m] << " " << (*eleEta)[m] << " " << (*eleLoose)[m] ;
				elesPt->push_back( (*elePt)[m] );
				numEle++;
			}
		}
	}
	///////////////////////////////////////////////////*/

	///// LHE content
	if ( !isData && !LHEcont ) {

		// all this section is based on https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L215-L281
		/////////// GEN weight
		Handle<GenEventInfoProduct> genEvtInfo; 
		iEvent.getByToken( generator_, genEvtInfo );
		genWeight = genEvtInfo->weight();	
		//vector<double> evtWeights = genEvtInfo->weights();

		/////////// LHE weight
		Handle<LHEEventProduct> lheEvtInfo;
		iEvent.getByToken( extLHEProducer_, lheEvtInfo );
		// This function will return all the weights from the lheEvtInfo (scaleWeights, pdfWeights, alphaWeights)
		getWeights( lheEvtInfo, lhaPdfId, scaleWeights, pdfWeights, alphaWeights );

		//// Calculating the RMS from the pdfWeights, per event
		double sum = accumulate(pdfWeights.begin(), pdfWeights.end(), 0.0);
		double mean = sum / pdfWeights.size();
		vector<double> diff(pdfWeights.size());
		transform(pdfWeights.begin(), pdfWeights.end(), diff.begin(), [mean](double x) { return x - mean; });
		double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		pdfWeight = sqrt(sq_sum / pdfWeights.size());
		////////////////////////////////////////////////////
	}
	///////////////////////////////////////////////////*/


	////////// Check trigger fired
	bool ORTriggers = false;
	if ( isData ) ORTriggers = checkORListOfTriggerBits( triggerNamesList, triggerBit, triggerPrescale, triggerPass, false );
	else ORTriggers = true;
	///////////////////////////////////////////////////*/
	

	// PU Reweight
	if ( isData ) puWeight = 1;
	else puWeight = PUWeight_.getPUWeight( *trueNInt, *bunchCross );
	//puWeight = 1;
	histos1D_[ "PUWeight" ]->Fill( puWeight );
	lumiWeight = scale;
	double totalWeight = puWeight * lumiWeight;
	///////////////////////////////////////////////////*/
	
	cutmap["Processed"] += 1;

	int numPV = *NPV;
	vector< myJet > JETS;
	vector< float > tmpTriggerMass;
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
		double sysJEC = 0;
		double sysJER = 1;
		if ( !isData ) {
			if ( systematics.Contains("JESUp") ){
				double JESUp = uncertainty( rawJet, jetCorrUnc, true );
				sysJEC = ( + JESUp );
			} else if  ( systematics.Contains("JESDown") ){
				double JESDown = uncertainty( rawJet, jetCorrUnc, false );
				sysJEC = ( - JESDown );
			}

			// Based on this: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCMSDataAnalysisSchool2015KoreaJetExercise#Resolution
			if ( systematics.Contains("JERUp") ){
				double JERUp = getJER( (*jetEta)[i], 1 );
				double corrJetPt = (*jetPt)[i] * JEC;
				double deltaPtUp = ( corrJetPt - (*jetGenPt)[i] ) * (JERUp-1.0);
				sysJER = max( 0.0, ( corrJetPt + deltaPtUp ) / corrJetPt );
			} else if  ( systematics.Contains("JERDown") ){
				double JERDown = getJER( (*jetEta)[i], -1 );
				double corrJetPt = (*jetPt)[i] * JEC;
				double deltaPtDown = ( corrJetPt - (*jetGenPt)[i] ) * (JERDown-1.0);
				sysJER = max( 0.0, ( corrJetPt + deltaPtDown ) / corrJetPt );
			}
		} 
		corrJet = rawJet* ( ( JEC * sysJER ) + sysJEC  );

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
	

	if ( ORTriggers ) {
		
		cutmap["Trigger"] += 1;
		histos1D_[ "HT_cutTrigger" ]->Fill( HT , totalWeight );
		histos1D_[ "jetNum_cutTrigger" ]->Fill( numJets, totalWeight );
		for (int i = 0; i < numJets; i++) {
			histos1D_[ "jetPt_cutTrigger" ]->Fill( JETS[i].p4.Pt(), totalWeight );
			histos1D_[ "jetEta_cutTrigger" ]->Fill( JETS[i].p4.Eta(), totalWeight );
		}
		event		= *ievent;
		run		= *Run;
		lumi		= *Lumi;

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

			if( ( numJets == 4 ) && ( HT > cutAK4HT ) && ( JETS[3].p4.Pt() > cutAK4jetPt ) ){
				
				myJet j1, j2, j3, j4;
				int tmpMin = - 999;
				bool maxdR1234 = 0, maxdR1324 = 0, maxdR1423 = 0;

				if ( massPairing ) {

					TLorentzVector j12, j34, j13, j24, j14, j23;
					vector<double> tmpmassPair;

					j12 = JETS[0].p4 + JETS[1].p4;
					j34 = JETS[2].p4 + JETS[3].p4;
					double massPair1234 = TMath::Abs( j12.M() - j34.M() )/ (j12.M() + j34.M());
					tmpmassPair.push_back( massPair1234 );

					j13 = JETS[0].p4 + JETS[2].p4;
					j24 = JETS[1].p4 + JETS[3].p4;
					double massPair1324 = TMath::Abs( j13.M() - j24.M() )/ (j13.M() + j24.M());
					tmpmassPair.push_back( massPair1324 );

					j14 = JETS[0].p4 + JETS[3].p4;
					j23 = JETS[1].p4 + JETS[2].p4;
					double massPair1423 = TMath::Abs( j14.M() - j23.M() )/ (j14.M() + j23.M());
					tmpmassPair.push_back( massPair1423 );
					
					tmpMin = min_element(tmpmassPair.begin(), tmpmassPair.end()) - tmpmassPair.begin();

				} else {
				
					vector<double> tmpDijetR;
					double dR12 = JETS[0].p4.DeltaR( JETS[1].p4 );
					double dR34 = JETS[2].p4.DeltaR( JETS[3].p4 );
					maxdR1234 = (dR12 > dR34);
					double deltaR1234 = abs( dR12 - 0.8 )  + abs( dR34 - 0.8 );
					tmpDijetR.push_back( deltaR1234 );

					double dR13 = JETS[0].p4.DeltaR( JETS[2].p4 );
					double dR24 = JETS[1].p4.DeltaR( JETS[3].p4 );
					maxdR1324 = (dR13 > dR24);
					double deltaR1324 = abs( dR13 - 0.8 )  + abs( dR24 - 0.8 );
					tmpDijetR.push_back( deltaR1324 );

					double dR14 = JETS[0].p4.DeltaR( JETS[3].p4 );
					double dR23 = JETS[1].p4.DeltaR( JETS[2].p4 );
					maxdR1423 = (dR14 > dR23);
					double deltaR1423 = abs( dR14 - 0.8 )  + abs( dR23 - 0.8 );
					tmpDijetR.push_back( deltaR1423 );

					tmpMin = min_element(tmpDijetR.begin(), tmpDijetR.end()) - tmpDijetR.begin();
				}

				if( tmpMin == 0 ){
					if ( maxdR1234 ){
						j1 = JETS[0];
						j2 = JETS[1];
						j3 = JETS[2];
						j4 = JETS[3];
					} else{
						j1 = JETS[2];
						j2 = JETS[3];
						j3 = JETS[0];
						j4 = JETS[1];
					}
				} else if ( tmpMin == 1 ) {
					if ( maxdR1324 ){
						j1 = JETS[0];
						j2 = JETS[2];
						j3 = JETS[1];
						j4 = JETS[3];
					} else {
						j1 = JETS[1];
						j2 = JETS[3];
						j3 = JETS[0];
						j4 = JETS[2];
					}
				} else if ( tmpMin == 2 ) {
					if ( maxdR1423 ){
						j1 = JETS[0];
						j2 = JETS[3];
						j3 = JETS[1];
						j4 = JETS[2];
					} else {
						j1 = JETS[1];
						j2 = JETS[2];
						j3 = JETS[0];
						j4 = JETS[3];
					}
				}

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

				/// btag scale factors
				double j1btagSF = 1, j2btagSF = 1, j3btagSF = 1, j4btagSF = 1;
				if ( !isData ) {
					string sysType;
					if ( systematics.Contains("BtagUp") ) sysType = "up";
					else if ( systematics.Contains("BtagDown") ) sysType = "down";
					else sysType = "central";
					string measurementType = "comb";

					j1btagSF = btagSF ( btagCSVFile, j1.p4.Pt(), j1.p4.Eta(), j1.hadronFlavour, sysType, measurementType );
					//LogWarning("btag") << j1btagSF << " " << j1.hadronFlavour;
					j2btagSF = btagSF ( btagCSVFile, j2.p4.Pt(), j2.p4.Eta(), j2.hadronFlavour, sysType, measurementType );
					j3btagSF = btagSF ( btagCSVFile, j3.p4.Pt(), j3.p4.Eta(), j3.hadronFlavour, sysType, measurementType );
					j4btagSF = btagSF ( btagCSVFile, j4.p4.Pt(), j4.p4.Eta(), j4.hadronFlavour, sysType, measurementType );
				}
				///////////////////////////////////////////
				
				if ( mkTree ) {
					jetsPt->push_back( j1.p4.Pt() );
					jetsPt->push_back( j2.p4.Pt() );
					jetsPt->push_back( j3.p4.Pt() );
					jetsPt->push_back( j4.p4.Pt() );
					jetsEta->push_back( j1.p4.Eta() );
					jetsEta->push_back( j2.p4.Eta() );
					jetsEta->push_back( j3.p4.Eta() );
					jetsEta->push_back( j4.p4.Eta() );
					jetsPhi->push_back( j1.p4.Phi() );
					jetsPhi->push_back( j2.p4.Phi() );
					jetsPhi->push_back( j3.p4.Phi() );
					jetsPhi->push_back( j4.p4.Phi() );
					jetsE->push_back( j1.p4.E() );
					jetsE->push_back( j2.p4.E() );
					jetsE->push_back( j3.p4.E() );
					jetsE->push_back( j4.p4.E() );
					jetsQGL->push_back( j1.qgl );
					jetsQGL->push_back( j2.qgl );
					jetsQGL->push_back( j3.qgl );
					jetsQGL->push_back( j4.qgl );
					jetsCSVv2->push_back( j1.btagCSVv2 );
					jetsCSVv2->push_back( j2.btagCSVv2 );
					jetsCSVv2->push_back( j3.btagCSVv2 );
					jetsCSVv2->push_back( j4.btagCSVv2 );
					jetsCSVv2SF->push_back( j1btagSF );
					jetsCSVv2SF->push_back( j2btagSF );
					jetsCSVv2SF->push_back( j3btagSF );
					jetsCSVv2SF->push_back( j4btagSF );
					jetsCMVAv2->push_back( j1.btagCMVAv2 );
					jetsCMVAv2->push_back( j2.btagCMVAv2 );
					jetsCMVAv2->push_back( j3.btagCMVAv2 );
					jetsCMVAv2->push_back( j4.btagCMVAv2 );
					jet1Pt = JETS[0].p4.Pt();
					jet2Pt = JETS[1].p4.Pt();
					jet3Pt = JETS[2].p4.Pt();
					jet4Pt = JETS[3].p4.Pt();

					RUNAtree->Fill();
				}

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

				for (int i = 0; i < numJets; i++) {
					histos1D_[ "neutralHadronEnergyFrac_cutBestPair" ]->Fill( JETS[i].nhf, totalWeight );
					histos1D_[ "neutralEmEnergyFrac_cutBestPair" ]->Fill( JETS[i].nEMf, totalWeight );
					histos1D_[ "chargedHadronEnergyFrac_cutBestPair" ]->Fill( JETS[i].chf, totalWeight );
					histos1D_[ "chargedEmEnergyFrac_cutBestPair" ]->Fill( JETS[i].cEMf, totalWeight );
					histos1D_[ "numConst_cutBestPair" ]->Fill( JETS[i].numConst, totalWeight );
					histos1D_[ "chargedMultiplicity_cutBestPair" ]->Fill( JETS[i].chm, totalWeight );
				}

				histos1D_[ "massAve_cutBestPair" ]->Fill( massAve, totalWeight );
				histos1D_[ "massAsym_cutBestPair" ]->Fill( massAsym, totalWeight );
				histos1D_[ "deltaEta_cutBestPair" ]->Fill( deltaEta, totalWeight );
				histos1D_[ "deltaR_cutBestPair" ]->Fill( deltaR , totalWeight );
				histos1D_[ "cosThetaStar1_cutBestPair" ]->Fill( cosThetaStar1 , totalWeight );
				histos1D_[ "cosThetaStar2_cutBestPair" ]->Fill( cosThetaStar2 , totalWeight );
				histos2D_[ "deltavsMassAve_cutBestPair" ]->Fill( massAve, delta1 , totalWeight );
				histos2D_[ "deltavsMassAve_cutBestPair" ]->Fill( massAve, delta2 , totalWeight );
				histos2D_[ "dijetsEta_cutBestPair" ]->Fill( eta1, eta2, totalWeight );

				if ( massPairing ) {
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
				} else {
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

			/*if ( numJets > 4 ){
				for (unsigned int ijet = 4; ijet < JETS.size(); ijet++) {
					jetsPt->push_back( JETS[ ijet ].p4.Pt() );
					jetsEta->push_back( JETS[ ijet ].p4.Eta() );
					jetsPhi->push_back( JETS[ ijet ].p4.Phi() );
					jetsE->push_back( JETS[ ijet ].p4.E() );
					jetsQGL->push_back( JETS[ ijet ].qgl );
				}
			}*/
		}
	}
	JETS.clear();
}


// ------------ method called once each job just before starting event loop  ------------
void RUNResolvedAnalysis::beginJob() {

	// Calculate PUWeight
	if ( !isData ) PUWeight_.generateWeights( dataPUFile );

	if ( mkTree ) {
		RUNAtree = fs_->make< TTree >("RUNATree", "RUNATree"); 
		RUNAtree->Branch( "run", &run, "run/I" );
		RUNAtree->Branch( "lumi", &lumi, "lumi/I" );
		RUNAtree->Branch( "event", &event, "event/I" );
		RUNAtree->Branch( "numJets", &numJets, "numJets/I" );
		RUNAtree->Branch( "numPV", &numPV, "numPV/I" );
		RUNAtree->Branch( "puWeight", &puWeight, "puWeight/F" );
		RUNAtree->Branch( "lumiWeight", &lumiWeight, "lumiWeight/F" );
		RUNAtree->Branch( "pdfWeight", &pdfWeight, "pdfWeight/F" );
		RUNAtree->Branch( "genWeight", &genWeight, "genWeight/F" );
		RUNAtree->Branch( "HT", &HT, "HT/F" );
		RUNAtree->Branch( "MET", &MET, "MET/F" );
		RUNAtree->Branch( "mass1", &mass1, "mass1/F" );
		RUNAtree->Branch( "mass2", &mass2, "mass2/F" );
		RUNAtree->Branch( "massAve", &massAve, "massAve/F" );
		RUNAtree->Branch( "delta1", &delta1, "delta1/F" );
		RUNAtree->Branch( "delta2", &delta2, "delta2/F" );
		RUNAtree->Branch( "massAsym", &massAsym, "massAsym/F" );
		RUNAtree->Branch( "eta1", &eta1, "eta1/F" );
		RUNAtree->Branch( "eta2", &eta2, "eta2/F" );
		RUNAtree->Branch( "deltaEta", &deltaEta, "deltaEta/F" );
		RUNAtree->Branch( "deltaR", &deltaR, "deltaR/F" );
		RUNAtree->Branch( "cosThetaStar1", &cosThetaStar1, "cosThetaStar1/F" );
		RUNAtree->Branch( "cosThetaStar2", &cosThetaStar2, "cosThetaStar2/F" );
		RUNAtree->Branch( "jetsPt", "vector<float>", &jetsPt);
		RUNAtree->Branch( "jetsEta", "vector<float>", &jetsEta);
		RUNAtree->Branch( "jetsPhi", "vector<float>", &jetsPhi);
		RUNAtree->Branch( "jetsE", "vector<float>", &jetsE);
		RUNAtree->Branch( "jetsQGL", "vector<float>", &jetsQGL);
		RUNAtree->Branch( "jetsCSVv2", "vector<float>", &jetsCSVv2);
		RUNAtree->Branch( "jetsCSVv2SF", "vector<float>", &jetsCSVv2SF);
		RUNAtree->Branch( "jetsCMVAv2", "vector<float>", &jetsCMVAv2);
		RUNAtree->Branch( "jet1Pt", &jet1Pt, "jet1Pt/F" );
		RUNAtree->Branch( "jet2Pt", &jet2Pt, "jet2Pt/F" );
		RUNAtree->Branch( "jet3Pt", &jet3Pt, "jet3Pt/F" );
		RUNAtree->Branch( "jet4Pt", &jet4Pt, "jet4Pt/F" );
		RUNAtree->Branch( "muonsPt", "vector<float>", &muonsPt);
		RUNAtree->Branch( "numMuon", &numMuon, "numMuon/I" );
		RUNAtree->Branch( "elesPt", "vector<float>", &elesPt);
		RUNAtree->Branch( "numEle", &numEle, "numEle/I" );
		//RUNAtree->Branch( "scaleWeights", &scaleWeights );
		//RUNAtree->Branch( "pdfWeights", &pdfWeights );
		//RUNAtree->Branch( "alphaWeights", &alphaWeights );
	}

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


	histos1D_[ "jetPt_cutTrigger" ] = fs_->make< TH1D >( "jetPt_cutTrigger", "jetPt_cutTrigger", 100, 0., 1000. );
	histos1D_[ "jetPt_cutTrigger" ]->Sumw2();
	histos1D_[ "jetEta_cutTrigger" ] = fs_->make< TH1D >( "jetEta_cutTrigger", "jetEta_cutTrigger", 100, -5., 5. );
	histos1D_[ "jetEta_cutTrigger" ]->Sumw2();
	histos1D_[ "jetNum_cutTrigger" ] = fs_->make< TH1D >( "jetNum_cutTrigger", "jetNum_cutTrigger", 10, 0., 10. );
	histos1D_[ "jetNum_cutTrigger" ]->Sumw2();
	histos1D_[ "HT_cutTrigger" ] = fs_->make< TH1D >( "HT_cutTrigger", "HT_cutTrigger", 300, 0., 3000. );
	histos1D_[ "HT_cutTrigger" ]->Sumw2();


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
	cutLabels.push_back("Trigger");
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
void RUNResolvedAnalysis::endJob() {

	int ibin = 1;
	for( const string &ivec : cutLabels ) {
		histos1D_["hcutflow"]->SetBinContent( ibin, cutmap[ ivec ] );
		histos1D_["hcutflow"]->GetXaxis()->SetBinLabel( ibin, ivec.c_str() );
		ibin++;
	}

}

void RUNResolvedAnalysis::clearVariables() {

	jetsPt->clear();
	jetsEta->clear();
	jetsPhi->clear();
	jetsE->clear();
	jetsQGL->clear();
	jetsCSVv2->clear();
	jetsCSVv2SF->clear();
	jetsCMVAv2->clear();

}

void RUNResolvedAnalysis::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<double>("cutAK4jetPt", 80);
	desc.add<double>("cutAK4HT", 800);
	desc.add<double>("cutAK4MassAsym", 0.20);
	desc.add<double>("cutDelta", 180);
	desc.add<double>("cutDeltaEtaDijetSyst", .75);
	desc.add<bool>("isData", false);
	desc.add<bool>("LHEcont", false);
	desc.add<bool>("mkTree", false);
	desc.add<bool>("massPairing", false);
	desc.add<double>("scale", 1);
	desc.add<string>("dataPUFile", "supportFiles/PileupData2015D_JSON_10-23-2015.root");
	desc.add<string>("jecVersion", "supportFiles/Summer15_25nsV6");
	desc.add<string>("btagCSVFile", "supportFiles/CSVv2_ichep.csv");
	desc.add<string>("systematics", "None");
	vector<string> HLTPass;
	HLTPass.push_back("HLT_PFHT800");
	desc.add<vector<string>>("triggerPass",	HLTPass);

	desc.add<InputTag>("Lumi", 	InputTag("eventInfo:evtInfoLumiBlock"));
	desc.add<InputTag>("Run", 	InputTag("eventInfo:evtInfoRunNumber"));
	desc.add<InputTag>("Event", 	InputTag("eventInfo:evtInfoEventNumber"));
	desc.add<InputTag>("generator", 	InputTag("generator"));
	desc.add<InputTag>("extLHEProducer", 	InputTag("externalLHEProducer"));
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
	// Trigger
	desc.add<InputTag>("triggerPrescale",		InputTag("TriggerUserData:triggerPrescaleTree"));
	desc.add<InputTag>("triggerBit",		InputTag("TriggerUserData:triggerBitTree"));
	desc.add<InputTag>("triggerName",		InputTag("TriggerUserData:triggerNameTree"));
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
      
void RUNResolvedAnalysis::beginRun(const Run& iRun, const EventSetup& iSetup){

	/* Weights from scale variations, PDFs etc. are stored in the relative product. 
	 * Notice that to be used they need to be renormalized to the central event weight
	 * at LHE level which may be different from genEvtInfo->weight()
	 */
	if (!isData && !LHEcont) {
		Handle<LHERunInfoProduct> lheRunInfo;
		iRun.getByLabel( "externalLHEProducer", lheRunInfo );

		if (lheRunInfo.isValid()) {
			// Check which PDF set was used
			lhaPdfId = lheRunInfo->heprup().PDFSUP.first;
			std::cout<<"LHE: LHA PDF ID = "<<lhaPdfId<< "\nLHE:   --> For more info about the sets, check: https://lhapdf.hepforge.org/pdfsets.html"<<std::endl;

			// Check headers
			LHERunInfoProduct lheRunInfoProduct = *(lheRunInfo.product());
			typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
			size_t iHead = 0;
			for (headers_const_iterator header=lheRunInfoProduct.headers_begin(); header!=lheRunInfoProduct.headers_end(); header++){
				if (header->tag()=="initrwgt") {
					//std::cout<<"LHE: "<<iHead<<" "<<header->tag()<<std::endl;
					for (auto line : header->lines()) {
						//std::cout<<"LHE: "<<line;
						// Fix buggy powheg samples
						if (lhaPdfId==-1 && line.find("weight id=\"2001\"")!=std::string::npos) {
						if (line.find("PDF set = 260001")!=std::string::npos) lhaPdfId = 260000;
						else if (line.find("PDF set = 260401")!=std::string::npos) lhaPdfId = 260400;
						}
					}
				}
				iHead++;
			}
		}
	}

	if (isData) {
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
		
}
      
void RUNResolvedAnalysis::endRun(const Run& iRun, const EventSetup& iSetup){
	triggerNamesList.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(RUNResolvedAnalysis);
