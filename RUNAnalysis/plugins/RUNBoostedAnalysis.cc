// -*- C++ -*-
//
// Package:    RUNA/RUNAnalysis
// Class:      RUNBoostedAnalysis
// Original Author:  Alejandro Gomez Espinosa
//         Created:  Tue, 14 Oct 2014 23:13:13 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
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
class RUNBoostedAnalysis : public EDAnalyzer {

	public:
		explicit RUNBoostedAnalysis(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~RUNBoostedAnalysis();

	private:
		virtual void beginJob() override;
		virtual void analyze(const Event&, const EventSetup&) override;
		virtual void endJob() override;
		virtual void beginRun(const Run&, const EventSetup&) override;
		virtual void endRun(Run const&, EventSetup const&) override;

		//virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
		//virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

		// ----------member data ---------------------------
		string PUMethod;
		PUReweighter PUWeight_;
		int lhaPdfId ;

		Service<TFileService> fs_;
		TTree *RUNAtree;
		map< string, TH1D* > histos1D_;
		map< string, TH2D* > histos2D_;
		vector< string > cutLabels;
		map< string, double > cutmap;

		bool isData;
		bool isTTbar;
		bool LHEcont;
		bool mkTree;
		bool sortInTau21;
		bool sortInMass;
		double scale;
		double cutAK8HT;
		double cutAK8jetPt;
		double cutAK8MassAsym;
		double cutTau21;
		double cutDeltaEtaDijet;
		string dataPUFile;
		string jecVersion;
		string btagCSVFile;
		string btagMVAFile;
		TString systematics;

		vector<string> triggerPass, triggerNamesList;
		vector<JetCorrectorParameters> jetPar;
		FactorizedJetCorrector * jetJECAK8;
		vector<JetCorrectorParameters> massPar;
		FactorizedJetCorrector * massJECAK8;
		JetCorrectionUncertainty *jetCorrUnc;

		vector<float> *muonsPt = new std::vector<float>();
		vector<float> *elesPt  = new std::vector<float>();
		ULong64_t event = 0;
		int numJets = 0, numPV = 0, numEle = 0, numMuon = 0;
		unsigned int lumi = 0, run=0;
		float AK4HT = 0, HT = 0, trimmedMass = -999, 
		      puWeight = -999, genWeight = -999, lumiWeight = -999, pdfWeight = -999, MET = -999,
		      jet1Pt = -999, jet1Eta = -999, jet1Phi = -999, jet1E = -999, jet1M = -999, jet1btagCSVv2 = -9999, jet1btagCMVAv2 = -9999, jet1btagDoubleB = -9999,
		      jet2Pt = -999, jet2Eta = -999, jet2Phi = -999, jet2E = -999, jet2M = -999, jet2btagCSVv2 = -9999, jet2btagCMVAv2 = -9999, jet2btagDoubleB = -9999,
		      jet1btagCSVv2SF = 1, jet2btagCSVv2SF = 1, jet1btagCMVAv2SF = 1, jet2btagCMVAv2SF = 1,
		      subjet11Pt = -999, subjet11Eta = -999, subjet11Phi = -999, subjet11E = -999, subjet11btagCSVv2 = -9999, subjet11btagCMVAv2 = -9999, 
		      subjet12Pt = -999, subjet12Eta = -999, subjet12Phi = -999, subjet12E = -999, subjet12btagCSVv2 = -9999, subjet12btagCMVAv2 = -9999, 
		      subjet21Pt = -999, subjet21Eta = -999, subjet21Phi = -999, subjet21E = -999, subjet21btagCSVv2 = -9999, subjet21btagCMVAv2 = -9999, 
		      subjet22Pt = -999, subjet22Eta = -999, subjet22Phi = -999, subjet22E = -999, subjet22btagCSVv2 = -9999, subjet22btagCMVAv2 = -9999,
		      genPartonPt1 = -999, genPartonMass1 = -999, genPartonDau11ID = -999, genPartonDau12ID = -999, 
		      genPartonPt2 = -999, genPartonMass2 = -999, genPartonDau21ID = -999, genPartonDau22ID = -999, 
		      //massAve = -9999, massAsym = -9999, 
		      jet1PrunedMass = -9999, jet2PrunedMass = -9999,
		      jet1SoftDropMass = -9999, jet2SoftDropMass = -9999,
		      //trimmedMassAve = -9999, trimmedMassAsym = -9999, 
		      prunedMassAve = -9999, prunedMassAsym = -9999, 
		      //filteredMassAve = -9999, filteredMassAsym = -9999, 
		      //softDropMassAve = -9999, softDropMassAsym = -9999, 
		      jet1CosThetaStar = -9999, //jet2CosThetaStar = -9999, 
		      deltaEtaDijet = -9999,
		      jet1Tau21 = -9999, jet1Tau31 = -9999, jet1Tau32 = -9999,
		      jet2Tau21 = -9999, jet2Tau31 = -9999, jet2Tau32 = -9999;
		      //jet1SubjetPtRatio = -999, jet2SubjetPtRatio = -999, jet1SubjetMass21Ratio = -999, jet1Subjet112MassRatio = -999, jet1Subjet1JetMassRatio = - 999, jet1Subjet212MassRatio = - 999, jet1Subjet2JetMassRatio = - 999,
		      //jet2SubjetMass21Ratio = -999, jet2Subjet112MassRatio = -999, jet2Subjet1JetMassRatio = - 999, jet2Subjet212MassRatio = - 999, jet2Subjet2JetMassRatio = - 999, 
		      //cosPhi13412 = -9999, cosPhi31234 = -9999,
		      //dalitzY1 = -9999, dalitzY2 = -9999, dalitzY3 = -9999, dalitzY4 = -9999, dalitzY5 = -9999, dalitzY6 = -9999, 
		      //dalitzX1 = -9999, dalitzX2 = -9999, dalitzX3 = -9999, dalitzX4 = -9999, dalitzX5 = -9999, dalitzX6 = -9999;
		vector<float> scaleWeights, pdfWeights, alphaWeights;


		/// Jets
		EDGetTokenT<vector<float>> jetAK4Pt_;
		EDGetTokenT<vector<float>> jetAK4Eta_;
		EDGetTokenT<vector<float>> jetAK4Phi_;
		EDGetTokenT<vector<float>> jetAK4E_;
		EDGetTokenT<vector<float>> jetPt_;
		EDGetTokenT<vector<float>> jetEta_;
		EDGetTokenT<vector<float>> jetPhi_;
		EDGetTokenT<vector<float>> jetE_;
		EDGetTokenT<vector<float>> jetTrimmedMass_;
		EDGetTokenT<vector<float>> jetPrunedMass_;
		EDGetTokenT<vector<float>> jetFilteredMass_;
		EDGetTokenT<vector<float>> jetSoftDropMass_;
		EDGetTokenT<vector<float>> jetTau1_;
		EDGetTokenT<vector<float>> jetTau2_;
		EDGetTokenT<vector<float>> jetTau3_;
		EDGetTokenT<vector<float>> jetNSubjets_;
		EDGetTokenT<vector<float>> jetSubjetIndex0_;
		EDGetTokenT<vector<float>> jetSubjetIndex1_;
		EDGetTokenT<vector<float>> jetSubjetIndex2_;
		EDGetTokenT<vector<float>> jetSubjetIndex3_;
		EDGetTokenT<vector<vector<int>>> jetKeys_;
		EDGetTokenT<vector<float>> jetCSVv2_;
		EDGetTokenT<vector<float>> jetCMVAv2_;
		EDGetTokenT<vector<float>> jetDoubleB_;
		EDGetTokenT<vector<float>> jetArea_;
		EDGetTokenT<vector<float>> jetGenPt_;
		EDGetTokenT<vector<float>> jetGenEta_;
		EDGetTokenT<vector<float>> jetGenPhi_;
		EDGetTokenT<vector<float>> jetGenE_;
		EDGetTokenT<vector<float>> jetHadronFlavour_;

		/// Event variables
		EDGetTokenT<vector<float>> metPt_;
		EDGetTokenT<int> NPV_;
		EDGetTokenT<int> trueNInt_;
		EDGetTokenT<vector<int>> bunchCross_;
		EDGetTokenT<vector<float>> rho_;
		EDGetTokenT<vector<int>> puNumInt_;
		EDGetTokenT<unsigned int> lumi_;
		EDGetTokenT<unsigned int> run_;
		EDGetTokenT<ULong64_t> event_;
		EDGetTokenT<GenEventInfoProduct> generator_;
		EDGetTokenT<GenParticleCollection> genParticles_;
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
		EDGetTokenT<vector<float>> chargedHadronMultiplicity_;
		EDGetTokenT<vector<float>> neutralMultiplicity_;
		EDGetTokenT<vector<float>> chargedMultiplicity_;
		EDGetTokenT<vector<float>> muonEnergyFrac_; 

		// Subjets
		EDGetTokenT<vector<float>> subjetPt_;
		EDGetTokenT<vector<float>> subjetEta_;
		EDGetTokenT<vector<float>> subjetPhi_;
		EDGetTokenT<vector<float>> subjetE_;
		EDGetTokenT<vector<float>> subjetMass_;
		EDGetTokenT<vector<float>> subjetCSVv2_;
		EDGetTokenT<vector<float>> subjetCMVAv2_;

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
RUNBoostedAnalysis::RUNBoostedAnalysis(const ParameterSet& iConfig):
	jetAK4Pt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetAK4Pt"))),
	jetAK4Eta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetAK4Eta"))),
	jetAK4Phi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetAK4Phi"))),
	jetAK4E_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetAK4E"))),
	jetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPt"))),
	jetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetEta"))),
	jetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPhi"))),
	jetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetE"))),
	jetTrimmedMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTrimmedMass"))),
	jetPrunedMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetPrunedMass"))),
	jetFilteredMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetFilteredMass"))),
	jetSoftDropMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSoftDropMass"))),
	jetTau1_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTau1"))),
	jetTau2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTau2"))),
	jetTau3_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetTau3"))),
	jetNSubjets_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetNSubjets"))),
	jetSubjetIndex0_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex0"))),
	jetSubjetIndex1_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex1"))),
	jetSubjetIndex2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex2"))),
	jetSubjetIndex3_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetSubjetIndex3"))),
	jetKeys_(consumes<vector<vector<int>>>(iConfig.getParameter<InputTag>("jetKeys"))),
	jetCSVv2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCSVv2"))),
	jetCMVAv2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetCMVAv2"))),
	jetDoubleB_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetDoubleB"))),
	jetArea_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetArea"))),
	jetGenPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenPt"))),
	jetGenEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenEta"))),
	jetGenPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenPhi"))),
	jetGenE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetGenE"))),
	jetHadronFlavour_(consumes<vector<float>>(iConfig.getParameter<InputTag>("jetHadronFlavour"))),
	metPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("metPt"))),
	NPV_(consumes<int>(iConfig.getParameter<InputTag>("NPV"))),
	trueNInt_(consumes<int>(iConfig.getParameter<InputTag>("trueNInt"))),
	bunchCross_(consumes<vector<int>>(iConfig.getParameter<InputTag>("bunchCross"))),
	rho_(consumes<vector<float>>(iConfig.getParameter<InputTag>("rho"))),
	puNumInt_(consumes<vector<int>>(iConfig.getParameter<InputTag>("puNumInt"))),
	lumi_(consumes<unsigned int>(iConfig.getParameter<InputTag>("Lumi"))),
	run_(consumes<unsigned int>(iConfig.getParameter<InputTag>("Run"))),
	event_(consumes<ULong64_t>(iConfig.getParameter<InputTag>("Event"))),
	generator_(consumes<GenEventInfoProduct>(iConfig.getParameter<InputTag>("generator"))),
	genParticles_(consumes<GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles"))),
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
	muonEnergyFrac_(consumes<vector<float>>(iConfig.getParameter<InputTag>("muonEnergyFrac"))),
	// Subjets
	subjetPt_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetPt"))),
	subjetEta_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetEta"))),
	subjetPhi_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetPhi"))),
	subjetE_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetE"))),
	subjetMass_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetMass"))),
	subjetCSVv2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetCSVv2"))),
	subjetCMVAv2_(consumes<vector<float>>(iConfig.getParameter<InputTag>("subjetCMVAv2"))),
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
	cutAK8HT 	= iConfig.getParameter<double>("cutAK8HT");
	cutAK8jetPt 	= iConfig.getParameter<double>("cutAK8jetPt");
	cutAK8MassAsym 	= iConfig.getParameter<double>("cutAK8MassAsym");
	cutTau21 	= iConfig.getParameter<double>("cutTau21");
	cutDeltaEtaDijet= iConfig.getParameter<double>("cutDeltaEtaDijet");
	isData 		= iConfig.getParameter<bool>("isData");
	isTTbar 	= iConfig.getParameter<bool>("isTTbar");
	LHEcont 	= iConfig.getParameter<bool>("LHEcont");
	mkTree 		= iConfig.getParameter<bool>("mkTree");
	sortInTau21	= iConfig.getParameter<bool>("sortInTau21");
	sortInMass	= iConfig.getParameter<bool>("sortInMass");
	dataPUFile 	= iConfig.getParameter<string>("dataPUFile");
	PUMethod 	= iConfig.getParameter<string>("PUMethod");
	jecVersion 	= iConfig.getParameter<string>("jecVersion");
	btagCSVFile 	= iConfig.getParameter<string>("btagCSVFile");
	btagMVAFile 	= iConfig.getParameter<string>("btagMVAFile");
	systematics 	= iConfig.getParameter<string>("systematics");
	triggerPass 	= iConfig.getParameter<vector<string>>("triggerPass");

	/////// JECs
	string prefix;
	if (isData) prefix = jecVersion + "_DATA_";
	else prefix = jecVersion + "_MC_";

	// all jet
	vector<string> jecAK8PayloadNames_;
	jecAK8PayloadNames_.push_back(prefix + "L1FastJet_AK8PF"+PUMethod+".txt");
	jecAK8PayloadNames_.push_back(prefix + "L2Relative_AK8PF"+PUMethod+".txt");
	jecAK8PayloadNames_.push_back(prefix + "L3Absolute_AK8PF"+PUMethod+".txt");
	if (isData) jecAK8PayloadNames_.push_back(prefix + "L2L3Residual_AK8PF"+PUMethod+".txt");

	for ( vector<string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		jetPar.push_back(pars);
	}
	jetJECAK8 = new FactorizedJetCorrector(jetPar);

	// jet mass
	vector<string> massjecAK8PayloadNames_;
	massjecAK8PayloadNames_.push_back(prefix + "L2Relative_AK8PF"+PUMethod+".txt");
	massjecAK8PayloadNames_.push_back(prefix + "L3Absolute_AK8PF"+PUMethod+".txt");
	if (isData) massjecAK8PayloadNames_.push_back(prefix + "L2L3Residual_AK8PF"+PUMethod+".txt");

	for ( vector<string>::const_iterator payloadBegin = massjecAK8PayloadNames_.begin(), payloadEnd = massjecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters massPars(*ipayload);
		massPar.push_back(massPars);
	}
	massJECAK8 = new FactorizedJetCorrector(massPar);

	// jec uncertainty
	JetCorrectorParameters jecUncParam( prefix + "Uncertainty_AK8PF"+PUMethod+".txt");
	jetCorrUnc  = new JetCorrectionUncertainty( jecUncParam );
	////////////////////////////////////////////////////
}


RUNBoostedAnalysis::~RUNBoostedAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void RUNBoostedAnalysis::analyze(const Event& iEvent, const EventSetup& iSetup) {


	Handle<vector<float> > jetAK4Pt;
	iEvent.getByToken(jetAK4Pt_, jetAK4Pt);

	Handle<vector<float> > jetAK4Eta;
	iEvent.getByToken(jetAK4Eta_, jetAK4Eta);

	Handle<vector<float> > jetAK4Phi;
	iEvent.getByToken(jetAK4Phi_, jetAK4Phi);

	Handle<vector<float> > jetAK4E;
	iEvent.getByToken(jetAK4E_, jetAK4E);

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

	Handle<vector<float> > jetFilteredMass;
	iEvent.getByToken(jetFilteredMass_, jetFilteredMass);

	Handle<vector<float> > jetSoftDropMass;
	iEvent.getByToken(jetSoftDropMass_, jetSoftDropMass);

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

	Handle<vector<float> > jetCSVv2;
	iEvent.getByToken(jetCSVv2_, jetCSVv2);

	Handle<vector<float> > jetCMVAv2;
	iEvent.getByToken(jetCMVAv2_, jetCMVAv2);

	Handle<vector<float> > jetDoubleB;
	iEvent.getByToken(jetDoubleB_, jetDoubleB);

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

	Handle<vector<float> > metPt;
	iEvent.getByToken(metPt_, metPt);

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

	Handle<vector<float> > muonEnergyFrac;
	iEvent.getByToken(muonEnergyFrac_, muonEnergyFrac);

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

	Handle<vector<float> > subjetCSVv2;
	iEvent.getByToken(subjetCSVv2_, subjetCSVv2);

	Handle<vector<float> > subjetCMVAv2;
	iEvent.getByToken(subjetCMVAv2_, subjetCMVAv2);

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

	Handle< GenParticleCollection > genParticles;
	iEvent.getByToken( genParticles_, genParticles );

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
	
	
	////// TTbar genInfo (all this piece of code is SPECIFIC for ttbar
	if (isTTbar) {

		for( const auto & p : *genParticles ) {  

			const reco::Candidate * mother = p.mother();
			if( p.status() == 22 ) {

				if( p.pdgId() == 6 ) { 
					genPartonPt1 = p.pt();
					genPartonMass1 = p.mass();
				}
				if( p.pdgId() == -6 ) { 
					genPartonPt2 = p.pt();
					genPartonMass2 = p.mass();
				}
				if( (mother->pdgId() == 6) && ( TMath::Abs(p.pdgId()) == 24 )){
					//LogWarning("pat") << p.pdgId();
					genPartonDau11ID = p.pdgId();
				}
				if( (mother->pdgId() == -6) && ( TMath::Abs(p.pdgId()) == 24 )) genPartonDau21ID = p.pdgId();

			}
			if( p.status() == 23 ) {
				if( (mother->pdgId() == 6) && ( TMath::Abs(p.pdgId()) == 5 )) genPartonDau12ID = p.pdgId();
				if( (mother->pdgId() == -6) && ( TMath::Abs(p.pdgId()) == 5 )) genPartonDau22ID = p.pdgId();
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

	////////// Check trigger fired
	bool ORTriggers = false;
	if ( isData ) ORTriggers = checkORListOfTriggerBits( triggerNamesList, triggerBit, triggerPrescale, triggerPass, false );
	else ORTriggers = true;
	///////////////////////////////////////////////////*/
	
	////////// PU Reweight
	if ( isData ) puWeight = 1;
	else puWeight = PUWeight_.getPUWeight( *trueNInt, *bunchCross );
	histos1D_[ "PUWeight" ]->Fill( puWeight );
	lumiWeight = scale;
	double totalWeight = puWeight * lumiWeight;
	////////////////////////////////////////////////////
	
	///////// AK4 jets to model PFHT trigger
	AK4HT = 0;
	for (size_t q = 0; q < jetAK4Pt->size(); q++) {

		if ( TMath::Abs( (*jetAK4Eta)[q] ) > 2.4 ) continue;
		if ( (*jetAK4Pt)[q] < 40.0 ) continue;
		AK4HT += (*jetAK4Pt)[q];
	}
	////////////////////////////////////////////////////

	/// Applying kinematic, trigger and jet ID
	cutmap["Processed"] += 1;
	vector< myJet > JETS;
	vector< float > tmpTriggerMass;
	bool cutHT = 0;
	int numberJets = 0;
	HT = 0;

	for (size_t i = 0; i < jetPt->size(); i++) {

		if( TMath::Abs( (*jetEta)[i] ) > 2.4 ) continue;

		string typeOfJetID = "tightLepVetoJetID";
		bool jetId = jetID( (*jetEta)[i], (*jetE)[i], (*jecFactor)[i], (*neutralHadronEnergyFrac)[i], (*neutralEmEnergyFrac)[i], (*chargedHadronEnergyFrac)[i], (*muonEnergyFrac)[i], (*chargedEmEnergyFrac)[i], (*chargedMultiplicity)[i], (*neutralMultiplicity)[i], typeOfJetID ); 

		if( (*jetPt)[i] < 50 ) continue; // just to reduce time

		TLorentzVector tmpJet, rawJet, corrJet, genJet, smearJet;
		tmpJet.SetPtEtaPhiE( (*jetPt)[i], (*jetEta)[i], (*jetPhi)[i], (*jetE)[i] );
		rawJet = tmpJet* (*jecFactor)[i] ;

		double JEC = corrections( rawJet, (*jetArea)[i], (*rho)[i] ,*NPV, jetJECAK8); 
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

		if( ( corrJet.Pt() > cutAK8jetPt ) && jetId ) { 

			HT += corrJet.Pt();
			tmpTriggerMass.push_back( (*jetTrimmedMass)[i] );
			++numberJets;

			double massJEC = corrections( rawJet, (*jetArea)[i], (*rho)[i] ,*NPV, massJECAK8); 
			double corrMass = corrJet.M() * ( ( massJEC * sysJER ) + sysJEC  );
			double corrTrimmedMass = (*jetTrimmedMass)[i] * ( ( massJEC * sysJER ) + sysJEC  );
			double corrPrunedMass = (*jetPrunedMass)[i] * ( ( massJEC * sysJER ) + sysJEC  );
			double corrSoftDropMass = (*jetSoftDropMass)[i] * ( ( massJEC * sysJER ) + sysJEC  );
			double corrFilteredMass = (*jetFilteredMass)[i] * ( ( massJEC * sysJER ) + sysJEC  );

			/// Vector of zeros
			TLorentzVector tmpSubjet0, tmpSubjet1, tmpZeros;
			double tmpSubjet0BtagCSVv2 = -999, tmpSubjet0BtagCMVAv2 = -999, tmpSubjet1BtagCSVv2 = -999, tmpSubjet1BtagCMVAv2 = -999;
			tmpZeros.SetPtEtaPhiE( 0, 0, 0, 0 );

			for (size_t j = 0; j < subjetPt->size(); j++) {
				if( j == (*jetSubjetIndex0)[i] ) {
					//LogWarning("subjets0") << j << " " << (*jetSubjetIndex0)[i] << " " <<  subjetPt->size() << " " << (*subjetPt)[j];
					tmpSubjet0.SetPtEtaPhiE( (*subjetPt)[j], (*subjetEta)[j], (*subjetPhi)[j], (*subjetE)[j] );
					tmpSubjet0BtagCSVv2 = (*subjetCSVv2)[j];
					tmpSubjet0BtagCMVAv2 = (*subjetCMVAv2)[j];
				} //else tmpSubjet0 = tmpZeros ; 
					
				if( j == (*jetSubjetIndex1)[i] ) {
					tmpSubjet1.SetPtEtaPhiE( (*subjetPt)[j], (*subjetEta)[j], (*subjetPhi)[j], (*subjetE)[j] );
					tmpSubjet1BtagCSVv2 = (*subjetCSVv2)[j];
					tmpSubjet1BtagCMVAv2 = (*subjetCMVAv2)[j];
				} //else tmpSubjet1 = tmpZeros ; 
			}

			//if ( (*jetCSVv2)[i] > 0.244 ) bTagCSVv2 = 1; 	// CSVv2L
			//if ( (*jetCSVv2)[i] > 0.679 ) bTagCSVv2 = 1; 	// CSVv2M
			//if ( (*jetCSVv2V1)[i] > 0.405 ) bTagCSVv2 = 1; 	// CSVv2V1L
			//if ( (*jetCSVv2V1)[i] > 0.783 ) bTagCSVv2 = 1; 	// CSVv2V1M
			
			double jec = 1. / ( (*jecFactor)[i] );
			histos1D_[ "oldJetPt" ]->Fill( (*jetPt)[i], totalWeight );
			histos1D_[ "jetPt" ]->Fill( corrJet.Pt(), totalWeight );
			histos1D_[ "rawJetPt" ]->Fill( rawJet.Pt(), totalWeight );
			histos1D_[ "oldJetEta" ]->Fill( (*jetEta)[i], totalWeight );
			histos1D_[ "jetEta" ]->Fill( corrJet.Eta(), totalWeight );
			histos1D_[ "rawJetEta" ]->Fill( rawJet.Eta(), totalWeight );
			histos1D_[ "oldJetMass" ]->Fill( rawJet.M(), totalWeight );
			histos1D_[ "jetMass" ]->Fill( corrMass, totalWeight );
			histos1D_[ "jetTrimmedMass" ]->Fill( corrTrimmedMass, totalWeight );
			histos1D_[ "jetPrunedMass" ]->Fill( corrPrunedMass, totalWeight );
			histos1D_[ "jetFilteredMass" ]->Fill( corrFilteredMass, totalWeight );
			histos1D_[ "jetSoftDropMass" ]->Fill( corrSoftDropMass, totalWeight );
			histos1D_[ "rawJetMass" ]->Fill( rawJet.M(), totalWeight );
			histos1D_[ "neutralHadronEnergyFrac" ]->Fill( (*neutralHadronEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "neutralEmEnergyFrac" ]->Fill( (*neutralEmEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "chargedHadronEnergyFrac" ]->Fill( (*chargedHadronEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "chargedEmEnergyFrac" ]->Fill( (*chargedEmEnergyFrac)[i] * jec, totalWeight );
			histos1D_[ "numConst" ]->Fill( (*chargedMultiplicity)[i] + (*neutralMultiplicity)[i], totalWeight );
			histos1D_[ "chargedMultiplicity" ]->Fill( (*chargedMultiplicity)[i], totalWeight );

			myJet tmpJET;
			tmpJET.p4 = corrJet;
			tmpJET.subjet0 = tmpSubjet0;
			tmpJET.subjet0BtagCSVv2 = tmpSubjet0BtagCSVv2;
			tmpJET.subjet0BtagCMVAv2 = tmpSubjet0BtagCMVAv2;
			tmpJET.subjet1 = tmpSubjet1;
			tmpJET.subjet1BtagCSVv2 = tmpSubjet1BtagCSVv2;
			tmpJET.subjet1BtagCMVAv2 = tmpSubjet1BtagCMVAv2;
			tmpJET.mass = corrMass;
			tmpJET.trimmedMass = corrTrimmedMass;
			tmpJET.prunedMass = corrPrunedMass;
			tmpJET.filteredMass = corrFilteredMass;
			tmpJET.softDropMass = corrSoftDropMass;
			tmpJET.tau1 = (*jetTau1)[i];
			tmpJET.tau2 = (*jetTau2)[i];
			tmpJET.tau3 = (*jetTau3)[i];
			tmpJET.btagCSVv2 = (*jetCSVv2)[i];
			tmpJET.btagCMVAv2 = (*jetCMVAv2)[i];
			tmpJET.btagDoubleB = (*jetDoubleB)[i];
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

	numPV = *NPV;
	MET = (*metPt)[0];
	numJets = numberJets;
	if ( sortInTau21 ) sort(JETS.begin(), JETS.end(), [](const myJet &p1, const myJet &p2) { float jet1tau21 = -999, jet2tau21 = -999; jet1tau21 = (p1.tau2/p1.tau1); jet2tau21 = (p2.tau2/p2.tau1); return jet1tau21 < jet2tau21; }); 
	if ( sortInMass ) sort(JETS.begin(), JETS.end(), [](const myJet &p1, const myJet &p2) { return p1.prunedMass > p2.prunedMass; }); 
	//sort(JETS.begin(), JETS.end(), [](const JETtype &p1, const JETtype &p2) { TLorentzVector tmpP1, tmpP2; tmpP1 = p1.p4; tmpP2 = p2.p4;  return tmpP1.M() > tmpP2.M(); }); 
	histos1D_[ "jetNum" ]->Fill( numJets, totalWeight );
	histos1D_[ "NPV_NOPUWeight" ]->Fill( numPV );
	histos1D_[ "NPV" ]->Fill( numPV, totalWeight );
	if ( HT > 0 ) histos1D_[ "HT" ]->Fill( HT, totalWeight );
	if ( HT > cutAK8HT ) cutHT = 1;

	sort(tmpTriggerMass.begin(), tmpTriggerMass.end(), [](const float p1, const float p2) { return p1 > p2; }); 
	if ( ( tmpTriggerMass.size()> 0 ) ) { //&& ( tmpTriggerMass[0] > cutTrimmedMassvalue) ){
		trimmedMass = tmpTriggerMass[0];
		if (HT > 0) histos2D_[ "jetTrimmedMassHT" ]->Fill( tmpTriggerMass[0], HT, totalWeight );
	}
	tmpTriggerMass.clear();

	if ( ORTriggers ) {
		//LogWarning("fired") << HT << " " << trimmedMass;
		cutmap["Trigger"] += 1; 

		histos1D_[ "HT_cutTrigger" ]->Fill( HT, totalWeight );
		histos1D_[ "MET_cutTrigger" ]->Fill( MET, totalWeight );
		histos1D_[ "METHT_cutTrigger" ]->Fill( MET/HT, totalWeight );
		histos1D_[ "NPV_cutTrigger" ]->Fill( numPV, totalWeight );
		histos1D_[ "jetNum_cutTrigger" ]->Fill( numJets, totalWeight );
		histos2D_[ "jetTrimmedMassHT_cutTrigger" ]->Fill( trimmedMass, HT, totalWeight );

		int kdum = 0;
		for (size_t k = 0; k < JETS.size(); k++) {
			histos1D_[ "jetPt_cutTrigger" ]->Fill( JETS[k].p4.Pt(), totalWeight );
			histos1D_[ "jetEta_cutTrigger" ]->Fill( JETS[k].p4.Eta(), totalWeight );
			histos1D_[ "jetMass_cutTrigger" ]->Fill( JETS[k].mass, totalWeight );

			if ( (++kdum) == 1 ) {
				histos1D_[ "jet1Pt_cutTrigger" ]->Fill( JETS[0].p4.Pt(), totalWeight );
				histos1D_[ "jet1Mass_cutTrigger" ]->Fill( JETS[0].mass, totalWeight );
				histos1D_[ "jet1Eta_cutTrigger" ]->Fill( JETS[0].p4.Eta(), totalWeight );
			}
		 }
		//////////////////////////////////////////////////////////////////////////////
		
		vector<double> dalitz1, dalitz2;
		vector<TLorentzVector> jet1SubjetsTLV, jet2SubjetsTLV;

		if ( JETS.size() > 1 ) {

			cutmap["Dijet"] += 1;
			
			histos1D_[ "HT_cutDijet" ]->Fill( HT, totalWeight );
			histos1D_[ "MET_cutDijet" ]->Fill( MET, totalWeight );
			histos1D_[ "METHT_cutDijet" ]->Fill( MET/HT, totalWeight );
			histos1D_[ "NPV_cutDijet" ]->Fill( numPV, totalWeight );
			histos1D_[ "jetNum_cutDijet" ]->Fill( numJets, totalWeight );
			histos1D_[ "jet1Pt_cutDijet" ]->Fill( JETS[0].p4.Pt(), totalWeight );
			histos1D_[ "jet1Eta_cutDijet" ]->Fill( JETS[0].p4.Eta(), totalWeight );
			histos1D_[ "jet1Mass_cutDijet" ]->Fill( JETS[0].mass, totalWeight );
			histos1D_[ "jet2Pt_cutDijet" ]->Fill( JETS[1].p4.Pt(), totalWeight );
			histos1D_[ "jet2Eta_cutDijet" ]->Fill( JETS[1].p4.Eta(), totalWeight );
			histos1D_[ "jet2Mass_cutDijet" ]->Fill( JETS[1].mass, totalWeight );
			for (size_t k = 0; k < JETS.size(); k++) {
				histos1D_[ "neutralHadronEnergyFrac_cutDijet" ]->Fill( JETS[k].nhf, totalWeight );
				histos1D_[ "neutralEmEnergyFrac_cutDijet" ]->Fill( JETS[k].nEMf, totalWeight );
				histos1D_[ "chargedHadronEnergyFrac_cutDijet" ]->Fill( JETS[k].chf, totalWeight );
				histos1D_[ "chargedEmEnergyFrac_cutDijet" ]->Fill( JETS[k].cEMf, totalWeight );
				histos1D_[ "numConst_cutDijet" ]->Fill( JETS[k].numConst, totalWeight );
				histos1D_[ "chargedMultiplicity_cutDijet" ]->Fill( JETS[k].chm, totalWeight );
			}


			// Mass average and asymmetry
			//massAve = massAverage( JETS[0].mass, JETS[1].mass );
			//massAsym = massAsymmetry( JETS[0].mass, JETS[1].mass );
			//trimmedMassAve = massAverage( JETS[0].trimmedMass, JETS[1].trimmedMass );
			//trimmedMassAsym = massAsymmetry( JETS[0].trimmedMass, JETS[1].trimmedMass );
			prunedMassAve = massAverage( JETS[0].prunedMass, JETS[1].prunedMass );
			prunedMassAsym = massAsymmetry( JETS[0].prunedMass, JETS[1].prunedMass );
			//filteredMassAve = massAverage( JETS[0].filteredMass, JETS[1].filteredMass );
			//filteredMassAsym = massAsymmetry( JETS[0].filteredMass, JETS[1].filteredMass );
			//softDropMassAve = massAverage( JETS[0].softDropMass, JETS[1].softDropMass );
			//softDropMassAsym = massAsymmetry( JETS[0].softDropMass, JETS[1].softDropMass );
			//////////////////////////////////////////////////////////////////////////
			
			// Btag
			jet1btagCSVv2 = JETS[0].btagCSVv2;
			jet2btagCSVv2 = JETS[1].btagCSVv2;
			jet1btagCMVAv2 = JETS[0].btagCMVAv2;
			jet2btagCMVAv2 = JETS[1].btagCMVAv2;
			jet1btagDoubleB = JETS[0].btagDoubleB;
			jet2btagDoubleB = JETS[1].btagDoubleB;

			/// btag scale factors
			if ( !isData ) {
				string sysType;
				if ( systematics.Contains("BtagUp") ) sysType = "up";
				else if ( systematics.Contains("BtagDown") ) sysType = "down";
				else sysType = "central";
				string measurementTypeCSV = "comb";
				string measurementTypeCMVA = "ttbar";

				jet1btagCSVv2SF = btagSF(btagCSVFile, JETS[0].p4.Pt(), JETS[0].p4.Eta(), JETS[0].hadronFlavour, sysType, measurementTypeCSV );
				jet2btagCSVv2SF = btagSF(btagCSVFile, JETS[1].p4.Pt(), JETS[1].p4.Eta(), JETS[1].hadronFlavour, sysType, measurementTypeCSV );
				jet1btagCMVAv2SF = btagSF(btagMVAFile, JETS[0].p4.Pt(), JETS[0].p4.Eta(), JETS[0].hadronFlavour, sysType, measurementTypeCMVA );
				jet2btagCMVAv2SF = btagSF(btagMVAFile, JETS[1].p4.Pt(), JETS[1].p4.Eta(), JETS[1].hadronFlavour, sysType, measurementTypeCMVA );
			}
			///////////////////////////////////////////

			// Dijet eta
			deltaEtaDijet = deltaValue( JETS[0].p4.Eta(), JETS[1].p4.Eta() );

			// Cos theta star
			jet1CosThetaStar = calculateCosThetaStar( JETS[0].p4, JETS[1].p4 ) ;
			//jet2CosThetaStar = calculateCosThetaStar( JETS[1].p4, JETS[0].p4 ) ;

			// Nsubjetiness
			jet1Tau21 = JETS[0].tau2 / JETS[0].tau1;
			jet1Tau31 = JETS[0].tau3 / JETS[0].tau1;
			jet1Tau32 = JETS[0].tau3 / JETS[0].tau2;
			jet2Tau21 = JETS[1].tau2 / JETS[1].tau1;
			jet2Tau31 = JETS[1].tau3 / JETS[1].tau1;
			jet2Tau32 = JETS[1].tau3 / JETS[1].tau2;
			////////////////////////////////////////////////////////////////////////////////



			// Cut Pt
			//if (( JETS[0].p4.Pt() > 500. ) && ( JETS[1].p4.Pt() > 450. ) ) cutJetPt = 1 ;
			//if ( cutHT && cutJetPt ) {
			if ( cutHT ) {

				cutmap["HT"] += 1;

				if ( mkTree ) {
					///// Variables for tree
					event		= *ievent;
					run		= *Run;
					lumi		= *Lumi;
					jet1Pt 		= JETS[0].p4.Pt();
					jet1Eta 	= JETS[0].p4.Eta();
					jet1Phi 	= JETS[0].p4.Phi();
					jet1E 		= JETS[0].p4.E();
					jet1M 		= JETS[0].p4.M();
					jet1PrunedMass	= JETS[0].prunedMass;
					jet1SoftDropMass	= JETS[0].softDropMass;
					jet2Pt 		= JETS[1].p4.Pt();
					jet2Eta 	= JETS[1].p4.Eta();
					jet2Phi 	= JETS[1].p4.Phi();
					jet2E 		= JETS[1].p4.E();
					jet2M 		= JETS[1].p4.M();
					jet2PrunedMass	= JETS[1].prunedMass;
					jet2SoftDropMass	= JETS[1].softDropMass;
					subjet11Pt	= JETS[0].subjet0.Pt();
					subjet11Eta	= JETS[0].subjet0.Eta();
					subjet11Phi	= JETS[0].subjet0.Phi();
					subjet11E	= JETS[0].subjet0.E();
					subjet11btagCSVv2 = JETS[0].subjet0BtagCSVv2;
					subjet11btagCMVAv2 = JETS[0].subjet0BtagCMVAv2;
					subjet12Pt	= JETS[0].subjet1.Pt();
					subjet12Eta	= JETS[0].subjet1.Eta();
					subjet12Phi	= JETS[0].subjet1.Phi();
					subjet12E	= JETS[0].subjet1.E();
					subjet12btagCSVv2 = JETS[0].subjet1BtagCSVv2;
					subjet12btagCMVAv2 = JETS[0].subjet1BtagCMVAv2;
					subjet21Pt	= JETS[1].subjet0.Pt();
					subjet21Eta	= JETS[1].subjet0.Eta();
					subjet21Phi	= JETS[1].subjet0.Phi();
					subjet21E	= JETS[1].subjet0.E();
					subjet21btagCSVv2 = JETS[1].subjet0BtagCSVv2;
					subjet21btagCMVAv2 = JETS[1].subjet0BtagCMVAv2;
					subjet22Pt	= JETS[1].subjet1.Pt();
					subjet22Eta	= JETS[1].subjet1.Eta();
					subjet22Phi	= JETS[1].subjet1.Phi();
					subjet22E	= JETS[1].subjet1.E();
					subjet22btagCSVv2 = JETS[1].subjet1BtagCSVv2;
					subjet22btagCMVAv2 = JETS[1].subjet1BtagCMVAv2;
					RUNAtree->Fill();  
				}

				histos1D_[ "HT_cutEffTrigger" ]->Fill( HT, totalWeight );
				histos1D_[ "MET_cutEffTrigger" ]->Fill( MET, totalWeight );
				histos1D_[ "METHT_cutEffTrigger" ]->Fill( MET/HT, totalWeight );
				histos1D_[ "NPV_cutEffTrigger" ]->Fill( numPV, totalWeight );
				histos1D_[ "jetNum_cutEffTrigger" ]->Fill( numJets, totalWeight );

				histos1D_[ "jet1NeutralHadronEnergyFrac_cutEffTrigger" ]->Fill( JETS[0].nhf, totalWeight );
				histos1D_[ "jet1NeutralEmEnergyFrac_cutEffTrigger" ]->Fill( JETS[0].nEMf, totalWeight );
				histos1D_[ "jet1ChargedHadronEnergyFrac_cutEffTrigger" ]->Fill( JETS[0].chf, totalWeight );
				histos1D_[ "jet1ChargedEmEnergyFrac_cutEffTrigger" ]->Fill( JETS[0].cEMf, totalWeight );
				histos1D_[ "jet1NumConst_cutEffTrigger" ]->Fill( JETS[0].numConst, totalWeight );
				histos1D_[ "jet1ChargedMultiplicity_cutEffTrigger" ]->Fill( JETS[0].chm, totalWeight );
				histos1D_[ "jet1Pt_cutEffTrigger" ]->Fill( JETS[0].p4.Pt(), totalWeight );
				histos1D_[ "jet1Eta_cutEffTrigger" ]->Fill( JETS[0].p4.Eta(), totalWeight );
				histos1D_[ "jet1Mass_cutEffTrigger" ]->Fill( JETS[0].mass, totalWeight );
				histos1D_[ "jet1PrunedMass_cutEffTrigger" ]->Fill( JETS[0].prunedMass, totalWeight );
				histos1D_[ "jet1TrimmedMass_cutEffTrigger" ]->Fill( JETS[0].trimmedMass, totalWeight );
				histos1D_[ "jet1SoftDropMass_cutEffTrigger" ]->Fill( JETS[0].softDropMass, totalWeight );
				histos1D_[ "jet1FilteredMass_cutEffTrigger" ]->Fill( JETS[0].filteredMass, totalWeight );

				histos1D_[ "jet2NeutralHadronEnergyFrac_cutEffTrigger" ]->Fill( JETS[1].nhf, totalWeight );
				histos1D_[ "jet2NeutralEmEnergyFrac_cutEffTrigger" ]->Fill( JETS[1].nEMf, totalWeight );
				histos1D_[ "jet2ChargedHadronEnergyFrac_cutEffTrigger" ]->Fill( JETS[1].chf, totalWeight );
				histos1D_[ "jet2ChargedEmEnergyFrac_cutEffTrigger" ]->Fill( JETS[1].cEMf, totalWeight );
				histos1D_[ "jet2NumConst_cutEffTrigger" ]->Fill( JETS[1].numConst, totalWeight );
				histos1D_[ "jet2ChargedMultiplicity_cutEffTrigger" ]->Fill( JETS[1].chm, totalWeight );
				histos1D_[ "jet2Pt_cutEffTrigger" ]->Fill( JETS[1].p4.Pt(), totalWeight );
				histos1D_[ "jet2Eta_cutEffTrigger" ]->Fill( JETS[1].p4.Eta(), totalWeight );
				histos1D_[ "jet2Mass_cutEffTrigger" ]->Fill( JETS[1].mass, totalWeight );
				histos1D_[ "jet2PrunedMass_cutEffTrigger" ]->Fill( JETS[1].prunedMass, totalWeight );
				histos1D_[ "jet2TrimmedMass_cutEffTrigger" ]->Fill( JETS[1].trimmedMass, totalWeight );
				histos1D_[ "jet2SoftDropMass_cutEffTrigger" ]->Fill( JETS[1].softDropMass, totalWeight );
				histos1D_[ "jet2FilteredMass_cutEffTrigger" ]->Fill( JETS[1].filteredMass, totalWeight );
				for (size_t k = 0; k < JETS.size(); k++) {
					histos1D_[ "neutralHadronEnergyFrac_cutEffTrigger" ]->Fill( JETS[k].nhf, totalWeight );
					histos1D_[ "neutralEmEnergyFrac_cutEffTrigger" ]->Fill( JETS[k].nEMf, totalWeight );
					histos1D_[ "chargedHadronEnergyFrac_cutEffTrigger" ]->Fill( JETS[k].chf, totalWeight );
					histos1D_[ "chargedEmEnergyFrac_cutEffTrigger" ]->Fill( JETS[k].cEMf, totalWeight );
					histos1D_[ "numConst_cutEffTrigger" ]->Fill( JETS[k].numConst, totalWeight );
					histos1D_[ "chargedMultiplicity_cutEffTrigger" ]->Fill( JETS[k].chm, totalWeight );
				}

				histos1D_[ "jet1Tau21_cutEffTrigger" ]->Fill( jet1Tau21, totalWeight );
				histos1D_[ "jet2Tau21_cutEffTrigger" ]->Fill( jet2Tau21, totalWeight );
				histos1D_[ "deltaEtaDijet_cutEffTrigger" ]->Fill( deltaEtaDijet, totalWeight );
				histos1D_[ "prunedMassAsym_cutEffTrigger" ]->Fill( prunedMassAsym, totalWeight );
				histos1D_[ "massAve_cutEffTrigger" ]->Fill( prunedMassAve, totalWeight );
				histos1D_[ "jet1btagCSVv2_cutEffTrigger" ]->Fill( jet1btagCSVv2, totalWeight );
				histos1D_[ "jet2btagCSVv2_cutEffTrigger" ]->Fill( jet2btagCSVv2, totalWeight );

				if ( ( jet1Tau21 < cutTau21 ) && ( jet2Tau21 < cutTau21 ) ) {

					cutmap["Tau21"] += 1;
					histos1D_[ "jet1Tau21_cutTau21" ]->Fill( jet1Tau21, totalWeight );
					histos1D_[ "jet2Tau21_cutTau21" ]->Fill( jet2Tau21, totalWeight );
					histos1D_[ "deltaEtaDijet_cutTau21" ]->Fill( deltaEtaDijet, totalWeight );
					histos1D_[ "prunedMassAsym_cutTau21" ]->Fill( prunedMassAsym, totalWeight );
					histos1D_[ "massAve_cutTau21" ]->Fill( prunedMassAve, totalWeight );

					if ( prunedMassAsym < cutAK8MassAsym ) {

						cutmap["MassAsym"] += 1;
						histos1D_[ "deltaEtaDijet_cutMassAsym" ]->Fill( deltaEtaDijet, totalWeight );
						histos1D_[ "prunedMassAsym_cutMassAsym" ]->Fill( prunedMassAsym, totalWeight );
						histos1D_[ "massAve_cutMassAsym" ]->Fill( prunedMassAve, totalWeight );

						if ( deltaEtaDijet < cutDeltaEtaDijet ) {

							cutmap["DeltaEtaDijet"] += 1;
							histos1D_[ "deltaEtaDijet_cutDeltaEtaDijet" ]->Fill( deltaEtaDijet, totalWeight );
							histos1D_[ "massAve_cutDeltaEtaDijet" ]->Fill( prunedMassAve, totalWeight );

							if ( ( jet1btagCSVv2 > 0.8 ) && ( jet2btagCSVv2 > 0.8 ) ) {
								cutmap["Btag"] += 1;
								histos1D_[ "massAve_cutBtag" ]->Fill( prunedMassAve, totalWeight );
							}
						
						}
					}
					
					///// Regular ABCD
					if ( ( prunedMassAsym < cutAK8MassAsym ) && ( deltaEtaDijet < cutDeltaEtaDijet ) ) {
						histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_A" ]->Fill( prunedMassAve, totalWeight );
						histos2D_[ "prunedMassAsymVsdeltaEtaDijet_A" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
					} else if ( ( prunedMassAsym < cutAK8MassAsym ) && ( deltaEtaDijet > cutDeltaEtaDijet ) ) {
						histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_B" ]->Fill( prunedMassAve, totalWeight );
						histos2D_[ "prunedMassAsymVsdeltaEtaDijet_B" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
					} else if ( ( prunedMassAsym > cutAK8MassAsym ) && ( deltaEtaDijet < cutDeltaEtaDijet ) ) {
						histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_C" ]->Fill( prunedMassAve, totalWeight );
						histos2D_[ "prunedMassAsymVsdeltaEtaDijet_C" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
					} else {
						histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_D" ]->Fill( prunedMassAve, totalWeight );
						histos2D_[ "prunedMassAsymVsdeltaEtaDijet_D" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
					}

					
					///// Regular ABCD
					if ( ( jet1btagCSVv2 > 0.8 ) && ( jet2btagCSVv2 > 0.8 ) ) {
						if ( ( prunedMassAsym < cutAK8MassAsym ) && ( deltaEtaDijet < cutDeltaEtaDijet ) ) {
							histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_A" ]->Fill( prunedMassAve, totalWeight );
							histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_A" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
						} else if ( ( prunedMassAsym < cutAK8MassAsym ) && ( deltaEtaDijet > cutDeltaEtaDijet ) ) {
							histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_B" ]->Fill( prunedMassAve, totalWeight );
							histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_B" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
						} else if ( ( prunedMassAsym > cutAK8MassAsym ) && ( deltaEtaDijet < cutDeltaEtaDijet ) ) {
							histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_C" ]->Fill( prunedMassAve, totalWeight );
							histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_C" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
						} else {
							histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_D" ]->Fill( prunedMassAve, totalWeight );
							histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_D" ]->Fill( prunedMassAsym, deltaEtaDijet, totalWeight );
						}
					}
				}

				if ( ( jet1Tau21 < cutTau21 ) && ( jet2Tau21 < cutTau21 ) && ( prunedMassAsym < cutAK8MassAsym ) )
					histos1D_[ "deltaEtaDijet_n-1" ]->Fill( deltaEtaDijet, totalWeight );
				if ( ( jet1Tau21 < cutTau21 ) && ( jet2Tau21 < cutTau21 ) && ( deltaEtaDijet < cutDeltaEtaDijet ) ) 
					histos1D_[ "prunedMassAsym_n-1" ]->Fill( prunedMassAsym, totalWeight );
				if ( ( prunedMassAsym < cutAK8MassAsym ) && ( deltaEtaDijet < cutDeltaEtaDijet ) ) {
					histos1D_[ "jet1Tau21_n-1" ]->Fill( jet1Tau21, totalWeight );
					histos1D_[ "jet2Tau21_n-1" ]->Fill( jet2Tau21, totalWeight );
				}
			}
		}
	}
	JETS.clear();
	scaleWeights.clear();
	pdfWeights.clear();
	alphaWeights.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void RUNBoostedAnalysis::beginJob() {

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
		//RUNAtree->Branch( "AK4HT", &AK4HT, "AK4HT/F" );
		//RUNAtree->Branch( "trimmedMass", &trimmedMass, "trimmedMass/F" );
		RUNAtree->Branch( "jet1Pt", &jet1Pt, "jet1Pt/F" );
		RUNAtree->Branch( "jet1Eta", &jet1Eta, "jet1Eta/F" );
		RUNAtree->Branch( "jet1Phi", &jet1Phi, "jet1Phi/F" );
		RUNAtree->Branch( "jet1E", &jet1E, "jet1E/F" );
		RUNAtree->Branch( "jet1M", &jet1M, "jet1M/F" );
		RUNAtree->Branch( "jet1PrunedMass", &jet1PrunedMass, "jet1PrunedMass/F" );
		RUNAtree->Branch( "jet1SoftDropMass", &jet1SoftDropMass, "jet1SoftDropMass/F" );
		RUNAtree->Branch( "jet1btagCSVv2", &jet1btagCSVv2, "jet1btagCSVv2/F" );
		RUNAtree->Branch( "jet1btagCMVAv2", &jet1btagCMVAv2, "jet1btagCMVAv2/F" );
		RUNAtree->Branch( "jet1btagDoubleB", &jet1btagDoubleB, "jet1btagDoubleB/F" );
		RUNAtree->Branch( "jet1btagCSVv2SF", &jet1btagCSVv2SF, "jet1btagCSVv2SF/F" );
		RUNAtree->Branch( "jet1btagCMVAv2SF", &jet1btagCMVAv2SF, "jet1btagCMVAv2SF/F" );
		RUNAtree->Branch( "jet2Pt", &jet2Pt, "jet2Pt/F" );
		RUNAtree->Branch( "jet2Eta", &jet2Eta, "jet2Eta/F" );
		RUNAtree->Branch( "jet2Phi", &jet2Phi, "jet2Phi/F" );
		RUNAtree->Branch( "jet2E", &jet2E, "jet2E/F" );
		RUNAtree->Branch( "jet2M", &jet2M, "jet2M/F" );
		RUNAtree->Branch( "jet2PrunedMass", &jet2PrunedMass, "jet2PrunedMass/F" );
		RUNAtree->Branch( "jet2SoftDropMass", &jet2SoftDropMass, "jet2SoftDropMass/F" );
		RUNAtree->Branch( "jet2btagCSVv2", &jet2btagCSVv2, "jet2btagCSVv2/F" );
		RUNAtree->Branch( "jet2btagCMVAv2", &jet2btagCMVAv2, "jet2btagCMVAv2/F" );
		RUNAtree->Branch( "jet2btagDoubleB", &jet2btagDoubleB, "jet2btagDoubleB/F" );
		RUNAtree->Branch( "jet2btagCSVv2SF", &jet2btagCSVv2SF, "jet2btagCSVv2SF/F" );
		RUNAtree->Branch( "jet2btagCMVAv2SF", &jet2btagCMVAv2SF, "jet2btagCMVAv2SF/F" );
		RUNAtree->Branch( "subjet11Pt", &subjet11Pt, "subjet11Pt/F" );
		RUNAtree->Branch( "subjet11Eta", &subjet11Eta, "subjet11Eta/F" );
		RUNAtree->Branch( "subjet11Phi", &subjet11Phi, "subjet11Phi/F" );
		RUNAtree->Branch( "subjet11E", &subjet11E, "subjet11E/F" );
		RUNAtree->Branch( "subjet11btagCSVv2", &subjet11btagCSVv2, "subjet11btagCSVv2/F" );
		RUNAtree->Branch( "subjet11btagCMVAv2", &subjet11btagCMVAv2, "subjet11btagCMVAv2/F" );
		RUNAtree->Branch( "subjet12Pt", &subjet12Pt, "subjet12Pt/F" );
		RUNAtree->Branch( "subjet12Eta", &subjet12Eta, "subjet12Eta/F" );
		RUNAtree->Branch( "subjet12Phi", &subjet12Phi, "subjet12Phi/F" );
		RUNAtree->Branch( "subjet12E", &subjet12E, "subjet12E/F" );
		RUNAtree->Branch( "subjet12btagCSVv2", &subjet12btagCSVv2, "subjet12btagCSVv2/F" );
		RUNAtree->Branch( "subjet12btagCMVAv2", &subjet12btagCMVAv2, "subjet12btagCMVAv2/F" );
		RUNAtree->Branch( "subjet21Pt", &subjet21Pt, "subjet21Pt/F" );
		RUNAtree->Branch( "subjet21Eta", &subjet21Eta, "subjet21Eta/F" );
		RUNAtree->Branch( "subjet21Phi", &subjet21Phi, "subjet21Phi/F" );
		RUNAtree->Branch( "subjet21E", &subjet21E, "subjet21E/F" );
		RUNAtree->Branch( "subjet21btagCSVv2", &subjet21btagCSVv2, "subjet21btagCSVv2/F" );
		RUNAtree->Branch( "subjet21btagCMVAv2", &subjet21btagCMVAv2, "subjet21btagCMVAv2/F" );
		RUNAtree->Branch( "subjet22Pt", &subjet22Pt, "subjet22Pt/F" );
		RUNAtree->Branch( "subjet22Eta", &subjet22Eta, "subjet22Eta/F" );
		RUNAtree->Branch( "subjet22Phi", &subjet22Phi, "subjet22Phi/F" );
		RUNAtree->Branch( "subjet22E", &subjet22E, "subjet22E/F" );
		RUNAtree->Branch( "subjet22btagCSVv2", &subjet22btagCSVv2, "subjet22btagCSVv2/F" );
		RUNAtree->Branch( "subjet22btagCMVAv2", &subjet22btagCMVAv2, "subjet22btagCMVAv2/F" );
		//RUNAtree->Branch( "massAve", &massAve, "massAve/F" );
		//RUNAtree->Branch( "massAsym", &massAsym, "massAsym/F" );
		//RUNAtree->Branch( "trimmedMassAve", &trimmedMassAve, "trimmedMassAve/F" );
		//RUNAtree->Branch( "trimmedMassAsym", &trimmedMassAsym, "trimmedMassAsym/F" );
		RUNAtree->Branch( "prunedMassAve", &prunedMassAve, "prunedMassAve/F" );
		RUNAtree->Branch( "prunedMassAsym", &prunedMassAsym, "prunedMassAsym/F" );
		//RUNAtree->Branch( "filteredMassAve", &filteredMassAve, "filteredMassAve/F" );
		//RUNAtree->Branch( "filteredMassAsym", &filteredMassAsym, "filteredMassAsym/F" );
		//RUNAtree->Branch( "softDropMassAve", &softDropMassAve, "softDropMassAve/F" );
		//RUNAtree->Branch( "softDropMassAsym", &softDropMassAsym, "softDropMassAsym/F" );
		RUNAtree->Branch( "deltaEtaDijet", &deltaEtaDijet, "deltaEtaDijet/F" );
		RUNAtree->Branch( "jet1CosThetaStar", &jet1CosThetaStar, "jet1CosThetaStar/F" );
		//RUNAtree->Branch( "jet2CosThetaStar", &jet2CosThetaStar, "jet2CosThetaStar/F" );
		RUNAtree->Branch( "jet1Tau21", &jet1Tau21, "jet1Tau21/F" );
		RUNAtree->Branch( "jet1Tau31", &jet1Tau31, "jet1Tau31/F" );
		RUNAtree->Branch( "jet1Tau32", &jet1Tau32, "jet1Tau32/F" );	
		RUNAtree->Branch( "jet2Tau21", &jet2Tau21, "jet2Tau21/F" );
		RUNAtree->Branch( "jet2Tau31", &jet2Tau31, "jet2Tau31/F" );
		RUNAtree->Branch( "jet2Tau32", &jet2Tau32, "jet2Tau32/F" );	
		//RUNAtree->Branch( "scaleWeights", &scaleWeights );
		//RUNAtree->Branch( "pdfWeights", &pdfWeights );
		//RUNAtree->Branch( "alphaWeights", &alphaWeights );
		RUNAtree->Branch( "genPartonPt1", &genPartonPt1, "genPartonPt1/F" );
		RUNAtree->Branch( "genPartonMass1", &genPartonMass1, "genPartonMass1/F" );
		RUNAtree->Branch( "genPartonDau11ID", &genPartonDau11ID, "genPartonDau11ID/F" );
		RUNAtree->Branch( "genPartonDau12ID", &genPartonDau12ID, "genPartonDau12ID/F" );
		RUNAtree->Branch( "genPartonPt2", &genPartonPt2, "genPartonPt2/F" );
		RUNAtree->Branch( "genPartonMass2", &genPartonMass2, "genPartonMass2/F" );
		RUNAtree->Branch( "genPartonDau21ID", &genPartonDau21ID, "genPartonDau21ID/F" );
		RUNAtree->Branch( "genPartonDau22ID", &genPartonDau22ID, "genPartonDau22ID/F" );
		RUNAtree->Branch( "muonsPt", "vector<float>", &muonsPt);
		RUNAtree->Branch( "numMuon", &numMuon, "numMuon/I" );
		RUNAtree->Branch( "elesPt", "vector<float>", &elesPt);
		RUNAtree->Branch( "numEle", &numEle, "numEle/I" );
	}


	histos1D_[ "oldJetPt" ] = fs_->make< TH1D >( "oldJetPt", "oldJetPt", 200, 0., 2000. );
	histos1D_[ "jetPt" ] = fs_->make< TH1D >( "jetPt", "jetPt", 200, 0., 2000. );
	histos1D_[ "rawJetPt" ] = fs_->make< TH1D >( "rawJetPt", "rawJetPt", 200, 0., 2000. );
	histos1D_[ "oldJetEta" ] = fs_->make< TH1D >( "oldJetEta", "oldJetEta", 100, -5., 5. );
	histos1D_[ "jetEta" ] = fs_->make< TH1D >( "jetEta", "jetEta", 100, -5., 5. );
	histos1D_[ "rawJetEta" ] = fs_->make< TH1D >( "rawJetEta", "rawJetEta", 100, -5., 5. );
	histos1D_[ "oldJetMass" ] = fs_->make< TH1D >( "oldJetMass", "oldJetMass", 600, 0., 600. );
	histos1D_[ "jetMass" ] = fs_->make< TH1D >( "jetMass", "jetMass", 600, 0., 600. );
	histos1D_[ "rawJetMass" ] = fs_->make< TH1D >( "rawJetMass", "rawJetMass", 600, 0., 600. );
	histos1D_[ "jetNum" ] = fs_->make< TH1D >( "jetNum", "jetNum", 10, 0., 10. );
	histos1D_[ "jetTrimmedMass" ] = fs_->make< TH1D >( "jetTrimmedMass", "jetTrimmedMass", 600, 0., 600. );
	histos1D_[ "jetPrunedMass" ] = fs_->make< TH1D >( "jetPrunedMass", "jetPrunedMass", 600, 0., 600. );
	histos1D_[ "jetFilteredMass" ] = fs_->make< TH1D >( "jetFilteredMass", "jetFilteredMass", 600, 0., 600. );
	histos1D_[ "jetSoftDropMass" ] = fs_->make< TH1D >( "jetSoftDropMass", "jetSoftDropMass", 600, 0., 600. );
	histos1D_[ "HT" ] = fs_->make< TH1D >( "HT", "HT", 500, 0., 5000. );
	histos1D_[ "NPV_NOPUWeight" ] = fs_->make< TH1D >( "NPV_NOPUWeight", "NPV_NOPUWeight", 80, 0., 80. );
	histos1D_[ "NPV" ] = fs_->make< TH1D >( "NPV", "NPV", 80, 0., 80. );
	histos1D_[ "PUWeight" ] = fs_->make< TH1D >( "PUWeight", "PUWeight", 50, 0., 5. );
	histos2D_[ "jetTrimmedMassHT" ] = fs_->make< TH2D >( "jetTrimmedMassHT", "jetTrimmedMassHT", 30, 0., 300., 500, 0., 5000. );
	histos1D_[ "neutralHadronEnergyFrac" ] = fs_->make< TH1D >( "neutralHadronEnergyFrac", "neutralHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "neutralEmEnergyFrac" ] = fs_->make< TH1D >( "neutralEmEnergyFrac", "neutralEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedEmEnergyFrac" ] = fs_->make< TH1D >( "chargedEmEnergyFrac", "chargedEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedHadronEnergyFrac" ] = fs_->make< TH1D >( "chargedHadronEnergyFrac", "chargedHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedMultiplicity" ] = fs_->make< TH1D >( "chargedMultiplicity", "chargedMultiplicity", 50, 0., 1. );
	histos1D_[ "numConst" ] = fs_->make< TH1D >( "numConst", "numConst", 200, 0., 200. );

	histos1D_[ "HT_cutTrigger" ] = fs_->make< TH1D >( "HT_cutTrigger", "HT_cutTrigger", 500, 0., 5000. );
	histos1D_[ "MET_cutTrigger" ] = fs_->make< TH1D >( "MET_cutTrigger", "MET_cutTrigger", 20, 0., 200. );
	histos1D_[ "METHT_cutTrigger" ] = fs_->make< TH1D >( "METHT_cutTrigger", "METHT_cutTrigger", 50, 0., 1. );
	histos1D_[ "NPV_cutTrigger" ] = fs_->make< TH1D >( "NPV_cutTrigger", "NPV_cutTrigger", 80, 0., 80. );
	histos1D_[ "jetPt_cutTrigger" ] = fs_->make< TH1D >( "jetPt_cutTrigger", "jetPt_cutTrigger", 200, 0., 2000. );
	histos1D_[ "jetEta_cutTrigger" ] = fs_->make< TH1D >( "jetEta_cutTrigger", "jetEta_cutTrigger", 100, -5., 5. );
	histos1D_[ "jetNum_cutTrigger" ] = fs_->make< TH1D >( "jetNum_cutTrigger", "jetNum_cutTrigger", 10, 0., 10. );
	histos1D_[ "jetMass_cutTrigger" ] = fs_->make< TH1D >( "jetMass_cutTrigger", "jetMass_cutTrigger", 600, 0., 600. );
	histos1D_[ "jet1Pt_cutTrigger" ] = fs_->make< TH1D >( "jet1Pt_cutTrigger", "jet1Pt_cutTrigger", 200, 0., 2000. );
	histos1D_[ "jet1Eta_cutTrigger" ] = fs_->make< TH1D >( "jet1Eta_cutTrigger", "jet1Eta_cutTrigger", 100, -5., 5. );
	histos1D_[ "jet1Mass_cutTrigger" ] = fs_->make< TH1D >( "jet1Mass_cutTrigger", "jet1Mass_cutTrigger", 600, 0., 600. );
	histos2D_[ "jetTrimmedMassHT_cutTrigger" ] = fs_->make< TH2D >( "jetTrimmedMassHT_cutTrigger", "jetTrimmedMassHT_cutTrigger", 30, 0., 300., 500, 0., 5000. );

	histos1D_[ "HT_cutDijet" ] = fs_->make< TH1D >( "HT_cutDijet", "HT_cutDijet", 500, 0., 5000. );
	histos1D_[ "NPV_cutDijet" ] = fs_->make< TH1D >( "NPV_cutDijet", "NPV_cutDijet", 80, 0., 80. );
	histos1D_[ "jetNum_cutDijet" ] = fs_->make< TH1D >( "jetNum_cutDijet", "jetNum_cutDijet", 10, 0., 10. );
	histos1D_[ "jet1Pt_cutDijet" ] = fs_->make< TH1D >( "jet1Pt_cutDijet", "jet1Pt_cutDijet", 200, 0., 2000. );
	histos1D_[ "jet1Eta_cutDijet" ] = fs_->make< TH1D >( "jet1Eta_cutDijet", "jet1Eta_cutDijet", 100, -5., 5. );
	histos1D_[ "jet1Mass_cutDijet" ] = fs_->make< TH1D >( "jet1Mass_cutDijet", "jet1Mass_cutDijet", 600, 0., 600. );
	histos1D_[ "jet2Pt_cutDijet" ] = fs_->make< TH1D >( "jet2Pt_cutDijet", "jet2Pt_cutDijet", 200, 0., 2000. );
	histos1D_[ "jet2Eta_cutDijet" ] = fs_->make< TH1D >( "jet2Eta_cutDijet", "jet2Eta_cutDijet", 100, -5., 5. );
	histos1D_[ "jet2Mass_cutDijet" ] = fs_->make< TH1D >( "jet2Mass_cutDijet", "jet2Mass_cutDijet", 600, 0., 600. );
	histos1D_[ "neutralHadronEnergyFrac_cutDijet" ] = fs_->make< TH1D >( "neutralHadronEnergyFrac_cutDijet", "neutralHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "neutralEmEnergyFrac_cutDijet" ] = fs_->make< TH1D >( "neutralEmEnergyFrac_cutDijet", "neutralEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedHadronEnergyFrac_cutDijet" ] = fs_->make< TH1D >( "chargedHadronEnergyFrac_cutDijet", "chargedHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedEmEnergyFrac_cutDijet" ] = fs_->make< TH1D >( "chargedEmEnergyFrac_cutDijet", "chargedEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedMultiplicity_cutDijet" ] = fs_->make< TH1D >( "chargedMultiplicity_cutDijet", "chargedMultiplicity", 50, 0., 1. );
	histos1D_[ "numConst_cutDijet" ] = fs_->make< TH1D >( "numConst_cutDijet", "numConst", 200, 0., 200. );
	histos1D_[ "MET_cutDijet" ] = fs_->make< TH1D >( "MET_cutDijet", "MET_cutDijet", 20, 0., 200. );
	histos1D_[ "METHT_cutDijet" ] = fs_->make< TH1D >( "METHT_cutDijet", "METHT_cutDijet", 50, 0., 1. );


	histos1D_[ "HT_cutEffTrigger" ] = fs_->make< TH1D >( "HT_cutEffTrigger", "HT_cutEffTrigger", 500, 0., 5000. );
	histos1D_[ "NPV_cutEffTrigger" ] = fs_->make< TH1D >( "NPV_cutEffTrigger", "NPV_cutEffTrigger", 80, 0., 80. );
	histos1D_[ "jetNum_cutEffTrigger" ] = fs_->make< TH1D >( "jetNum_cutEffTrigger", "jetNum_cutEffTrigger", 10, 0., 10. );

	histos1D_[ "jet1NeutralHadronEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet1NeutralHadronEnergyFrac_cutEffTrigger", "jet1NeutralHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet1NeutralEmEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet1NeutralEmEnergyFrac_cutEffTrigger", "jet1NeutralEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet1ChargedHadronEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet1ChargedHadronEnergyFrac_cutEffTrigger", "jet1ChargedHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet1ChargedEmEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet1ChargedEmEnergyFrac_cutEffTrigger", "jet1ChargedEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet1ChargedMultiplicity_cutEffTrigger" ] = fs_->make< TH1D >( "jet1ChargedMultiplicity_cutEffTrigger", "jet1ChargedMultiplicity", 50, 0., 1. );
	histos1D_[ "jet1NumConst_cutEffTrigger" ] = fs_->make< TH1D >( "jet1NumConst_cutEffTrigger", "jet1NumConst", 200, 0., 200. );
	histos1D_[ "jet1Pt_cutEffTrigger" ] = fs_->make< TH1D >( "jet1Pt_cutEffTrigger", "jet1Pt_cutEffTrigger", 200, 0., 2000. );
	histos1D_[ "jet1Eta_cutEffTrigger" ] = fs_->make< TH1D >( "jet1Eta_cutEffTrigger", "jet1Eta_cutEffTrigger", 100, -5., 5. );
	histos1D_[ "jet1Mass_cutEffTrigger" ] = fs_->make< TH1D >( "jet1Mass_cutEffTrigger", "jet1Mass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet1PrunedMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet1PrunedMass_cutEffTrigger", "jet1PrunedMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet1TrimmedMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet1TrimmedMass_cutEffTrigger", "jet1TrimmedMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet1FilteredMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet1FilteredMass_cutEffTrigger", "jet1FilteredMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet1SoftDropMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet1SoftDropMass_cutEffTrigger", "jet1SoftDropMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet2NeutralHadronEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet2NeutralHadronEnergyFrac_cutEffTrigger", "jet2NeutralHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet2NeutralEmEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet2NeutralEmEnergyFrac_cutEffTrigger", "jet2NeutralEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet2ChargedHadronEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet2ChargedHadronEnergyFrac_cutEffTrigger", "jet2ChargedHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet2ChargedEmEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "jet2ChargedEmEnergyFrac_cutEffTrigger", "jet2ChargedEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "jet2ChargedMultiplicity_cutEffTrigger" ] = fs_->make< TH1D >( "jet2ChargedMultiplicity_cutEffTrigger", "jet2ChargedMultiplicity", 50, 0., 1. );
	histos1D_[ "jet2NumConst_cutEffTrigger" ] = fs_->make< TH1D >( "jet2NumConst_cutEffTrigger", "jet2NumConst", 200, 0., 200. );
	histos1D_[ "jet2Pt_cutEffTrigger" ] = fs_->make< TH1D >( "jet2Pt_cutEffTrigger", "jet2Pt_cutEffTrigger", 200, 0., 2000. );
	histos1D_[ "jet2Eta_cutEffTrigger" ] = fs_->make< TH1D >( "jet2Eta_cutEffTrigger", "jet2Eta_cutEffTrigger", 100, -5., 5. );
	histos1D_[ "jet2Mass_cutEffTrigger" ] = fs_->make< TH1D >( "jet2Mass_cutEffTrigger", "jet2Mass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet2PrunedMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet2PrunedMass_cutEffTrigger", "jet2PrunedMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet2TrimmedMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet2TrimmedMass_cutEffTrigger", "jet2TrimmedMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet2SoftDropMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet2SoftDropMass_cutEffTrigger", "jet2SoftDropMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "jet2FilteredMass_cutEffTrigger" ] = fs_->make< TH1D >( "jet2FilteredMass_cutEffTrigger", "jet2FilteredMass_cutEffTrigger", 600, 0., 600. );
	histos1D_[ "MET_cutEffTrigger" ] = fs_->make< TH1D >( "MET_cutEffTrigger", "MET_cutEffTrigger", 20, 0., 200. );
	histos1D_[ "METHT_cutEffTrigger" ] = fs_->make< TH1D >( "METHT_cutEffTrigger", "METHT_cutEffTrigger", 50, 0., 1. );
	histos1D_[ "neutralHadronEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "neutralHadronEnergyFrac_cutEffTrigger", "neutralHadronEnergyFrac_cutEffTrigger", 50, 0., 1. );
	histos1D_[ "neutralEmEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "neutralEmEnergyFrac_cutEffTrigger", "neutralEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedHadronEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "chargedHadronEnergyFrac_cutEffTrigger", "chargedHadronEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedEmEnergyFrac_cutEffTrigger" ] = fs_->make< TH1D >( "chargedEmEnergyFrac_cutEffTrigger", "chargedEmEnergyFrac", 50, 0., 1. );
	histos1D_[ "chargedMultiplicity_cutEffTrigger" ] = fs_->make< TH1D >( "chargedMultiplicity_cutEffTrigger", "chargedMultiplicity", 50, 0., 1. );
	histos1D_[ "numConst_cutEffTrigger" ] = fs_->make< TH1D >( "numConst_cutEffTrigger", "numConst", 200, 0., 200. );

	histos1D_[ "jet1Tau21_cutEffTrigger" ] = fs_->make< TH1D >( "jet1Tau21_cutEffTrigger", "jet1Tau21", 20, 0., 1. );
	histos1D_[ "jet2Tau21_cutEffTrigger" ] = fs_->make< TH1D >( "jet2Tau21_cutEffTrigger", "jet2Tau21", 20, 0., 1. );
	histos1D_[ "prunedMassAsym_cutEffTrigger" ] = fs_->make< TH1D >( "prunedMassAsym_cutEffTrigger", "prunedMassAsym", 20, 0., 1. );
	histos1D_[ "deltaEtaDijet_cutEffTrigger" ] = fs_->make< TH1D >( "deltaEtaDijet_cutEffTrigger", "deltaEtaDijet", 100, 0., 5. );
	histos1D_[ "massAve_cutEffTrigger" ] = fs_->make< TH1D >( "massAve_cutEffTrigger", "massAve", 500, 0., 500. );
	histos1D_[ "jet1btagCSVv2_cutEffTrigger" ] = fs_->make< TH1D >( "jet1btagCSVv2_cutEffTrigger", "jet1btagCSVv2", 20, 0., 1. );
	histos1D_[ "jet2btagCSVv2_cutEffTrigger" ] = fs_->make< TH1D >( "jet2btagCSVv2_cutEffTrigger", "jet2btagCSVv2", 20, 0., 1. );

	histos1D_[ "jet1Tau21_cutTau21" ] = fs_->make< TH1D >( "jet1Tau21_cutTau21", "jet1Tau21", 20, 0., 1. );
	histos1D_[ "jet2Tau21_cutTau21" ] = fs_->make< TH1D >( "jet2Tau21_cutTau21", "jet2Tau21", 20, 0., 1. );
	histos1D_[ "prunedMassAsym_cutTau21" ] = fs_->make< TH1D >( "prunedMassAsym_cutTau21", "prunedMassAsym", 20, 0., 1. );
	histos1D_[ "deltaEtaDijet_cutTau21" ] = fs_->make< TH1D >( "deltaEtaDijet_cutTau21", "deltaEtaDijet", 100, 0., 5. );
	histos1D_[ "massAve_cutTau21" ] = fs_->make< TH1D >( "massAve_cutTau21", "massAve", 500, 0., 500. );

	histos1D_[ "prunedMassAsym_cutMassAsym" ] = fs_->make< TH1D >( "prunedMassAsym_cutMassAsym", "prunedMassAsym", 20, 0., 1. );
	histos1D_[ "deltaEtaDijet_cutMassAsym" ] = fs_->make< TH1D >( "deltaEtaDijet_cutMassAsym", "deltaEtaDijet", 100, 0., 5. );
	histos1D_[ "massAve_cutMassAsym" ] = fs_->make< TH1D >( "massAve_cutMassAsym", "massAve", 500, 0., 500. );

	histos1D_[ "deltaEtaDijet_cutDeltaEtaDijet" ] = fs_->make< TH1D >( "deltaEtaDijet_cutDeltaEtaDijet", "deltaEtaDijet", 100, 0., 5. );
	histos1D_[ "massAve_cutDeltaEtaDijet" ] = fs_->make< TH1D >( "massAve_cutDeltaEtaDijet", "massAve", 500, 0., 500. );

	histos1D_[ "massAve_cutBtag" ] = fs_->make< TH1D >( "massAve_cutBtag", "massAve", 500, 0., 500. );

	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_A" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_A", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_A" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_A", "prunedMassAsymVsdeltaEtaDijet_A", 20, 0., 1., 100, 0., 5. );
	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_B" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_B", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_B" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_B", "prunedMassAsymVsdeltaEtaDijet_B", 20, 0., 1., 100, 0., 5. );
	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_C" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_C", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_C" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_C", "prunedMassAsymVsdeltaEtaDijet_C", 20, 0., 1., 100, 0., 5. );
	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_D" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_D", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_D" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_D", "prunedMassAsymVsdeltaEtaDijet_D", 20, 0., 1., 100, 0., 5. );

	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_A" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_btag_A", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_A" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_btag_A", "prunedMassAsymVsdeltaEtaDijet_btag_A", 20, 0., 1., 100, 0., 5. );
	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_B" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_btag_B", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_B" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_btag_B", "prunedMassAsymVsdeltaEtaDijet_btag_B", 20, 0., 1., 100, 0., 5. );
	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_C" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_btag_C", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_C" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_btag_C", "prunedMassAsymVsdeltaEtaDijet_btag_C", 20, 0., 1., 100, 0., 5. );
	histos1D_[ "massAve_prunedMassAsymVsdeltaEtaDijet_btag_D" ] = fs_->make< TH1D >( "massAve_prunedMassAsymVsdeltaEtaDijet_btag_D", "massAve", 500, 0., 500. );
	histos2D_[ "prunedMassAsymVsdeltaEtaDijet_btag_D" ] = fs_->make< TH2D >( "prunedMassAsymVsdeltaEtaDijet_btag_D", "prunedMassAsymVsdeltaEtaDijet_btag_D", 20, 0., 1., 100, 0., 5. );

	histos1D_[ "jet1Tau21_n-1" ] = fs_->make< TH1D >( "jet1Tau21_n-1", "jet1Tau21", 20, 0., 1. );
	histos1D_[ "jet2Tau21_n-1" ] = fs_->make< TH1D >( "jet2Tau21_n-1", "jet2Tau21", 20, 0., 1. );
	histos1D_[ "prunedMassAsym_n-1" ] = fs_->make< TH1D >( "prunedMassAsym_n-1", "prunedMassAsym", 20, 0., 1. );
	histos1D_[ "deltaEtaDijet_n-1" ] = fs_->make< TH1D >( "deltaEtaDijet_n-1", "deltaEtaDijet", 100, 0., 5. );

	///// Sumw2 all the histos
	for( auto const& histo : histos1D_ ) histos1D_[ histo.first ]->Sumw2();
	for( auto const& histo : histos2D_ ) histos2D_[ histo.first ]->Sumw2();


	cutLabels.push_back("Processed");
	cutLabels.push_back("Trigger");
	cutLabels.push_back("Dijet");
	cutLabels.push_back("HT");
	cutLabels.push_back("Tau21");
	cutLabels.push_back("MassAsym");
	cutLabels.push_back("DeltaEtaDijet");
	cutLabels.push_back("Btag");
	histos1D_[ "hcutflow" ] = fs_->make< TH1D >("cutflow","cut flow", cutLabels.size(), 0.5, cutLabels.size() +0.5 );
	histos1D_[ "hcutflow" ]->Sumw2();
	for( const string &ivec : cutLabels ) cutmap[ ivec ] = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void RUNBoostedAnalysis::endJob() {

	int ibin = 1;
	for( const string &ivec : cutLabels ) {
		histos1D_["hcutflow"]->SetBinContent( ibin, cutmap[ ivec ] );
		histos1D_["hcutflow"]->GetXaxis()->SetBinLabel( ibin, ivec.c_str() );
		ibin++;
	}

}

void RUNBoostedAnalysis::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<bool>("isData", false);
	desc.add<bool>("isTTbar", false);
	desc.add<bool>("LHEcont", false);
	desc.add<bool>("mkTree", false);
	desc.add<bool>("sortInTau21", false);
	desc.add<bool>("sortInMass", false);
	desc.add<string>("dataPUFile", "supportFiles/PileupData2015D_JSON_latest.root");
	desc.add<string>("jecVersion", "supportFiles/Fall15_25nsV2");
	desc.add<string>("btagCSVFile", "supportFiles/CSVv2_Moriond17_B_H.csv");
	desc.add<string>("btagMVAFile", "supportFiles/cMVAv2_Moriond17_B_H.csv");
	desc.add<string>("systematics", "None");
	desc.add<string>("PUMethod", "chs");
	desc.add<double>("scale", 1);
	desc.add<double>("cutAK8HT", 900.);
	desc.add<double>("cutAK8jetPt", 150.);
	desc.add<double>("cutAK8MassAsym", 0.1);
	desc.add<double>("cutTau21", 0.45);
	desc.add<double>("cutDeltaEtaDijet", 1.5);
	vector<string> HLTPass;
	HLTPass.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50");
	desc.add<vector<string>>("triggerPass",	HLTPass);

	desc.add<InputTag>("Lumi", 	InputTag("eventInfo:evtInfoLumiBlock"));
	desc.add<InputTag>("Run", 	InputTag("eventInfo:evtInfoRunNumber"));
	desc.add<InputTag>("Event", 	InputTag("eventInfo:evtInfoEventNumber"));
	desc.add<InputTag>("generator", InputTag("generator"));
	desc.add<InputTag>("extLHEProducer", 	InputTag("externalLHEProducer"));
	desc.add<InputTag>("genParticles", 	InputTag("filteredPrunedGenParticles"));
	desc.add<InputTag>("bunchCross", 	InputTag("eventUserData:puBX"));
	desc.add<InputTag>("rho", 	InputTag("vertexInfo:rho"));
	desc.add<InputTag>("puNumInt", 	InputTag("eventUserData:puNInt"));
	desc.add<InputTag>("trueNInt", 	InputTag("eventUserData:puNtrueInt"));
	desc.add<InputTag>("NPV", 	InputTag("vertexInfo:npv"));
	desc.add<InputTag>("jetAK4Pt", 	InputTag("jetsAK4CHS:jetAK4CHSPt"));
	desc.add<InputTag>("jetAK4Eta", 	InputTag("jetsAK4CHS:jetAK4CHSEta"));
	desc.add<InputTag>("jetAK4Phi", 	InputTag("jetsAK4CHS:jetAK4CHSPhi"));
	desc.add<InputTag>("jetAK4E", 	InputTag("jetsAK4CHS:jetAK4CHSE"));
	desc.add<InputTag>("jetPt", 	InputTag("jetsAK8CHS:jetAK8CHSPt"));
	desc.add<InputTag>("jetEta", 	InputTag("jetsAK8CHS:jetAK8CHSEta"));
	desc.add<InputTag>("jetPhi", 	InputTag("jetsAK8CHS:jetAK8CHSPhi"));
	desc.add<InputTag>("jetE", 	InputTag("jetsAK8CHS:jetAK8CHSE"));
	desc.add<InputTag>("jetTrimmedMass", 	InputTag("jetsAK8CHS:jetAK8CHStrimmedMassCHS"));
	desc.add<InputTag>("jetPrunedMass", 	InputTag("jetsAK8CHS:jetAK8CHSprunedMassCHS"));
	desc.add<InputTag>("jetFilteredMass", 	InputTag("jetsAK8CHS:jetAK8CHSfilteredMassCHS"));
	desc.add<InputTag>("jetSoftDropMass", 	InputTag("jetsAK8CHS:jetAK8CHSsoftDropMassCHS"));
	desc.add<InputTag>("jetTau1", 	InputTag("jetsAK8CHS:jetAK8CHStau1CHS"));
	desc.add<InputTag>("jetTau2", 	InputTag("jetsAK8CHS:jetAK8CHStau2CHS"));
	desc.add<InputTag>("jetTau3", 	InputTag("jetsAK8CHS:jetAK8CHStau3CHS"));
	desc.add<InputTag>("jetNSubjets", 	InputTag("jetsAK8CHS:jetAK8CHSnSubjets"));
	desc.add<InputTag>("jetSubjetIndex0", 	InputTag("jetsAK8CHS:jetAK8CHSvSubjetIndex0"));
	desc.add<InputTag>("jetSubjetIndex1", 	InputTag("jetsAK8CHS:jetAK8CHSvSubjetIndex1"));
	desc.add<InputTag>("jetSubjetIndex2", 	InputTag("jetsAK8CHS:jetAK8CHSvSubjetIndex2"));
	desc.add<InputTag>("jetSubjetIndex3", 	InputTag("jetsAK8CHS:jetAK8CHSvSubjetIndex3"));
	desc.add<InputTag>("jetKeys", 	InputTag("jetKeysAK8CHS"));
	desc.add<InputTag>("jetCSVv2", 	InputTag("jetsAK8CHS:jetAK8CHSCSVv2"));
	desc.add<InputTag>("jetCMVAv2", 	InputTag("jetsAK8CHS:jetAK8CHSCMVAv2"));
	desc.add<InputTag>("jetDoubleB", 	InputTag("jetsAK8CHS:jetAK8CHSDoubleBAK8"));
	desc.add<InputTag>("jetArea", 	InputTag("jetsAK8CHS:jetAK8CHSjetArea"));
	desc.add<InputTag>("jetGenPt", 	InputTag("jetsAK8CHS:jetAK8CHSGenJetPt"));
	desc.add<InputTag>("jetGenEta", 	InputTag("jetsAK8CHS:jetAK8CHSGenJetEta"));
	desc.add<InputTag>("jetGenPhi", 	InputTag("jetsAK8CHS:jetAK8CHSGenJetPhi"));
	desc.add<InputTag>("jetGenE", 	InputTag("jetsAK8CHS:jetAK8CHSGenJetE"));
	desc.add<InputTag>("jetHadronFlavour", 	InputTag("jetsAK8CHS:jetAK8CHSHadronFlavour"));
	desc.add<InputTag>("metPt", 	InputTag("metFull:metFullPt"));
	// JetID
	desc.add<InputTag>("jecFactor", 		InputTag("jetsAK8CHS:jetAK8CHSjecFactor0"));
	desc.add<InputTag>("neutralHadronEnergyFrac", 	InputTag("jetsAK8CHS:jetAK8CHSneutralHadronEnergyFrac"));
	desc.add<InputTag>("neutralEmEnergyFrac", 		InputTag("jetsAK8CHS:jetAK8CHSneutralEmEnergyFrac"));
	desc.add<InputTag>("chargedEmEnergyFrac", 		InputTag("jetsAK8CHS:jetAK8CHSchargedEmEnergyFrac"));
	desc.add<InputTag>("muonEnergyFrac", 		InputTag("jetsAK8CHS:jetAK8CHSMuonEnergy"));
	desc.add<InputTag>("chargedHadronEnergyFrac", 	InputTag("jetsAK8CHS:jetAK8CHSchargedHadronEnergyFrac"));
	desc.add<InputTag>("neutralMultiplicity",	InputTag("jetsAK8CHS:jetAK8CHSneutralMultiplicity"));
	desc.add<InputTag>("chargedMultiplicity", 	InputTag("jetsAK8CHS:jetAK8CHSchargedMultiplicity"));
	// Trigger
	desc.add<InputTag>("triggerPrescale",		InputTag("TriggerUserData:triggerPrescaleTree"));
	desc.add<InputTag>("triggerBit",		InputTag("TriggerUserData:triggerBitTree"));
	desc.add<InputTag>("triggerName",		InputTag("TriggerUserData:triggerNameTree"));
	// Subjets
	desc.add<InputTag>("subjetPt", 	InputTag("subjetsAK8CHS:subjetAK8CHSPt"));
	desc.add<InputTag>("subjetEta", 	InputTag("subjetsAK8CHS:subjetAK8CHSEta"));
	desc.add<InputTag>("subjetPhi", 	InputTag("subjetsAK8CHS:subjetAK8CHSPhi"));
	desc.add<InputTag>("subjetE", 	InputTag("subjetsAK8CHS:subjetAK8CHSE"));
	desc.add<InputTag>("subjetMass", 	InputTag("subjetsAK8CHS:subjetAK8CHSMass"));
	desc.add<InputTag>("subjetCSVv2", 	InputTag("subjetsAK8CHS:subjetAK8CHSCSVv2"));
	desc.add<InputTag>("subjetCMVAv2", 	InputTag("subjetsAK8CHS:subjetAK8CHSCMVAv2"));
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
      
void RUNBoostedAnalysis::beginRun(const Run& iRun, const EventSetup& iSetup){

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

void RUNBoostedAnalysis::endRun(const Run& iRun, const EventSetup& iSetup){
	triggerNamesList.clear();
}



//define this as a plug-in
DEFINE_FWK_MODULE(RUNBoostedAnalysis);
