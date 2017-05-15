#include "RUNA/RUNAnalysis/interface/PUReweighter.h"


//==============================================================================================
// Get weight factor dependent on number of added PU interactions
float PUReweighter::getPUWeight(const int trueNInt, std::vector<int> bunchCrossing ){

	//std::vector<float> puWeigths_;
	//generateWeights( nameOfDataDistribution );

	float trueNPV = -1;
	for (unsigned int i = 0; i < bunchCrossing.size(); i++) {
		if( bunchCrossing[i] == 0 ) trueNPV = trueNInt; 
	}
	float w = LumiWeights_.weight( trueNPV );
	return w;
}

//==============================================================================================
// Generate weights for given data PU distribution
// Code adapted from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
// weights for Moriond17 are taken from https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py
void PUReweighter::generateWeights(const std::string& nameOfDataDistribution) {

	// Get data distribution from file
	TFile file(nameOfDataDistribution.c_str(), "READ");
	TH1* data_npu_estimated = NULL;
	file.GetObject("pileup",data_npu_estimated);
	if( data_npu_estimated == NULL ) {
		std::cerr << "\n\nERROR in PUReweighter: Histogram 'pileup' does not exist in file '" << nameOfDataDistribution << "'\n.";
		throw std::exception();
	}
	data_npu_estimated->SetDirectory(0);
	file.Close();

	// Store probabilites for each pu bin
	unsigned int trueNIntMax = 0;

	trueNIntMax = 75;
	float npuWinter15_25ns[ trueNIntMax ] = {
			0.0000178653,
			0.0000256602,
			0.0000527857,
			0.0000888954,
			0.000109362,
			0.000140973,
			0.000240998,
			0.00071209,
			0.00130121,
			0.00245255,
			0.00502589,
			0.00919534,
			0.0146697,
			0.0204126,
			0.0267586,
			0.0337697,
			0.0401478,
			0.0450159,
			0.0490577,
			0.0524855,
			0.0548159,
			0.0559937,
			0.0554468,
			0.0537687,
			0.0512055,
			0.0476713,
			0.0435312,
			0.0393107,
			0.0349812,
			0.0307413,
			0.0272425,
			0.0237115,
			0.0208329,
			0.0182459,
			0.0160712,
			0.0142498,
			0.012804,
			0.011571,
			0.010547,
			0.00959489,
			0.00891718,
			0.00829292,
			0.0076195,
			0.0069806,
			0.0062025,
			0.00546581,
			0.00484127,
			0.00407168,
			0.00337681,
			0.00269893,
			0.00212473,
			0.00160208,
			0.00117884,
			0.000859662,
			0.000569085,
			0.000365431,
			0.000243565,
			0.00015688,
			0.0000988128,
			0.0000653783,
			0.0000373924,
			0.0000261382,
			0.000020307,
			0.0000173032,
			0.00001435,
			0.0000136486,
			0.0000135555,
			0.0000137491,
			0.0000134255,
			0.0000133987,
			0.0000134061,
			0.0000134211,
			0.0000134177,
			0.0000132959,
			0.0000133287
		/*
		0.000829312873542,
 		0.00124276120498,
 		0.00339329181587,
 		0.00408224735376,
 		0.00383036590008,
		0.00659159288946,
 		0.00816022734493,
 		0.00943640833116,
 		0.0137777376066,
 		0.017059392038,
 		0.0213193035468,
 		0.0247343174676,
 		0.0280848773878,
 		0.0323308476564,
 		0.0370394341409,
 		0.0456917721191,
 		0.0558762890594,
 		0.0576956187107,
 		0.0625325287017,
 		0.0591603758776,
 		0.0656650815128,
 		0.0678329011676,
 		0.0625142146389,
 		0.0548068448797,
 		0.0503893295063,
 		0.040209818868,
 		0.0374446988111,
 		0.0299661572042,
 		0.0272024759921,
 		0.0219328403791,
 		0.0179586571619,
 		0.0142926728247,
 		0.00839941654725,
 		0.00522366397213,
 		0.00224457976761,
 		0.000779274977993,
 		0.000197066585944,
 		7.16031761328e-05,
 		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
 		0.0,
 		0.0,
		0.0 */
	};

	// Check that binning of data-profile matches MC scenario
	if( trueNIntMax != static_cast<unsigned int>(data_npu_estimated->GetNbinsX()) ) {
		std::cerr << "\n\nERROR number of bins (" << data_npu_estimated->GetNbinsX() << ") in data PU-profile does not match number of bins (" << trueNIntMax << ")" << std::endl;
		throw std::exception();
	}

	std::vector<float> result(trueNIntMax,0.);
	std::vector<float> puMC;
	std::vector<float> puData;
	for(unsigned int npu = 0; npu < trueNIntMax; ++npu) {
		puMC.push_back( npuWinter15_25ns[npu] );
		puData.push_back( data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu)) );
	}

	LumiWeights_ = edm::LumiReWeighting(puMC, puData);
}

