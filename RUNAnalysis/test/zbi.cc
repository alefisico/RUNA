#include <TRandom3.h>
#include <TMath.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>

#include <vector>
#include <sstream>
#include <cassert>

double pValueToZScore(double pvalue)
{
  if(pvalue<=0.0 && numeric_limits<double>::has_infinity) {
    return numeric_limits<double>::infinity();
  } else if(pvalue<=0.0) {
    return numeric_limits<double>::max();
  } else if(pvalue>=1.0 && numeric_limits<double>::has_infinity) {
    return -numeric_limits<double>::infinity();
  } else if(pvalue>=1.0) {
    return -numeric_limits<double>::max();
  } else if(pvalue==0.5) {
    return 0.0;
  } else {
    return TMath::Sqrt(2.0)*TMath::ErfcInverse(2*pvalue);
  }
}

double zScoreToPValue(double zscore)
{
  if(numeric_limits<double>::has_infinity && zscore==numeric_limits<double>::infinity()) {
    return 0.0;
  } else if(!numeric_limits<double>::has_infinity && zscore==numeric_limits<double>::max()) {
    return 0.0;
  } else if(numeric_limits<double>::has_infinity && zscore==-numeric_limits<double>::infinity()) {
    return 1.0;
  } else if(!numeric_limits<double>::has_infinity && zscore==-numeric_limits<double>::max()) {
    return 1.0;
  } else {
    return TMath::Erfc(zscore/TMath::Sqrt(2.0))/2.0;
  }
}


void getQuantiles(std::vector<double>& limits, double &median_, std::pair<double, double>& onesigma_, std::pair<double, double>& twosigma_) {
  unsigned int nit=limits.size();
  if(nit==0) return;

  // sort the vector with limits
  std::sort(limits.begin(), limits.end());

  // median for the expected limit
  median_ = TMath::Median(nit, &limits[0]);

  // quantiles for the expected limit bands
  double prob[4]; // array with quantile boundaries
  prob[0] = 0.021;
  prob[1] = 0.159;
  prob[2] = 0.841;
  prob[3] = 0.979;

  // array for the results
  double quantiles[4];
  TMath::Quantiles(nit, 4, &limits[0], quantiles, prob); // evaluate quantiles
  onesigma_.first=quantiles[1];
  onesigma_.second=quantiles[2];
  twosigma_.first=quantiles[0];
  twosigma_.second=quantiles[3];
    
  return;
}


// lognormal distribution where par[0] is the median and par[1] is the variance
double lognormal(double *x, double *par)
{
  if(par[0]<0.0) {
    std::cout << "par[0] = " << par[0] << std::endl;
    assert(0);
  }
  if(par[1]<0.0) {
    std::cout << "par[1] = " << par[1] << std::endl;
    assert(0);
  }

  double m0=par[0];
  double k=par[1]/par[0]+1.;
  double s=TMath::Log(k);
  return TMath::LogNormal(x[0], s, 0.0, m0);
}

// get a random lognormal value given a value for the median and variance
class RandomLognormal
{
public:
  RandomLognormal(double median, double variance);
  RandomLognormal(double median, double variance, double min, double max);
  virtual ~RandomLognormal();

  double getRandom(void) const;

private:

  TF1* lognormfcn_;
  static int counter_;
};

RandomLognormal::RandomLognormal(double median, double variance)
{
  double min = std::max(0.0, median-5.0*variance);
  double max = (median+5.0*variance);

  // create function
  std::ostringstream oss;
  oss << "_Random_Lognormal__lognormfcn_" << (counter_++);
  lognormfcn_ = new TF1(oss.str().c_str(), lognormal, min, max, 2);
  lognormfcn_->SetParameter(0, median);
  lognormfcn_->SetParameter(1, variance);
}

RandomLognormal::RandomLognormal(double median, double variance, double min, double max)
{
  // create function
  std::ostringstream oss;
  oss << "_Random_Lognormal__lognormfcn_" << (counter_++);
  lognormfcn_ = new TF1(oss.str().c_str(), lognormal, min, max, 2);
  lognormfcn_->SetParameter(0, median);
  lognormfcn_->SetParameter(1, variance);
}

RandomLognormal::~RandomLognormal()
{
  delete lognormfcn_;
}

double RandomLognormal::getRandom(void) const
{
  return lognormfcn_->GetRandom();
}

int RandomLognormal::counter_=0;


// obs = number of _observed_ events (i.e. signal+background)
// back = expected number of background events
// backe = uncertainty on the expected number of background events
double pbi(double obs, double bkg, double bkge)
{
  assert(obs>=0);
  assert(bkg>0);
  double P_BI;
  if(obs>0) {
    double tau = bkg/(bkge*bkge);
    double n_off = tau*bkg;
    P_BI = TMath::BetaIncomplete(1./(1.+tau), obs, n_off+1);
  } else {
    P_BI = 1.0;
  }

  // require that the z-score cannot exceed the number of events
  double Z_BI = pValueToZScore(P_BI);
  if(Z_BI>=obs) Z_BI=obs;
  return zScoreToPValue(Z_BI);
}

double zbi(double obs, double bkg, double bkge)
{
  return pValueToZScore(pbi(obs, bkg, bkge));
}


// sig = expected number of signal events
// bkg = expected number of background events
// bkge = uncertainty on the expected number of background events
// ntoys = how many toys to throw
double expectedZbi(double sig, double bkg, double bkge, int ntoys=10000)
{
  // setup random stuff
  static TRandom3 *rand=0;
  if(!rand) rand=new TRandom3(31415);
  RandomLognormal randlogn(bkg, bkge);
  
  TH1D* hPBIs=new TH1D("hPBIs","P_BIs",102,-0.01,1.01);
  std::vector<double> pbis;
  double average=0;
  for(int toy=0; toy<ntoys; toy++) {
    double sig_obs=rand->PoissonD(sig);		// Tote changed "Poisson" to "PoissonD", which returns a double (see TRandom ROOT docs).
    double thing = randlogn.getRandom();
    double bkg_obs=rand->PoissonD(thing);
    double nobs=sig_obs+bkg_obs;
    if (nobs < 0) {
    	cout << sig << "   " << sig_obs << "   " << bkg << "   " << bkg_obs << "   " << nobs << "   " << thing << endl;
    }
    double pbival=pbi(nobs, bkg, bkge);
    pbis.push_back(pbival);
    hPBIs->Fill(pbival);
    average+=pbival;
  }

  hPBIs->Draw();

  average/=ntoys;
  double median;
  std::pair<double, double> onesigma;
  std::pair<double, double> twosigma;
  
  getQuantiles(pbis, median, onesigma, twosigma);
  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << "observed pbi: " << pbi(sig+bkg, bkg, bkge) << std::endl;
  std::cout << "median expected pbi: " << median << std::endl;
  std::cout << "average expected pbi: " << average << std::endl;
  std::cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << std::endl;
  std::cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << "observed zbi: " << pValueToZScore(pbi(sig+bkg, bkg, bkge)) << std::endl;
  std::cout << "median expected zbi: " << pValueToZScore(median) << std::endl;
  std::cout << "average expected zbi: " << pValueToZScore(average) << std::endl;
  std::cout << "+/-1 sigma band: [ " << pValueToZScore(onesigma.second) << " , " << pValueToZScore(onesigma.first) << " ] " << std::endl;
  std::cout << "+/-2 sigma band: [ " << pValueToZScore(twosigma.second) << " , " << pValueToZScore(twosigma.first) << " ] " << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;

  return pValueToZScore(median);
}


