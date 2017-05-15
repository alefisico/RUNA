export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148
cd /afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNAnalysis/test/

python RUNMiniResolvedAnalyzer.py -s QCDPtAll -q Pt -m 300 -p run -v v08
#python RUNMiniResolvedAnalyzer.py -s QCDHTAll -q HT -m 300 -p run -v v08

#python RUNMiniBoostedAnalyzer.py -s QCDPtAll -q Pt -m 100 -v v08
#python RUNMiniBoostedAnalyzer.py -s QCDHTAll -q HT -m 100 -v v08
#python RUNMiniBoostedAnalyzer.py -s TT -m 100 -v v08
#python RUNMiniBoostedAnalyzer.py -s WJetsToQQ -m 100 -v v08

#python RUNMiniBoostedAnalyzer.py -m 190 -s DATA -r high -g prunedPuppi
#python RUNBkgEstimation.py -m 110 -p DATA -r low -g pruned 
