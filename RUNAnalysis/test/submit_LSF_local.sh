export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148
cd /afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNAnalysis/test/
#cmsRun RUNFullAnalysis_cfg.py local=1 PROC=JetHT_Run2016B jecVersion=supportFiles/Spring16_25nsV8BCD
#cmsRun RUNFullAnalysis_cfg.py local=1 PROC=JetHT_Run2016C jecVersion=supportFiles/Spring16_25nsV8BCD
#cmsRun RUNFullAnalysis_cfg.py local=1 PROC=JetHT_Run2016D jecVersion=supportFiles/Spring16_25nsV8BCD
#cmsRun RUNFullAnalysis_cfg.py local=1 PROC=JetHT_Run2016E jecVersion=supportFiles/Spring16_25nsV8E
#cmsRun RUNFullAnalysis_cfg.py local=1 PROC=JetHT_Run2016F jecVersion=supportFiles/Spring16_25nsV8F
#cmsRun RUNFullAnalysis_cfg.py local=1 PROC=JetHT_Run2016G jecVersion=supportFiles/Spring16_25nsV8p2
#cmsRun RUNFullAnalysis_cfg.py local=1 PROC=JetHT_Run2016H jecVersion=supportFiles/Spring16_25nsV8p2

cmsRun RUNFullAnalysis_cfg.py local=1 PROC=QCD_HT500to700 jecVersion=supportFiles/Spring16_25nsV8BCD systematics=0 namePUFile=supportFiles/PileupData2016BCDEFGH_271036-283685_69200.root  

#cmsRun simpleMatching.py
#python RUNMiniAnalyzer.py -s QCDPt
#python RUNMiniBoostedAnalyzer.py -m 110 -s DATA -r low -g pruned 
#python RUNMiniBoostedAnalyzer.py -m 190 -s DATA -r high -g pruned 
#python RUNMiniBoostedAnalyzer.py -m 110 -s DATA -r low -g softDropPuppi
#python RUNMiniBoostedAnalyzer.py -m 190 -s DATA -r high -g softDropPuppi
#python RUNMiniBoostedAnalyzer.py -m 110 -s DATA -r low -g softDrop
#python RUNMiniBoostedAnalyzer.py -m 110 -s DATA -r low -g prunedPuppi
#python RUNMiniBoostedAnalyzer.py -m 190 -s DATA -r high -g softDrop
#python RUNMiniBoostedAnalyzer.py -m 190 -s DATA -r high -g prunedPuppi

#python RUNBkgEstimation.py -m 110 -p DATA -r low -g pruned 
#python RUNBkgEstimation.py -m 110 -p DATA -r low -g softDrop
#python RUNBkgEstimation.py -m 110 -p DATA -r low -g softDropPuppi
#python RUNBkgEstimation.py -m 110 -p DATA -r low -g prunedPuppi
#python RUNBkgEstimation.py -m 190 -p DATA -r high -g softDrop
#python RUNBkgEstimation.py -m 190 -p DATA -r high -g pruned 
#python RUNBkgEstimation.py -m 190 -p DATA -r high -g softDropPuppi
#python RUNBkgEstimation.py -m 190 -p DATA -r high -g prunedPuppi
