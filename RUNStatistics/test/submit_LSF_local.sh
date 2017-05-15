export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148
cd /afs/cern.ch/work/a/algomez/RPVStops/CMSSW_7_4_7/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/a/algomez/RPVStops/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/
#source runCombine.sh all fullCLs
#source runCombine.sh all Bias delta UDD312
source runCombine.sh 300 Bias delta UDD312
#source runCombine.sh all Bias delta_2CSVv2L UDD323
