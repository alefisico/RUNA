#! /bin/bash
export CUR_DIR=$PWD

export SCRAM_ARCH="slc6_amd64_gcc493"
export VO_CMS_SW_DIR="/cms/base/cmssoft"
export COIN_FULL_INDIRECT_RENDERING=1
source /cms/base/cmssoft/cmsset_default.sh

# Change to your CMSSW software version
# Shown for c shell
# Also change 'username' to your username
#export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148
cd /cms/gomez/Substructure/CMSSW_8_1_0/src/
eval `scramv1 runtime -sh`
cd /cms/gomez/Substructure/CMSSW_8_0_20/src/RUNA/RUNStatistics/test/
source runCombine.sh all fullCLs jet1Tau32_BinReso UDD312 v09p15
#source runCombine.sh all Bias delta UDD312 v09p1 1000
#source runCombine.sh all Bias delta_2CSVv2L UDD323 v09p1 1000
