#! /bin/bash

# If you have additional custom enhancements to your shell 
# environment, you may need to add them here

export CUR_DIR=$PWD

export SCRAM_ARCH="slc6_amd64_gcc493"
export VO_CMS_SW_DIR="/cms/base/cmssoft"
export COIN_FULL_INDIRECT_RENDERING=1
source /cms/base/cmssoft/cmsset_default.sh

# Change to your CMSSW software version
# Shown for c shell
# Also change 'username' to your username
cd /cms/gomez/Substructure/CMSSW_8_0_20/src/
eval `scramv1 runtime -sh`
date

# Switch to your working directory below
cd  $CUR_DIR
#for mass in 80 90 100 110 120 130 140 150 170 180 190 210 220 230 240 300 350;
#do
#	python RUNOptimization.py -p calcROC -m ${mass} -b True
#done
#python RUNMiniBoostedAnalyzer.py -b -m 100 -s QCDPtAll -q Pt -v v09
python RUNMiniResolvedAnalyzer.py -b -m 100 -s QCDPtAll -q Pt -v v09
source tmpRunJobs.sh UDD312
#python RUNMiniBoostedAnalyzer.py -b -m 100 -s QCDHTAll -q HT -v v08
