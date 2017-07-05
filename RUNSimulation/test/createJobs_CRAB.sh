#!/bin/bash 

###########################################################################
###
###  Simple bash script to create python files for simulation.
###
###  Alejandro Gomez Espinosa
###  gomez@physics.rutgers.edu
###
###  How to run: 
###  source createJobs_CRAB.sh
###  (If it is an executable (chmod +x createJobs.sh):
###  ./createJobs_CRAB.sh 
###
###  You must change the PARAMETERS according with your needs. 
###  Initially is the only part that you should modify.
###
###########################################################################

######################################
### PARAMETERS
#####################################

totalNumberEvents=${1}
stop1=${2}		## You can use this parameters later to make everything simpler. 
coupling=${3}	 	##### UDD312 is for stop to jj, UDD323 for stop to bj

Name=RPVStopStopToJets_${coupling}_M-${stop1}_TuneCUETP8M1_13TeV-madgraph-pythia8      
#Name=AtoZ


#####################################################
#### Here is where the code starts.. 
#### Initially you shouldn't modify this part
#####################################################
empty=""
shortName=${Name/_TuneCUETP8M1_13TeV-madgraph-pythia8/$dummy}
echo " Creating directories..."
Main_Dir=$PWD 
####### Working directory
Working_Dir=${Main_Dir}/${shortName}	
if [ -d $Working_Dir ]; then
	rm -rf $Working_Dir
	mkdir -p $Working_Dir
else
	mkdir -p $Working_Dir
fi

cd $Working_Dir/

#user=`echo $USER`
#sed -i 's/algomez/'"${user}"'/' *

##############################################
##### Create the python file for Ntuples
##############################################

echo " Creating python file for RAWSIM (different PU scenarios).. "
step1PythonFile="step1_${shortName}_RAWSIM.py"
cp ${Main_Dir}/step1_RAW_v8020.py ${step1PythonFile}
cp ${Main_Dir}/MinBias_TuneCUETP8M1_13TeV-pythia8_cfi.py  . 
sed -i 's/outputFile/'"${shortName}"'_RAWSIM/' ${step1PythonFile}

step2PythonFile="step2_${shortName}_AODSIM.py"
cp ${Main_Dir}/step2_AODSIM_v8020.py ${step2PythonFile}
sed -i 's/outputFile/'"${shortName}"'_AODSIM/' ${step2PythonFile}

step3PythonFile="step3_${shortName}_MiniAOD.py"
cp ${Main_Dir}/step3_MiniAOD_v8020.py ${step3PythonFile}
sed -i 's/outputFile/'"${shortName}"'_MINIAODSIM/' ${step3PythonFile}

step23PythonFile="step23_${shortName}_AOD_MiniAOD.py"
cp ${Main_Dir}/step23_AOD_MINIAOD_v8020.py ${step23PythonFile}
sed -i 's/outputFile/'"${shortName}"'_MINIAODSIM/' ${step23PythonFile}

########################################################
######### Small file with the commands for condor
########################################################
echo " Creating crab files .... "

crabFileStep1=crab3_${shortName}_RAWSIM_step1.py
crabFileStep2=crab3_${shortName}_AODSIM_step2.py
crabFileStep3=crab3_${shortName}_MiniAOD_step3.py
crabFileStep23=crab3_${shortName}_AOD_MINIAODSIM_step23.py

cp ${Main_Dir}/crab3_RAW.py  ${crabFileStep1}
sed -i 's/test/'"${step1PythonFile}"'/' ${crabFileStep1}
sed -i 's/NAME/'"${shortName}"'_RAWSIM/' ${crabFileStep1}
sed -i 's/PROC/Moriond17_RAWSIM/' ${crabFileStep1}
sed -i 's/None/ADD_YOUR_DATASET_HERE/' ${crabFileStep1}

cp ${Main_Dir}/crab3_RAW.py  ${crabFileStep2}
sed -i 's/test/'"${step2PythonFile}"'/' ${crabFileStep2}
sed -i 's/NAME/'"${shortName}"'_AODSIM/' ${crabFileStep2}
sed -i 's/PROC/Moriond17_AODSIM/' ${crabFileStep2}
sed -i 's/None/ADD_YOUR_DATASET_HERE/' ${crabFileStep2}

cp ${Main_Dir}/crab3_RAW.py  ${crabFileStep3}
sed -i 's/test/'"${step3PythonFile}"'/' ${crabFileStep3}
sed -i 's/NAME/'"${shortName}"'_MINIAODSIM/' ${crabFileStep3}
sed -i 's/PROC/Moriond17_MINIAODSIM/' ${crabFileStep3}
sed -i 's/None/ADD_YOUR_DATASET_HERE/' ${crabFileStep3}

cp ${Main_Dir}/crab3_RAW.py  ${crabFileStep23}
sed -i 's/test/'"${step23PythonFile}"'/' ${crabFileStep23}
sed -i 's/NAME/'"${shortName}"'_MINIAODSIM/' ${crabFileStep23}
sed -i 's/PROC/Moriond17_MINIAODSIM/' ${crabFileStep23}
sed -i 's/None/ADD_YOUR_DATASET_HERE/' ${crabFileStep23}


#################################
##### To make it run
#################################
echo 'To make it run: 
First load the libraries (only once per session):
source /cvmfs/cms.cern.ch/crab3/crab.sh

Create and submit your jobs (Example for step0):
cd '${shortName}'
crab submit '${crabFileStep0}' 

To check the status:
crab status '${shortName}'

Once you have the dataset from step0 or step1, for example:
sed -i "s/ADD_YOUR_DATASET_HERE/\/RPVSt100tojj_13TeV_pythia8_GENSIM\/algomez-RPVSt100tojj_13TeV_pythia8_GENSIM-62459d50bdc5c4568f334137235e3bfc\/USER/g" crab3*RAWSIM*
' >> README

echo 'This script creates CRAB3 config files by default. Crab3 instructions in the README file'
echo 'Have a nice day :D '

