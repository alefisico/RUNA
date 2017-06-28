if [ "$1" == "all" ]
then
	if [ $2 == "Resolved" ] || [ $2 == "Bias" ]
	then
		if [ $4 == "UDD312" ]
		then
			masses="200 220 240 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100"
		else
			masses="200 220 240 280 300 350 450 500 550 600 650 700 750 800 850 950 1000 1100"
		fi
	else
		#masses="80 90 100 110 120 130 140 150 170 180 190 210 220 230 240 300"
		masses="80 100 120 140 160 180 200 220 240 300 350"
	fi

else
	masses="$1"
fi
decay=$4
version=$5

for mass in $masses
do
	echo "======= Running datacard_RPVStopStopToJets_${decay}_M-${mass}"

	if [ $2 == "Resolved" ]
	then
		combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}.txt -n ${decay}RPVSt_M-${mass}_Resolved_${3}_${version}

	elif [ $2 == "Bias" ]
	then
		numTests=${6}
		combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M GenerateOnly --setPhysicsModelParameters pdf_index=0 --toysFrequentist -t ${numTests} --expectSignal 1 --saveToys -n ${decay}RPVSt_M-${mass}_Resolved_${3}_${version}_Index0_signal1_10k --freezeNuisances pdf_index
		combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M GenerateOnly --setPhysicsModelParameters pdf_index=0 --toysFrequentist -t ${numTests} --expectSignal 0 --saveToys -n ${decay}RPVSt_M-${mass}_Resolved_${3}_${version}_Index0_signal0_10k --freezeNuisances pdf_index
		for ind in 0 1 2 3
		do
			echo "======= Running Index ${ind} for mass ${mass}"
			combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M MaxLikelihoodFit  --setPhysicsModelParameters pdf_index=${ind} --toysFile higgsCombine${decay}RPVSt_M-${mass}_Resolved_${3}_${version}_Index0_signal1_10k.GenerateOnly.mH120.123456.root  -t ${numTests} --rMin -10 --rMax 10 --freezeNuisances pdf_index -n _RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_signal1_${version}_Index0ToIndex${ind}

			combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M MaxLikelihoodFit  --setPhysicsModelParameters pdf_index=${ind} --toysFile higgsCombine${decay}RPVSt_M-${mass}_Resolved_${3}_${version}_Index0_signal0_10k.GenerateOnly.mH120.123456.root  -t ${numTests} --rMin -10 --rMax 10 --freezeNuisances pdf_index  -n _RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_signal0_${version}_Index0ToIndex${ind}
		done

	elif [ $2 == "fullCLs" ]
	then
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_withMC_altBkg_NOSys_v08p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n ${decay}RPVSt_M-${mass}_Boosted_NOSys_v02
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_withMC_altBkg_NOSys_v08p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n ${decay}RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.025
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_withMC_altBkg_NOSys_v08p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n ${decay}RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.16
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_withMC_altBkg_NOSys_v08p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n ${decay}RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.5
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_withMC_altBkg_NOSys_v08p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n ${decay}RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.84
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_withMC_altBkg_NOSys_v08p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n ${decay}RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.975
		hadd higgsCombine${decay}RPVSt_M-${mass}_Boosted_NOSys_v02.HybridNewAll.mH120.root higgsCombine${decay}RPVSt_M-${mass}_Boosted_NOSys_v02.HybridNew*root 

	elif [ $2 == "pseudo" ]
	then
		for ((i=0;i<=$2;i++)); do
			echo "+++++++ Running pseudoExperiment "${i}
			combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_altBkg_signalInjectionTest${i}_Bin5_v05p3_bins.txt -n ${decay}RPVSt_M-${mass}_altBkg_signalInjectionTest${i}_Bin5_v05p3 
		done
	else
		combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_${2}_${3}_${version}_bins.txt -n ${decay}RPVSt_M-${mass}_${2}_Boosted_${3}_${version}

	fi
done
