if [ "$1" == "all" ]
then
	if [ $2 == "Resolved" ] || [ $2 == "Bias" ]
	then
		masses="300 350 400 450 500 550 600 650"
	else
		#masses="80 90 100 110 120 130 140 150 170 180 190 210 220 230 240 300"
		masses="80 100 120 140 170 180 190 230 240 "
	fi

else
	masses="$1"
fi

for mass in $masses
do
	echo "======= Running datacard_RPVStopStopToJets_UDD312_M-${mass}"

	if [ $2 == "Resolved" ]
	then
		combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_Resolved_${3}_v02p1.txt -n UDD312RPVSt_M-${mass}_Resolved_${3}_v02p1
		#combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_Resolved_${3}_v02.txt -n UDD312RPVSt_M-${mass}_Resolved_${3}_v02

	elif [ $2 == "Bias" ]
	then
		combine Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_Resolved_delta_BiasTest_v02p1.txt -M GenerateOnly --setPhysicsModelParameters pdf_index=0 --toysFrequentist -t 10000 --expectSignal 1 --saveToys -n UDD312RPVSt_M-${mass}_Resolved_delta_v02p1_Index0_10k --freezeNuisances pdf_index
		for ind in 1 2 3
		do
			echo "======= Running Index ${ind} for mass ${mass}"
			combine Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_Resolved_delta_BiasTest_v02p1.txt -M MaxLikelihoodFit  --setPhysicsModelParameters pdf_index=${ind} --toysFile higgsCombineUDD312RPVSt_M-${mass}_Resolved_delta_v02p1_Index0_10k.GenerateOnly.mH120.123456.root  -t 10000 --rMin -10 --rMax 10 --freezeNuisances pdf_index
			mv mlfit.root mlfit_RPVStopStopToJets_UDD312_M-${mass}_Resolved_delta_BiasTest_v02p1_Index0ToIndex${ind}.root 
		done

	elif [ $2 == "fullCLs" ]
	then
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_withMC_altBkg_NOSys_v02p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n UDD312RPVSt_M-${mass}_Boosted_NOSys_v02
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_withMC_altBkg_NOSys_v02p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n UDD312RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.025
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_withMC_altBkg_NOSys_v02p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n UDD312RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.16
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_withMC_altBkg_NOSys_v02p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n UDD312RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.5
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_withMC_altBkg_NOSys_v02p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n UDD312RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.84
		combine -M HybridNew --testStat=LHC --frequentist Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_withMC_altBkg_NOSys_v02p1_bins.txt -T 2000 -H ProfileLikelihood --fork 4 -n UDD312RPVSt_M-${mass}_Boosted_NOSys_v02 --expectedFromGrid 0.975
		hadd higgsCombineUDD312RPVSt_M-${mass}_Boosted_NOSys_v02.HybridNewAll.mH120.root higgsCombineUDD312RPVSt_M-${mass}_Boosted_NOSys_v02.HybridNew*root 

	elif [ $2 == "pseudo" ]
	then
		for ((i=0;i<=$2;i++)); do
			echo "+++++++ Running pseudoExperiment "${i}
			combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_altBkg_signalInjectionTest${i}_Bin5_v05p3_bins.txt -n UDD312RPVSt_M-${mass}_altBkg_signalInjectionTest${i}_Bin5_v05p3 
		done
	else
		combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_UDD312_M-${mass}_${2}_withMC_altBkg_NOSys_v02p1_bins.txt -n UDD312RPVSt_M-${mass}_${2}_Boosted_NOSys_v02

	fi
done
