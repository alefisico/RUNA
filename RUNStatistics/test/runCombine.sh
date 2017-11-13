if [ "$1" == "all" ]
then
	if [ $2 == "Resolved" ] || [ $2 == "Bias" ]
	then
		if [ $2 == "Bias" ]
		then 
			masses="200 300 400 500 600 700 800 900 1000"
		else
			if [ $4 == "UDD312" ]
			then
				masses="200 220 240 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100 1200 1300"
			else
				masses="200 220 240 260 280 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000" # 1100 1200"
			fi
		fi
	elif [ $2 == "final" ] 
	then
		if [ $4 == "UDD312" ]
		then
			masses="80 100 120 140 160 180 200 220 240 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100 1200"
		else
			masses="80 100 120 140 160 180 200 220 240 260 280 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100 1200"
		fi
	else 
		#masses="80 100 120 140 160 180 200 2"
		if [ $4 == "UDD312" ]
		then
			masses="80 100 120 140 160 180 200 220 240 300 350 400 450 500 550"
		else
			masses="80 100 120 140 160 180 200 220 240 260 280 300 350"
		fi
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
		combine -M AsymptoticLimits Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}.txt -n _RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}

	elif [ $2 == "Bias" ]
	then
		numTests=${6}
		combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M GenerateOnly --setPhysicsModelParameters pdf_index=0 --toysFrequentist -t ${numTests} --expectSignal 1 --saveToys -n _RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}_Index0_signal1_${numTests} --freezeNuisances pdf_index
		#combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M GenerateOnly --setPhysicsModelParameters pdf_index=0 --toysFrequentist -t ${numTests} --expectSignal 0 --saveToys -n _RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}_Index0_signal0_${numTests} --freezeNuisances pdf_index
		for ind in 0 1 2 3
		do
			echo "======= Running Index ${ind} for mass ${mass}"
			combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M MaxLikelihoodFit  --setPhysicsModelParameters pdf_index=${ind} --toysFile higgsCombine_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}_Index0_signal1_${numTests}.GenerateOnly.mH120.123456.root  -t ${numTests} --rMin -10 --rMax 10 --freezeNuisances pdf_index -n _RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_signal1_${version}_Index0ToIndex${ind}

			#combine Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_${version}.txt -M MaxLikelihoodFit  --setPhysicsModelParameters pdf_index=${ind} --toysFile higgsCombine_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}_Index0_signal0_${numTests}.GenerateOnly.mH120.123456.root  -t ${numTests} --rMin -10 --rMax 10 --freezeNuisances pdf_index  -n _RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_BiasTest_signal0_${version}_Index0ToIndex${ind}
		done

	elif [ $2 == "fullCLs" ]
	then
		NAME="_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}"
		combine Datacards/datacard${NAME}.txt -M HybridNew --LHCmode LHC-limits -n ${NAME} 
		combine Datacards/datacard${NAME}.txt -M HybridNew --LHCmode LHC-limits -n ${NAME} --expectedFromGrid 0.025
		combine Datacards/datacard${NAME}.txt -M HybridNew --LHCmode LHC-limits -n ${NAME} --expectedFromGrid 0.16
		combine Datacards/datacard${NAME}.txt -M HybridNew --LHCmode LHC-limits -n ${NAME} --expectedFromGrid 0.50
		combine Datacards/datacard${NAME}.txt -M HybridNew --LHCmode LHC-limits -n ${NAME} --expectedFromGrid 0.84
		combine Datacards/datacard${NAME}.txt -M HybridNew --LHCmode LHC-limits -n ${NAME} --expectedFromGrid 0.975
		hadd higgsCombine${NAME}.HybridNewAll.mH120.root higgsCombine${NAME}.HybridNew*root 

	elif [ $2 == "pseudo" ]
	then
		for ((i=0;i<=$2;i++)); do
			echo "+++++++ Running pseudoExperiment "${i}
			combine -M Asymptotic Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_altBkg_signalInjectionTest${i}_Bin5_v05p3_bins.txt -n ${decay}RPVSt_M-${mass}_altBkg_signalInjectionTest${i}_Bin5_v05p3 
		done
	elif [ $2 == "Boosted" ]
	then
		sufix="_RPVStopStopToJets_${decay}_M-${mass}_${2}_${3}_${version}"
		nameDatacard="Datacards/datacard${sufix}_bins.txt"
		if [ -f "${nameDatacard}" ]
		then
			rm ${nameDatacard} 
		fi
		combineCards.py Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_pruned_${3}_${version}_bin*.txt > ${nameDatacard}
		combine -M AsymptoticLimits ${nameDatacard} -n ${sufix} 
	
	elif [ $2 == "final" ] 
	then
		if [ "${mass}" -le 360 -a "${mass}" -ge 190 ]
		then
			echo "combining boosted and resolved cards"
			combineCards.py -S Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_pruned_jet1Tau32_Bin5_v09p2_bins.txt Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}.txt > Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_final_${version}.txt
			sed 's/Datacards\/\//\//g' -i  Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_final_${version}.txt

		elif [ "${mass}" -le 190 ]
		then
			cp  Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_pruned_jet1Tau32_Bin5_v09p2_bins.txt Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_final_${version}.txt

		elif [ "${mass}" -ge 360 ]
		then
			cp Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_Resolved_${3}_${version}.txt Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_final_${version}.txt
		fi

		combine -M AsymptoticLimits Datacards/datacard_RPVStopStopToJets_${decay}_M-${mass}_final_${version}.txt -n _RPVStopStopToJets_${decay}_M-${mass}_final_${version}

	fi
done
