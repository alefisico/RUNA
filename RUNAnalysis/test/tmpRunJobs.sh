#./RUNOptimization.py -m 80 -p calcROC -g Puppi -q Pt
#./RUNOptimization.py -m 100 -p calcROC -g Puppi -q Pt
#./RUNOptimization.py -m 170 -p calcROC -g Puppi -q Pt
#./RUNOptimization.py -m 240 -p calcROC -g Puppi -q Pt
#./RUNOptimization.py -m 80 -p calcROC -q Pt
#./RUNOptimization.py -m 100 -p calcROC -q Pt
#./RUNOptimization.py -m 170 -p calcROC -q Pt
#./RUNOptimization.py -m 240 -p calcROC -q Pt
#./RUNOptimization.py -m 80 -p plotROC -g Puppi -q Pt
#./RUNOptimization.py -m 100 -p plotROC -g Puppi -q Pt
#./RUNOptimization.py -m 170 -p plotROC -g Puppi -q Pt
#./RUNOptimization.py -m 240 -p plotROC -g Puppi -q Pt
#./RUNOptimization.py -m 80 -p plotROC -q Pt
#./RUNOptimization.py -m 100 -p plotROC -q Pt
#./RUNOptimization.py -m 170 -p plotROC -q Pt
#./RUNOptimization.py -m 240 -p plotROC -q Pt

#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s RPV -v v06
#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s JetHT_Run2016C -v v06
#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s TTJets -v v06
#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s WJetsToQQ -v v06
#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s ZJetsToQQ -v v06
#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s Dibosons -v v06
#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s QCDPtAll -v v06 -q Pt
#./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s QCDPtAll -v v06 -q HT
#
#./RUNMiniResolvedAnalyzer.py -m 300 -s QCDHTAll -v v06 -q HT
#./RUNMiniResolvedAnalyzer.py -m 300 -s TT -v v06


if [[ $1 == *312* ]]; then
	massList="200 220 240 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100 1200 1300 1500"
else
	massList="200 220 240 280 300 350 450 500 550 600 650 700 750 800 850 950 1000 1100 1200 1300 1500"
fi
for mass in $massList
do
	echo "Running ${1} mass " $mass
	#./RUNOptimization.py -p calcROC -b Resolved -v v08 -m $mass -d $1 -l 36000 -q Pt
	#./RUNOptimization.py -p plotROC -b Resolved -v v08 -m $mass -d $1 -l 36000 -q Pt 
	./RUNMiniResolvedAnalyzer.py -v v08 -s RPV -m $mass -d $1
	#./RUNFitter.py -p Bias -v v02p1 -m $mass -l 2646 -C delta -t
done
