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

#./RUNMiniBoostedAnalyzer.py -m 100 -s TT -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s WJetsToQQ -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s ZJetsToQQ -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s JetHT_Run2016B -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s JetHT_Run2016C -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s JetHT_Run2016D -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s JetHT_Run2016E -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s JetHT_Run2016F -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s JetHT_Run2016G -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s JetHT_Run2016H -v v08 -b
#./RUNMiniBoostedAnalyzer.py -m 100 -s Dibosons -v v08 -b
##./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s QCDPtAll -v v06 -q Pt
##./RUNMiniBoostedAnalyzer.py -m 100 -p single  -s QCDPtAll -v v06 -q HT


if [[ $1 == *312* ]]; then
#	#massList="200 220 240 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100 1200 1300 1500"
	massList="80 100 120 140 160 180 200 220 240 300 350"
else
	massList="200 220 240 280 300 350 450 500 550 600 650 700 750 800 850 950 1000 1100 1200 1300 1500"
fi
for mass in $massList
do
	echo "Running ${1} mass " $mass
#	#./RUNOptimization.py -p calcROC -b Resolved -v v08 -m $mass -d $1 -l 36000 -q Pt
#	#./RUNOptimization.py -p plotROC -b Resolved -v v08 -m $mass -d $1 -l 36000 -q Pt 
#	#./RUNMiniResolvedAnalyzer.py -v v08 -s RPV -m $mass -d $1
	./RUNMiniBoostedAnalyzer.py -v v08 -s RPV -m $mass -d $1 -b
#	#./RUNFitter.py -p Bias -v v02p1 -m $mass -l 2646 -C delta -t
done
