mkdir -p raw

out="raw/2.a.PeakMemParty.RAW.csv"
outStats="raw/2.a.PeakMemParty.csv"

rm -f $out
rm -f $outStats 

for dataFile in $(ls data/party_*) 
do
peakMem=`cat $dataFile | grep 'Peak Memory' | awk '{ print $4}'`
echo "$dataFile,$peakMem" >> $out
done


maxVal=`cat $out | awk -F',' '{ print $2 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print(max(x))"`
echo "max,$maxVal" >> $outStats
minVal=`cat $out | awk -F',' '{ print $2 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print(min(x))"`
echo "min,$minVal" >> $outStats
avg=`cat $out | awk -F',' '{ print $2 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print(sum([int(i) for i in x]) / len(x))"`
echo "avg,$avg" >> $outStats
stdev=`cat $out | awk -F',' '{ print $2 }' | python3 -c "import sys; import statistics; x=sys.stdin.read().split('\n')[:-1]; print(statistics.stdev([int(i) for i in x]))"`
echo "std,$stdev" >> $outStats
