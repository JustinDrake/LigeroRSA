# serializing
outSerializing="raw/3.a.Serializing.csv"
rm -f $outSerializing
cat data/coordinator.log | grep 'Serializing Data' | awk -F',' '{ print $2","$7 }' | awk '{ print $3 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print('\n'.join([i.split(',')[0] + ',' + str(float(i.split(',')[1].split(':')[0]) * 60 + float(i.split(',')[1].split(':')[1])) for i in x]))" >> $outSerializing
outDeseril="raw/3.a.Deserializing.csv"
rm -f $outDeseril
cat data/coordinator.log | grep '3\.a\.' | awk -F',' '{ print $7 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print('\n'.join([str(float(i.split(':')[0]) * 60 + float(i.split(':')[1])) for i in x]))" >> $outDeseril
