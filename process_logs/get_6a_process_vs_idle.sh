mkdir -p raw

out="raw/6.a.TimeToVerifyZkProof.csv"
rm -f $out

echo "filename,idle,compute" >> $out
for dataFile in $(ls data/distributed_*) 
do
idle=`cat $dataFile | grep '6.a. verification' | grep idle | tail -n +2  | awk -F',' '{ print $9 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print(sum([float(i.split(':')[0]) * 60 + float(i.split(':')[1]) for i in x]))"`
compute=`cat $dataFile | grep '6.a. verification' | grep compute | awk -F',' '{ print $9 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print(sum([float(i.split(':')[0]) * 60 + float(i.split(':')[1]) for i in x]))"`
echo "$dataFile,$idle,$compute" >> $out
done
