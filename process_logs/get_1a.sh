out="raw/1a.Overall.csv"
rm -f $out 
cat data/coordinator.log | grep "Overall" | awk -F',' '{ print $2","$7}' >> $out
