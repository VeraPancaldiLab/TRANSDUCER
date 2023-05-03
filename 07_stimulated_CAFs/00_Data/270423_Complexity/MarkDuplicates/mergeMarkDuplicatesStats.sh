isfirst=1
for file in metrics_*.bam.txt
do
	while [ $isfirst == 1 ]
	do
		awk -v F=$file 'FNR==7{print "File\t"$0}; FNR==8{print F"\t"$0}' $file
		isfirst=0
	done
	awk -v F=$file 'FNR==8{print F"\t"$0}' "$file"
	
done

