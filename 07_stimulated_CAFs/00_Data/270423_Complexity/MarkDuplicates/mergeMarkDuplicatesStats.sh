isfirst=1
for file in metrics_*.bam.txt
do
	file_clean=${file#metrics_sorted_}

	while [ $isfirst == 1 ]
	do
		awk -v F=$file_clean 'FNR==7{print "FILE\t"$0}; FNR==8{print F"\t"$0}' $file
		isfirst=0
	done
	awk -v F=$file_clean 'FNR==8{print F"\t"$0}' "$file"
	
done
