$isfirst = 1
for file in metrics_*.bam.txt
do
	while [ $isfirst == 1 ]
	do
		awk '/## METRICS CLASS/,/## HISTOGRAM/' "$file" | grep -v "##"
		$isfirst = 0
	done
	awk '/## METRICS CLASS/,/## HISTOGRAM/' "$file" | grep -v "##"
	
done

