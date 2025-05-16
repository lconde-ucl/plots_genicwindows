#!/bin/bash


#----------------------------------------------------------------
#- NOTE : I need to do this in an interactive shell, otherwise
#-	some of the bigwigs are too big and the bigWigAverageOverBed
#-	step gets killed
#-
#-	qrsh -l tmem=64G,h_vmem=64G
#----------------------------------------------------------------

canonical="knownCanonical.txt"
species="hs"
chromsizes="hg38.chrom.sizes"

for file in bedgraphs_nopseudo_hs/*hs.inputnormalised.bg bedgraphs_hs/*hs.inputnormalised.bg; do

	echo $file

	dir=""
	if [[ "$file" == *"nopseudo"* ]]; then
	    dir="_nopseudo"
	fi

	echo $dir


	#---------------------------------------
	#- bedgraphs have windows collapsed, so 
	#- I need to separate them in 100bp windows for the plots
	#---------------------------------------
	output1="$(basename $file .bg).bg.bed"
	output1b="$(basename $file .bg).bg.names.bed"
	output2="$(basename $file .bg).bg.bw"
	output3="$(basename $file .bg).100bp.bed"
	output3b="$(basename $file .bg).100bp.bg"
	output="$(basename $file .merged.mapped.${species}.inputnormalised.bg).100bp.genic.bg"

	#- make windows (add a name to the windows or bigWigA verageOverBed complains)
	/share/apps/bedtools2-2.30.0/bin/bedtools makewindows -b $file -w 100 > bedgraphs${dir}_${species}/$output1
	awk -F'\t' '{print $0"\t"$1":"$2"-"$3}' bedgraphs${dir}_${species}/$output1 > bedgraphs${dir}_${species}/$output1b

	#- convert bedgraph to bigwig
	bedGraphToBigWig $file $chromsizes bedgraphs${dir}_${species}/$output2

	#- get average signal over bed
	echo "bigWigAverageOverBed bedgraphs${dir}_${species}/$output2 bedgraphs${dir}_${species}/$output1b bedgraphs${dir}_${species}/$output3b -bedOut=bedgraphs${dir}_${species}/$output3"
	bigWigAverageOverBed bedgraphs${dir}_${species}/$output2 bedgraphs${dir}_${species}/$output1b bedgraphs${dir}_${species}/$output3b -bedOut=bedgraphs${dir}_${species}/$output3
 
	#-----------------------------------------
	#- keep the ones that overlap genes
	#-----------------------------------------
	#- keep only those that overlap genes
	/share/apps/bedtools2-2.30.0/bin/bedtools intersect -a bedgraphs${dir}_${species}/$output3 -b $canonical -wa > bedgraphs${dir}_${species}/$output

	#- remove intermediate files
	rm -rf bedgraphs${dir}_${species}/$output1 bedgraphs${dir}_${species}/$output1b bedgraphs${dir}_${species}/$output2 bedgraphs${dir}_${species}/$output3 bedgraphs${dir}_${species}/$output3b

done


canonical="knownCanonical.Mm.txt"
species="mm"
chromsizes="mm10.chrom.sizes"

for file in bedgraphs_nopseudo_mm/*mm.inputnormalised.bg bedgraphs_mm/*mm.inputnormalised.bg; do

	echo $file

	dir=""
	if [[ "$file" == *"nopseudo"* ]]; then
	    dir="_nopseudo"
	fi

	echo $dir

	#---------------------------------------
	#- bedgraphs have windows collapsed, so 
	#- I need to separate them in 100bp windows for the plots
	#---------------------------------------
	output1="$(basename $file .bg).bg.bed"
	output1b="$(basename $file .bg).bg.names.bed"
	output2="$(basename $file .bg).bg.bw"
	output3="$(basename $file .bg).100bp.bed"
	output3b="$(basename $file .bg).100bp.bg"
	output="$(basename $file .merged.mapped.${species}.inputnormalised.bg).100bp.genic.bg"

	#- make windows (add a name to the windows or bigWigA verageOverBed complains)
	/share/apps/bedtools2-2.30.0/bin/bedtools makewindows -b $file -w 100 > bedgraphs${dir}_${species}/$output1
	awk -F'\t' '{print $0"\t"$1":"$2"-"$3}' bedgraphs${dir}_${species}/$output1 > bedgraphs${dir}_${species}/$output1b

	#- convert bedgraph to bigwig
	bedGraphToBigWig $file $chromsizes bedgraphs${dir}_${species}/$output2

	#- get average signal over bed
	echo "bigWigAverageOverBed bedgraphs${dir}_${species}/$output2 bedgraphs${dir}_${species}/$output1b bedgraphs${dir}_${species}/$output3b -bedOut=bedgraphs${dir}_${species}/$output3"
	bigWigAverageOverBed bedgraphs${dir}_${species}/$output2 bedgraphs${dir}_${species}/$output1b bedgraphs${dir}_${species}/$output3b -bedOut=bedgraphs${dir}_${species}/$output3
 
	#-----------------------------------------
	#- keep the ones that overlap genes
	#-----------------------------------------
	#- keep only those that overlap genes
	/share/apps/bedtools2-2.30.0/bin/bedtools intersect -a bedgraphs${dir}_${species}/$output3 -b $canonical -wa > bedgraphs${dir}_${species}/$output

	#- remove intermediate files
	rm -rf bedgraphs${dir}_${species}/$output1 bedgraphs${dir}_${species}/$output1b bedgraphs${dir}_${species}/$output2 bedgraphs${dir}_${species}/$output3 bedgraphs${dir}_${species}/$output3b

done

