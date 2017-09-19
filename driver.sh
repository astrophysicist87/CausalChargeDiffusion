#! /usr/bin/env bash

echo '************************************************'
echo Started everything at `date`
echo '************************************************'
echo

#make clean && make

#set up results directory
resultsDir="results"
#if [ -d "$resultsDir" ]
#then
#	rm -rf $resultsDir
#fi
#mkdir $resultsDir

#get snapshots first (for pions; white and colored noise)
#for frac in $(seq 0.05 0.05 0.5)
#do
#	#white noise
#	./run_CCD 1 1 100 $resultsDir $frac &
#	#colored noise
#	./run_CCD 1 1 0.3333333 $resultsDir $frac &
#done

#wait

#get rest of snapshots
#for frac in $(seq 0.55 0.05 1)
#do
#	#white noise
#	./run_CCD 1 1 100 $resultsDir $frac &
#	#colored noise
#	./run_CCD 1 1 0.3333333 $resultsDir $frac &
#done

#wait

#note that vQ2=0.333 and vQ2=100(white) already done
for vQ2 in 0.1 0.25 0.5 0.75 1.0 2.0 5.0 10.0
do
	./run_CCD 1 1 $vQ2 $resultsDir &
done

wait

for vQ2 in 0.1 0.25 0.5 0.75 1.0 2.0 5.0 10.0
do
	./run_CCD 2 2 $vQ2 $resultsDir &
done

wait

for vQ2 in 0.1 0.25 0.5 0.75 1.0 2.0 5.0 10.0
do
	./run_CCD 3 3 $vQ2 $resultsDir &
done

wait

echo
echo '************************************************'
echo Finished everything at `date`
echo '************************************************'
