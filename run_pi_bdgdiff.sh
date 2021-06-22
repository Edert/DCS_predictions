
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL

files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}"
  
  NAME1=$(basename $4 -rep1_mm.bam)
  NAME2=$(basename $7 -rep1_mm.bam)
  LOG="log/$TOOL/$SET.log"
  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  #export PYTHONPATH=$PYTHONPATH:/apps/MACS2-2.1.0.20140616/lib/python2.7/site-packages/
  #active virtual environment for MACS2
  source /proj/chipseq_norm_diffbind_062017/analysis/03_db_analysis/python-virtual-environments/macs2_env/bin/activate

  STARTTIME=`date +%s.%N`


  #sharp 200
  #merge both replicates and use -B to get bedgraph files
  macs2 callpeak -B -t $4 $5 -c $6 -f BAM -n $NAME1 --nomodel --extsize 200 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG
  macs2 callpeak -B -t $7 $8 -c $9 -f BAM -n $NAME2 --nomodel --extsize 200 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG
  
  #egrep "tags after filtering in treatment|tags after filtering in control" cond2_peaks.xls
  D1=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME1"_peaks.xls" | cut -d' ' -f7)
  D2=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME2"_peaks.xls" | cut -d' ' -f7)
  
  PREPDONE=`date +%s.%N`
  
  #bgdiff 
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $D1 --d2 $D2 -C 1 --o-prefix diff_c1_vs_c2_1 --outdir results/$TOOL/$SET 2>> $LOG

  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" >  results/$TOOL/$SET/time.txt
  echo "1 $TIMEDIFF" >>  results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "1 $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
      
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating

  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_1.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_1.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_1.bed"

  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_1_c1*.bed results/$TOOL/$SET/mem.txt


  #bgdiff with scaling from spike-in 
  SF1=$(echo "print int(((${10}+${11})/2)*1000000)" | python2) #mean spike-in scaling-factor of the replicates
  SF2=$(echo "print int(((${12}+${13})/2)*1000000)" | python2)
   
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $SF1 --d2 $SF2 -C 1 --o-prefix diff_c1_vs_c2_si1 --outdir results/$TOOL/$SET 2>> $LOG
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "1si $TIMEDIFF" >>  results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "1si $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
  
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #cat results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_*.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"1"\t"$5}' > results/$TOOL/$SET/$OUT_NAME"_si1.csv"
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_1si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_1si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_1si.bed"
 
  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_si1_c1*.bed results/$TOOL/$SET/mem.txt
  rm results/$TOOL/$SET/$NAME1*peaks.narrowPeak results/$TOOL/$SET/$NAME1*summits.bed
  rm results/$TOOL/$SET/$NAME2*peaks.narrowPeak results/$TOOL/$SET/$NAME2*summits.bed
  rm results/$TOOL/$SET/$NAME1*xls results/$TOOL/$SET/$NAME2*xls
  rm results/$TOOL/$SET/$NAME1*control_lambda.bdg results/$TOOL/$SET/$NAME1*treat_pileup.bdg
  rm results/$TOOL/$SET/$NAME2*control_lambda.bdg results/$TOOL/$SET/$NAME2*treat_pileup.bdg
  
 
  #sharp model
  #merge both replicates and use -B to get bedgraph files
  macs2 callpeak -B -t $4 $5 -c $6 -f BAM -n $NAME1 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG
  macs2 callpeak -B -t $7 $8 -c $9 -f BAM -n $NAME2 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG
  
  #egrep "tags after filtering in treatment|tags after filtering in control" cond2_peaks.xls
  D1=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME1"_peaks.xls" | cut -d' ' -f7)
  D2=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME2"_peaks.xls" | cut -d' ' -f7)
  
  STARTTIME=`date +%s.%N`
  
  #bgdiff 
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $D1 --d2 $D2 -C 1 --o-prefix diff_c1_vs_c2_1 --outdir results/$TOOL/$SET 2>> $LOG
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >>  results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "2 $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
  
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_2.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_2.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_2.bed"

  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_1_c1*.bed results/$TOOL/$SET/mem.txt


  #bgdiff with scaling from spike-in 
  SF1=$(echo "print int(((${10}+${11})/2)*1000000)" | python2) #mean spike-in scaling-factor of the replicates
  SF2=$(echo "print int(((${12}+${13})/2)*1000000)" | python2)
   
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $SF1 --d2 $SF2 -C 1 --o-prefix diff_c1_vs_c2_si1 --outdir results/$TOOL/$SET 2>> $LOG
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2si $TIMEDIFF" >> results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "2si $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
  
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #cat results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_*.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"1"\t"$5}' > results/$TOOL/$SET/$OUT_NAME"_si1.csv"
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_2si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_2si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_2si.bed"
 
  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_si1_c1*.bed results/$TOOL/$SET/*model.r results/$TOOL/$SET/mem.txt
  rm results/$TOOL/$SET/$NAME1*peaks.narrowPeak results/$TOOL/$SET/$NAME1*summits.bed
  rm results/$TOOL/$SET/$NAME2*peaks.narrowPeak results/$TOOL/$SET/$NAME2*summits.bed
  rm results/$TOOL/$SET/$NAME1*xls results/$TOOL/$SET/$NAME2*xls
  rm results/$TOOL/$SET/$NAME1*control_lambda.bdg results/$TOOL/$SET/$NAME1*treat_pileup.bdg
  rm results/$TOOL/$SET/$NAME2*control_lambda.bdg results/$TOOL/$SET/$NAME2*treat_pileup.bdg
  
  
  #broad 200
  macs2 callpeak --broad -B -t $4 $5 -c $6 -f BAM -n $NAME1 --nomodel --extsize 200 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG
  macs2 callpeak --broad -B -t $7 $8 -c $9 -f BAM -n $NAME2 --nomodel --extsize 200 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG

  #egrep "tags after filtering in treatment|tags after filtering in control" cond2_peaks.xls
  D1=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME1"_peaks.xls" | cut -d' ' -f7)
  D2=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME2"_peaks.xls" | cut -d' ' -f7)

  STARTTIME=`date +%s.%N`
  
  #bgdiff broad
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $D1 --d2 $D2 -C 1 --o-prefix diff_c1_vs_c2_1 --outdir results/$TOOL/$SET 2>> $LOG
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "3 $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
  
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_3.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_3.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_3.bed"

  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_1_c1*.bed results/$TOOL/$SET/mem.txt


  #bgdiff with scaling from spike-in  broad
  SF1=$(echo "print int(((${10}+${11})/2)*1000000)" | python2) #mean spike-in scaling-factor of the replicates
  SF2=$(echo "print int(((${12}+${13})/2)*1000000)" | python2) 
  
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $SF1 --d2 $SF2 -C 1 --o-prefix diff_c1_vs_c2_si1 --outdir results/$TOOL/$SET 2>> $LOG
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3si $TIMEDIFF" >>  results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "3si $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
  
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #cat results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_*.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"1"\t"$5}' > results/$TOOL/$SET/$OUT_NAME"_si1.csv"
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_3si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_3si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_3si.bed"
  
  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_si1_c1*.bed results/$TOOL/$SET/mem.txt
  rm results/$TOOL/$SET/$NAME1*peaks.broadPeak results/$TOOL/$SET/$NAME1*gappedPeak
  rm results/$TOOL/$SET/$NAME2*peaks.broadPeak results/$TOOL/$SET/$NAME2*gappedPeak
  rm results/$TOOL/$SET/$NAME1*xls results/$TOOL/$SET/$NAME2*xls
  rm results/$TOOL/$SET/$NAME1*control_lambda.bdg results/$TOOL/$SET/$NAME1*treat_pileup.bdg
  rm results/$TOOL/$SET/$NAME2*control_lambda.bdg results/$TOOL/$SET/$NAME2*treat_pileup.bdg
  
  
  
  #broad model
  macs2 callpeak --broad -B -t $4 $5 -c $6 -f BAM -n $NAME1 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG
  macs2 callpeak --broad -B -t $7 $8 -c $9 -f BAM -n $NAME2 --outdir results/$TOOL/$SET -q 0.01 2>> $LOG

  #egrep "tags after filtering in treatment|tags after filtering in control" cond2_peaks.xls
  D1=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME1"_peaks.xls" | cut -d' ' -f7)
  D2=$(egrep "tags after filtering in control" results/$TOOL/$SET/$NAME2"_peaks.xls" | cut -d' ' -f7)

  STARTTIME=`date +%s.%N`
  
  #bgdiff broad
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $D1 --d2 $D2 -C 1 --o-prefix diff_c1_vs_c2_1 --outdir results/$TOOL/$SET 2>> $LOG
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4 $TIMEDIFF" >>  results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "4 $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
  
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_4.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_4.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_4.bed"

  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_1_c1*.bed results/$TOOL/$SET/mem.txt


  #bgdiff with scaling from spike-in  broad
  SF1=$(echo "print int(((${10}+${11})/2)*1000000)" | python2) #mean spike-in scaling-factor of the replicates
  SF2=$(echo "print int(((${12}+${13})/2)*1000000)" | python2) 
  
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o results/$TOOL/$SET/mem.txt -f "%K %M" macs2 bdgdiff \
  --t1 results/$TOOL/$SET/$NAME1"_treat_pileup.bdg" --c1 results/$TOOL/$SET/$NAME1"_control_lambda.bdg" \
  --t2 results/$TOOL/$SET/$NAME2"_treat_pileup.bdg" --c2 results/$TOOL/$SET/$NAME2"_control_lambda.bdg" \
  --d1 $SF1 --d2 $SF2 -C 1 --o-prefix diff_c1_vs_c2_si1 --outdir results/$TOOL/$SET 2>> $LOG
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4si $TIMEDIFF" >>  results/$TOOL/$SET/time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' results/$TOOL/$SET/mem.txt )
  echo "4si $MEMUSAGE" >> results/$TOOL/$SET/memory.txt
  
  #save result, add 1 or 2 for up in sample 1 or sample 2 and save as one csv
  #cat results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_*.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"1"\t"$5}' > results/$TOOL/$SET/$OUT_NAME"_si1.csv"
  #add log2 fold-change of -1.3 and 1.3 to get over the threshold of 0.7 we were simulating
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond1.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1.3"}'  | sort -k1,1 -k2,2n >  results/$TOOL/$SET/$OUT_NAME"_4si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_cond2.bed | grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-1.3"}' | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_4si.bed"
  cut -f1,2,3,5 results/$TOOL/$SET/diff_c1_vs_c2_si1_c1.0_common.bed| grep -v "^track" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}'    | sort -k1,1 -k2,2n >> results/$TOOL/$SET/$OUT_NAME"_4si.bed"
  
  #clean up
  rm results/$TOOL/$SET/diff_c1_vs_c2_si1_c1*.bed results/$TOOL/$SET/*model.r results/$TOOL/$SET/mem.txt
  rm results/$TOOL/$SET/$NAME1*peaks.broadPeak results/$TOOL/$SET/$NAME1*gappedPeak
  rm results/$TOOL/$SET/$NAME2*peaks.broadPeak results/$TOOL/$SET/$NAME2*gappedPeak
  rm results/$TOOL/$SET/$NAME1*xls results/$TOOL/$SET/$NAME2*xls
  rm results/$TOOL/$SET/$NAME1*control_lambda.bdg results/$TOOL/$SET/$NAME1*treat_pileup.bdg
  rm results/$TOOL/$SET/$NAME2*control_lambda.bdg results/$TOOL/$SET/$NAME2*treat_pileup.bdg
  
  
  #leave virtual environment
  deactivate
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
