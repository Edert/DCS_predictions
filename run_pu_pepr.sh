
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15"

  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  PREPDONE=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" PePr -c $4,$5 -i $6 --chip2 $7,$8 --input2 $9 -f bam --diff -n $OUT_NAME --peaktype sharp 1>> $LOG 2>&1 #default sharp
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e $OUT_NAME"__PePr_chip1_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip1_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  fi
  if [ -e $OUT_NAME"__PePr_chip2_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip2_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t-"$7}' | sort -k1,1 -k2,2n >> $OUT_NAME"_1.bed"
  fi
  touch $OUT_NAME"_1.bed"
  
  #clean up
  rm -f *debug.log *parameters.txt $OUT_NAME"__PePr_chip1_peaks.bed" $OUT_NAME"__PePr_chip2_peaks.bed" mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" PePr -c $4,$5 -i $6 --chip2 $7,$8 --input2 $9 -f bam --diff -n $OUT_NAME --normalization inter-group --num-processors 1 >> $LOG 2>&1 #default broad
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e $OUT_NAME"__PePr_chip1_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip1_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  fi
  if [ -e $OUT_NAME"__PePr_chip2_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip2_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t-"$7}' | sort -k1,1 -k2,2n >> $OUT_NAME"_2.bed"
  fi
  touch $OUT_NAME"_2.bed"
  
  #clean up
  rm -f *debug.log *parameters.txt $OUT_NAME"__PePr_chip1_peaks.bed" $OUT_NAME"__PePr_chip2_peaks.bed" mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" PePr -c $4,$5 -i $6 --chip2 $7,$8 --input2 $9 -f bam --diff -n $OUT_NAME --normalization intra-group --num-processors 1  >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e $OUT_NAME"__PePr_chip1_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip1_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME"_3.bed"
  fi
  if [ -e $OUT_NAME"__PePr_chip2_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip2_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t-"$7}' | sort -k1,1 -k2,2n >> $OUT_NAME"_3.bed"
  fi
  touch $OUT_NAME"_3.bed"
  
  #clean up
  rm -f *debug.log *parameters.txt $OUT_NAME"__PePr_chip1_peaks.bed" $OUT_NAME"__PePr_chip2_peaks.bed" mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" PePr -c $4,$5 -i $6 --chip2 $7,$8 --input2 $9 -f bam --diff -n $OUT_NAME --normalization scale --num-processors 1  >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "4 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e $OUT_NAME"__PePr_chip1_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip1_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME"_4.bed"
  fi
  if [ -e $OUT_NAME"__PePr_chip2_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip2_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t-"$7}' | sort -k1,1 -k2,2n >> $OUT_NAME"_4.bed"
  fi
  touch $OUT_NAME"_4.bed"
  
  #clean up
  rm -f *debug.log *parameters.txt $OUT_NAME"__PePr_chip1_peaks.bed" $OUT_NAME"__PePr_chip2_peaks.bed" mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" PePr -c $4,$5 -i $6 --chip2 $7,$8 --input2 $9 -f bam --diff -n $OUT_NAME --normalization no --num-processors 1  >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "5 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "5 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e $OUT_NAME"__PePr_chip1_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip1_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME"_5.bed"
  fi
  if [ -e $OUT_NAME"__PePr_chip2_peaks.bed" ]; then
     cat $OUT_NAME"__PePr_chip2_peaks.bed" | awk '{print $1"\t"$2"\t"$3"\t"$9"\t-"$7}' | sort -k1,1 -k2,2n >> $OUT_NAME"_5.bed"
  fi
  touch $OUT_NAME"_5.bed"
  
  #clean up
  rm -f *debug.log *parameters.txt $OUT_NAME"__PePr_chip1_peaks.bed" $OUT_NAME"__PePr_chip2_peaks.bed" mem.txt
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
