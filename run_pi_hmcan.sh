
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15"

  LOG="../../../log/$TOOL/$SET.log"
  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  #create files
  echo $4 > C1_files.txt
  echo $5 >> C1_files.txt
  echo $6 > C1_control.txt
  echo $7 > C2_files.txt
  echo $8 >> C2_files.txt
  echo $9 > C2_control.txt
  
  PREPDONE=`date +%s.%N`
  
  #run 
  /usr/bin/time -o mem.txt -f "%K %M" /apps/hmcan-diff/HMCan-diff --name hmcan_out --C1_label S1 --C2_label S2 --C1_ChIP C1_files.txt --C2_ChIP C2_files.txt --C1_Control C1_control.txt --C2_Control C2_control.txt --format BAM --genomePath /ssd/references/THOR/mm10/ --GCProfile /ssd/references/THOR/mm10/chr19.cnp --C1_minLength 160 --C1_medLength 200 --C1_maxLength 240 --C2_minLength 160 --C2_medLength 200 --C2_maxLength 240 --largeBin 50000 --blackListFile /ssd/references/THOR/mm10/chr19.blacklist.bed --fold_change 0.7 --pvalue 1.0   >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e hmcan_out_peaks.bed ]; then
     #reformat for eval
     cat hmcan_out_peaks.bed | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$8}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi
  
  #clean up
  rm -f hmcan_out_peaks.bed hmcan_out_regions.bed mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run 
  /usr/bin/time -o mem.txt -f "%K %M" /apps/hmcan-diff/HMCan-diff --name hmcan_out --C1_label S1 --C2_label S2 --C1_ChIP C1_files.txt --C2_ChIP C2_files.txt --C1_Control C1_control.txt --C2_Control C2_control.txt --format BAM --genomePath /ssd/references/THOR/mm10/ --GCProfile /ssd/references/THOR/mm10/chr19.cnp --C1_minLength 160 --C1_medLength 200 --C1_maxLength 240 --C2_minLength 160 --C2_medLength 200 --C2_maxLength 240 --largeBin 50000 --blackListFile /ssd/references/THOR/mm10/chr19.blacklist.bed --fold_change 0.7 --pvalue 1.0 --negativeBinomial >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e hmcan_out_peaks.bed ]; then
     #reformat for eval
     cat hmcan_out_peaks.bed | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$8}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi
  #clean up
  rm -f hmcan_out_peaks.bed hmcan_out_regions.bed mem.txt
  
  #with sacling factors
  SF11=$(echo "print 1/($10/1000000.)" | python)
  SF12=$(echo "print 1/($11/1000000.)" | python)
  SF21=$(echo "print 1/($12/1000000.)" | python)
  SF22=$(echo "print 1/($13/1000000.)" | python)
  
  STARTTIME=`date +%s.%N`
  
  #run 
  /usr/bin/time -o mem.txt -f "%K %M" /apps/hmcan-diff/HMCan-diff --name hmcan_out --C1_label S1 --C2_label S2 --C1_ChIP C1_files.txt --C2_ChIP C2_files.txt --C1_Control C1_control.txt --C2_Control C2_control.txt --format BAM --genomePath /ssd/references/THOR/mm10/ --GCProfile /ssd/references/THOR/mm10/chr19.cnp --C1_minLength 160 --C1_medLength 200 --C1_maxLength 240 --C2_minLength 160 --C2_medLength 200 --C2_maxLength 240 --largeBin 50000 --blackListFile /ssd/references/THOR/mm10/chr19.blacklist.bed --fold_change 0.7 --pvalue 1.0 --C1_spikeInReadCounts $SF11,$SF12 --C2_spikeInReadCounts $SF21,$SF22 >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  if [ -e hmcan_out_peaks.bed ]; then
     #reformat for eval
     cat hmcan_out_peaks.bed | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$8}' | sort -k1,1 -k2,2n > $OUT_NAME"_3.bed"
  else
     #create empty file
     touch $OUT_NAME"_3.bed"
  fi
  #clean up
  rm -f hmcan_out_peaks.bed hmcan_out_regions.bed mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run 
  /usr/bin/time -o mem.txt -f "%K %M" /apps/hmcan-diff/HMCan-diff --name hmcan_out --C1_label S1 --C2_label S2 --C1_ChIP C1_files.txt --C2_ChIP C2_files.txt --C1_Control C1_control.txt --C2_Control C2_control.txt --format BAM --genomePath /ssd/references/THOR/mm10/ --GCProfile /ssd/references/THOR/mm10/chr19.cnp --C1_minLength 160 --C1_medLength 200 --C1_maxLength 240 --C2_minLength 160 --C2_medLength 200 --C2_maxLength 240 --largeBin 50000 --blackListFile /ssd/references/THOR/mm10/chr19.blacklist.bed --fold_change 0.7 --pvalue 1.0 --negativeBinomial --C1_spikeInReadCounts $SF11,$SF12 --C2_spikeInReadCounts $SF21,$SF22 >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "4 $MEMUSAGE" >> memory.txt
  
  if [ -e hmcan_out_peaks.bed ]; then
     #reformat for eval
     cat hmcan_out_peaks.bed | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$8}' | sort -k1,1 -k2,2n > $OUT_NAME"_4.bed"
  else
     #create empty file
     touch $OUT_NAME"_4.bed"
  fi
  #clean up
  rm -f C1_files.txt C2_files.txt C1_control.txt C2_control.txt hmcan_out_peaks.bed hmcan_out_regions.bed mem.txt
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
