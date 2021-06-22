
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
  
  bamToBed -i $4 > S1.bed
  bamToBed -i $5 >> S1.bed
  bamToBed -i $7 > S2.bed
  bamToBed -i $8 >> S2.bed
  
  cp /apps/ChIPDiff/config.txt config.txt
  
  PREPDONE=`date +%s.%N`
  
  #run 
  /usr/bin/time -o mem.txt -f "%K %M" /apps/ChIPDiff/ChIPDiff S1.bed S2.bed /ssd/references/THOR/mm10/chr19.chrom.size config.txt $OUT_NAME >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME".bin" ]; then
     #reformat for eval
     cat $OUT_NAME".bin" | awk '{lfold=(($5==0) ? "inf" : log($4/$5)/log(2)); p=(($6 > $8) ? $6 : $8); if($9!=0) print $1"\t"$3"\t"($3+1000)"\t"p"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi
  
  #clean up
  rm -f config.txt $OUT_NAME".bin" $OUT_NAME".hmm" $OUT_NAME".region" mem.txt


  #run with non default values 1
  echo "maxIterationNum	500" > config.txt
  echo "minP	0.95" >> config.txt
  echo "maxTrainingSeqNum	10000" >> config.txt
  echo "minFoldChange	2.0" >> config.txt
  echo "minRegionDist	100" >> config.txt
  
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o mem.txt -f "%K %M" /apps/ChIPDiff/ChIPDiff S1.bed S2.bed /ssd/references/THOR/mm10/chr19.chrom.size config.txt $OUT_NAME >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME".bin" ]; then
     #reformat for eval
     cat $OUT_NAME".bin" | awk '{lfold=(($5==0) ? "inf" : log($4/$5)/log(2)); p=(($6 > $8) ? $6 : $8); if($9!=0) print $1"\t"$3"\t"($3+1000)"\t"p"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi
  
  #clean up
  rm -f config.txt $OUT_NAME".bin" $OUT_NAME".hmm" $OUT_NAME".region" mem.txt
  
  
  #run with non default values 2
  echo "maxIterationNum	500" > config.txt
  echo "minP	0.95" >> config.txt
  echo "maxTrainingSeqNum	10000" >> config.txt
  echo "minFoldChange	0.7" >> config.txt
  echo "minRegionDist	100" >> config.txt
  
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o mem.txt -f "%K %M" /apps/ChIPDiff/ChIPDiff S1.bed S2.bed /ssd/references/THOR/mm10/chr19.chrom.size config.txt $OUT_NAME >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME".bin" ]; then
     #reformat for eval
     cat $OUT_NAME".bin" | awk '{lfold=(($5==0) ? "inf" : log($4/$5)/log(2)); p=(($6 > $8) ? $6 : $8); if($9!=0) print $1"\t"$3"\t"($3+1000)"\t"p"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_3.bed"
  else
     #create empty file
     touch $OUT_NAME"_3.bed"
  fi
  
  #clean up
  rm -f config.txt S1.bed S2.bed $OUT_NAME".bin" $OUT_NAME".hmm" $OUT_NAME".region" mem.txt

else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
