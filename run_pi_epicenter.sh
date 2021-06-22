
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
  
  samtools merge s1.bam $4 $5
  samtools merge s2.bam $7 $8
  
  WIN=100
  
  PREPDONE=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" /apps/epicenter-1.7.0.8/EpiCenter -t 0 -w $WIN -i bam s1.bam s2.bam >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e "s1_s2-"$WIN".tscan" ]; then
     #reformat for eval
     cat "s1_s2-"$WIN".tscan" | grep -v "^chr\s" | awk '{print $1"\t"$2"\t"$2+$3"\t"$12"\t"$10}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi
  
  #clean up
  rm -rf *.ratio *.tscan mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" /apps/epicenter-1.7.0.8/EpiCenter -t 1 -w $WIN -i bam s1.bam s2.bam >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e "s1_s2-"$WIN".tscan" ]; then
     #reformat for eval
     cat "s1_s2-"$WIN".tscan" | grep -v "^chr\s" | awk '{print $1"\t"$2"\t"$2+$3"\t"$12"\t"$10}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi
  
  #clean up
  rm -rf *.ratio *.tscan mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" /apps/epicenter-1.7.0.8/EpiCenter -t 2 -w $WIN -i bam s1.bam s2.bam >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  if [ -e "s1_s2-"$WIN".fscan" ]; then
     #reformat for eval
     cat "s1_s2-"$WIN".fscan" | grep -v "^chr\s" | awk '{print $1"\t"$2"\t"$2+$3"\t"$12"\t"$10}' | sort -k1,1 -k2,2n > $OUT_NAME"_3.bed"
  else
     #create empty file
     touch $OUT_NAME"_3.bed"
  fi
  
  #clean up
  rm -rf *.ratio *.fscan mem.txt
  
  
  
  WIN=1000
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" /apps/epicenter-1.7.0.8/EpiCenter -t 0 -w $WIN -i bam s1.bam s2.bam >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "4 $MEMUSAGE" >> memory.txt
  
  if [ -e "s1_s2-"$WIN".tscan" ]; then
     #reformat for eval
     cat "s1_s2-"$WIN".tscan" | grep -v "^chr\s" | awk '{print $1"\t"$2"\t"$2+$3"\t"$12"\t"$10}' | sort -k1,1 -k2,2n > $OUT_NAME"_4.bed"
  else
     #create empty file
     touch $OUT_NAME"_4.bed"
  fi
  
  #clean up
  rm -rf *.ratio *.tscan mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" /apps/epicenter-1.7.0.8/EpiCenter -t 1 -w $WIN -i bam s1.bam s2.bam >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "5 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "5 $MEMUSAGE" >> memory.txt
  
  if [ -e "s1_s2-"$WIN".tscan" ]; then
     #reformat for eval
     cat "s1_s2-"$WIN".tscan" | grep -v "^chr\s" | awk '{print $1"\t"$2"\t"$2+$3"\t"$12"\t"$10}' | sort -k1,1 -k2,2n > $OUT_NAME"_5.bed"
  else
     #create empty file
     touch $OUT_NAME"_5.bed"
  fi
  
  #clean up
  rm -rf *.ratio *.tscan mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run
  /usr/bin/time -o mem.txt -f "%K %M" /apps/epicenter-1.7.0.8/EpiCenter -t 2 -w $WIN -i bam s1.bam s2.bam >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "6 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "6 $MEMUSAGE" >> memory.txt
  
  if [ -e "s1_s2-"$WIN".fscan" ]; then
     #reformat for eval
     cat "s1_s2-"$WIN".fscan" | grep -v "^chr\s" | awk '{print $1"\t"$2"\t"$2+$3"\t"$12"\t"$10}' | sort -k1,1 -k2,2n > $OUT_NAME"_6.bed"
  else
     #create empty file
     touch $OUT_NAME"_6.bed"
  fi
  
  #clean up
  rm -rf *.ratio *.fscan mem.txt
  rm -rf s1.bam s2.bam
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
