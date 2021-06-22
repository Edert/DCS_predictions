
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
  
  #creat tags
  /apps/homer_4.11/bin/makeTagDirectory S11 $4 2>> $LOG
  /apps/homer_4.11/bin/makeTagDirectory S12 $5 2>> $LOG
  /apps/homer_4.11/bin/makeTagDirectory S21 $7 2>> $LOG
  /apps/homer_4.11/bin/makeTagDirectory S22 $8 2>> $LOG
  /apps/homer_4.11/bin/makeTagDirectory IN1 $6 2>> $LOG
  /apps/homer_4.11/bin/makeTagDirectory IN2 $9 2>> $LOG

  PREPDONE=`date +%s.%N`
  
  #run  homer
  /usr/bin/time -o mem.txt -f "%K %M" /apps/homer_4.11/bin/getDifferentialPeaksReplicates.pl -t S11/ S12/ -b S21/ S22/ -i IN1/ IN2/ -f 0.7 -q 1 -style factor > results.txt 2>> $LOG 
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "fa $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "fa $MEMUSAGE" >> memory.txt
  
  lines=$(wc -l results.txt | awk '{print $1}')
  if [ $lines -ge 2 ]; then
     #reformat for eval
     cat results.txt | awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$26"\t"$24}' | grep -v "^Chr" | sort -k1,1 -k2,2n > $OUT_NAME"_fa.bed"
  else
     #create empty file
     touch $OUT_NAME"_fa.bed"
  fi
  
  rm -rf results.txt mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run  homer
  /usr/bin/time -o mem.txt -f "%K %M" /apps/homer_4.11/bin/getDifferentialPeaksReplicates.pl -t S11/ S12/ -b S21/ S22/ -i IN1/ IN2/ -f 0.7 -q 1 -style histone > results.txt 2>> $LOG 
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "hi $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "hi $MEMUSAGE" >> memory.txt
  
  lines=$(wc -l results.txt | awk '{print $1}')
  if [ $lines -ge 2 ]; then
     #reformat for eval
     cat results.txt | awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$26"\t"$24}' | grep -v "^Chr" | sort -k1,1 -k2,2n > $OUT_NAME"_hi.bed"
  else
     #create empty file
     touch $OUT_NAME"_hi.bed"
  fi
  
  rm -rf results.txt mem.txt
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run  homer
  /usr/bin/time -o mem.txt -f "%K %M" /apps/homer_4.11/bin/getDifferentialPeaksReplicates.pl -t S11/ S12/ -b S21/ S22/ -i IN1/ IN2/ -f 0.7 -q 1 -style tss > results.txt 2>> $LOG 
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "ts $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "ts $MEMUSAGE" >> memory.txt
  
  lines=$(wc -l results.txt | awk '{print $1}')
  if [ $lines -ge 2 ]; then
     #reformat for eval
     cat results.txt | awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$26"\t"$24}' | grep -v "^Chr" | sort -k1,1 -k2,2n > $OUT_NAME"_ts.bed"
  else
     #create empty file
     touch $OUT_NAME"_ts.bed"
  fi
  
  rm -rf results.txt mem.txt
  
  #clean up
  rm -rf results.txt S11 S12 S21 S22 IN1 IN2
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
