
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13"
  
  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  bamToBed -i $4 > S11.bed
  bamToBed -i $5 > S12.bed
  bamToBed -i $6 > IN1.bed
  bamToBed -i $7 > S21.bed
  bamToBed -i $8 > S22.bed
  bamToBed -i $9 > IN2.bed
  
  PREPDONE=`date +%s.%N`
  
  #default run
  /usr/bin/time -o mem.txt -f "%K %M" diffReps.pl --report report --treatment S11.bed S12.bed --control S21.bed S22.bed --btr IN1.bed --bco IN2.bed --chrlen /ssd/references/THOR/mm10/chr19.chrom.size >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e report ]; then
    #move results
    cat report | grep -v "^#" | grep -v "^Chrom" | awk '{print $1"\t"$2"\t"$3"\t"$14"\t"$12}' | sort -k1,1 -k2,2n >  $OUT_NAME"_1.bed"
  else
    #create empty file
    touch $OUT_NAME"_1.bed"
  fi
  
  #clean up
  rm -f report report.hotspot norm.txt mem.txt



  #with scaling factors
  echo "treatment $10 $11" > norm.txt
  echo "control $12 $13" >> norm.txt
  
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o mem.txt -f "%K %M" diffReps.pl --report report --treatment S11.bed S12.bed --control S21.bed S22.bed --btr IN1.bed --bco IN2.bed --chrlen /ssd/references/THOR/mm10/chr19.chrom.size --norm norm.txt >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e report ]; then
    #move results
    cat report | grep -v "^#" | grep -v "^Chrom" | awk '{print $1"\t"$2"\t"$3"\t"$14"\t"$12}' | sort -k1,1 -k2,2n >  $OUT_NAME"_2.bed"
  else
    #create empty file
    touch $OUT_NAME"_2.bed"
  fi
  
  #clean up
  rm -f S11.bed S12.bed S21.bed S22.bed IN1.bed IN2.bed report report.hotspot norm.txt mem.txt

else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
