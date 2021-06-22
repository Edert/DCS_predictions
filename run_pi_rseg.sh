
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
  
  export LC_ALL=C
  bamToBed -i $4 > tmp.bed
  bamToBed -i $5 >> tmp.bed
  cat tmp.bed | sort -k1,1 -k2,2n -k3,3n -k6,6 > S1.bed
  bamToBed -i $7 > tmp.bed
  bamToBed -i $8 >> tmp.bed
  cat tmp.bed | sort -k1,1 -k2,2n -k3,3n -k6,6 > S2.bed
  rm tmp.bed
  
  PREPDONE=`date +%s.%N`

  #run 
  /usr/bin/time -o mem.txt -f "%K %M" /apps/rseg-0.4.9/bin/rseg-diff -c /ssd/references/THOR/mm10/chr19.chromsize.bed -out domains.bed -d /ssd/references/THOR/mm10/deadzones-mm10-chr19-k50_sub.bed -v -mode 3 S1.bed S2.bed >> $LOG 2>&1

  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e domains.bed ]; then
     #reformat for eval
     cat domains.bed | awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$5}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi

  #clean up
  rm -f S1.bed S2.bed domains.bed mem.txt
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
