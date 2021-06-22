
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
  
  bamToBed -i $4 > S1.bed
  bamToBed -i $5 >> S1.bed
  bamToBed -i $6 > IN1.bed
  bamToBed -i $7 > S2.bed
  bamToBed -i $8 >> S2.bed
  bamToBed -i $9 > IN2.bed
 
  #sharp
  WINDOW=50
  GAP=100 
  
  PREPDONE=`date +%s.%N`
  
  #run it
  /usr/bin/time -o mem.txt -f "%K %M" sicer_df -t S2.bed S1.bed -c IN2.bed IN1.bed -s mm10 -w $WINDOW -g $GAP -fdr 0.01 -fdr_df 1 >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  #move results
  #reformat for eval if p-value is exactly 1 for one direction use the FDR and log2 fold-change from other direction
  cat S2-and-S1-W$WINDOW-G$GAP-summary | grep -v "^#chrom" | awk '{if($9 == 1){$10=$13; $8=$11} print $1"\t"$2"\t"$3"\t"$10"\t"log($8)/log(2)}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"

  #clean up
  rm -f S2-W$WINDOW-G$GAP* S1-W$WINDOW-G$GAP* mem.txt
  rm -f S1-W$WINDOW-normalized.wig S2-W$WINDOW-normalized.wig
  rm -f S2-and-S1-W$WINDOW-G$GAP-summary S2-vs-S1-W$WINDOW-G$GAP-E1000-union.island
  
  
  #default
  WINDOW=100
  GAP=200 
  
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o mem.txt -f "%K %M" sicer_df -t S2.bed S1.bed -c IN2.bed IN1.bed -s mm10 -w $WINDOW -g $GAP -fdr 0.01 -fdr_df 1 >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  #move results
  #reformat for eval if p-value is exactly 1 for one direction use the FDR and log2 fold-change from other direction
  cat S2-and-S1-W$WINDOW-G$GAP-summary | grep -v "^#chrom" | awk '{if($9 == 1){$10=$13; $8=$11} print $1"\t"$2"\t"$3"\t"$10"\t"log($8)/log(2)}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"

  #clean up
  rm -f S2-W$WINDOW-G$GAP* S1-W$WINDOW-G$GAP* mem.txt
  rm -f S1-W$WINDOW-normalized.wig S2-W$WINDOW-normalized.wig
  rm -f S2-and-S1-W$WINDOW-G$GAP-summary S2-vs-S1-W$WINDOW-G$GAP-E1000-union.island


  #broader
  WINDOW=200
  GAP=400
  
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o mem.txt -f "%K %M" sicer_df -t S2.bed S1.bed -c IN2.bed IN1.bed -s mm10 -w $WINDOW -g $GAP -fdr 0.01 -fdr_df 1 >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  #move results
  #reformat for eval if p-value is exactly 1 for one direction use the FDR and log2 fold-change from other direction
  cat S2-and-S1-W$WINDOW-G$GAP-summary | grep -v "^#chrom" | awk '{if($9 == 1){$10=$13; $8=$11} print $1"\t"$2"\t"$3"\t"$10"\t"log($8)/log(2)}' | sort -k1,1 -k2,2n > $OUT_NAME"_3.bed"

  #clean up
  rm -f S2-W$WINDOW-G$GAP* S1-W$WINDOW-G$GAP* mem.txt
  rm -f S1-W$WINDOW-normalized.wig S2-W$WINDOW-normalized.wig
  rm -f S2-and-S1-W$WINDOW-G$GAP-summary S2-vs-S1-W$WINDOW-G$GAP-E1000-union.island

  rm S1.bed S2.bed IN1.bed IN2.bed
  
else
  echo "results/$TOOL/$SET/bed\" already exists exiting..."
fi
