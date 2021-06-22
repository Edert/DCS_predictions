
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
 
  /apps/chipnorm/chipnorminput S1.bed /ssd/references/THOR/mm10/chr19.chrom.size >> $LOG 2>&1
  /apps/chipnorm/chipnorminput S2.bed /ssd/references/THOR/mm10/chr19.chrom.size >> $LOG 2>&1
  /apps/chipnorm/chipnorminput IN1.bed /ssd/references/THOR/mm10/chr19.chrom.size >> $LOG 2>&1
  /apps/chipnorm/chipnorminput IN2.bed /ssd/references/THOR/mm10/chr19.chrom.size >> $LOG 2>&1
  
  cp /apps/chipnorm/ChIPnorm_to_send_mod.m .
  cp /ssd/references/THOR/mm10/chr19.chrom.size .
  
  #clean up
  rm -rf S1.wig S2.wig IN1.wig IN2.wig
  
  PREPDONE=`date +%s.%N`
  
  #default run win100
  /usr/bin/time -o mem.txt -f "%K %M" /apps/matlab/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -nojvm -r "ChIPnorm_to_send_mod('S1_matlabInput.txt', 'IN1_matlabInput.txt', 'S2_matlabInput.txt', 'IN2_matlabInput.txt', 1, 'chr19.chrom.size', 'b', 100, 't', 0.7);exit"  >> $LOG 2>&1
 
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e diff_reg_S1_S2_thresh_0.bed ]; then
    #move results
    cat diff_reg_S1_S2_thresh_0.bed | grep -v "non_differential" |  awk '{print $1"\t"$2"\t"$3"\t0\t1"}' | sort -k1,1 -k2,2n >  $OUT_NAME"_1.bed"
  else
    #create empty file
    touch $OUT_NAME"_1.bed"
  fi
  
  #clean up
  rm -f diff_reg_S1_S2_thresh_0.bed S1_S2_thresh_0.7.bed mem.txt
  
  
  #default run win1000
  STARTTIME=`date +%s.%N`
  
  /usr/bin/time -o mem.txt -f "%K %M" /apps/matlab/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -nojvm -r "ChIPnorm_to_send_mod('S1_matlabInput.txt', 'IN1_matlabInput.txt', 'S2_matlabInput.txt', 'IN2_matlabInput.txt', 1, 'chr19.chrom.size', 'b', 1000, 't', 0.7);exit"  >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e diff_reg_S1_S2_thresh_0.bed ]; then
    #move results
    cat diff_reg_S1_S2_thresh_0.bed | grep -v "non_differential" |  awk '{print $1"\t"$2"\t"$3"\t0\t1"}' | sort -k1,1 -k2,2n >  $OUT_NAME"_2.bed"
  else
    #create empty file
    touch $OUT_NAME"_2.bed"
  fi
  
  #clean up
  rm -f ChIPnorm_to_send_mod.m chr19.chrom.size diff_reg_S1_S2_thresh_0.bed S1_S2_thresh_0.7.bed mem.txt
  rm -f S1_matlabInput.txt S2_matlabInput.txt IN1_matlabInput.txt IN2_matlabInput.txt
  rm -rf S1.bed S2.bed IN1.bed IN2.bed ChIPnorm_to_send_mod.m chr19.chrom.size 
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
