
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13"

  NAME11=$(basename $4 _mm.bam)
  NAME12=$(basename $5 _mm.bam)
  NAME21=$(basename $7 _mm.bam)
  NAME22=$(basename $8 _mm.bam)
  
  SCRIPT="script.R"
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  bamToBed -i $4 | cut -f1,2,3,6 > S1.bed
  bamToBed -i $5 | cut -f1,2,3,6 >> S1.bed
  #bamToBed -i $6 > IN1.bed
  bamToBed -i $7 | cut -f1,2,3,6 > S2.bed
  bamToBed -i $8 | cut -f1,2,3,6 >> S2.bed
  #bamToBed -i $9 > IN2.bed
  
  #Copy the modified sh and r scripts to the directory where the bed files are stored.
  cp /apps/MAnorm/MAnorm.r .
  cp /apps/MAnorm/MAnorm.sh .
  
  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt

  for PCALLER in ../../../results_peaks/*; do
  #for PCALLER in "../../../results_peaks/sicer"; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"

      cat $PMODE/s11_peaks.bed $PMODE/s12_peaks.bed | sort -k1,1 -k2,2n > s1_peaks.bed
      cat $PMODE/s21_peaks.bed $PMODE/s22_peaks.bed | sort -k1,1 -k2,2n > s2_peaks.bed
      
      bedtools merge -i s1_peaks.bed | sort -k1,1 -k2,2n > s1m_peaks.bed
      bedtools merge -i s2_peaks.bed | sort -k1,1 -k2,2n > s2m_peaks.bed
      
      #shift 200
      STARTTIME=`date +%s.%N`
      
      #sh MAnorm.sh s1m_peaks.bed s2m_peaks.bed S1.bed S2.bed 200 200 >> $LOG 2>&1
      /usr/bin/time -o mem.txt -f "%K %M" sh MAnorm.sh s1m_peaks.bed s2m_peaks.bed S1.bed S2.bed 200 200 >> $LOG 2>&1
      
      #save log
      cat Rcommand.out >> $LOG
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_1.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_1 $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_1 $MEMUSAGE" >> memory.txt
      
      #reformat for eval if p-value is exactly 1 for one direction use the FDR and log2 fold-change from other direction
      if [ -e MAnorm_result_commonPeak_merged.xls ]; then
         #reformat for eval
         cat MAnorm_result_commonPeak_merged.xls | grep -v "^chr\sstart" | awk '{ print $1"\t"$2"\t"$3"\t"10^-$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
      
      #clean up
      rm -f MAnorm_result.xls MAnorm_result_commonPeak_merged.xls MAplot_after_rescaling.png MAplot_before_rescaling.png mem.txt
      
      
      #shift 400
      STARTTIME=`date +%s.%N`
      
      #sh MAnorm.sh s1m_peaks.bed s2m_peaks.bed S1.bed S2.bed 400 400 >> $LOG 2>&1
      /usr/bin/time -o mem.txt -f "%K %M" sh MAnorm.sh s1m_peaks.bed s2m_peaks.bed S1.bed S2.bed 400 400 >> $LOG 2>&1
      
      #save log
      cat Rcommand.out >> $LOG
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_2.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_2 $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_2 $MEMUSAGE" >> memory.txt
      
      #reformat for eval if p-value is exactly 1 for one direction use the FDR and log2 fold-change from other direction
      if [ -e MAnorm_result_commonPeak_merged.xls ]; then
         #reformat for eval
         cat MAnorm_result_commonPeak_merged.xls | grep -v "^chr\sstart" | awk '{ print $1"\t"$2"\t"$3"\t"10^-$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
      
      #clean up
      rm -f MAnorm_result.xls MAnorm_result_commonPeak_merged.xls MAplot_after_rescaling.png MAplot_before_rescaling.png mem.txt
      
      
      #shift 1000
      STARTTIME=`date +%s.%N`
      
      #sh MAnorm.sh s1m_peaks.bed s2m_peaks.bed S1.bed S2.bed 1000 1000 >> $LOG 2>&1
      /usr/bin/time -o mem.txt -f "%K %M" sh MAnorm.sh s1m_peaks.bed s2m_peaks.bed S1.bed S2.bed 1000 1000 >> $LOG 2>&1
      
      #save log
      cat Rcommand.out >> $LOG
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_3.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_3 $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_3 $MEMUSAGE" >> memory.txt
      
      #reformat for eval if p-value is exactly 1 for one direction use the FDR and log2 fold-change from other direction
      if [ -e MAnorm_result_commonPeak_merged.xls ]; then
         #reformat for eval
         cat MAnorm_result_commonPeak_merged.xls | grep -v "^chr\sstart" | awk '{ print $1"\t"$2"\t"$3"\t"10^-$9"\t"$7}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
      
      #clean up
      rm -f .RData Rcommand.out sample1_peaks.wig sample2_peaks.wig s1_peaks.bed s2_peaks.bed s1m_peaks.bed s2m_peaks.bed mem.txt
      rm -f MAnorm_result.xls MAnorm_result_commonPeak_merged.xls MAplot_after_rescaling.png MAplot_before_rescaling.png
    done
  done
  
  rm -f MAnorm.r MAnorm.sh S1.bed S2.bed 
      
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
