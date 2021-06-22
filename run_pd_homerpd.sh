
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
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  
  for PCALLER in ../../../results_peaks/*; do
  #for PCALLER in "../../../results_peaks/sicer"; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
      
      cat $PMODE/s11_peaks.bed $PMODE/s12_peaks.bed $PMODE/s21_peaks.bed $PMODE/s22_peaks.bed | sort -k1,1 -k2,2n > peaks.bed
      bedtools merge -i peaks.bed | sort -k1,1 -k2,2n > m_peaks.bed
      
      STARTTIME=`date +%s.%N`
      
      #/apps/homer_4.11/bin/getDifferentialPeaksReplicates.pl -t S11/ S12/ -b S21/ S22/ -i IN1/ IN2/ -f 0.7 -q 1 -p m_peaks.bed > results.txt 2>> $LOG 
      /usr/bin/time -o mem.txt -f "%K %M" /apps/homer_4.11/bin/getDifferentialPeaksReplicates.pl -t S11/ S12/ -b S21/ S22/ -i IN1/ IN2/ -f 0.7 -q 1 -p m_peaks.bed > results.txt 2>> $LOG 
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT".bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT" $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT" $MEMUSAGE" >> memory.txt
      
      lines=$(wc -l results.txt | awk '{print $1}')
      if [ $lines -ge 2 ]; then
         #reformat for eval
         cat results.txt | awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$26"\t"$24}' | grep -v "^Chr" | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
  
      #clean up
      rm -rf m_peaks.bed peaks.bed mem.txt
    done
  done
  
  rm -rf results.txt S11 S12 S21 S22 IN1 IN2
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
