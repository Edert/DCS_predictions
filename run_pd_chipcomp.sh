
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
  
  bamToBed -i $4 | cut -f1,2,3 > S11.bed
  bamToBed -i $5 | cut -f1,2,3 > S12.bed
  bamToBed -i $6 | cut -f1,2,3 > IN1.bed
  bamToBed -i $7 | cut -f1,2,3 > S21.bed
  bamToBed -i $8 | cut -f1,2,3 > S22.bed
  bamToBed -i $9 | cut -f1,2,3 > IN2.bed
  
  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt

  for PCALLER in ../../../results_peaks/*; do
  #for PCALLER in "../../../results_peaks/sicer"; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
      
      cut -f1,2,3 $PMODE/s11_peaks.bed > s11_peaks.bed
      cut -f1,2,3 $PMODE/s12_peaks.bed > s12_peaks.bed
      cut -f1,2,3 $PMODE/s21_peaks.bed > s21_peaks.bed
      cut -f1,2,3 $PMODE/s22_peaks.bed > s22_peaks.bed
      
      # in R diffbind
      echo "library(ChIPComp)" > $SCRIPT
      echo "conf=data.frame(" >> $SCRIPT
      echo "SampleID=1:4," >> $SCRIPT
      echo "condition=c(\"s1\",\"s1\",\"s2\",\"s2\")," >> $SCRIPT
      echo "factor=c(\"test\",\"test\",\"test\",\"test\")," >> $SCRIPT
      echo "ipReads=c(\"S11.bed\",\"S12.bed\",\"S21.bed\",\"S22.bed\")," >> $SCRIPT
      echo "ctReads=c(\"IN1.bed\",\"IN1.bed\",\"IN2.bed\",\"IN2.bed\")," >> $SCRIPT
      echo "peaks=c(\"s11_peaks.bed\",\"s12_peaks.bed\",\"s21_peaks.bed\",\"s22_peaks.bed\")" >> $SCRIPT
      echo ")" >> $SCRIPT
      echo "design=as.data.frame(lapply(conf[,c(\"condition\",\"factor\")],as.numeric))-1" >> $SCRIPT
      echo "design=as.data.frame(model.matrix(~condition,design))" >> $SCRIPT
      echo "countSet=makeCountSet(conf,design, filetype=\"bed\",species=\"other\",binsize=10)" >> $SCRIPT
      echo "countSet=ChIPComp(countSet)" >> $SCRIPT
      
      #save as csv..
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file=\"results.csv\", as.data.frame(countSet\$db),quote = FALSE)" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #get average memory or maximum in kbytes
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      #R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT".bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT" $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT" $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
        #reformat for eval
        cat results.csv | grep -v "^chr"  | awk '{lfold=(($7+$8)==0) ? "inf" : log((($5+$6)/2)/(($7+$8)/2))/log(2); print $2"\t"$3"\t"$4"\t"$14"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
        #create empty file
        touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv s11_peaks.bed s12_peaks.bed s21_peaks.bed s22_peaks.bed mem.txt
      
    done
  done
  
  rm -f S11.bed S12.bed IN1.bed S21.bed S22.bed IN2.bed

else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
