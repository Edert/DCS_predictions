NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL

files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13"
  
  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  SCRIPT="script.R"
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  bamToBed -i $4 > S1.bed
  bamToBed -i $5 >> S1.bed
  bamToBed -i $6 > IN1.bed
  bamToBed -i $7 > S2.bed
  bamToBed -i $8 >> S2.bed
  bamToBed -i $9 > IN2.bed
  
  # in R norm1 bin 100
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='100',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='100',species2='mouse',norm.method='1',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT

  PREPDONE=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R

  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_1.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_1.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_1.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_1.bed"
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  
  
  # in R norm2 bin 100
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='100',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='100',species2='mouse',norm.method='2',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_2.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_2.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_2.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_2.bed"
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  
  
  # in R norm3 bin 100
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='100',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='100',species2='mouse',norm.method='3',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_3.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_3.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_3.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_3.bed" 
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  
  
  # in R norm4 bin 100
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='100',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='100',species2='mouse',norm.method='4',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "4 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_4.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_4.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_4.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_4.bed"
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  


  # in R norm1 bin 1000
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='1000',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='1000',species2='mouse',norm.method='1',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "5 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "5 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_5.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_5.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_5.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_5.bed"
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  
  
  # in R norm2 bin 1000
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='1000',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='1000',species2='mouse',norm.method='2',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "6 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "6 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_6.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_6.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_6.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_6.bed"
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  
  
  # in R norm3 bin 1000
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='1000',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='1000',species2='mouse',norm.method='3',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "7 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "7 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_7.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_7.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_7.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_7.bed"
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  
  
  # in R norm4 bin 1000
  echo "library(QChIPat)" > $SCRIPT
  echo "peak.compare(peak.method='b',inpstr1='S1.bed',fmatstr1='b',ctl1=T,ctlfile1='IN1.bed',bin1='1000',species1='mouse',inpstr2='S2.bed',fmatstr2='b',ctl2=T,ctlfile2='IN2.bed',bin2='1000',species2='mouse',norm.method='4',wil.pair=F,p.value='0.05',ratio='0.7')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "8 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "8 $MEMUSAGE" >> memory.txt
  
  #reformat for eval
  if [ -e Unchanged.xls ]; then
    cut -f4,6,7,15,16 Unchanged.xls | grep -v "^Chromosome" > $OUT_NAME"_8.bed"
  fi
  if [ -e BindingPattern.xls ]; then
    cut -f4,6,7,15,16 BindingPattern.xls | grep -v "^Chromosome" >> $OUT_NAME"_8.bed"
  fi
  if [ -e GeneDesert.xls ]; then
    cut -f3,5,6,14,15 GeneDesert.xls | grep -v "^Chromosome" >> $OUT_NAME"_8.bed"
  fi
  #in case no result files were created...
  touch $OUT_NAME"_8.bed"
  
  #save log 
  cat script.Rout >> $LOG
  #clean up
  rm -f BindingPattern.xls Unchanged.xls GeneDesert.xls $SCRIPT script.Rout mem.txt
  rm -f SummaryInformation.txt Enriched_S1.bed_desert.wig Enriched_S1.bed_pattern.wig Enriched_S2.bed_desert.wig Enriched_S2.bed_pattern.wig
  

  #clean up
  rm -f .RData S1.bed S2.bed IN1.bed IN2.bed
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi

