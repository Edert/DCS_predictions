
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
  
  bamToBed -i $4 | cut -f1,2,3,6 > S11.bed
  bamToBed -i $5 | cut -f1,2,3,6 > S12.bed
  #bamToBed -i $6 > IN1.bed
  bamToBed -i $7 | cut -f1,2,3,6 > S21.bed
  bamToBed -i $8 | cut -f1,2,3,6 > S22.bed
  #bamToBed -i $9 > IN2.bed
  
  
  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt

  for PCALLER in ../../../results_peaks/*; do
  #for PCALLER in "../../../results_peaks/sicer"; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
      
      cat $PMODE/s11_peaks.bed | sort -k1,1 -k2,2n > s11_peaks.bed
      cat $PMODE/s12_peaks.bed | sort -k1,1 -k2,2n > s12_peaks.bed
      cat $PMODE/s21_peaks.bed | sort -k1,1 -k2,2n > s21_peaks.bed
      cat $PMODE/s22_peaks.bed | sort -k1,1 -k2,2n > s22_peaks.bed
      
      profile_bins --peaks=s11_peaks.bed,s12_peaks.bed,s21_peaks.bed,s22_peaks.bed \
      --reads=S11.bed,S12.bed,S21.bed,S22.bed --labs=s1,s1,s2,s2 -n MAnorm2_result  >> $LOG 2>&1
      
      
      # in R #norm1 #fit1
      echo "library(MAnorm2)" > $SCRIPT
      echo "data <- read.table('MAnorm2_result_profile_bins.xls',header = T)" >> $SCRIPT
      echo "norm <- normalize(data, count = 4:7, occupancy = 8:11)" >> $SCRIPT
      echo "conds <- list(s1 = bioCond(norm[4:5], norm[8:9], name = 's1'), s2 = bioCond(norm[6:7], norm[10:11], name = 's2'))" >> $SCRIPT
      echo "conds <- normBioCond(conds)" >> $SCRIPT
      echo "conds <- fitMeanVarCurve(conds, method = 'parametric', occupy.only = TRUE, init.coef = c(0.1, 10))" >> $SCRIPT
      echo "res <- diffTest(conds[[1]], conds[[2]])" >> $SCRIPT
      echo "out <- cbind(as.character(norm\$chrom),norm\$start,norm\$end,res\$padj,res\$Mval)" >> $SCRIPT
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(out, file=\"results.csv\", sep='\t', col.names=FALSE, row.names=FALSE, quote = FALSE)" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT

      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_1.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_1 $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_1 $MEMUSAGE" >> memory.txt
      
      #reformat for eval 
      if [ -e results.csv ]; then
         #reformat for eval
         cat results.csv | grep -v "^chr\sstart" | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout results.csv mem.txt
      
      
      # in R #norm1 #fit2
      echo "library(MAnorm2)" > $SCRIPT
      echo "data <- read.table('MAnorm2_result_profile_bins.xls',header = T)" >> $SCRIPT
      echo "norm <- normalize(data, count = 4:7, occupancy = 8:11)" >> $SCRIPT
      echo "conds <- list(s1 = bioCond(norm[4:5], norm[8:9], name = 's1'), s2 = bioCond(norm[6:7], norm[10:11], name = 's2'))" >> $SCRIPT
      echo "conds <- fitMeanVarCurve(conds, method = 'local', occupy.only = FALSE)" >> $SCRIPT
      echo "conds <- estimatePriorDf(conds, occupy.only = TRUE)" >> $SCRIPT
      echo "res <- diffTest(conds[[1]], conds[[2]])" >> $SCRIPT
      echo "out <- cbind(as.character(norm\$chrom),norm\$start,norm\$end,res\$padj,res\$Mval)" >> $SCRIPT
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(out, file=\"results.csv\", sep='\t', col.names=FALSE, row.names=FALSE, quote = FALSE)" >> $SCRIPT
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_2.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_2 $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_2 $MEMUSAGE" >> memory.txt
      
      #reformat for eval 
      if [ -e results.csv ]; then
         #reformat for eval
         cat results.csv | grep -v "^chr\sstart" | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout results.csv mem.txt
      
      
      # in R  #norm2 #fit1
      echo "library(MAnorm2)" > $SCRIPT
      echo "data <- read.table('MAnorm2_result_profile_bins.xls',header = T)" >> $SCRIPT
      echo "norm <- normalize(data, count = 4:5, occupancy = 8:9)" >> $SCRIPT
      echo "norm <- normalize(norm, count = 6:7, occupancy = 10:11)" >> $SCRIPT
      echo "conds <- list(s1 = bioCond(norm[4:5], norm[8:9], name = 's1'), s2 = bioCond(norm[6:7], norm[10:11], name = 's2'))" >> $SCRIPT
      echo "conds <- normBioCond(conds)" >> $SCRIPT
      echo "conds <- fitMeanVarCurve(conds, method = 'parametric', occupy.only = TRUE, init.coef = c(0.1, 10))" >> $SCRIPT
      echo "res <- diffTest(conds[[1]], conds[[2]])" >> $SCRIPT
      echo "out <- cbind(as.character(norm\$chrom),norm\$start,norm\$end,res\$padj,res\$Mval)" >> $SCRIPT
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(out, file=\"results.csv\", sep='\t', col.names=FALSE, row.names=FALSE, quote = FALSE)" >> $SCRIPT
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_3.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_3 $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_3 $MEMUSAGE" >> memory.txt
      
      #reformat for eval 
      if [ -e results.csv ]; then
         #reformat for eval
         cat results.csv | grep -v "^chr\sstart" | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout results.csv mem.txt
      
      # in R #norm2 #fit2
      echo "library(MAnorm2)" > $SCRIPT
      echo "data <- read.table('MAnorm2_result_profile_bins.xls',header = T)" >> $SCRIPT
      echo "norm <- normalize(data, count = 4:5, occupancy = 8:9)" >> $SCRIPT
      echo "norm <- normalize(norm, count = 6:7, occupancy = 10:11)" >> $SCRIPT
      echo "conds <- list(s1 = bioCond(norm[4:5], norm[8:9], name = 's1'), s2 = bioCond(norm[6:7], norm[10:11], name = 's2'))" >> $SCRIPT
      echo "conds <- fitMeanVarCurve(conds, method = 'local', occupy.only = FALSE)" >> $SCRIPT
      echo "conds <- estimatePriorDf(conds, occupy.only = TRUE)" >> $SCRIPT
      echo "res <- diffTest(conds[[1]], conds[[2]])" >> $SCRIPT
      echo "out <- cbind(as.character(norm\$chrom),norm\$start,norm\$end,res\$padj,res\$Mval)" >> $SCRIPT
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(out, file=\"results.csv\", sep='\t', col.names=FALSE, row.names=FALSE, quote = FALSE)" >> $SCRIPT
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_4.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_4 $TIMEDIFF" >> time.txt
      
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_4 $MEMUSAGE" >> memory.txt
      
      #reformat for eval 
      if [ -e results.csv ]; then
         #reformat for eval
         cat results.csv | grep -v "^chr\sstart" | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
         #create empty file
         touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout results.csv mem.txt
      
      #clean up peak files
      rm -f s11_peaks.bed s12_peaks.bed s21_peaks.bed s22_peaks.bed MAnorm2_result_profile_bins.xls MAnorm2_result_profile_bins_log.txt
    done
  done
  
  rm -f S11.bed S12.bed S21.bed S22.bed 
      
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
