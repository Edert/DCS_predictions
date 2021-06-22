
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15"

  SCRIPT="script.R"
  LOG="../../../log/$TOOL/$SET.log"
  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  samtools merge s1.bam $4 $5 
  samtools merge s2.bam $7 $8 
  samtools index s1.bam
  samtools index s2.bam
  
  # in R 
  echo "library(normr)" > $SCRIPT
  echo "genome <-data.frame(UCSC_seqlevel=c(\"chr19\"),UCSC_seqlength=c(59128983))" >> $SCRIPT

  echo "countConfigSE <- countConfigSingleEnd(binsize = 200, mapq = 20,filteredFlag = 1024,shift = 0)" >> $SCRIPT
  echo "de <- diffR(treatment = 's1.bam',control = 's2.bam',genome = genome,countConfig= countConfigSE)" >> $SCRIPT
  echo "res <- data.frame(getQvalues(de),getRanges(de),getEnrichment(de))" >> $SCRIPT
  echo "res <- res[complete.cases(res), ]" >> $SCRIPT
  #save as file
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"results.csv\", res,quote = FALSE,row.names = F,col.names = F,sep='\t')" >> $SCRIPT
  
  PREPDONE=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e results.csv ]; then
     #reformat for eval
     cat results.csv | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$8}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi
  #save log     
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT script.Rout results.csv mem.txt
  
  
  echo "library(normr)" > $SCRIPT
  echo "genome <-data.frame(UCSC_seqlevel=c(\"chr19\"),UCSC_seqlength=c(59128983))" >> $SCRIPT

  echo "countConfigSE <- countConfigSingleEnd(binsize = 2000, mapq = 20,filteredFlag = 1024,shift = 0)" >> $SCRIPT
  echo "de <- diffR(treatment = 's1.bam',control = 's2.bam',genome = genome,countConfig= countConfigSE)" >> $SCRIPT
  echo "res <- data.frame(getQvalues(de),getRanges(de),getEnrichment(de))" >> $SCRIPT
  echo "res <- res[complete.cases(res), ]" >> $SCRIPT
  #save as file
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"results.csv\", res,quote = FALSE,row.names = F,col.names = F,sep='\t')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e results.csv ]; then
     #reformat for eval
     cat results.csv | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$8}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi
  #save log     
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT script.Rout results.csv mem.txt
  rm -rf s1.bam s2.bam s1.bam.bai s2.bam.bai
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
