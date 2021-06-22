
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

  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt

  for PCALLER in ../../../results_peaks/*; do
  #for PCALLER in "../../../results_peaks/sicer"; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
      
      #merge all peaks to one file
      cut -f1,2,3 $PMODE/s11_peaks.bed > peaks.bed
      cut -f1,2,3 $PMODE/s12_peaks.bed >> peaks.bed
      cut -f1,2,3 $PMODE/s21_peaks.bed >> peaks.bed
      cut -f1,2,3 $PMODE/s22_peaks.bed >> peaks.bed

      bedtools sort -i peaks.bed > peaks_sorted.bed
      bedtools merge -i peaks_sorted.bed > peaks_merged.bed
      
      #get coverage per peak from the merged file
      bedtools multicov -bams $4 $5 $7 $8 -bed peaks_merged.bed > coverage.txt

      # in R edgeR
      echo "library(edgeR)" > $SCRIPT
      echo "cov <- read.table('coverage.txt')" >> $SCRIPT
      echo "cov\$id <- with(cov, paste0(V1,'_',V2,'_',V3))" >> $SCRIPT
      echo "group <- c('S1','S1','S2','S2')" >> $SCRIPT
      echo "cnts <- cov[,4:7]" >> $SCRIPT
      echo "y <- DGEList(counts=cnts, group=group)" >> $SCRIPT
      echo "y <- calcNormFactors(y)" >> $SCRIPT
      echo "design <- model.matrix(~group)" >> $SCRIPT
      echo "y <- estimateDisp(y, design)" >> $SCRIPT
      echo "fit <- glmQLFit(y, design)" >> $SCRIPT
      echo "lrt <- glmLRT(fit, coef = 2)" >> $SCRIPT
      echo "results_edgeR <- as.data.frame(topTags(lrt, n = nrow(cnts), sort.by = 'none'))" >> $SCRIPT
      echo "res <- cbind(as.character(cov\$V1),cov\$V2,cov\$V3,results_edgeR\$FDR,results_edgeR\$logFC)" >> $SCRIPT
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file='results.csv', res,quote = FALSE,row.names=F,col.names=F, sep ='\t')" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
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
        cat results.csv | sort -k1,1 -k2,2n > $OUT_NAME
      else
        #create empty file
        touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv peaks.bed peaks_sorted.bed peaks_merged.bed coverage.txt mem.txt
      
    done
  done
  
  rm -f S11.bed S12.bed IN1.bed S21.bed S22.bed IN2.bed

else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
