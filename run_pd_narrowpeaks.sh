
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
  
  SF1=$(echo $10 | awk '{printf "%.10f", $1-7}')
  SF2=$(echo $11 | awk '{printf "%.10f", $1-7}')
  SF3=$(echo $12 | awk '{printf "%.10f", $1-7}')
  SF4=$(echo $13 | awk '{printf "%.10f", $1-7}')
  
  SCRIPT="script.R"
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  #bigwig files, no norm
  bamCoverage -b $4 -o S11.bw --binSize 1 -p 1 >> $LOG 2>&1
  bamCoverage -b $5 -o S12.bw --binSize 1 -p 1 >> $LOG 2>&1
  bamCoverage -b $7 -o S21.bw --binSize 1 -p 1 >> $LOG 2>&1
  bamCoverage -b $8 -o S22.bw --binSize 1 -p 1 >> $LOG 2>&1
  
  PREPDONE=`date +%s.%N` #stop recodring here as we prepare the data for all parameter sets at once...
  
  #run with different normalization: RPKM
  bamCoverage -b $4 -o S11r.bw --binSize 1 -p 1 --normalizeUsing RPKM >> $LOG 2>&1
  bamCoverage -b $5 -o S12r.bw --binSize 1 -p 1 --normalizeUsing RPKM >> $LOG 2>&1
  bamCoverage -b $7 -o S21r.bw --binSize 1 -p 1 --normalizeUsing RPKM >> $LOG 2>&1
  bamCoverage -b $8 -o S22r.bw --binSize 1 -p 1 --normalizeUsing RPKM >> $LOG 2>&1
  #run with different normalization: CPM
  bamCoverage -b $4 -o S11c.bw --binSize 1 -p 1 --normalizeUsing CPM >> $LOG 2>&1
  bamCoverage -b $5 -o S12c.bw --binSize 1 -p 1 --normalizeUsing CPM >> $LOG 2>&1
  bamCoverage -b $7 -o S21c.bw --binSize 1 -p 1 --normalizeUsing CPM >> $LOG 2>&1
  bamCoverage -b $8 -o S22c.bw --binSize 1 -p 1 --normalizeUsing CPM >> $LOG 2>&1
  #run with different normalization: BPM
  bamCoverage -b $4 -o S11b.bw --binSize 1 -p 1 --normalizeUsing BPM >> $LOG 2>&1
  bamCoverage -b $5 -o S12b.bw --binSize 1 -p 1 --normalizeUsing BPM >> $LOG 2>&1
  bamCoverage -b $7 -o S21b.bw --binSize 1 -p 1 --normalizeUsing BPM >> $LOG 2>&1
  bamCoverage -b $8 -o S22b.bw --binSize 1 -p 1 --normalizeUsing BPM >> $LOG 2>&1
  #run with different normalization: RPGC
  #all bases minus N is effective genome size
  #grep -v ">" ../01_simulating/data/mm10_chr19.fasta | wc | awk '{print $3-$1}'
  #cat ../01_simulating/data/mm10_chr19.fasta | grep -v "^>" | tr -cd N | wc -c
  bamCoverage -b $4 -o S11g.bw --binSize 1 -p 1 --normalizeUsing RPGC --effectiveGenomeSize 58205856 >> $LOG 2>&1
  bamCoverage -b $5 -o S12g.bw --binSize 1 -p 1 --normalizeUsing RPGC --effectiveGenomeSize 58205856 >> $LOG 2>&1
  bamCoverage -b $7 -o S21g.bw --binSize 1 -p 1 --normalizeUsing RPGC --effectiveGenomeSize 58205856 >> $LOG 2>&1
  bamCoverage -b $8 -o S22g.bw --binSize 1 -p 1 --normalizeUsing RPGC --effectiveGenomeSize 58205856 >> $LOG 2>&1
  #run with scaling factor
  bamCoverage -b $4 -o S11s.bw --binSize 1 -p 1 --scaleFactor $SF1 >> $LOG 2>&1
  bamCoverage -b $5 -o S12s.bw --binSize 1 -p 1 --scaleFactor $SF2 >> $LOG 2>&1
  bamCoverage -b $7 -o S21s.bw --binSize 1 -p 1 --scaleFactor $SF3 >> $LOG 2>&1
  bamCoverage -b $8 -o S22s.bw --binSize 1 -p 1 --scaleFactor $SF4 >> $LOG 2>&1
  
  #get script
  cp /proj/chipseq_norm_diffbind_062017/analysis/03_db_analysis/bin/narrowpeaks/diffNGS.R .
  
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
      #bed file only 3 columns!
      bedtools merge -i peaks_sorted.bed | cut -f1,2,3 > peaks_merged.bed
      
      #bigwig files, no norm
      # in R 
      echo "library(NarrowPeaks)" > $SCRIPT
      echo "library(GenomicRanges)" >> $SCRIPT
      echo "require('fda')" >> $SCRIPT
      echo "require('ICSNP')" >> $SCRIPT
      echo "require('genomation')" >> $SCRIPT
      echo "require('GenomicRanges')" >> $SCRIPT
      echo "require('RColorBrewer')" >> $SCRIPT
      echo "setwd('.')" >> $SCRIPT
      echo "CNDS     <- c('S1','S1','S2','S2')" >> $SCRIPT
      echo "bws      <- c('S11.bw','S12.bw','S21.bw','S22.bw')" >> $SCRIPT
      echo "peaks    <-  'peaks_merged.bed'" >> $SCRIPT
      echo "Fl <- 0" >> $SCRIPT         # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
      echo "Nbasis <- 10" >> $SCRIPT    # number of B-spline basis used in the functiona PCA analysis
      echo "Bins <- 31" >> $SCRIPT      # Number of bins in genomation
      echo "bed <- readGeneric(peaks, keep.all.metadata = FALSE)" >> $SCRIPT 
      echo "results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse='_vs_'), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )" >> $SCRIPT 
      echo "source('diffNGS.R')" >> $SCRIPT 
      echo "x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis=Nbasis, NB=Bins)" >> $SCRIPT 
      echo "for (j in 1:length(x\$p.values)  ){" >> $SCRIPT 
      echo "  results\$pval[j] <- x\$p.values[[j]][1,2]" >> $SCRIPT 
      echo "  results\$avg.C2[j] <- mean( max(x\$fdaprofiles[[3]][j,]) , max(x\$fdaprofiles[[4]][j,]) ) " >> $SCRIPT 
      echo "  results\$avg.C1[j] <- mean( max(x\$fdaprofiles[[1]][j,]) , max(x\$fdaprofiles[[2]][j,]) ) " >> $SCRIPT 
      echo "  if ( results\$avg.C2[j] >= results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C2[j] + 0.001) / ( results\$avg.C1[j] + 0.001  )  )  }" >> $SCRIPT     # Increase
      echo "  if ( results\$avg.C2[j] <  results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C1[j] + 0.001) / ( results\$avg.C2[j]  + 0.001 )  )  }" >> $SCRIPT     # Decrease
      echo "}" >> $SCRIPT 
      echo "zeroP <- which(results\$pval == 0.0)" >> $SCRIPT 
      echo "results\$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])" >> $SCRIPT 
      echo "results\$fdr <- p.adjust(results\$pval, method = 'fdr')" >> $SCRIPT 
      echo "results_sorted  <- results[order(results\$fdr, decreasing = FALSE),] " >> $SCRIPT 
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file='results.csv', results_sorted,quote = FALSE,row.names=F,col.names=F, sep ='\t')" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_1 $TIMEDIFF" >> time.txt
      
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_1 $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
        #move results
        cat results.csv | grep -v "region.chr" |  awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9}' | sort -k1,1 -k2,2n >  $OUT_NAME"_1.bed"
      else
        #create empty file
        touch $OUT_NAME"_1.bed"
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv peaks.bed peaks_sorted.bed S2_vs_S1.mean.pdf mem.txt
      
      
      
      #run with different normalization: RPKM
      # in R 
      echo "library(NarrowPeaks)" > $SCRIPT
      echo "library(GenomicRanges)" >> $SCRIPT
      echo "require('fda')" >> $SCRIPT
      echo "require('ICSNP')" >> $SCRIPT
      echo "require('genomation')" >> $SCRIPT
      echo "require('GenomicRanges')" >> $SCRIPT
      echo "require('RColorBrewer')" >> $SCRIPT
      echo "setwd('.')" >> $SCRIPT
      echo "CNDS     <- c('S1','S1','S2','S2')" >> $SCRIPT
      echo "bws      <- c('S11r.bw','S12r.bw','S21r.bw','S22r.bw')" >> $SCRIPT
      echo "peaks    <-  'peaks_merged.bed'" >> $SCRIPT
      echo "Fl <- 0" >> $SCRIPT         # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
      echo "Nbasis <- 10" >> $SCRIPT    # number of B-spline basis used in the functiona PCA analysis
      echo "Bins <- 31" >> $SCRIPT      # Number of bins in genomation
      echo "bed <- readGeneric(peaks, keep.all.metadata = FALSE)" >> $SCRIPT 
      echo "results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse='_vs_'), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )" >> $SCRIPT 
      echo "source('diffNGS.R')" >> $SCRIPT 
      echo "x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis=Nbasis, NB=Bins)" >> $SCRIPT 
      echo "for (j in 1:length(x\$p.values)  ){" >> $SCRIPT 
      echo "  results\$pval[j] <- x\$p.values[[j]][1,2]" >> $SCRIPT 
      echo "  results\$avg.C2[j] <- mean( max(x\$fdaprofiles[[3]][j,]) , max(x\$fdaprofiles[[4]][j,]) ) " >> $SCRIPT 
      echo "  results\$avg.C1[j] <- mean( max(x\$fdaprofiles[[1]][j,]) , max(x\$fdaprofiles[[2]][j,]) ) " >> $SCRIPT 
      echo "  if ( results\$avg.C2[j] >= results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C2[j] + 0.001) / ( results\$avg.C1[j] + 0.001  )  )  }" >> $SCRIPT     # Increase
      echo "  if ( results\$avg.C2[j] <  results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C1[j] + 0.001) / ( results\$avg.C2[j]  + 0.001 )  )  }" >> $SCRIPT     # Decrease
      echo "}" >> $SCRIPT 
      echo "zeroP <- which(results\$pval == 0.0)" >> $SCRIPT 
      echo "results\$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])" >> $SCRIPT 
      echo "results\$fdr <- p.adjust(results\$pval, method = 'fdr')" >> $SCRIPT 
      echo "results_sorted  <- results[order(results\$fdr, decreasing = FALSE),] " >> $SCRIPT 
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file='results.csv', results_sorted,quote = FALSE,row.names=F,col.names=F, sep ='\t')" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_2 $TIMEDIFF" >> time.txt
      
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_2 $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
        #move results
        cat results.csv | grep -v "region.chr" |  awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9}' | sort -k1,1 -k2,2n >  $OUT_NAME"_2.bed"
      else
        #create empty file
        touch $OUT_NAME"_2.bed"
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv peaks.bed peaks_sorted.bed S2_vs_S1.mean.pdf mem.txt
      
      
      
      #run with different normalization: CPM
      # in R 
      echo "library(NarrowPeaks)" > $SCRIPT
      echo "library(GenomicRanges)" >> $SCRIPT
      echo "require('fda')" >> $SCRIPT
      echo "require('ICSNP')" >> $SCRIPT
      echo "require('genomation')" >> $SCRIPT
      echo "require('GenomicRanges')" >> $SCRIPT
      echo "require('RColorBrewer')" >> $SCRIPT
      echo "setwd('.')" >> $SCRIPT
      echo "CNDS     <- c('S1','S1','S2','S2')" >> $SCRIPT
      echo "bws      <- c('S11c.bw','S12c.bw','S21c.bw','S22c.bw')" >> $SCRIPT
      echo "peaks    <-  'peaks_merged.bed'" >> $SCRIPT
      echo "Fl <- 0" >> $SCRIPT         # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
      echo "Nbasis <- 10" >> $SCRIPT    # number of B-spline basis used in the functiona PCA analysis
      echo "Bins <- 31" >> $SCRIPT      # Number of bins in genomation
      echo "bed <- readGeneric(peaks, keep.all.metadata = FALSE)" >> $SCRIPT 
      echo "results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse='_vs_'), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )" >> $SCRIPT 
      echo "source('diffNGS.R')" >> $SCRIPT 
      echo "x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis=Nbasis, NB=Bins)" >> $SCRIPT 
      echo "for (j in 1:length(x\$p.values)  ){" >> $SCRIPT 
      echo "  results\$pval[j] <- x\$p.values[[j]][1,2]" >> $SCRIPT 
      echo "  results\$avg.C2[j] <- mean( max(x\$fdaprofiles[[3]][j,]) , max(x\$fdaprofiles[[4]][j,]) ) " >> $SCRIPT 
      echo "  results\$avg.C1[j] <- mean( max(x\$fdaprofiles[[1]][j,]) , max(x\$fdaprofiles[[2]][j,]) ) " >> $SCRIPT 
      echo "  if ( results\$avg.C2[j] >= results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C2[j] + 0.001) / ( results\$avg.C1[j] + 0.001  )  )  }" >> $SCRIPT     # Increase
      echo "  if ( results\$avg.C2[j] <  results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C1[j] + 0.001) / ( results\$avg.C2[j]  + 0.001 )  )  }" >> $SCRIPT     # Decrease
      echo "}" >> $SCRIPT 
      echo "zeroP <- which(results\$pval == 0.0)" >> $SCRIPT 
      echo "results\$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])" >> $SCRIPT 
      echo "results\$fdr <- p.adjust(results\$pval, method = 'fdr')" >> $SCRIPT 
      echo "results_sorted  <- results[order(results\$fdr, decreasing = FALSE),] " >> $SCRIPT 
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file='results.csv', results_sorted,quote = FALSE,row.names=F,col.names=F, sep ='\t')" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_3 $TIMEDIFF" >> time.txt
      
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_3 $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
        #move results
        cat results.csv | grep -v "region.chr" |  awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9}' | sort -k1,1 -k2,2n >  $OUT_NAME"_3.bed"
      else
        #create empty file
        touch $OUT_NAME"_3.bed"
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv peaks.bed peaks_sorted.bed S2_vs_S1.mean.pdf mem.txt
      
      
      
      #run with different normalization: BPM
      # in R 
      echo "library(NarrowPeaks)" > $SCRIPT
      echo "library(GenomicRanges)" >> $SCRIPT
      echo "require('fda')" >> $SCRIPT
      echo "require('ICSNP')" >> $SCRIPT
      echo "require('genomation')" >> $SCRIPT
      echo "require('GenomicRanges')" >> $SCRIPT
      echo "require('RColorBrewer')" >> $SCRIPT
      echo "setwd('.')" >> $SCRIPT
      echo "CNDS     <- c('S1','S1','S2','S2')" >> $SCRIPT
      echo "bws      <- c('S11b.bw','S12b.bw','S21b.bw','S22b.bw')" >> $SCRIPT
      echo "peaks    <-  'peaks_merged.bed'" >> $SCRIPT
      echo "Fl <- 0" >> $SCRIPT         # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
      echo "Nbasis <- 10" >> $SCRIPT    # number of B-spline basis used in the functiona PCA analysis
      echo "Bins <- 31" >> $SCRIPT      # Number of bins in genomation
      echo "bed <- readGeneric(peaks, keep.all.metadata = FALSE)" >> $SCRIPT 
      echo "results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse='_vs_'), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )" >> $SCRIPT 
      echo "source('diffNGS.R')" >> $SCRIPT 
      echo "x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis=Nbasis, NB=Bins)" >> $SCRIPT 
      echo "for (j in 1:length(x\$p.values)  ){" >> $SCRIPT 
      echo "  results\$pval[j] <- x\$p.values[[j]][1,2]" >> $SCRIPT 
      echo "  results\$avg.C2[j] <- mean( max(x\$fdaprofiles[[3]][j,]) , max(x\$fdaprofiles[[4]][j,]) ) " >> $SCRIPT 
      echo "  results\$avg.C1[j] <- mean( max(x\$fdaprofiles[[1]][j,]) , max(x\$fdaprofiles[[2]][j,]) ) " >> $SCRIPT 
      echo "  if ( results\$avg.C2[j] >= results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C2[j] + 0.001) / ( results\$avg.C1[j] + 0.001  )  )  }" >> $SCRIPT     # Increase
      echo "  if ( results\$avg.C2[j] <  results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C1[j] + 0.001) / ( results\$avg.C2[j]  + 0.001 )  )  }" >> $SCRIPT     # Decrease
      echo "}" >> $SCRIPT 
      echo "zeroP <- which(results\$pval == 0.0)" >> $SCRIPT 
      echo "results\$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])" >> $SCRIPT 
      echo "results\$fdr <- p.adjust(results\$pval, method = 'fdr')" >> $SCRIPT 
      echo "results_sorted  <- results[order(results\$fdr, decreasing = FALSE),] " >> $SCRIPT 
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file='results.csv', results_sorted,quote = FALSE,row.names=F,col.names=F, sep ='\t')" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_4 $TIMEDIFF" >> time.txt
      
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_4 $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
        #move results
        cat results.csv | grep -v "region.chr" |  awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9}' | sort -k1,1 -k2,2n >  $OUT_NAME"_4.bed"
      else
        #create empty file
        touch $OUT_NAME"_4.bed"
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv peaks.bed peaks_sorted.bed S2_vs_S1.mean.pdf mem.txt
      
      
      
      #run with different normalization: RPGC
      # in R 
      echo "library(NarrowPeaks)" > $SCRIPT
      echo "library(GenomicRanges)" >> $SCRIPT
      echo "require('fda')" >> $SCRIPT
      echo "require('ICSNP')" >> $SCRIPT
      echo "require('genomation')" >> $SCRIPT
      echo "require('GenomicRanges')" >> $SCRIPT
      echo "require('RColorBrewer')" >> $SCRIPT
      echo "setwd('.')" >> $SCRIPT
      echo "CNDS     <- c('S1','S1','S2','S2')" >> $SCRIPT
      echo "bws      <- c('S11g.bw','S12g.bw','S21g.bw','S22g.bw')" >> $SCRIPT
      echo "peaks    <-  'peaks_merged.bed'" >> $SCRIPT
      echo "Fl <- 0" >> $SCRIPT         # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
      echo "Nbasis <- 10" >> $SCRIPT    # number of B-spline basis used in the functiona PCA analysis
      echo "Bins <- 31" >> $SCRIPT      # Number of bins in genomation
      echo "bed <- readGeneric(peaks, keep.all.metadata = FALSE)" >> $SCRIPT 
      echo "results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse='_vs_'), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )" >> $SCRIPT 
      echo "source('diffNGS.R')" >> $SCRIPT 
      echo "x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis=Nbasis, NB=Bins)" >> $SCRIPT 
      echo "for (j in 1:length(x\$p.values)  ){" >> $SCRIPT 
      echo "  results\$pval[j] <- x\$p.values[[j]][1,2]" >> $SCRIPT 
      echo "  results\$avg.C2[j] <- mean( max(x\$fdaprofiles[[3]][j,]) , max(x\$fdaprofiles[[4]][j,]) ) " >> $SCRIPT 
      echo "  results\$avg.C1[j] <- mean( max(x\$fdaprofiles[[1]][j,]) , max(x\$fdaprofiles[[2]][j,]) ) " >> $SCRIPT 
      echo "  if ( results\$avg.C2[j] >= results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C2[j] + 0.001) / ( results\$avg.C1[j] + 0.001  )  )  }" >> $SCRIPT     # Increase
      echo "  if ( results\$avg.C2[j] <  results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C1[j] + 0.001) / ( results\$avg.C2[j]  + 0.001 )  )  }" >> $SCRIPT     # Decrease
      echo "}" >> $SCRIPT 
      echo "zeroP <- which(results\$pval == 0.0)" >> $SCRIPT 
      echo "results\$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])" >> $SCRIPT 
      echo "results\$fdr <- p.adjust(results\$pval, method = 'fdr')" >> $SCRIPT 
      echo "results_sorted  <- results[order(results\$fdr, decreasing = FALSE),] " >> $SCRIPT 
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file='results.csv', results_sorted,quote = FALSE,row.names=F,col.names=F, sep ='\t')" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_5 $TIMEDIFF" >> time.txt
      
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_5 $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
        #move results
        cat results.csv | grep -v "region.chr" |  awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9}' | sort -k1,1 -k2,2n >  $OUT_NAME"_5.bed"
      else
        #create empty file
        touch $OUT_NAME"_5.bed"
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv peaks.bed peaks_sorted.bed S2_vs_S1.mean.pdf mem.txt
      
      
      
      #run with scaling factor bw´s
      # in R 
      echo "library(NarrowPeaks)" > $SCRIPT
      echo "library(GenomicRanges)" >> $SCRIPT
      echo "require('fda')" >> $SCRIPT
      echo "require('ICSNP')" >> $SCRIPT
      echo "require('genomation')" >> $SCRIPT
      echo "require('GenomicRanges')" >> $SCRIPT
      echo "require('RColorBrewer')" >> $SCRIPT
      echo "setwd('.')" >> $SCRIPT
      echo "CNDS     <- c('S1','S1','S2','S2')" >> $SCRIPT
      echo "bws      <- c('S11s.bw','S12s.bw','S21s.bw','S22s.bw')" >> $SCRIPT
      echo "peaks    <-  'peaks_merged.bed'" >> $SCRIPT
      echo "Fl <- 0" >> $SCRIPT         # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
      echo "Nbasis <- 10" >> $SCRIPT    # number of B-spline basis used in the functiona PCA analysis
      echo "Bins <- 31" >> $SCRIPT      # Number of bins in genomation
      echo "bed <- readGeneric(peaks, keep.all.metadata = FALSE)" >> $SCRIPT 
      echo "results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse='_vs_'), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )" >> $SCRIPT 
      echo "source('diffNGS.R')" >> $SCRIPT 
      echo "x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis=Nbasis, NB=Bins)" >> $SCRIPT 
      echo "for (j in 1:length(x\$p.values)  ){" >> $SCRIPT 
      echo "  results\$pval[j] <- x\$p.values[[j]][1,2]" >> $SCRIPT 
      echo "  results\$avg.C2[j] <- mean( max(x\$fdaprofiles[[3]][j,]) , max(x\$fdaprofiles[[4]][j,]) ) " >> $SCRIPT 
      echo "  results\$avg.C1[j] <- mean( max(x\$fdaprofiles[[1]][j,]) , max(x\$fdaprofiles[[2]][j,]) ) " >> $SCRIPT 
      echo "  if ( results\$avg.C2[j] >= results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C2[j] + 0.001) / ( results\$avg.C1[j] + 0.001  )  )  }" >> $SCRIPT     # Increase
      echo "  if ( results\$avg.C2[j] <  results\$avg.C1[j]) {  results\$log2fc[j]   <-   log2(  ( results\$avg.C1[j] + 0.001) / ( results\$avg.C2[j]  + 0.001 )  )  }" >> $SCRIPT     # Decrease
      echo "}" >> $SCRIPT 
      echo "zeroP <- which(results\$pval == 0.0)" >> $SCRIPT 
      echo "results\$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])" >> $SCRIPT 
      echo "results\$fdr <- p.adjust(results\$pval, method = 'fdr')" >> $SCRIPT 
      echo "results_sorted  <- results[order(results\$fdr, decreasing = FALSE),] " >> $SCRIPT 
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file='results.csv', results_sorted,quote = FALSE,row.names=F,col.names=F, sep ='\t')" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_6 $TIMEDIFF" >> time.txt
      
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_6 $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
        #move results
        cat results.csv | grep -v "region.chr" |  awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9}' | sort -k1,1 -k2,2n >  $OUT_NAME"_6.bed"
      else
        #create empty file
        touch $OUT_NAME"_6.bed"
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv peaks.bed peaks_sorted.bed S2_vs_S1.mean.pdf mem.txt
      
      
      #clean up for next peak caller
      rm -rf peaks_merged.bed
      
    done
  done
  
  #clean up
  rm -rf diffNGS.R 
  rm S*.bw 
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
