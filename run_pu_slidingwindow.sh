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
  
  # in R TF edgeR
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "bam.files <- c('$4','$5','$7','$8')" >> $SCRIPT
  echo "standard.chr <- paste0('chr', c(19))" >> $SCRIPT
  echo "param <- readParam(minq=50, restrict=standard.chr)" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=10, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)" >> $SCRIPT
  echo "win.data <- normFactors(bins, se.out=win.data)" >> $SCRIPT #now normFactors instead of normOffsets if we want TMM
  echo "normfacs <- win.data\$norm.factors" >> $SCRIPT
  echo "y.bin <- asDGEList(bins)" >> $SCRIPT
  echo "bin.ab <- aveLogCPM(asDGEList(bins))" >> $SCRIPT
  echo "adjc <- cpm(y.bin, log=TRUE)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type='global')" >> $SCRIPT
  echo "min.fc <- 3" >> $SCRIPT
  echo "keep <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- subset(win.data,keep)" >> $SCRIPT
  echo "genotype <- factor(c('S1','S1','S2','S2'))" >> $SCRIPT
  echo "design <- model.matrix(~0+genotype)" >> $SCRIPT
  echo "colnames(design) <- levels(genotype)" >> $SCRIPT
  echo "y <- asDGEList(filtered.data, norm.factors=normfacs)" >> $SCRIPT
  echo "y <- estimateDisp(y, design)" >> $SCRIPT
  echo "fit <- glmQLFit(y, design, robust=TRUE)" >> $SCRIPT
  echo "contrast <- makeContrasts(S2-S1, levels=design)" >> $SCRIPT #Loess
  echo "res <- glmQLFTest(fit, contrast=contrast)" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res\$table)" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res\$table)" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(filtered.data[tabbest\$best]))),best.logFC=tabbest\$logFC,best.FDR=tabbest\$FDR)" >> $SCRIPT
  echo "outres <- data.frame(out.ranges)[,c('seqnames','start','end','best.FDR','best.logFC')]" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"results.csv\", outres, quote = FALSE, row.names = F, col.names = F, sep='\t')" >> $SCRIPT

  PREPDONE=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "edger_tf $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "edger_tf $MEMUSAGE" >> memory.txt
  
  #save result
  #reformat for eval 
  if [ -e results.csv ]; then
    cat results.csv | sort -k1,1 -k2,2n > $OUT_NAME"_edger_tf.bed"
  else
    #create empty file
    touch $OUT_NAME"_edger_tf.bed"
  fi
  
  #save log 
  cat script.Rout >> $LOG

  #clean up
  rm -f $SCRIPT .RData script.Rout results.csv mem.txt
  
 
  # in R histone edgeR
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "bam.files <- c('$4','$5','$7','$8')" >> $SCRIPT
  echo "standard.chr <- paste0('chr', c(19))" >> $SCRIPT
  echo "param <- readParam(minq=50, restrict=standard.chr)" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=150, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type='global')" >> $SCRIPT
  echo "min.fc <- 3" >> $SCRIPT
  echo "keep <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- subset(win.data,keep)" >> $SCRIPT
  echo "filtered.data <- normOffsets(filtered.data, se.out=TRUE)" >> $SCRIPT #loess  (former [, type='loess'] was used... in normOffsets)
  echo "celltype <- factor(c('S1','S1','S2','S2'))" >> $SCRIPT
  echo "design <- model.matrix(~0+celltype)" >> $SCRIPT
  echo "colnames(design) <- levels(celltype)" >> $SCRIPT
  echo "y <- asDGEList(filtered.data)" >> $SCRIPT
  echo "y <- estimateDisp(y, design)" >> $SCRIPT
  echo "fit <- glmQLFit(y, design, robust=TRUE)" >> $SCRIPT
  echo "contrast <- makeContrasts(S2-S1, levels=design)" >> $SCRIPT
  echo "res <- glmQLFTest(fit, contrast=contrast)" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res\$table)" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res\$table)" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(filtered.data[tabbest\$best]))),best.logFC=tabbest\$logFC,best.FDR=tabbest\$FDR)" >> $SCRIPT
  echo "outres <- data.frame(out.ranges)[,c('seqnames','start','end','best.FDR','best.logFC')]" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"results.csv\", outres, quote = FALSE, row.names = F, col.names = F, sep='\t')" >> $SCRIPT

  STARTTIME=`date +%s.%N`

  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "edger_hist $TIMEDIFF" >> time.txt 
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "edger_hist $MEMUSAGE" >> memory.txt
  
  if [ -e results.csv ]; then
     #reformat for eval
     cat results.csv | sort -k1,1 -k2,2n > $OUT_NAME"_edger_hist.bed"
  else
     #create empty file
     touch $OUT_NAME"_edger_hist.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.csv mem.txt
  
  
  # in R tf deseq2
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "library(DESeq2)" >> $SCRIPT
  echo "bam.files <- c('$4','$5','$7','$8')" >> $SCRIPT
  echo "standard.chr <- paste0('chr', c(19))" >> $SCRIPT
  echo "param <- readParam(minq=50, restrict=standard.chr)" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=10, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)" >> $SCRIPT
  echo "win.data <- normFactors(bins, se.out=win.data)" >> $SCRIPT #now normFactors instead of normOffsets if we want TMM
  echo "normfacs <- win.data\$norm.factors" >> $SCRIPT
  echo "y.bin <- asDGEList(bins)" >> $SCRIPT
  echo "bin.ab <- aveLogCPM(asDGEList(bins))" >> $SCRIPT
  echo "adjc <- cpm(y.bin, log=TRUE)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type='global')" >> $SCRIPT
  echo "min.fc <- 3" >> $SCRIPT
  echo "keep <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- subset(win.data,keep)" >> $SCRIPT
  echo "se <- filtered.data" >> $SCRIPT
  echo "cond <- factor(c('S1','S1','S2','S2'))" >> $SCRIPT
  echo "colData(se) <- cbind(colData(se),condition=cond)" >> $SCRIPT
  echo "ddsSE <- DESeqDataSet(se, design = ~condition)" >> $SCRIPT
  echo "ddsSE <- DESeq(ddsSE)" >> $SCRIPT
  echo "res <- data.frame(results(ddsSE, contrast=c('condition','S2','S1')))" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(se), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res, fc.col=c(2),pval.col=c(6))" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res,pval.col=c(6))" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(se[tabbest\$best]))),best.logFC=tabbest\$log2FoldChange,best.padj=tabbest\$padj)" >> $SCRIPT
  echo "outres <- data.frame(out.ranges)[,c('seqnames','start','end','best.padj','best.logFC')]" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"results.csv\", outres, quote = FALSE, row.names = F, col.names = F, sep='\t')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`

  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "deseq2_tf $TIMEDIFF" >> time.txt
   
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "deseq2_tf $MEMUSAGE" >> memory.txt
  
  if [ -e results.csv ]; then
     #reformat for eval
     cat results.csv | sort -k1,1 -k2,2n > $OUT_NAME"_deseq2_tf.bed"
  else
     #create empty file
     touch $OUT_NAME"_deseq2_tf.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.csv mem.txt


  # in R histone deseq2
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "library(DESeq2)" >> $SCRIPT
  echo "bam.files <- c('$4','$5','$7','$8')" >> $SCRIPT
  echo "standard.chr <- paste0('chr', c(19))" >> $SCRIPT
  echo "param <- readParam(minq=50, restrict=standard.chr)" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=150, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type='global')" >> $SCRIPT
  echo "min.fc <- 3" >> $SCRIPT
  echo "keep <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- subset(win.data,keep)" >> $SCRIPT
  echo "filtered.data <- normOffsets(filtered.data, se.out=TRUE)" >> $SCRIPT #loess  (former [, type='loess'] was used... in normOffsets)
  echo "se <- filtered.data" >> $SCRIPT
  echo "cond <- factor(c('S1','S1','S2','S2'))" >> $SCRIPT
  echo "colData(se) <- cbind(colData(se),condition=cond)" >> $SCRIPT
  echo "ddsSE <- DESeqDataSet(se, design = ~condition)" >> $SCRIPT
  echo "ddsSE <- DESeq(ddsSE)" >> $SCRIPT
  echo "res <- data.frame(results(ddsSE, contrast=c('condition','S2','S1')))" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(se), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res, fc.col=c(2),pval.col=c(6))" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res,pval.col=c(6))" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(se[tabbest\$best]))),best.logFC=tabbest\$log2FoldChange,best.padj=tabbest\$padj)" >> $SCRIPT
  echo "outres <- data.frame(out.ranges)[,c('seqnames','start','end','best.padj','best.logFC')]" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"results.csv\", outres, quote = FALSE, row.names = F, col.names = F, sep='\t')" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`

  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "deseq2_hist $TIMEDIFF" >> time.txt 
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "deseq2_hist $MEMUSAGE" >> memory.txt
  
  if [ -e results.csv ]; then
     #reformat for eval
     cat results.csv | sort -k1,1 -k2,2n > $OUT_NAME"_deseq2_hist.bed"
  else
     #create empty file
     touch $OUT_NAME"_deseq2_hist.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.csv mem.txt
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
