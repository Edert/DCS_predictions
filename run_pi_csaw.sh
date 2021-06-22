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
  
  # in R diffbind TMM
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "bam.files<-c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
  echo "standard.chr <- paste0(\"chr\", c(19))" >> $SCRIPT
  echo "param <- readParam(minq=1, restrict=standard.chr)" >> $SCRIPT
  echo "celltype <- c(\"sample1\",\"sample1\",\"sample2\",\"sample2\")" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "min.fc <- 0.7" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=100, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type=\"global\")" >> $SCRIPT
  echo "keep1 <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- win.data[keep1,]" >> $SCRIPT
  echo "celltype <- factor(celltype)" >> $SCRIPT
  echo "design <- model.matrix(~0+celltype)" >> $SCRIPT
  echo "colnames(design) <- levels(celltype)" >> $SCRIPT 
  echo "binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)" >> $SCRIPT
  echo "filtered.data <- normOffsets(binned, se.out=filtered.data)" >> $SCRIPT
  echo "y <- asDGEList(filtered.data, norm.factors=filtered.data\$norm.factors)" >> $SCRIPT
  echo "y <- estimateDisp(y, design)" >> $SCRIPT
  echo "fit <- glmQLFit(y, design, robust=TRUE)" >> $SCRIPT
  echo "contrast <- makeContrasts(sample1-sample2, levels=design)" >> $SCRIPT
  echo "res <- glmQLFTest(fit, contrast=contrast)" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res\$table)" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res\$table)" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(filtered.data[tabbest\$best]))),best.logFC=tabbest\$logFC)" >> $SCRIPT
  echo "write.table(file=\"csaw_TMM.csv\", as.data.frame(out.ranges), quote = FALSE, sep = \"\t\")" >> $SCRIPT
  
  PREPDONE=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1_TMM $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1_TMM $MEMUSAGE" >> memory.txt
  
  #save result
  #reformat for eval take FDR and best log2-foldchange 
  if [ -e csaw_TMM.csv ]; then
    cut -f2,3,4,11,14 csaw_TMM.csv | grep -v "^start" | sort -k1,1 -k2,2n > $OUT_NAME"_1_TMM.bed"
  else
    #create empty file
    touch $OUT_NAME"_1_TMM.bed"
  fi
  
  #save log 
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT .RData script.Rout csaw_TMM.csv mem.txt
  
  
  # in R diffbind SPIKE-IN
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "bam.files<-c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
  echo "standard.chr <- paste0(\"chr\", c(19))" >> $SCRIPT
  echo "param <- readParam(minq=1, restrict=standard.chr)" >> $SCRIPT
  echo "celltype <- c(\"sample1\",\"sample1\",\"sample2\",\"sample2\")" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "min.fc <- 0.7" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=100, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type=\"global\")" >> $SCRIPT
  echo "keep1 <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- win.data[keep1,]" >> $SCRIPT
  echo "celltype <- factor(celltype)" >> $SCRIPT
  echo "design <- model.matrix(~0+celltype)" >> $SCRIPT
  echo "colnames(design) <- levels(celltype)" >> $SCRIPT
  echo "spike.facs <- c($10,$11,$12,$13)" >> $SCRIPT
  echo "y <- asDGEList(filtered.data, norm.factors=spike.facs)" >> $SCRIPT
  echo "y <- estimateDisp(y, design)" >> $SCRIPT
  echo "fit <- glmQLFit(y, design, robust=TRUE)" >> $SCRIPT
  echo "contrast <- makeContrasts(sample1-sample2, levels=design)" >> $SCRIPT
  echo "res <- glmQLFTest(fit, contrast=contrast)" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res\$table)" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res\$table)" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(filtered.data[tabbest\$best]))),best.logFC=tabbest\$logFC)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"csaw_si.csv\", as.data.frame(out.ranges), quote = FALSE, sep = \"\t\")" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "1_si $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1_si $MEMUSAGE" >> memory.txt
  
  #save result
  #reformat for eval take FDR and best log2-foldchange 
  if [ -e csaw_si.csv ]; then
    cut -f2,3,4,11,14 csaw_si.csv | grep -v "^start" | sort -k1,1 -k2,2n > $OUT_NAME"_1_si.bed"
  else
    #create empty file
    touch $OUT_NAME"_1_si.bed"
  fi
  
  #save log 
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT .RData script.Rout csaw_si.csv mem.txt
  
 
  # in R diffbind TMM, w = 400
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "bam.files<-c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
  echo "standard.chr <- paste0(\"chr\", c(19))" >> $SCRIPT
  echo "param <- readParam(minq=1, restrict=standard.chr)" >> $SCRIPT
  echo "celltype <- c(\"sample1\",\"sample1\",\"sample2\",\"sample2\")" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "min.fc <- 0.7" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=400, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type=\"global\")" >> $SCRIPT
  echo "keep1 <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- win.data[keep1,]" >> $SCRIPT
  echo "celltype <- factor(celltype)" >> $SCRIPT
  echo "design <- model.matrix(~0+celltype)" >> $SCRIPT
  echo "colnames(design) <- levels(celltype)" >> $SCRIPT 
  echo "binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)" >> $SCRIPT
  echo "filtered.data <- normOffsets(binned, se.out=filtered.data)" >> $SCRIPT
  echo "y <- asDGEList(filtered.data, norm.factors=filtered.data\$norm.factors)" >> $SCRIPT
  echo "y <- estimateDisp(y, design)" >> $SCRIPT
  echo "fit <- glmQLFit(y, design, robust=TRUE)" >> $SCRIPT
  echo "contrast <- makeContrasts(sample1-sample2, levels=design)" >> $SCRIPT
  echo "res <- glmQLFTest(fit, contrast=contrast)" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res\$table)" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res\$table)" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(filtered.data[tabbest\$best]))),best.logFC=tabbest\$logFC)" >> $SCRIPT
  echo "write.table(file=\"csaw_TMM.csv\", as.data.frame(out.ranges), quote = FALSE, sep = \"\t\")" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2_TMM $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2_TMM $MEMUSAGE" >> memory.txt
  
  #save result
  #reformat for eval take FDR and best log2-foldchange 
  if [ -e csaw_TMM.csv ]; then
    cut -f2,3,4,11,14 csaw_TMM.csv | grep -v "^start" | sort -k1,1 -k2,2n > $OUT_NAME"_2_TMM.bed"
  else
    #create empty file
    touch $OUT_NAME"_2_TMM.bed"
  fi
  
  #save log 
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT .RData script.Rout csaw_TMM.csv mem.txt
  
  
  
  # in R diffbind SPIKE-IN. w 400
  echo "library(rtracklayer)" > $SCRIPT
  echo "library(csaw)" >> $SCRIPT
  echo "library(edgeR)" >> $SCRIPT
  echo "bam.files<-c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
  echo "standard.chr <- paste0(\"chr\", c(19))" >> $SCRIPT
  echo "param <- readParam(minq=1, restrict=standard.chr)" >> $SCRIPT
  echo "celltype <- c(\"sample1\",\"sample1\",\"sample2\",\"sample2\")" >> $SCRIPT
  echo "x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))" >> $SCRIPT
  echo "frag.len <- maximizeCcf(x)" >> $SCRIPT
  echo "min.fc <- 0.7" >> $SCRIPT
  echo "win.data <- windowCounts(bam.files, param=param, width=400, ext=frag.len)" >> $SCRIPT
  echo "bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)" >> $SCRIPT
  echo "filter.stat <- filterWindows(win.data, bins, type=\"global\")" >> $SCRIPT
  echo "keep1 <- filter.stat\$filter > log2(min.fc)" >> $SCRIPT
  echo "filtered.data <- win.data[keep1,]" >> $SCRIPT
  echo "celltype <- factor(celltype)" >> $SCRIPT
  echo "design <- model.matrix(~0+celltype)" >> $SCRIPT
  echo "colnames(design) <- levels(celltype)" >> $SCRIPT
  echo "spike.facs <- c($10,$11,$12,$13)" >> $SCRIPT
  echo "y <- asDGEList(filtered.data, norm.factors=spike.facs)" >> $SCRIPT
  echo "y <- estimateDisp(y, design)" >> $SCRIPT
  echo "fit <- glmQLFit(y, design, robust=TRUE)" >> $SCRIPT
  echo "contrast <- makeContrasts(sample1-sample2, levels=design)" >> $SCRIPT
  echo "res <- glmQLFTest(fit, contrast=contrast)" >> $SCRIPT
  echo "merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)" >> $SCRIPT
  echo "tabcom <- combineTests(merged\$id, res\$table)" >> $SCRIPT
  echo "tabbest <- getBestTest(merged\$id, res\$table)" >> $SCRIPT
  echo "out.ranges <- merged\$region" >> $SCRIPT
  echo "elementMetadata(out.ranges) <- data.frame(tabcom,best.pos=mid(ranges(rowRanges(filtered.data[tabbest\$best]))),best.logFC=tabbest\$logFC)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file=\"csaw_si.csv\", as.data.frame(out.ranges), quote = FALSE, sep = \"\t\")" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2_si $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2_si $MEMUSAGE" >> memory.txt
  
  #save result
  #reformat for eval take FDR and best log2-foldchange 
  if [ -e csaw_si.csv ]; then
    cut -f2,3,4,11,14 csaw_si.csv | grep -v "^start" | sort -k1,1 -k2,2n > $OUT_NAME"_2_si.bed"
  else
    #create empty file
    touch $OUT_NAME"_2_si.bed"
  fi
  
  #save log 
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT .RData script.Rout csaw_si.csv mem.txt
   
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi

