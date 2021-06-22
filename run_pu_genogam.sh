NAME=$1
SET=$2
TOOL=$3
LENchr19=61431566
WORKER=5

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
  
  # in R 
  echo "start_time <- Sys.time()" > $SCRIPT
  echo "library(GenoGAM)" >> $SCRIPT
  #echo "BiocParallel::register(BiocParallel::SnowParam(worker = 5))" >> $SCRIPT  
  echo "BiocParallel::register(BiocParallel::MulticoreParam(worker = $WORKER))" >> $SCRIPT  
  echo "ID<-c('S11','S12','S21','S22')" >> $SCRIPT
  echo "file<-c('$4','$5','$7','$8')" >> $SCRIPT
  echo "paired <- c(FALSE,FALSE,FALSE,FALSE)" >> $SCRIPT
  echo "S2vsS1<-c(0,0,1,1)" >> $SCRIPT
  echo "S1vsS2<-c(1,1,0,0)" >> $SCRIPT
  echo "expDesign <- data.frame(ID,file,paired,S2vsS1,S1vsS2)" >> $SCRIPT
  echo "chunkSize <- 75000" >> $SCRIPT
  echo "overhangSize <- 1000" >> $SCRIPT
  echo "design = ~s(x, by = S2vsS1) + s(x, by = S1vsS2)" >> $SCRIPT
  echo "settings <- GenoGAMSettings(chromosomeList = 'chr19')" >> $SCRIPT
  #echo "ggd <- GenoGAMDataSet(expDesign, directory = '/', chunkSize = chunkSize, overhangSize = overhangSize, design=design, settings=settings, hdf5 = TRUE)" >> $SCRIPT
  echo "ggd <- GenoGAMDataSet(expDesign, directory = '/', chunkSize = chunkSize, overhangSize = overhangSize, design=design, settings=settings, hdf5 = F, split=TRUE)" >> $SCRIPT
  echo "ggd <- computeSizeFactors(ggd)" >> $SCRIPT
  echo "result <- genogam(ggd)" >> $SCRIPT
  #echo "computeSignificance(result)" >> $SCRIPT
  echo "end_time <- Sys.time()" >> $SCRIPT
  #S1
  echo "peaks <- callPeaks(result)" >> $SCRIPT
  echo "peaks1 <- data.frame(peaks\$\`s(x):S2vsS1\`\$chromosome,peaks\$\`s(x):S2vsS1\`\$pos-100,peaks\$\`s(x):S2vsS1\`\$pos+100,peaks\$\`s(x):S2vsS1\`\$fdr)" >> $SCRIPT
  echo "colnames(peaks1) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "gpeaks1 <- makeGRangesFromDataFrame(peaks1)" >> $SCRIPT
  echo "peaks2 <- data.frame(peaks\$\`s(x):S1vsS2\`\$chromosome,peaks\$\`s(x):S1vsS2\`\$pos-100,peaks\$\`s(x):S1vsS2\`\$pos+100,peaks\$\`s(x):S1vsS2\`\$fdr)" >> $SCRIPT
  echo "colnames(peaks2) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "gpeaks2 <- makeGRangesFromDataFrame(peaks2)" >> $SCRIPT
  echo "gpeaks <- c(gpeaks1,gpeaks2)" >> $SCRIPT
  echo "rgpeaks <- reduce(gpeaks)" >> $SCRIPT
  echo "res1 <- computeRegionSignificance(result,rgpeaks)" >> $SCRIPT
  echo "resout <- data.frame(data.frame(res1\$\`s(x):S2vsS1\`)\$seqnames,data.frame(res1\$\`s(x):S2vsS1\`)\$start,data.frame(res1\$\`s(x):S2vsS1\`)\$end,pmin(data.frame(res1\$\`s(x):S2vsS1\`)\$FDR, data.frame(res1\$\`s(x):S1vsS2\`)\$FDR))" >> $SCRIPT
  echo "colnames(resout) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='mypeaks_s1.txt', resout,quote = FALSE, sep ='	',row.names = F,col.names = F)" >> $SCRIPT
  echo "end_times1 <- Sys.time()" >> $SCRIPT
  #S2
  echo "peaks <- callPeaks(result,peakType = 'broad',maxgap=100)" >> $SCRIPT
  echo "peaks1 <- data.frame(peaks\$\`s(x):S2vsS1\`\$seqnames,peaks\$\`s(x):S2vsS1\`\$start, peaks\$\`s(x):S2vsS1\`\$end,peaks\$\`s(x):S2vsS1\`\$fdr)" >> $SCRIPT
  echo "colnames(peaks1) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "gpeaks1 <- makeGRangesFromDataFrame(peaks1)" >> $SCRIPT
  echo "peaks2 <- data.frame(peaks\$\`s(x):S1vsS2\`\$seqnames,peaks\$\`s(x):S1vsS2\`\$start, peaks\$\`s(x):S1vsS2\`\$end,peaks\$\`s(x):S1vsS2\`\$fdr)" >> $SCRIPT
  echo "colnames(peaks2) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "gpeaks2 <- makeGRangesFromDataFrame(peaks2)" >> $SCRIPT
  echo "gpeaks <- c(gpeaks1,gpeaks2)" >> $SCRIPT
  echo "rgpeaks <- reduce(gpeaks)" >> $SCRIPT
  echo "res2 <- computeRegionSignificance(result,rgpeaks)" >> $SCRIPT
  echo "resout <- data.frame(data.frame(res2\$\`s(x):S2vsS1\`)\$seqnames,data.frame(res2\$\`s(x):S2vsS1\`)\$start,data.frame(res2\$\`s(x):S2vsS1\`)\$end,pmin(data.frame(res2\$\`s(x):S2vsS1\`)\$FDR, data.frame(res2\$\`s(x):S1vsS2\`)\$FDR))" >> $SCRIPT
  echo "colnames(resout) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='mypeaks_s2.txt', resout,quote = FALSE, sep ='	',row.names = F,col.names = F)" >> $SCRIPT
  echo "end_times2 <- Sys.time()" >> $SCRIPT
  #S3
  echo "peaks <- callPeaks(result,peakType = 'broad', maxgap=10000)" >> $SCRIPT
  echo "peaks1 <- data.frame(peaks\$\`s(x):S2vsS1\`\$seqnames,peaks\$\`s(x):S2vsS1\`\$start, peaks\$\`s(x):S2vsS1\`\$end,peaks\$\`s(x):S2vsS1\`\$fdr)" >> $SCRIPT
  echo "colnames(peaks1) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "gpeaks1 <- makeGRangesFromDataFrame(peaks1)" >> $SCRIPT
  echo "peaks2 <- data.frame(peaks\$\`s(x):S1vsS2\`\$seqnames,peaks\$\`s(x):S1vsS2\`\$start, peaks\$\`s(x):S1vsS2\`\$end,peaks\$\`s(x):S1vsS2\`\$fdr)" >> $SCRIPT
  echo "colnames(peaks2) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "gpeaks2 <- makeGRangesFromDataFrame(peaks2)" >> $SCRIPT
  echo "gpeaks <- c(gpeaks1,gpeaks2)" >> $SCRIPT
  echo "rgpeaks <- reduce(gpeaks)" >> $SCRIPT
  echo "res3 <- computeRegionSignificance(result,rgpeaks)" >> $SCRIPT
  echo "resout <- data.frame(data.frame(res3\$\`s(x):S2vsS1\`)\$seqnames, data.frame(res3\$\`s(x):S2vsS1\`)\$start,data.frame(res3\$\`s(x):S2vsS1\`)\$end,pmin(data.frame(res3\$\`s(x):S2vsS1\`)\$FDR, data.frame(res3\$\`s(x):S1vsS2\`)\$FDR))" >> $SCRIPT
  echo "colnames(resout) <- c('chromosome','start','end','fdr')" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='mypeaks_s3.txt', resout,quote = FALSE, sep ='	',row.names = F,col.names = F)" >> $SCRIPT
  echo "end_times3 <- Sys.time()" >> $SCRIPT
  #calc timings
  echo "prep <- difftime(end_time, start_time, units = 's')" >> $SCRIPT
  echo "s1 <- prep + difftime(end_times1, end_time, units = 's')" >> $SCRIPT
  echo "s2 <- prep + difftime(end_times2, end_times1, units = 's')" >> $SCRIPT
  echo "s3 <- prep + difftime(end_times3, end_times2, units = 's')" >> $SCRIPT
  echo "write.table(file='time.txt',append = T,rbind(cbind('s1',s1),cbind('s2',s2),cbind('s3',s3)),sep = ' ',quote = F,row.names = F,col.names = F)" >> $SCRIPT
  
  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  #the remainign timings will be done in R
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla script.R
  
  #save result s1
  if [ -e mypeaks_s1.txt ]; then
    #move results 
    cat mypeaks_s1.txt | grep -v "^start" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"1}' | sort -k1,1 -k2,2n > $OUT_NAME"_s1.bed"
  else
    #create empty file
    touch $OUT_NAME"_s1.bed"
  fi

  #save result s2
  if [ -e mypeaks_s2.txt ]; then
    #move results
    cat mypeaks_s2.txt | grep -v "^start" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"1}' | sort -k1,1 -k2,2n > $OUT_NAME"_s2.bed"
  else
    #create empty file
    touch $OUT_NAME"_s2.bed"
  fi
  
  #save result s3
  if [ -e mypeaks_s3.txt ]; then
    #move results
    cat mypeaks_s3.txt | grep -v "^start" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"1}' | sort -k1,1 -k2,2n > $OUT_NAME"_s3.bed"
  else
    #create empty file
    touch $OUT_NAME"_s3.bed"
  fi
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "s1 $MEMUSAGE" >> memory.txt
  echo "s2 $MEMUSAGE" >> memory.txt
  echo "s3 $MEMUSAGE" >> memory.txt
  
  #save log 
  cat script.Rout >> $LOG
 
  #clean up
  rm -f $SCRIPT .RData script.Rout mypeaks_s1.txt mypeaks_s2.txt mypeaks_s3.txt mem.txt
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
