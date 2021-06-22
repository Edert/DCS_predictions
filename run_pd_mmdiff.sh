
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
  
  #link bam files
  mkdir -p reads
  ln -s $4 reads/S11.bam
  ln -s $4.bai reads/S11.bam.bai
  ln -s $5 reads/S12.bam
  ln -s $5.bai reads/S12.bam.bai
  ln -s $6 reads/I1.bam
  ln -s $6.bai reads/I1.bam.bai
  ln -s $7 reads/S21.bam  
  ln -s $7.bai reads/S21.bam.bai 
  ln -s $8 reads/S22.bam
  ln -s $8.bai reads/S22.bam.bai
  ln -s $9 reads/I2.bam
  ln -s $9.bai reads/I2.bam.bai
  
  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt

  for PCALLER in ../../../results_peaks/*; do
  #for PCALLER in "../../../results_peaks/sicer"; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
      
      mkdir -p peaks 
      cat $PMODE/s11_peaks.bed > peaks/s11_peaks.bed
      cat $PMODE/s12_peaks.bed > peaks/s12_peaks.bed
      cat $PMODE/s21_peaks.bed > peaks/s21_peaks.bed
      cat $PMODE/s22_peaks.bed > peaks/s22_peaks.bed
     
      # in R diffbind
      echo "library(MMDiff)" > $SCRIPT
      echo "SampleID <- c(\"S11\",\"S12\",\"S21\",\"S22\")" >> $SCRIPT
      echo "Tissue <- c(\"NA\",\"NA\",\"NA\",\"NA\")" >> $SCRIPT
      echo "Factor <- c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
      echo "Condition <- c(1,1,2,2)" >> $SCRIPT
      echo "Replicate <- c(1,2,1,2)" >> $SCRIPT
      echo "bamReads <- c(\"reads/S11.bam\",\"reads/S12.bam\",\"reads/S21.bam\",\"reads/S22.bam\")" >> $SCRIPT
      echo "bamControl <- c(\"reads/I1.bam\",\"reads/I1.bam\",\"reads/I2.bam\",\"reads/I2.bam\")" >> $SCRIPT
      echo "Peaks <- c(\"peaks/s11_peaks.bed\",\"peaks/s12_peaks.bed\",\"peaks/s21_peaks.bed\",\"peaks/s22_peaks.bed\")" >> $SCRIPT
      echo "PeakCaller <- c(\"raw\",\"raw\",\"raw\",\"raw\")" >> $SCRIPT
      echo "myinputdata <- data.frame(SampleID,Tissue,Factor,Condition,Replicate,bamReads,bamControl,Peaks,PeakCaller,stringsAsFactors=FALSE)" >> $SCRIPT
      echo "mydata <- dba(sampleSheet=myinputdata, minOverlap=2)" >> $SCRIPT
      echo "Peaks <- dba.peakset(mydata,bRetrieve=TRUE)" >> $SCRIPT
      #echo "Profiles <- getPeakProfiles(mydata, Peaks,bin.length=50, save.files=FALSE, draw.on=FALSE, run.parallel=FALSE)" >> $SCRIPT
      echo "Profiles <- getPeakProfiles(mydata, Peaks,bin.length=50, save.files=FALSE, draw.on=FALSE)" >> $SCRIPT
      echo "Norm <- getNormFactors(Profiles, method = 'DESeq', SampleIDs=SampleID)" >> $SCRIPT
      #echo "Dists <- compHistDists(Norm, method='MMD', overWrite=TRUE, NormMethod='DESeq', run.parallel=FALSE, save.file=FALSE)" >> $SCRIPT
      #MMD failed for larger samples eg. set3_1 
      echo "Dists <- compHistDists(Norm, method='GMD', overWrite=TRUE, NormMethod='DESeq', save.file=FALSE)" >> $SCRIPT
      echo "group1 <- c(\"S11\",\"S12\")" >> $SCRIPT
      echo "group2 <- c(\"S21\",\"S22\")" >> $SCRIPT
      #echo "Pvals <- detPeakPvals(Dists, group1=group1, group2=group2, name1='s1_s2', name2='Null', method='MMD')" >> $SCRIPT
      echo "Pvals <- detPeakPvals(Dists, group1=group1, group2=group2, name1='s1_s2', name2='Null', method='GMD')" >> $SCRIPT
      #echo "out <- Pvals\$MD\$Pvals\$MMD[,7,drop=FALSE]" >> $SCRIPT
      echo "out <- Pvals\$MD\$Pvals\$GMD[,7,drop=FALSE]" >> $SCRIPT
      echo "out <- out[complete.cases(out), ,drop=FALSE]" >> $SCRIPT
      echo "final <- out[out < 0.05,,drop=FALSE]" >> $SCRIPT
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(final, file=\"results.csv\", sep='\t', col.names=FALSE, row.names=TRUE, quote = FALSE)" >> $SCRIPT
      
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
        cat results.csv | awk '{split($1,chr,":");split(chr[2],pos,"-"); print chr[1]"\t"pos[1]"\t"pos[2]"\t"$2"\t1"}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
        #create empty file
        touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout results.csv .Dists.Rdata .cluster.log mem.txt
      rm -rf peaks

    done
  done
  
  #clean up
  rm -rf reads 
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
