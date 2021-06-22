
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
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
      
      mkdir -p peaks 
      cat $PMODE/s11_peaks.bed > peaks/s11_peaks.bed
      cat $PMODE/s12_peaks.bed > peaks/s12_peaks.bed
      cat $PMODE/s21_peaks.bed > peaks/s21_peaks.bed
      cat $PMODE/s22_peaks.bed > peaks/s22_peaks.bed
      
      #create config
      echo "SampleID,ControlID,Tissue,Factor,Condition,Replicate,bamReads,bamControl,Peaks,PeakCaller" > mmdiff2.csv
      echo "" >> mmdiff2.csv
      echo "S11,I1,NA,sim,1,1,reads/S11.bam,reads/I1.bam,peaks/s11_peaks.bed,raw" >> mmdiff2.csv
      echo "S12,I1,NA,sim,1,2,reads/S12.bam,reads/I1.bam,peaks/s12_peaks.bed,raw" >> mmdiff2.csv
      echo "S21,I2,NA,sim,2,1,reads/S21.bam,reads/I2.bam,peaks/s21_peaks.bed,raw" >> mmdiff2.csv
      echo "S22,I2,NA,sim,2,2,reads/S22.bam,reads/I2.bam,peaks/s22_peaks.bed,raw" >> mmdiff2.csv
      
      # in R 
      echo "library(MMDiff2)" > $SCRIPT
      echo "ExperimentData <- list(genome='BSgenome.Mmusculus.UCSC.mm10',dataDir='.',sampleSheet='mmdiff2.csv')" >> $SCRIPT
      echo "MetaData <- list('ExpData' = ExperimentData)" >> $SCRIPT
      echo "MMD <- DBAmmd(MetaData)" >> $SCRIPT
      echo "library('DiffBind')" >> $SCRIPT
      echo "DBA <- dba(sampleSheet='mmdiff2.csv', minOverlap=2)" >> $SCRIPT
      echo "Peaks <- dba.peakset(DBA, bRetrieve = TRUE)" >> $SCRIPT
      echo "MMD <- setRegions(MMD,Peaks)" >> $SCRIPT
      echo "MMD <- getPeakReads(MMD, pairedEnd = FALSE, run.parallel = FALSE)" >> $SCRIPT
      echo "MMD <- estimateFragmentCenters(MMD)" >> $SCRIPT
      echo "MMD <- compDists(MMD)" >> $SCRIPT
      echo "MMD <- setContrast(MMD,contrast='byCondition')" >> $SCRIPT
      echo "MMD <- compPvals(MMD)" >> $SCRIPT
      echo "res <- reportResults(MMD)" >> $SCRIPT
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(res, file=\"results.csv\", sep='\t', col.names=FALSE, row.names=TRUE, quote = FALSE)" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
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
        cat results.csv | awk '{m1=($7+$8)/2; m2=($9+$10)/2; lfold=((m2==0) ? "inf" : log(m1/m2)/log(2));print $2"\t"$3"\t"$4"\t"$11"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME

      else
        #create empty file
        touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG

      #clean up
      rm -f $SCRIPT script.Rout results.csv mmdiff2.csv mem.txt
      rm -rf peaks 

    done
  done
  
  #clean up
  rm -rf reads 
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
