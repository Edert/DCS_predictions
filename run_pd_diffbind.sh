
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
  
  DBMETHODS="DBA_EDGER DBA_DESEQ2"
  #DBA_DESEQ DBA_EDGER_CLASSIC DBA_EDGER_GLM DBA_DESEQ_CLASSIC DBA_DESEQ_GLM DBA_EDGER_BLOCK DBA_DESEQ_BLOCK DBA_DESEQ2_BLOCK
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt

  for PCALLER in ../../../results_peaks/*; do
  #for PCALLER in "../../../results_peaks/sicer"; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
           
      cat $PMODE/s11_peaks.bed > s11_peaks.bed
      cat $PMODE/s12_peaks.bed > s12_peaks.bed
      cat $PMODE/s21_peaks.bed > s21_peaks.bed
      cat $PMODE/s22_peaks.bed > s22_peaks.bed

      for DBMETHOD in $DBMETHODS; do #loop through DB methods
        #echo "working with method: $DBMETHOD"
      
        # in R: bFullLibrarySize=F bSubControl=F
        echo "library(DiffBind)" > $SCRIPT
        echo "library(locfit)" >> $SCRIPT
        echo "SampleID=c(\"1_1\",\"1_2\",\"2_1\",\"2_2\")" >> $SCRIPT
        echo "Tissue=c(\"NA\",\"NA\",\"NA\",\"NA\")" >> $SCRIPT
        echo "Factor=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Condition=c(1,1,2,2)" >> $SCRIPT
        echo "Treatment=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Replicate=c(1,2,1,2)" >> $SCRIPT
        echo "bamReads=c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
        echo "ControlID=c(\"c1\",\"c1\",\"c2\",\"c2\")" >> $SCRIPT
        echo "bamControl=c(\"$6\",\"$6\",\"$8\",\"$8\")" >> $SCRIPT
        echo "Peaks=c(\"s11_peaks.bed\",\"s12_peaks.bed\",\"s21_peaks.bed\",\"s22_peaks.bed\")" >> $SCRIPT
        echo "PeakCaller=c(\"raw\",\"raw\",\"raw\",\"raw\")" >> $SCRIPT
        echo "myinputdata <- data.frame(SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller)" >> $SCRIPT
        echo "mydata = dba(sampleSheet=myinputdata)" >> $SCRIPT
        echo "mydata = dba.count(mydata, minOverlap=1)" >> $SCRIPT
        echo "mydata = dba.contrast(mydata, categories=DBA_CONDITION,minMembers=2)" >> $SCRIPT 
        echo "mydata = dba.analyze(mydata,method=$DBMETHOD,bFullLibrarySize=F,bSubControl=F)" >> $SCRIPT
        
        #save as csv..
        echo "options(scipen = 999)" >> $SCRIPT
        echo "rep = dba.report(mydata, bCounts=T, bCalled=T, bCalledDetail=T, th=1, method=$DBMETHOD, file=\"diff_peaks\")" >> $SCRIPT
        
        STARTTIME=`date +%s.%N`
        
        #run it...
        #R CMD BATCH --vanilla $SCRIPT
        /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
        
        #save result
        PSHORT=$(basename $PCALLER)
        MSHORT=$(basename $PMODE)
        OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_"$DBMETHOD"_1.bed"
        
        ENDTIME=`date +%s.%N`
        TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_1 $TIMEDIFF" >> time.txt

        MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_1 $MEMUSAGE" >> memory.txt
      
        if [ -e DBA_diff_peaks.csv ]; then
            #reformat for eval
            cut -f1,2,3,7,9 -d, DBA_diff_peaks.csv | grep -v "^\"Chr\",\"Start\"" | sed "s/\"//g" | sed "s/,/\t/g" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | sort -k1,1 -k2,2n > $OUT_NAME
        else
          #create empty file
          touch $OUT_NAME
        fi
        
        #save log     
        cat script.Rout >> $LOG
        
        #clean up
        rm -f $SCRIPT script.Rout .RData DBA_diff_peaks.csv mem.txt
        
        
        
        # bFullLibrarySize=T bSubControl=F
        echo "library(DiffBind)" > $SCRIPT
        echo "library(locfit)" >> $SCRIPT
        echo "SampleID=c(\"1_1\",\"1_2\",\"2_1\",\"2_2\")" >> $SCRIPT
        echo "Tissue=c(\"NA\",\"NA\",\"NA\",\"NA\")" >> $SCRIPT
        echo "Factor=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Condition=c(1,1,2,2)" >> $SCRIPT
        echo "Treatment=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Replicate=c(1,2,1,2)" >> $SCRIPT
        echo "bamReads=c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
        echo "ControlID=c(\"c1\",\"c1\",\"c2\",\"c2\")" >> $SCRIPT
        echo "bamControl=c(\"$6\",\"$6\",\"$8\",\"$8\")" >> $SCRIPT
        echo "Peaks=c(\"s11_peaks.bed\",\"s12_peaks.bed\",\"s21_peaks.bed\",\"s22_peaks.bed\")" >> $SCRIPT
        echo "PeakCaller=c(\"raw\",\"raw\",\"raw\",\"raw\")" >> $SCRIPT
        echo "myinputdata <- data.frame(SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller)" >> $SCRIPT
        echo "mydata = dba(sampleSheet=myinputdata)" >> $SCRIPT
        echo "mydata = dba.count(mydata, minOverlap=1)" >> $SCRIPT
        echo "mydata = dba.contrast(mydata, categories=DBA_CONDITION,minMembers=2)" >> $SCRIPT 
        echo "mydata = dba.analyze(mydata,method=$DBMETHOD,bFullLibrarySize=T,bSubControl=F)" >> $SCRIPT
        
        #save as csv..
        echo "options(scipen = 999)" >> $SCRIPT
        echo "rep = dba.report(mydata, bCounts=T, bCalled=T, bCalledDetail=T, th=1, method=$DBMETHOD, file=\"diff_peaks\")" >> $SCRIPT
        
        STARTTIME=`date +%s.%N`
        
        #run it...
        #R CMD BATCH --vanilla $SCRIPT
        /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
        
        #save result
        PSHORT=$(basename $PCALLER)
        MSHORT=$(basename $PMODE)
        OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_"$DBMETHOD"_2.bed"
        
        ENDTIME=`date +%s.%N`
        TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_2 $TIMEDIFF" >> time.txt
        
        MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_2 $MEMUSAGE" >> memory.txt
        
        if [ -e DBA_diff_peaks.csv ]; then
            #reformat for eval
            cut -f1,2,3,7,9 -d, DBA_diff_peaks.csv | grep -v "^\"Chr\",\"Start\"" | sed "s/\"//g" | sed "s/,/\t/g" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | sort -k1,1 -k2,2n > $OUT_NAME
        else
          #create empty file
          touch $OUT_NAME
        fi
        
        #save log     
        cat script.Rout >> $LOG
        
        #clean up
        rm -f $SCRIPT script.Rout .RData DBA_diff_peaks.csv mem.txt
        
        
        
        # bFullLibrarySize=F bSubControl=T
        echo "library(DiffBind)" > $SCRIPT
        echo "library(locfit)" >> $SCRIPT
        echo "SampleID=c(\"1_1\",\"1_2\",\"2_1\",\"2_2\")" >> $SCRIPT
        echo "Tissue=c(\"NA\",\"NA\",\"NA\",\"NA\")" >> $SCRIPT
        echo "Factor=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Condition=c(1,1,2,2)" >> $SCRIPT
        echo "Treatment=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Replicate=c(1,2,1,2)" >> $SCRIPT
        echo "bamReads=c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
        echo "ControlID=c(\"c1\",\"c1\",\"c2\",\"c2\")" >> $SCRIPT
        echo "bamControl=c(\"$6\",\"$6\",\"$8\",\"$8\")" >> $SCRIPT
        echo "Peaks=c(\"s11_peaks.bed\",\"s12_peaks.bed\",\"s21_peaks.bed\",\"s22_peaks.bed\")" >> $SCRIPT
        echo "PeakCaller=c(\"raw\",\"raw\",\"raw\",\"raw\")" >> $SCRIPT
        echo "myinputdata <- data.frame(SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller)" >> $SCRIPT
        echo "mydata = dba(sampleSheet=myinputdata)" >> $SCRIPT
        echo "mydata = dba.count(mydata, minOverlap=1)" >> $SCRIPT
        echo "mydata = dba.contrast(mydata, categories=DBA_CONDITION,minMembers=2)" >> $SCRIPT 
        echo "mydata = dba.analyze(mydata,method=$DBMETHOD,bFullLibrarySize=F,bSubControl=T)" >> $SCRIPT
        
        #save as csv..
        echo "options(scipen = 999)" >> $SCRIPT
        echo "rep = dba.report(mydata, bCounts=T, bCalled=T, bCalledDetail=T, th=1, method=$DBMETHOD, file=\"diff_peaks\")" >> $SCRIPT
        
        STARTTIME=`date +%s.%N`
        
        #run it...
        #R CMD BATCH --vanilla $SCRIPT
        /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
        
        #save result
        PSHORT=$(basename $PCALLER)
        MSHORT=$(basename $PMODE)
        OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_"$DBMETHOD"_3.bed"
        
        ENDTIME=`date +%s.%N`
        TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_3 $TIMEDIFF" >> time.txt
        
        MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_3 $MEMUSAGE" >> memory.txt
        
        if [ -e DBA_diff_peaks.csv ]; then
            #reformat for eval
            cut -f1,2,3,7,9 -d, DBA_diff_peaks.csv | grep -v "^\"Chr\",\"Start\"" | sed "s/\"//g" | sed "s/,/\t/g" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | sort -k1,1 -k2,2n > $OUT_NAME
        else
          #create empty file
          touch $OUT_NAME
        fi
        
        #save log     
        cat script.Rout >> $LOG
        
        #clean up
        rm -f $SCRIPT script.Rout .RData DBA_diff_peaks.csv mem.txt
        
        
        
        # bFullLibrarySize=T bSubControl=T
        echo "library(DiffBind)" > $SCRIPT
        echo "library(locfit)" >> $SCRIPT
        echo "SampleID=c(\"1_1\",\"1_2\",\"2_1\",\"2_2\")" >> $SCRIPT
        echo "Tissue=c(\"NA\",\"NA\",\"NA\",\"NA\")" >> $SCRIPT
        echo "Factor=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Condition=c(1,1,2,2)" >> $SCRIPT
        echo "Treatment=c(\"sim\",\"sim\",\"sim\",\"sim\")" >> $SCRIPT
        echo "Replicate=c(1,2,1,2)" >> $SCRIPT
        echo "bamReads=c(\"$4\",\"$5\",\"$7\",\"$8\")" >> $SCRIPT
        echo "ControlID=c(\"c1\",\"c1\",\"c2\",\"c2\")" >> $SCRIPT
        echo "bamControl=c(\"$6\",\"$6\",\"$8\",\"$8\")" >> $SCRIPT
        echo "Peaks=c(\"s11_peaks.bed\",\"s12_peaks.bed\",\"s21_peaks.bed\",\"s22_peaks.bed\")" >> $SCRIPT
        echo "PeakCaller=c(\"raw\",\"raw\",\"raw\",\"raw\")" >> $SCRIPT
        echo "myinputdata <- data.frame(SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller)" >> $SCRIPT
        echo "mydata = dba(sampleSheet=myinputdata)" >> $SCRIPT
        echo "mydata = dba.count(mydata, minOverlap=1)" >> $SCRIPT
        echo "mydata = dba.contrast(mydata, categories=DBA_CONDITION,minMembers=2)" >> $SCRIPT 
        echo "mydata = dba.analyze(mydata,method=$DBMETHOD,bFullLibrarySize=T,bSubControl=T)" >> $SCRIPT
        
        #save as csv..
        echo "options(scipen = 999)" >> $SCRIPT
        echo "rep = dba.report(mydata, bCounts=T, bCalled=T, bCalledDetail=T, th=1, method=$DBMETHOD, file=\"diff_peaks\")" >> $SCRIPT
        
        STARTTIME=`date +%s.%N`
        
        #run it...
        #R CMD BATCH --vanilla $SCRIPT
        /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
        
        #save result
        PSHORT=$(basename $PCALLER)
        MSHORT=$(basename $PMODE)
        OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_"$DBMETHOD"_4.bed"
        
        ENDTIME=`date +%s.%N`
        TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_4 $TIMEDIFF" >> time.txt
        
        MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
        echo $PSHORT"_"$MSHORT"_"$DBMETHOD"_4 $MEMUSAGE" >> memory.txt
        
        if [ -e DBA_diff_peaks.csv ]; then
            #reformat for eval
            cut -f1,2,3,7,9 -d, DBA_diff_peaks.csv | grep -v "^\"Chr\",\"Start\"" | sed "s/\"//g" | sed "s/,/\t/g" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | sort -k1,1 -k2,2n > $OUT_NAME
        else
          #create empty file
          touch $OUT_NAME
        fi
        
        #save log     
        cat script.Rout >> $LOG
        
        #clean up
        rm -f $SCRIPT script.Rout .RData DBA_diff_peaks.csv 
      done
      
      #clean up
      rm -rf s11_peaks.bed s12_peaks.bed s21_peaks.bed s22_peaks.bed mem.txt
       
    done
  done

else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
