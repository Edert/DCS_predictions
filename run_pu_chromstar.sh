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
  
  CURDIR=$(pwd)
  INDIR=$(dirname $4)
  S11=$(basename $4)
  S12=$(basename $5)
  I1=$(basename $6)
  S21=$(basename $7)
  S22=$(basename $8)
  I2=$(basename $9)
  
  # in R bin 100
  echo "library(chromstaR)" > $SCRIPT
  echo "file<-c('$S11','$S12','$S21','$S22')" >> $SCRIPT
  echo "mark<-c('sim','sim','sim','sim')" >> $SCRIPT
  echo "condition<-c('S1','S1','S2','S2')" >> $SCRIPT
  echo "replicate<-c('1','2','1','2')" >> $SCRIPT
  echo "pairedEndReads <- c(FALSE,FALSE,FALSE,FALSE)" >> $SCRIPT
  echo "controlFiles<-c('$I1','$I1','$I2','$I2')" >> $SCRIPT
  echo "experiment_table <- data.frame(file,mark,condition,replicate,pairedEndReads,controlFiles)" >> $SCRIPT
  echo "inputfolder <- '$INDIR'" >> $SCRIPT
  echo "outputfolder <- '$CURDIR'" >> $SCRIPT
  echo "chrominfo <- data.frame(chromosome=c('chr19'),length=c(61431566))" >> $SCRIPT
  echo "Chromstar(inputfolder, experiment.table=experiment_table,outputfolder=outputfolder, numCPU=1, binsize=100, stepsize=50,assembly=chrominfo, prefit.on.chr='chr19', chromosomes='chr19',mode='differential')" >> $SCRIPT 
  echo "model <- get(load(file.path(outputfolder,'multivariate','multivariate_mode-differential_mark-sim_binsize100_stepsize50.RData')))" >> $SCRIPT
  echo "diff.states <- stateBrewer(experiment_table, mode='differential',differential.states=TRUE)" >> $SCRIPT
  echo "exportCombinations(model, filename='results', include.states=diff.states,header=F)" >> $SCRIPT
  
  PREPDONE=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M"  R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e results_combinations.bed.gz ]; then
     gunzip results_combinations.bed.gz
     #reformat for eval
     cat results_combinations.bed | awk '{fold=(($4 == "[S1]") ? 0.7 : -0.7); print $1"\t"$2"\t"$3"\t"0"\t"fold}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT script.Rout results_combinations.bed chrominfo.tsv chromstaR.config experiment_table.tsv README.txt mem.txt
  rm -rf binned BROWSERFILES combined multivariate PLOTS univariate
  
  
  
  # in R bin 1000
  echo "library(chromstaR)" > $SCRIPT
  echo "file<-c('$S11','$S12','$S21','$S22')" >> $SCRIPT
  echo "mark<-c('sim','sim','sim','sim')" >> $SCRIPT
  echo "condition<-c('S1','S1','S2','S2')" >> $SCRIPT
  echo "replicate<-c('1','2','1','2')" >> $SCRIPT
  echo "pairedEndReads <- c(FALSE,FALSE,FALSE,FALSE)" >> $SCRIPT
  echo "controlFiles<-c('$I1','$I1','$I2','$I2')" >> $SCRIPT
  echo "experiment_table <- data.frame(file,mark,condition,replicate,pairedEndReads,controlFiles)" >> $SCRIPT
  echo "inputfolder <- '$INDIR'" >> $SCRIPT
  echo "outputfolder <- '$CURDIR'" >> $SCRIPT
  echo "chrominfo <- data.frame(chromosome=c('chr19'),length=c(61431566))" >> $SCRIPT
  echo "Chromstar(inputfolder, experiment.table=experiment_table,outputfolder=outputfolder, numCPU=1, binsize=1000, stepsize=500,assembly=chrominfo, prefit.on.chr='chr19', chromosomes='chr19',mode='differential')" >> $SCRIPT 
  echo "model <- get(load(file.path(outputfolder,'multivariate','multivariate_mode-differential_mark-sim_binsize1000_stepsize500.RData')))" >> $SCRIPT
  echo "diff.states <- stateBrewer(experiment_table, mode='differential',differential.states=TRUE)" >> $SCRIPT
  echo "exportCombinations(model, filename='results', include.states=diff.states,header=F)" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e results_combinations.bed.gz ]; then
     gunzip results_combinations.bed.gz
     #reformat for eval
     cat results_combinations.bed | awk '{fold=(($4 == "[S1]") ? 0.7 : -0.7); print $1"\t"$2"\t"$3"\t"0"\t"fold}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  
  #clean up
  rm -f $SCRIPT script.Rout results_combinations.bed chrominfo.tsv chromstaR.config experiment_table.tsv README.txt mem.txt
  rm -rf binned BROWSERFILES combined multivariate PLOTS univariate
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi

