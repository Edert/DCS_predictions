
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
  
  bamToBed -i $4 > S11.bed
  bamToBed -i $5 > S12.bed
  bamToBed -i $6 > IN1.bed
  bamToBed -i $7 > S21.bed
  bamToBed -i $8 > S22.bed
  bamToBed -i $9 > IN2.bed
  
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
      
      # in R diffbind
      echo "library(DBChIP)" > $SCRIPT
      echo "library(ShortRead)" >> $SCRIPT
      
      echo "treatments.to.run = c(\"S11.bed\",\"S12.bed\",\"S21.bed\",\"S22.bed\")" >> $SCRIPT
      echo "input.to.run = c(\"IN1.bed\",\"IN1.bed\",\"IN2.bed\",\"IN2.bed\")" >> $SCRIPT
      
      echo "treatment.label = c('s11','s12','s21','s22')" >> $SCRIPT
      echo "output.bed.file=c(\"s11_peaks.bed\",\"s12_peaks.bed\",\"s21_peaks.bed\",\"s22_peaks.bed\")" >> $SCRIPT
      echo "peaks.bed.list = list()" >> $SCRIPT
      echo "for (i in 1:length(treatment.label)) { #get chr and midpoint" >> $SCRIPT
      echo "    tmp = read.table(file = output.bed.file[i], header = F)" >> $SCRIPT
      echo "    tmp = data.frame(tmp[, 1], floor(tmp[, 2]+(tmp[, 3]-tmp[, 2])/2))" >> $SCRIPT
      echo "    colnames(tmp) = c('chr', 'pos')" >> $SCRIPT
      echo "    peaks.bed.list[[treatment.label[i]]] = tmp" >> $SCRIPT
      echo "}" >> $SCRIPT
      
      echo "bed.data.list = list()" >> $SCRIPT
      echo "for (i in 1:length(treatment.label)) {" >> $SCRIPT
      echo "    bed.data.list[[treatment.label[i]]] = treatments.to.run[i]" >> $SCRIPT
      echo "}" >> $SCRIPT
      
      echo "input.data.list = list()" >> $SCRIPT
      echo "for (i in 1:length(treatment.label)) {" >> $SCRIPT
      echo "    input.data.list[[treatment.label[i]]] = input.to.run[i]" >> $SCRIPT
      echo "}" >> $SCRIPT
      
      echo "ana.conds = factor(c(rep(\"S1\", 2), rep(\"S2\", 2)), levels = c(\"S1\", \"S2\"))" >> $SCRIPT
      echo "ana.conds.label = \"S1.vs.S2\"" >> $SCRIPT
      
      echo "bs.list = read.binding.site.list(peaks.bed.list)" >> $SCRIPT
      echo "consensus.site = site.merge(bs.list)" >> $SCRIPT
      echo "dat = load.data(chip.data.list = bed.data.list, input.data.list=input.data.list, conds = ana.conds, consensus.site = consensus.site,  data.type = \"BED\")" >> $SCRIPT
      
      echo "data.count = get.site.count(dat)" >> $SCRIPT
      echo "test.diff = test.diff.binding(data.count)" >> $SCRIPT
      echo "rept = report.peak(test.diff,FDR=1)" >> $SCRIPT
      
      #save as csv..
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file=\"results.csv\", as.data.frame(rept),quote = FALSE)" >> $SCRIPT
      
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
          cat results.csv | grep -v "^chr"  | awk '{print $2"\t"$3-100"\t"$3+100"\t"$9"\t"log($7)/log(2)}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
        #create empty file
        touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      
      #run with scaling factors
      echo "library(DBChIP)" > $SCRIPT
      echo "library(ShortRead)" >> $SCRIPT
      
      echo "treatments.to.run = c(\"S11.bed\",\"S12.bed\",\"S21.bed\",\"S22.bed\")" >> $SCRIPT
      echo "input.to.run = c(\"IN1.bed\",\"IN1.bed\",\"IN2.bed\",\"IN2.bed\")" >> $SCRIPT
      
      echo "treatment.label = c('s11','s12','s21','s22')" >> $SCRIPT
      echo "output.bed.file=c(\"s11_peaks.bed\",\"s12_peaks.bed\",\"s21_peaks.bed\",\"s22_peaks.bed\")" >> $SCRIPT
      echo "peaks.bed.list = list()" >> $SCRIPT
      echo "for (i in 1:length(treatment.label)) { #get chr and midpoint" >> $SCRIPT
      echo "    tmp = read.table(file = output.bed.file[i], header = F)" >> $SCRIPT
      echo "    tmp = data.frame(tmp[, 1], floor(tmp[, 2]+(tmp[, 3]-tmp[, 2])/2))" >> $SCRIPT
      echo "    colnames(tmp) = c('chr', 'pos')" >> $SCRIPT
      echo "    peaks.bed.list[[treatment.label[i]]] = tmp" >> $SCRIPT
      echo "}" >> $SCRIPT
      
      echo "bed.data.list = list()" >> $SCRIPT
      echo "for (i in 1:length(treatment.label)) {" >> $SCRIPT
      echo "    bed.data.list[[treatment.label[i]]] = treatments.to.run[i]" >> $SCRIPT
      echo "}" >> $SCRIPT
      
      echo "input.data.list = list()" >> $SCRIPT
      echo "for (i in 1:length(treatment.label)) {" >> $SCRIPT
      echo "    input.data.list[[treatment.label[i]]] = input.to.run[i]" >> $SCRIPT
      echo "}" >> $SCRIPT
      
      echo "ana.conds = factor(c(rep(\"S1\", 2), rep(\"S2\", 2)), levels = c(\"S1\", \"S2\"))" >> $SCRIPT
      echo "ana.conds.label = \"S1.vs.S2\"" >> $SCRIPT
      
      echo "bs.list = read.binding.site.list(peaks.bed.list)" >> $SCRIPT
      echo "consensus.site = site.merge(bs.list)" >> $SCRIPT
      echo "nvec <- c($10,$11,$12,$13)" >> $SCRIPT
      echo "dat = load.data(chip.data.list = bed.data.list, input.data.list=input.data.list, conds = ana.conds, consensus.site = consensus.site, norm.factor.vec= nvec, data.type = \"BED\")" >> $SCRIPT
      
      echo "data.count = get.site.count(dat)" >> $SCRIPT
      echo "test.diff = test.diff.binding(data.count)" >> $SCRIPT
      echo "rept = report.peak(test.diff,FDR=1)" >> $SCRIPT
      
      #save as csv..
      echo "options(scipen = 999)" >> $SCRIPT
      echo "write.table(file=\"results.csv\", as.data.frame(rept),quote = FALSE)" >> $SCRIPT
      
      STARTTIME=`date +%s.%N`
      
      #run it...
      #R CMD BATCH --vanilla $SCRIPT
      /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_sf.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_sf $TIMEDIFF" >> time.txt

      MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
      echo $PSHORT"_"$MSHORT"_sf $MEMUSAGE" >> memory.txt
      
      if [ -e results.csv ]; then
          #reformat for eval
          cat results.csv | grep -v "^chr"  | awk '{print $2"\t"$3-100"\t"$3+100"\t"$9"\t"log($7)/log(2)}' | sort -k1,1 -k2,2n > $OUT_NAME
      else
        #create empty file
        touch $OUT_NAME
      fi
      
      #save log     
      cat script.Rout >> $LOG
      
      #clean up
      rm -f $SCRIPT script.Rout .RData results.csv s11_peaks.bed s12_peaks.bed s21_peaks.bed s22_peaks.bed mem.txt
      
    done
  done
  
  rm -f S11.bed S12.bed IN1.bed S21.bed S22.bed IN2.bed
 
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
