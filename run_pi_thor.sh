
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}"

  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/

  STARTTIME=`date +%s.%N`
  
  #create config
  echo "#rep1" > THOR.config
  echo "$4" >> THOR.config
  echo "$5" >> THOR.config
  echo "#rep2" >> THOR.config
  echo "$7" >> THOR.config
  echo "$8" >> THOR.config
  echo "#chrom_sizes" >> THOR.config
  echo "/ssd/references/THOR/mm10/chr19.chrom.size" >> THOR.config
  echo "#genome" >> THOR.config
  echo "/ssd/references/THOR/mm10/chr19.fa" >> THOR.config
  echo "#inputs1" >> THOR.config
  echo "$6" >> THOR.config
  echo "$6" >> THOR.config
  echo "#inputs2" >> THOR.config
  echo "$9" >> THOR.config
  echo "$9" >> THOR.config
  
  source /proj/chipseq_norm_diffbind_062017/analysis/03_db_analysis/python-virtual-environments/thor_env/bin/activate
  
  PREPDONE=`date +%s.%N`
  
  #run thor
  /usr/bin/time -o mem.txt -f "%K %M" rgt-THOR -n $NAME THOR.config > $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e $NAME"-diffpeaks.bed" ]; then
     #reformat for eval
     cat $NAME"-diffpeaks.bed" | awk '{split($11,a,";"); split(a[1],b,":"); split(a[2],c,":"); m1=(b[1]+b[2])/2; m2=(c[1]+c[2])/2; lfold=((m2==0) ? "inf" : log(m1/m2)/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME".bed"
  else
     #create empty file
     touch $OUT_NAME".bed"
  fi
  
  #clean up
  rm -f *.bw *.info mem.txt
  rm -f $NAME"-diffpeaks.bed" $NAME"-diffpeaks.narrowPeak" $NAME"-uncor-diffpeaks.bed" $NAME"-uncor-diffpeaks.narrowPeak"
  
  
  
  IN1=$(echo ${14} | awk '{printf "%.10f", $1-7}') #10 decimal places
  IN2=$(echo ${15} | awk '{printf "%.10f", $1-7}') #this is just to avoid an error caused by thor...
  SF1=$(echo ${10} | awk '{printf "%.10f", $1-7}')
  SF2=$(echo ${11} | awk '{printf "%.10f", $1-7}')
  SF3=$(echo ${12} | awk '{printf "%.10f", $1-7}')
  SF4=$(echo ${13} | awk '{printf "%.10f", $1-7}')
  
  STARTTIME=`date +%s.%N`

  #run with si norm factors
  /usr/bin/time -o mem.txt -f "%K %M" rgt-THOR -n $NAME"_si" --scaling-factors $SF1,$SF2,$SF3,$SF4 --factors-inputs $IN1,$IN1,$IN2,$IN2 THOR.config >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "si $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "si $MEMUSAGE" >> memory.txt
  
  if [ -e $NAME"_si-diffpeaks.bed" ]; then
     #reformat for eval
     cat $NAME"_si-diffpeaks.bed" | awk '{split($11,a,";"); split(a[1],b,":"); split(a[2],c,":"); m1=(b[1]+b[2])/2; m2=(c[1]+c[2])/2; lfold=((m2==0) ? "inf" : log(m1/m2)/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_si.bed"
  else
     #create empty file
     touch $OUT_NAME"_si.bed"
  fi
  
  deactivate
  
  #clean up
  rm -f THOR.config *.bw *.info mem.txt
  rm -f $NAME"_si-diffpeaks.bed" $NAME"_si-diffpeaks.narrowPeak" $NAME"_si-uncor-diffpeaks.bed" $NAME"_si-uncor-diffpeaks.narrowPeak"

else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
