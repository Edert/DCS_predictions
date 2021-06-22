
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
  ODIN_PATH="/proj/chipseq_norm_diffbind_062017/analysis/03_db_analysis/python-virtual-environments/reg-gen-ODIN-0.4.1-release/rgt/ODIN/ODIN.py"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`

  samtools merge s1.bam $4 $5 
  samtools merge s2.bam $7 $8 
  samtools index s1.bam
  samtools index s2.bam

  #add dm_chr2L to files..
  cat /ssd/references/THOR/mm10/chr19.fa > chr19_plus.fa
  cat /proj/chipseq_norm_diffbind_062017/analysis/01_simulating/data/dm6_chr2L.fasta | sed "s/chr2L/dm_chr2L/" >> chr19_plus.fa
  
  cat /ssd/references/THOR/mm10/chr19.chrom.size > chr19_plus.chrom.size
  echo -e "dm_chr2L\t23513712" >> chr19_plus.chrom.size
  
  #active virtual environment for ODIN
  source /proj/chipseq_norm_diffbind_062017/analysis/03_db_analysis/python-virtual-environments/odin_env/bin/activate

  PREPDONE=`date +%s.%N`
  
  #run odin binomial
  /usr/bin/time -o mem.txt -f "%K %M" python $ODIN_PATH -n $OUT_NAME -f 0.7 -b 100 -p 1.0 --input-1=$6 --input-2=$9 s1.bam s2.bam chr19_plus.fa chr19_plus.chrom.size >> $LOG 2>&1

  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME"-diffpeaks.bed" ]; then
     #reformat for eval
     cat $OUT_NAME"-diffpeaks.bed" | awk '{split($11,a,","); lfold=((a[2]==0) ? "inf" : log(a[1]/a[2])/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi
  
  #clean up
  rm -f *.bw *.info  mem.txt
  rm -f $OUT_NAME"-diffpeaks.bed" $OUT_NAME"-diffpeaks.narrowPeak" $OUT_NAME"-uncor-diffpeaks.bed" $OUT_NAME"-uncor-diffpeaks.narrowPeak"
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run odin poisson
  /usr/bin/time -o mem.txt -f "%K %M" python $ODIN_PATH -n $OUT_NAME --dist poisson -f 0.7 -b 100 -p 1.0 --input-1=$6 --input-2=$9 s1.bam s2.bam chr19_plus.fa chr19_plus.chrom.size >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME"-diffpeaks.bed" ]; then
     #reformat for eval
     cat $OUT_NAME"-diffpeaks.bed" | awk '{split($11,a,","); lfold=((a[2]==0) ? "inf" : log(a[1]/a[2])/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi
  
  #clean up
  rm -f *.bw *.info  mem.txt
  rm -f $OUT_NAME"-diffpeaks.bed" $OUT_NAME"-diffpeaks.narrowPeak" $OUT_NAME"-uncor-diffpeaks.bed" $OUT_NAME"-uncor-diffpeaks.narrowPeak"
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run odin poisson-c
  /usr/bin/time -o mem.txt -f "%K %M" python $ODIN_PATH -n $OUT_NAME --dist poisson-c -f 0.7 -b 100 -p 1.0 --input-1=$6 --input-2=$9 s1.bam s2.bam chr19_plus.fa chr19_plus.chrom.size >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME"-diffpeaks.bed" ]; then
     #reformat for eval
     cat $OUT_NAME"-diffpeaks.bed" | awk '{split($11,a,","); lfold=((a[2]==0) ? "inf" : log(a[1]/a[2])/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_3.bed"
  else
     #create empty file
     touch $OUT_NAME"_3.bed"
  fi
  
  #clean up
  rm -f *.bw *.info  mem.txt
  rm -f $OUT_NAME"-diffpeaks.bed" $OUT_NAME"-diffpeaks.narrowPeak" $OUT_NAME"-uncor-diffpeaks.bed" $OUT_NAME"-uncor-diffpeaks.narrowPeak"
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run odin binomial 1000 bp bin size
  /usr/bin/time -o mem.txt -f "%K %M" python $ODIN_PATH -n $OUT_NAME -f 0.7 -b 1000 -p 1.0 --input-1=$6 --input-2=$9 s1.bam s2.bam chr19_plus.fa chr19_plus.chrom.size >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "4 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME"-diffpeaks.bed" ]; then
     #reformat for eval
     cat $OUT_NAME"-diffpeaks.bed" | awk '{split($11,a,","); lfold=((a[2]==0) ? "inf" : log(a[1]/a[2])/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_4.bed"
  else
     #create empty file
     touch $OUT_NAME"_4.bed"
  fi
  
  #clean up
  rm -f *.bw *.info  mem.txt
  rm -f $OUT_NAME"-diffpeaks.bed" $OUT_NAME"-diffpeaks.narrowPeak" $OUT_NAME"-uncor-diffpeaks.bed" $OUT_NAME"-uncor-diffpeaks.narrowPeak"
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run odin constraint 1000 bp bin size
  /usr/bin/time -o mem.txt -f "%K %M" python $ODIN_PATH -n $OUT_NAME --dist poisson -f 0.7 -b 1000 -p 1.0 --input-1=$6 --input-2=$9 s1.bam s2.bam chr19_plus.fa chr19_plus.chrom.size >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "5 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "5 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME"-diffpeaks.bed" ]; then
     #reformat for eval
     cat $OUT_NAME"-diffpeaks.bed" | awk '{split($11,a,","); lfold=((a[2]==0) ? "inf" : log(a[1]/a[2])/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_5.bed"
  else
     #create empty file
     touch $OUT_NAME"_5.bed"
  fi
  
  #clean up
  rm -f *.bw *.info  mem.txt
  rm -f $OUT_NAME"-diffpeaks.bed" $OUT_NAME"-diffpeaks.narrowPeak" $OUT_NAME"-uncor-diffpeaks.bed" $OUT_NAME"-uncor-diffpeaks.narrowPeak"
  
  
  
  STARTTIME=`date +%s.%N`
  
  #run odin poisson-c 1000 bp bin size
  /usr/bin/time -o mem.txt -f "%K %M" python $ODIN_PATH -n $OUT_NAME --dist poisson-c -f 0.7 -b 1000 -p 1.0 --input-1=$6 --input-2=$9 s1.bam s2.bam chr19_plus.fa chr19_plus.chrom.size >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "6 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "6 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME"-diffpeaks.bed" ]; then
     #reformat for eval
     cat $OUT_NAME"-diffpeaks.bed" | awk '{split($11,a,","); lfold=((a[2]==0) ? "inf" : log(a[1]/a[2])/log(2)); print $1"\t"$2"\t"$3"\t"a[3]"\t"lfold}' | sort -k1,1 -k2,2n > $OUT_NAME"_6.bed"
  else
     #create empty file
     touch $OUT_NAME"_6.bed"
  fi
  
  #clean up
  rm -f *.bw *.info mem.txt
  rm -f $OUT_NAME"-diffpeaks.bed" $OUT_NAME"-diffpeaks.narrowPeak" $OUT_NAME"-uncor-diffpeaks.bed" $OUT_NAME"-uncor-diffpeaks.narrowPeak"
  
  #leave virtual environment
  deactivate
  
  #clean up
  rm chr19_plus.fa chr19_plus.fa.fai chr19_plus.chrom.size 
  rm -rf s1.bam s2.bam s1.bam.bai s2.bam.bai
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
