
NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15"

  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  samtools merge s1.bam $4 $5 
  samtools merge s2.bam $7 $8 
  samtools index s1.bam
  samtools index s2.bam
  
  PREPDONE=`date +%s.%N`
  
  #bin 1000
  /usr/bin/time -o mem.txt -f "%K %M" /home/eder/R/x86_64-pc-linux-gnu-library/3.5/histoneHMM/bin/histoneHMM_call_regions.R -b 1000 -c /ssd/references/THOR/mm10/chr19.chrom.size -o S1 s1.bam >> $LOG 2>&1
  /usr/bin/time -a -o mem.txt -f "%K %M" /home/eder/R/x86_64-pc-linux-gnu-library/3.5/histoneHMM/bin/histoneHMM_call_regions.R -b 1000 -c /ssd/references/THOR/mm10/chr19.chrom.size -o S2 s2.bam >> $LOG 2>&1
  
  /usr/bin/time -a -o mem.txt -f "%K %M" /home/eder/R/x86_64-pc-linux-gnu-library/3.5/histoneHMM/bin/histoneHMM_call_differential.R S1.txt S2.txt >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  #get mean of mean and max of max mem usage
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt | awk 'BEGIN { max = "NaN" }{ sum += $1; max = (NR==1 || $2>max ? $2 : max) } END { if (NR > 0) print sum1/NR " " max}' )
  echo "1 $MEMUSAGE" >> memory.txt
      
  if [ -e histoneHMM_differential/sample1-vs-sample2.txt ]; then
     #reformat for eval
     cat histoneHMM_differential/sample1-vs-sample2.txt | grep -v "^chrom" | grep "sample1" | awk '{print $1"\t"$2"\t"$3"\t"1-$10"\t1"}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
     cat histoneHMM_differential/sample1-vs-sample2.txt | grep -v "^chrom" | grep "sample2" | awk '{print $1"\t"$2"\t"$3"\t"1-$11"\t-1"}' | sort -k1,1 -k2,2n >> $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi

  #clean up
  rm -f S1* S2* mem.txt
  rm -rf histoneHMM_differential
  
  
  STARTTIME=`date +%s.%N`
  
  #bin 100
  /usr/bin/time -o mem.txt -f "%K %M" /home/eder/R/x86_64-pc-linux-gnu-library/3.5/histoneHMM/bin/histoneHMM_call_regions.R -b 100 -c /ssd/references/THOR/mm10/chr19.chrom.size -o S1 s1.bam >> $LOG 2>&1
  /usr/bin/time -a -o mem.txt -f "%K %M" /home/eder/R/x86_64-pc-linux-gnu-library/3.5/histoneHMM/bin/histoneHMM_call_regions.R -b 100 -c /ssd/references/THOR/mm10/chr19.chrom.size -o S2 s2.bam >> $LOG 2>&1
   
  /usr/bin/time -a -o mem.txt -f "%K %M" /home/eder/R/x86_64-pc-linux-gnu-library/3.5/histoneHMM/bin/histoneHMM_call_differential.R S1.txt S2.txt >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt
  
  #get mean of mean and max of max mem usage
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt | awk 'BEGIN { max = "NaN" }{ sum += $1; max = (NR==1 || $2>max ? $2 : max) } END { if (NR > 0) print sum1/NR " " max}' )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e histoneHMM_differential/sample1-vs-sample2.txt ]; then
     #reformat for eval
     cat histoneHMM_differential/sample1-vs-sample2.txt | grep -v "^chrom" | grep "sample1" | awk '{print $1"\t"$2"\t"$3"\t"1-$10"\t1"}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
     cat histoneHMM_differential/sample1-vs-sample2.txt | grep -v "^chrom" | grep "sample2" | awk '{print $1"\t"$2"\t"$3"\t"1-$11"\t-1"}' | sort -k1,1 -k2,2n >> $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi

  #clean up
  rm -f S1* S2* mem.txt
  rm -rf histoneHMM_differential
  
  rm -rf s1.bam s2.bam s1.bam.bai s2.bam.bai
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
