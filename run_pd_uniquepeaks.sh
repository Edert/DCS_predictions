
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
  
  PREPDONE=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  
  for PCALLER in ../../../results_peaks/*; do
  
    for PMODE in $PCALLER/$SET/*; do
    
      echo "using $(basename $PCALLER) mode: $(basename $PMODE)"
      
      cut -f1,2,3 $PMODE/s11_peaks.bed > s1_peaks.bed
      cut -f1,2,3 $PMODE/s12_peaks.bed >> s1_peaks.bed
      cut -f1,2,3 $PMODE/s21_peaks.bed > s2_peaks.bed
      cut -f1,2,3 $PMODE/s22_peaks.bed >> s2_peaks.bed
      
      #merge
      cat s1_peaks.bed | sort -k1,1 -k2,2n > tmp1
      cat s2_peaks.bed | sort -k1,1 -k2,2n > tmp2
      
      bedtools merge -i tmp1 > s1_peaks.bed
      bedtools merge -i tmp2 > s2_peaks.bed
      
      rm tmp1 tmp2
      
      STARTTIME=`date +%s.%N`
      
      #compare
      /usr/bin/time -o mem.txt -f "%K %M" bedtools subtract -a s1_peaks.bed -b s2_peaks.bed > results1.bed
      /usr/bin/time -a -o mem.txt -f "%K %M" bedtools subtract -b s1_peaks.bed -a s2_peaks.bed > results2.bed
      
      #save result
      PSHORT=$(basename $PCALLER)
      MSHORT=$(basename $PMODE)
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT".bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT" $TIMEDIFF" >> time.txt
      
      #get mean of mean and max of max mem usage
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt | awk 'BEGIN { max = "NaN" }{ sum += $1; max = (NR==1 || $2>max ? $2 : max) } END { if (NR > 0) print sum1/NR " " max}' )
      echo $PSHORT"_"$MSHORT" $MEMUSAGE" >> memory.txt
      
      if [ -e results1.bed ]; then
        #reformat for eval
        cat results1.bed | awk '{print $1"\t"$2"\t"$3"\t0\t0.7"}' | sort -k1,1 -k2,2n > $OUT_NAME
        cat results2.bed | awk '{print $1"\t"$2"\t"$3"\t0\t-0.7"}' | sort -k1,1 -k2,2n >> $OUT_NAME
      else
        #create empty file
        touch $OUT_NAME
      fi
      
      
      STARTTIME=`date +%s.%N`
      
      #compare stringend (remove complete entry if there is at least 1 bp overlap)
      /usr/bin/time -o mem.txt -f "%K %M" bedtools subtract -a s1_peaks.bed -b s2_peaks.bed -A > results1.bed
      /usr/bin/time -a -o mem.txt -f "%K %M" bedtools subtract -b s1_peaks.bed -a s2_peaks.bed -A > results2.bed
      
      OUT_NAME=$(basename $NAME _sample1-rep1_mm)"_"$PSHORT"_"$MSHORT"_s.bed"
      
      ENDTIME=`date +%s.%N`
      TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
      echo $PSHORT"_"$MSHORT"_s $TIMEDIFF" >> time.txt
      
      #get mean of mean and max mem usage
      MEMUSAGE=$(sed '/non-zero status/d' mem.txt | awk '{ sum1 += $1; sum2 += $2 } END { if (NR > 0) print sum1/NR " " sum2/NR}' )
      echo $PSHORT"_"$MSHORT"_s $MEMUSAGE" >> memory.txt
      
      if [ -e results1.bed ]; then
        #reformat for eval
        cat results1.bed | awk '{print $1"\t"$2"\t"$3"\t0\t0.7"}' | sort -k1,1 -k2,2n > $OUT_NAME
        cat results2.bed | awk '{print $1"\t"$2"\t"$3"\t0\t-0.7"}' | sort -k1,1 -k2,2n >> $OUT_NAME
      else
        #create empty file
        touch $OUT_NAME
      fi
      
      #clean up
      rm -f results1.bed results2.bed s1_peaks.bed s2_peaks.bed mem.txt
      
    done
  done
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
