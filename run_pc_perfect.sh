
NAME=$1
SET=$2
TOOL=$3

mkdir -p results_peaks/$TOOL results_peaks/$TOOL/$SET
mkdir -p log/pc log/pc/$TOOL

files=$(ls results_peaks/$TOOL/$SET/*/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then
  
  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13"
  
  LOG="log/pc/$TOOL/$SET.log"
  RESULTSDIR="results_peaks/$TOOL/$SET"
  
  NAME=$(basename $4 _sample1-rep1_mm.bam)
  NAME11=$(basename $4 -rep1_mm.bam)
  NAME21=$(basename $7 -rep1_mm.bam)
  
  #instead of peak calling get total list of peaks from simulation
  PEAKbed="/proj/chipseq_norm_diffbind_062017/analysis/01_simulating/results/$SET/"$NAME"_all_peaks.bed"
  PEAKbed1="/proj/chipseq_norm_diffbind_062017/analysis/01_simulating/results/$SET/$NAME11"-peaks.bed
  PEAKbed2="/proj/chipseq_norm_diffbind_062017/analysis/01_simulating/results/$SET/$NAME21"-peaks.bed
  
  #all peaks
  PCRUN=$RESULTSDIR"/all_peaks"
  mkdir -p $PCRUN
  
  cat $PEAKbed | grep -v "^#" | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s11_peaks.bed #chr, start, stop, weight=1
  cat $PEAKbed | grep -v "^#" | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s12_peaks.bed
  cat $PEAKbed | grep -v "^#" | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s21_peaks.bed
  cat $PEAKbed | grep -v "^#" | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s22_peaks.bed
  
  
  #only db peaks
  PCRUN=$RESULTSDIR"/db_peaks"
  mkdir -p $PCRUN
  
  cat $PEAKbed1 | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s11_peaks.bed
  cat $PEAKbed1 | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s12_peaks.bed
  cat $PEAKbed2 | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s21_peaks.bed
  cat $PEAKbed2 | awk '{print $1"\t"$2"\t"$3"\t"1}' > $PCRUN/s22_peaks.bed
  
else
  echo "results_peaks/$TOOL/$SET/bed already exists exiting..."
fi
