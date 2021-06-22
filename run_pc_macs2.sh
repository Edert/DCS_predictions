
NAME=$1
SET=$2
TOOL=$3

MODE=200
  
mkdir -p results_peaks/$TOOL results_peaks/$TOOL/$SET
mkdir -p log/pc log/pc/$TOOL

files=$(ls results_peaks/$TOOL/$SET/*/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running peak calling tool: $TOOL with: $1 $2 $3 " #$4 $5 $6 $7 $8 $9"

  LOG="log/pc/$TOOL/$SET.log"
  RESULTSDIR="results_peaks/$TOOL/$SET"
  
  #export PYTHONPATH=$PYTHONPATH:/apps/MACS2-2.1.0.20140616/lib/python2.7/site-packages/
  #active virtual environment for MACS2
  source "/proj/chipseq_norm_diffbind_062017/analysis/03_db_analysis/python-virtual-environments/macs2_env/bin/activate"

  #genome size chr19: 61431566  bed file: 26969981  --> diff = 34461585
  #cat data/mm10_chr19_repeatmasker.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
  NAME11=$(basename $4 _mm.bam)
  NAME12=$(basename $5 _mm.bam)
  NAME21=$(basename $7 _mm.bam)
  NAME22=$(basename $8 _mm.bam)
  
  #sharp model
  PCRUN=$RESULTSDIR"/sharp_m"
  mkdir -p $PCRUN
  
  echo "macs2 model for $1 $2"
  macs2 callpeak -t $4 -c $6 -f BAM -n $NAME11 --outdir $PCRUN -q 0.05  2>> $LOG
  macs2 callpeak -t $5 -c $6 -f BAM -n $NAME12 --outdir $PCRUN -q 0.05  2>> $LOG
  macs2 callpeak -t $7 -c $9 -f BAM -n $NAME21 --outdir $PCRUN -q 0.05  2>> $LOG
  macs2 callpeak -t $8 -c $9 -f BAM -n $NAME22 --outdir $PCRUN -q 0.05  2>> $LOG
  
  #-g 3.4e7
  #-g 3.4e7
  #-g 3.4e7
  #-g 3.4e7
  
  #move results
  cat $PCRUN/$NAME11"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s11_peaks.bed #chr start stop score
  cat $PCRUN/$NAME12"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s12_peaks.bed
  cat $PCRUN/$NAME21"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s21_peaks.bed
  cat $PCRUN/$NAME22"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s22_peaks.bed
  
  rm -rf $PCRUN/*_peaks.narrowPeak $PCRUN/*peaks.xls $PCRUN/*summits.bed $PCRUN/*_model.r


  #sharp no model
  PCRUN=$RESULTSDIR"/sharp_f"
  mkdir -p $PCRUN

  echo "macs2 nomodel for $1 $2"
  macs2 callpeak -t $4 -c $6 -f BAM -n $NAME11 --outdir $PCRUN -q 0.05  --nomodel --extsize $MODE 2>> $LOG
  macs2 callpeak -t $5 -c $6 -f BAM -n $NAME12 --outdir $PCRUN -q 0.05  --nomodel --extsize $MODE 2>> $LOG
  macs2 callpeak -t $7 -c $9 -f BAM -n $NAME21 --outdir $PCRUN -q 0.05  --nomodel --extsize $MODE 2>> $LOG
  macs2 callpeak -t $8 -c $9 -f BAM -n $NAME22 --outdir $PCRUN -q 0.05  --nomodel --extsize $MODE 2>> $LOG
  
  #move results
  cat $PCRUN/$NAME11"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s11_peaks.bed #chr start stop score
  cat $PCRUN/$NAME12"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s12_peaks.bed
  cat $PCRUN/$NAME21"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s21_peaks.bed
  cat $PCRUN/$NAME22"_peaks.narrowPeak" | cut -f1,2,3,5 > $PCRUN/s22_peaks.bed
  
  rm -rf $PCRUN/*_peaks.narrowPeak $PCRUN/*peaks.xls $PCRUN/*summits.bed $PCRUN/*_model.r
  
  
  #broad model
  PCRUN=$RESULTSDIR"/broad_m"
  mkdir -p $PCRUN
  
  echo "macs2 broad model for $1 $2"
  macs2 callpeak --broad -t $4 -c $6 -f BAM -n $NAME11 --outdir $PCRUN -q 0.05  2>> $LOG
  macs2 callpeak --broad -t $5 -c $6 -f BAM -n $NAME12 --outdir $PCRUN -q 0.05  2>> $LOG
  macs2 callpeak --broad -t $7 -c $9 -f BAM -n $NAME21 --outdir $PCRUN -q 0.05  2>> $LOG
  macs2 callpeak --broad -t $8 -c $9 -f BAM -n $NAME22 --outdir $PCRUN -q 0.05  2>> $LOG

  #move results
  cat $PCRUN/$NAME11"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s11_peaks.bed
  cat $PCRUN/$NAME12"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s12_peaks.bed
  cat $PCRUN/$NAME21"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s21_peaks.bed
  cat $PCRUN/$NAME22"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s22_peaks.bed
  
  rm -rf $PCRUN/*_peaks.broadPeak $PCRUN/*peaks.xls $PCRUN/*gappedPeak $PCRUN/*_model.r
  
  
  #broad no model
  PCRUN=$RESULTSDIR"/broad_f"
  mkdir -p $PCRUN
  
  echo "macs2 broad nomodel for $1 $2"
  macs2 callpeak --broad -t $4 -c $6 -f BAM -n $NAME11 --outdir $PCRUN -q 0.05 --nomodel --extsize $MODE 2>> $LOG
  macs2 callpeak --broad -t $5 -c $6 -f BAM -n $NAME12 --outdir $PCRUN -q 0.05 --nomodel --extsize $MODE 2>> $LOG
  macs2 callpeak --broad -t $7 -c $9 -f BAM -n $NAME21 --outdir $PCRUN -q 0.05 --nomodel --extsize $MODE 2>> $LOG
  macs2 callpeak --broad -t $8 -c $9 -f BAM -n $NAME22 --outdir $PCRUN -q 0.05 --nomodel --extsize $MODE 2>> $LOG

  #move results
  cat $PCRUN/$NAME11"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s11_peaks.bed
  cat $PCRUN/$NAME12"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s12_peaks.bed
  cat $PCRUN/$NAME21"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s21_peaks.bed
  cat $PCRUN/$NAME22"_peaks.broadPeak" | cut -f1,2,3,5 > $PCRUN/s22_peaks.bed
  
  rm -rf $PCRUN/*_peaks.broadPeak $PCRUN/*peaks.xls $PCRUN/*gappedPeak $PCRUN/*_model.r
  
  #leave virtual environment
  deactivate
  
else
  echo "results_peaks/$TOOL/$SET/bed already exists exiting..."
fi
